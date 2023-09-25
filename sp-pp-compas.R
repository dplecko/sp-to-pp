root <- rprojroot::find_root(rprojroot::is_git_root)

library(fairadapt)
library(faircause)
library(ranger)
library(ggplot2)
library(latex2exp)
library(data.table)

# helpers
ts_boot <- function(x, y) {
  
  nrep <- 10000
  x_t <- rnorm(nrep, mean = mean(x), sd = sd(x))
  y_t <- rnorm(nrep, mean = mean(y), sd = sd(y))

  min(mean(x_t < y_t), 1 - mean(x_t < y_t))
}

compute_ippm <- function(pred, out, x, x0, name, nboot = 100) {

  breaks <- quantile(pred, probs = seq(0.1, 0.9, 0.1))
  bins <- .bincode(pred, c(-Inf, breaks, Inf))

  idx <- x == x0

  bin_diffs <- tapply(out[!idx], bins[!idx], mean) -
    tapply(out[idx], bins[idx], mean)

  ippm <- vapply(
    seq_len(nboot),
    function(i) {
      boot <- sample.int(length(pred), replace = TRUE)

      out_b <- out[boot]
      idx_b <- idx[boot]
      bins_b <- bins[boot]

      bin_diffs <- tapply(out_b[!idx_b], bins_b[!idx_b], mean) -
        tapply(out_b[idx_b], bins_b[idx_b], mean)

      mean(bin_diffs)
    }, numeric(1L)
  )

  data.frame(value = mean(ippm), sd = sd(ippm), measure = "ippm",
             outcome = name)
}

compute_roc <- function(probs, labs, nboot = 100) {
  
  aucs <- vapply(
    seq_len(nboot),
    function(i) {
      
      bts <- sample.int(length(probs), replace = TRUE)
      auc <- PRROC::roc.curve(scores.class0 = probs[bts], 
                              weights.class0 = labs[bts])$auc
      ifelse(auc >= 0.5, auc, 1 - auc)
    }, numeric(1L)
  )
  
  c(mean(aucs), sd(aucs))
}

compute_tv <- function(out, x, x0, name, nboot = 100) {
  
  idx <- x == x0
  tv <- vapply(
    seq_len(nboot),
    function(i) {
      
      boot <- sample.int(length(out), replace = TRUE)
      out_b <- out[boot]
      idx_b <- idx[boot]
      mean(out_b[idx_b]) - mean(out_b[!idx_b])
    }, numeric(1L)
  )
  data.frame(value = mean(tv), sd = sd(tv), measure = "tv", outcome = name)
}

set.seed(23)
data <- read.csv(file.path(root, "compas-scores-two-years.csv"))
col.keep <- which(
  names(data) %in% c("age", "sex", "juv_fel_count",
                     "juv_misd_count", "juv_other_count", "priors_count",
                     "c_charge_degree", "race", "two_year_recid", "decile_score")
)
data <- data[, col.keep]
data$race <- factor(data$race)
levels(data$race) <- c("Minority", "Minority", "Majority", "Minority",
                       "Minority", "Minority")
data$race <- relevel(data$race, "Majority")
cumsum(table(data$decile_score)) / sum(table(data$decile_score))
## decile_score > 4 represents high risk (approximately)
data$northpointe <- as.integer(data$decile_score > 4)
np_ippm <- compute_ippm(data$decile_score / 10, data$two_year_recid,
                         x = data$race, x0 = "Majority", "northpointe")
data$decile_score <- NULL
names(data) <- gsub("_count", "", names(data))
names(data)[which(names(data) == "c_charge_degree")] <- "charge"

# decomposing the true outcome
X <- "race"
Z <- c("age", "sex")
W <- c("juv_fel", "juv_misd", "juv_other", "priors", "charge")
Y <- c("two_year_recid")
two_year <- fairness_cookbook(data, X = X, W = W, Z = Z, Y = Y,
                              x0 = "Majority", x1 = "Minority", nboot1 = 5)

# decomposing the Northpointe predictions
Yhat <- "northpointe"
northpointe <- fairness_cookbook(data, X = X, W = W, Z = Z, Y = Yhat,
                                 x0 = "Majority", x1 = "Minority", nboot1 = 5)

# obtain fair predictions

# get adjacency matrix
vars <- c("age", "sex", "juv_fel","juv_misd",
          "juv_other", "priors",
          "charge", "race", "two_year_recid")

adj.mat <- array(0, dim = c(length(vars), length(vars)))
colnames(adj.mat) <- vars
rownames(adj.mat) <- colnames(adj.mat)

# adding the edges to the matrix
adj.mat[
  c("race", "sex", "age"), c("juv_fel", "juv_misd",
                             "juv_other", "priors",
                             "charge", "two_year_recid")] <- 1
adj.mat[c("juv_fel", "juv_misd", "juv_other"),
        c("priors", "charge", "two_year_recid")] <- 1
adj.mat["priors", c("charge", "two_year_recid")] <- 1
adj.mat["charge", "two_year_recid"] <- 1

cfd.mat <- adj.mat
cfd.mat[, ] <- 0L
cfd.mat[c("sex", "age"), "race"] <- cfd.mat["race", c("sex", "age")] <- 1L

fdp <- fairadapt::fairadapt(two_year_recid ~ ., prot.attr = "race",
                            train.data = data[, vars], adj.mat = adj.mat)

# obtain the adapted data
ad_dat <- fairadapt:::adaptedData(fdp)
ad_dat$race <- data$race

# obtain predictions based on the adapted data
adapt_rf <- ranger(two_year_recid ~ ., ad_dat,
                   probability = TRUE)
adapt_oob_prob <- adapt_rf$predictions[, 2]
fp_ippm <- compute_ippm(adapt_oob_prob, data$two_year_recid,
                        data$race, "Majority", "fairadapt")

# match the prevalence with NP and true
adapt_rf_class <- ranger(two_year_recid ~ ., ad_dat,
                         classification = TRUE)
adapt_oob <- adapt_rf_class$predictions

# recover probabilities too
ad_dat$two_year_recid <- adapt_oob

# decomposing the fair predictor
fairadapt_tvd <- fairness_cookbook(
  ad_dat, X = X, W = W, Z = Z, Y = "two_year_recid",
  x0 = "Majority", x1 = "Minority", nboot1 = 5
)

# add SP and PP predictors

# fit an unconstrained predictor using random forest
rf <- ranger(two_year_recid ~ ., data = data[, vars],
             probability = TRUE)

p_oob <- rf$predictions[, 2]

# SP predictor
th_min <- quantile(p_oob[data$race == "Minority"], 
                   probs = 1 - mean(data$northpointe))
th_maj <- quantile(p_oob[data$race == "Majority"], 
                   probs = 1 - mean(data$northpointe))
y_sp <- ifelse(data$race == "Minority", p_oob > th_min, p_oob > th_maj)
data$y_sp <- as.integer(y_sp)

sp_tvd <- fairness_cookbook(
  data, X = X, W = W, Z = Z, Y = "y_sp",
  x0 = "Majority", x1 = "Minority", nboot1 = 5
)
sp_ippm <- compute_ippm(data$y_sp, data$two_year_recid, data$race,
                        "Majority", "SP")

# PP predictor
th <- quantile(p_oob, probs = 1 - mean(data$northpointe))
data$y_pp <- as.integer(p_oob > th)

pp_tvd <- fairness_cookbook(
  data, X = X, W = W, Z = Z, Y = "y_pp",
  x0 = "Majority", x1 = "Minority", nboot1 = 5
)
pp_ippm <- compute_ippm(data$y_pp, data$two_year_recid, data$race,
                        "Majority", "PP")

# hypothesis testing for differences

# NP vs. true
ts_boot(
  subset(two_year$measures, measure == "ctfse")$value,
  subset(northpointe$measures, measure == "ctfse")$value
)

# FP vs. true
ts_boot(
  subset(two_year$measures, measure == "ctfse")$value,
  subset(fairadapt_tvd$measures, measure == "ctfse")$value
)

# side-by-side plot of three decompositions
res <- rbind(
  cbind(summary(two_year)$measures, outcome = "true"),
  cbind(summary(northpointe)$measures, outcome = "northpointe"),
  cbind(summary(fairadapt_tvd)$measures, outcome = "fairadapt"),
  cbind(summary(sp_tvd)$measures, outcome = "SP"),
  cbind(summary(pp_tvd)$measures, outcome = "PP"),
  np_ippm,
  fp_ippm,
  sp_ippm,
  pp_ippm
)


res <- res[res$measure %in% c("ctfde", "ctfie", "ctfse", "tv", "ippm"), ]

ggplot(
  res,
  aes(x = factor(measure, levels = c("ctfde", "ctfie", "ctfse", "tv", "ippm")),
      y = value,
      fill = factor(
        outcome, 
        levels = c("true", "northpointe", "fairadapt", "SP", "PP")
      )
     )
  ) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd),
                width=.2, position=position_dodge(.9)) +
  scale_fill_discrete(
    name = "Outcome",
    labels = c(TeX("$Y^{true}$"), TeX("$\\hat{Y}^{NP}$"),
               TeX("$\\hat{Y}^{FP}$"), TeX("$\\hat{Y}^{SP}$"),
               TeX("$\\hat{Y}^{PP}$"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 14),
    axis.title =  element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    title = element_text(size = 16)
  ) +
  xlab("Fairness Measure") + ylab("Value") +
  scale_x_discrete(labels = c(TeX("Ctf-DE"),
                              TeX("Ctf-IE"),
                              TeX("Ctf-SE"),
                              TeX("SPM"),
                              TeX("iPPM"))) +
  scale_y_continuous(labels = scales::percent) # +
  # ggtitle("Fairness Measures on the COMPAS dataset")
ggsave("compas-algorithm-1-revised.png", width = 7, height = 4)


# Pareto plots
cases <- list(empty = NULL, Z = c(Z), W = c(W), ZW = c(Z, W))
res <- NULL
nrep <- 10
for (seed in seq_len(nrep)) {
  
  set.seed(seed)
  
  for (i in seq_along(cases)) {
    
    res.vars <- cases[[i]] 
    adj.curr <- adj.mat
    if (!is.element("age", res.vars)) adj.curr["sex", Z] <- 1
    
    fdp <- fairadapt::fairadapt(two_year_recid ~ ., prot.attr = "race",
                                train.data = data[, vars], adj.mat = adj.mat,
                                res.vars = res.vars)
    
    # obtain the adapted data
    ad_dat <- fairadapt:::adaptedData(fdp)
    ad_dat$race <- data$race
    
    # obtain predictions based on the adapted data
    adapt_rf <- ranger(two_year_recid ~ ., ad_dat,
                       probability = TRUE)
    adapt_oob_prob <- adapt_rf$predictions[, 2]
    fp_ippm <- compute_ippm(adapt_oob_prob, data$two_year_recid,
                            data$race, "Majority", names(cases)[i])
    fp_ippm$seed <- seed
    
    auc_vals <- compute_roc(adapt_oob_prob, data$two_year_recid)
    fp_auc <- data.frame(value = auc_vals[1], sd = auc_vals[2], measure = "auc", 
                         outcome = names(cases)[i], seed = seed)
    
    # match the prevalence with NP and true
    adapt_rf_class <- ranger(two_year_recid ~ ., ad_dat,
                             classification = TRUE)
    adapt_oob <- as.integer(adapt_oob_prob > quantile(adapt_oob_prob, probs = 0.68)) # adapt_rf_class$predictions
    cat(names(cases)[i], ": prevalence =", mean(adapt_oob), "\n")
    fp_tv <- compute_tv(adapt_oob, data$race, "Minority", names(cases)[i])
    fp_tv$seed <- seed
    
    res <- rbind(res, fp_tv, fp_ippm, fp_auc)
  }
}

combine_musd <- function(mus, sds) {
  
  mixture <- rnorm(10000, mean = mus, sd = sds)
  list(value = mean(mixture), sd = sd(mixture))
}

res <- as.data.table(res)
res <- res[, combine_musd(value, sd), by = c("measure", "outcome")]

paretoA <- merge(res[res$measure == "tv", ], res[res$measure == "ippm", ],
                 by = "outcome")
pp_vs_sp <- ggplot(paretoA, aes(x = value.x, y = value.y, label = outcome, 
                                color = outcome)) +
  geom_errorbar(aes(ymin = value.y - sd.y, ymax = value.y + sd.y),
                linewidth = 1.5, width = 0.01) +
  geom_errorbarh(aes(xmin = value.x - sd.x, xmax = value.x + sd.x),
                 linewidth = 1.5, height = 0.006) +
  geom_point() + theme_bw() + xlab("Statistical Parity Measure (SPM)") +
  ylab("Integrated PPM") +
  scale_color_discrete(
    name = "BN set",
    labels = c(TeX("$\\oslash$"),
               TeX("$W$"),
               TeX("$Z$"),
               TeX("$\\{ Z, W \\}$")),
  ) +
  theme(
    legend.position = "bottom", 
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )


paretoB <- merge(res[res$measure == "tv", ], res[res$measure == "auc", ],
                 by = "outcome")
auc_vs_sp <- ggplot(paretoB, aes(x = value.x, y = value.y, label = outcome, 
                                 color = outcome)) +
  geom_errorbar(aes(ymin = value.y - sd.y, ymax = value.y + sd.y),
                linewidth = 1.5, width = 0.01) +
  geom_errorbarh(aes(xmin = value.x - sd.x, xmax = value.x + sd.x),
                 linewidth = 1.5, height = 0.0015) +
  geom_point() + theme_bw() + xlab("Statistical Parity Measure (SPM)") +
  ylab("Area Under ROC") +
  scale_color_discrete(
    name = "BN set",
    labels = c(TeX("$\\oslash$"),
               TeX("$W$"),
               TeX("$Z$"),
               TeX("$\\{ Z, W \\}$")),
  ) +
  theme(
    legend.position = "bottom", 
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )

# Remove x-axis title and legend from pp_vs_sp
pp_vs_sp <- pp_vs_sp + theme(axis.title.x = element_blank(), legend.position = "none")

# Remove x-axis title from auc_vs_sp
auc_vs_sp <- auc_vs_sp + theme(axis.title.x = element_blank())

# Combine plots one below the other
combined_plot <- cowplot::plot_grid(pp_vs_sp, auc_vs_sp, ncol = 1)

# Create drawing canvas
final_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(combined_plot, y = 0.05, height = 0.95) +
  cowplot::draw_label("Statistical Parity Measure (SPM)", x = 0.53, y = 0.025)

# Display final plot
final_plot

ggsave("pareto-plot-with-auc.png", width = 5 * 1.6, height = 4 * 1.6)

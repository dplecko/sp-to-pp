
library(ranger)
library(ggplot2)
library(fairadapt)
library(data.table)
library(xgboost)

decomp_ppm_byn <- function(n, type = "A1", seed = 23) {
  
  set.seed(seed)
  
  x <- rbinom(n, 1, 0.5)
  u1 <- rnorm(n)
  u2 <- rnorm(n)
  uy <- rnorm(n)
  
  eval_scm <- function(x, u1, u2, uy) {
    
    alpha <- 1
    beta <- 0.5
    
    w1 <- alpha * x + u1
    w2 <- w1^2 / 4 + u2 - x / 3
    y <- 1 / 6 * w1 * w2 + w1 + beta * x + uy
    
    data.frame(X = x, W1 = w1, W2 = w2, Y = y)
  }
  
  data <- eval_scm(x, u1, u2, uy)
  xgbcv <- xgb.cv(params = list(eta = 0.1), data = as.matrix(data[, 1:3]), 
                  label = data$Y, nrounds = 100, early_stopping_rounds = 3,
                  nfold = 10, verbose = FALSE)
  nrounds <- xgbcv$best_iteration
  xgb <- xgboost(params = list(eta = 0.1), data = as.matrix(data[, 1:3]), 
                 label = data$Y, nrounds = nrounds, verbose = FALSE)
  
  # construct Yhat (xgboost)
  data$Yhat <- predict(xgb, as.matrix(data[, 1:3]))
  
  # fix a decile of Yhat
  width <- 1 / (n / 50)
  d7 <- quantile(data$Yhat, probs = c(0.6 - width, 0.6 + width))
  idx1_d7 <- which(data$Yhat < d7[2] & data$Yhat >= d7[1] & data$X == 1)
  idx0_d7 <- which(data$Yhat < d7[2] & data$Yhat >= d7[1] & data$X == 0)
  
  # compute Y_{x0}
  # get adjacency matrix
  vars <- c("X", "W1", "W2", "Y")
  
  adj.mat <- array(0, dim = c(length(vars), length(vars)))
  colnames(adj.mat) <- vars
  rownames(adj.mat) <- colnames(adj.mat)
  
  # adding the edges to the matrix
  adj.mat[c("X"), c("W1", "W2", "Y")] <- adj.mat[c("W1"), c("W2", "Y")] <- 
    adj.mat[c("W2"), c("Y")] <- 1
  
  fa <- fairadapt(Y ~ ., prot.attr = "X", adj.mat = adj.mat,
                  train.data = data[, vars])
  Y_x0 <- adaptedData(fa)$Y
  Ys_x0 <- eval_scm(0, u1, u2, uy)$Y
  Y <- data$Y
  
  termI <- mean(Y[idx1_d7] - Y_x0[idx1_d7])
  termII <- mean(Y_x0[idx1_d7]) - mean(Y_x0[idx0_d7])
  ppm <- termI + termII
  
  terms <- list()
  if (type == "A2") bins <- 20 else bins <- 20
  for (i in seq_len(bins)) {
    
    d7 <- quantile(data$Yhat, probs = c(i-1, i) / bins)
    idx1_d7 <- which(data$Yhat < d7[2] & data$Yhat >= d7[1] & data$X == 1)
    idx0_d7 <- which(data$Yhat < d7[2] & data$Yhat >= d7[1] & data$X == 0)
    
    # gtruth
    termI_si <- mean(Y[idx1_d7] - Ys_x0[idx1_d7])
    
    # estimated
    termI_i <- mean(Y[idx1_d7] - Y_x0[idx1_d7])
    
    assertthat::assert_that(all(Y_x0[idx0_d7] == Y[idx0_d7]))
    termII_i <- mean(Y_x0[idx1_d7]) - mean(Y_x0[idx0_d7])
    ppm_i <- termI_i + termII_i
    
    # save and plot
    terms[[i]] <- c(termI_i, termI_si, termII_i, ppm_i)
  }
  
  if (type == "A1") {
    
    res <- as.data.frame(do.call(rbind, terms))
    names(res) <- c("termI", "termIs", "termII", "ppm")
    res$dec <- seq_len(nrow(res))
    res$seed <- seed
    res$n <- n
    return(res)
  } else if (type == "A2") {
    
    res <- as.data.frame(do.call(rbind, terms))[, c(1, 2)]
    res$dec <- seq_len(nrow(res))
    res$seed <- seed
    return(res)
  }
}

# experiment A1
n_grid <- c(1000, 2000, 3500, 5000, 7500, 10000)
nrep <- 20
a1res <-  NULL
for (n in n_grid) {
  
  for (i in seq_len(nrep)) {
    
    a1res <- rbind(a1res, decomp_ppm_byn(n, "A1", seed = i))
  }
}
a1res <- as.data.table(a1res)

p1 <- ggplot(
  melt(a1res, id.vars = c("n", "dec", "seed"))[, 
    mean(value, na.rm = TRUE), by = c("n", "variable", "seed")][ 
         variable != "termIs", list(mean(V1), sd(V1)), by = c("n", "variable")],
  aes(x = log(n), y = V1, color = variable)
) +
  geom_point() + geom_line() +
  geom_ribbon(aes(ymin = V1 - V2, ymax = V1 + V2, fill = variable), alpha = 0.4) +
  theme_bw() + ylab("Value") +
  scale_color_discrete(name = "Quantity", labels = c("iTerm I", "iTerm II", "iPPM")) +
  scale_fill_discrete(name = "Quantity", labels = c("iTerm I", "iTerm II", "iPPM")) +
  theme(legend.position = "bottom")

ggsave("pp-decomp-T12.png", plot = p1, width = 6, height = 4)
  
# experiment A2
a2res <- NULL
for (seed in seq_len(50)) {
  
  a2res <- rbind(a2res, decomp_ppm_byn(5000, "A2", seed = seed))
}

a2res <- as.data.table(a2res)

p2 <- ggplot(
  a2res[complete.cases(a2res), 
        list(termI_hat = mean(V1), sdev = sd(V1), termI_s = mean(V2)), 
        by = c("dec")],
  aes(x = dec, y = termI_hat)
) +
  geom_line(aes(color = "Estimated")) + 
  geom_ribbon(aes(ymin = termI_hat - sdev, ymax = termI_hat + sdev), alpha = 0.4) +
  geom_point(aes(x = dec, y = termI_s, color = "True"), size = 2) + theme_bw() +
  xlab("Decile") + ylab("Term (I) value") +
  scale_x_continuous(breaks = 1:20) + 
  scale_color_manual(name='Value',
                     breaks=c("Estimated", "True"),
                     values=c("Estimated" = "black", "True" = "red")) +
  theme(
    legend.position = "bottom"
  )

ggsave("pp-decomp-T1truth.png", plot = p2, width = 4 / 3 * 6, height = 4 / 3 * 4)

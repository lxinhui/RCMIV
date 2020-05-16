#' IV approach when X&Y have independent measurement error with continuous IV
#' @name rcmiv
#' @param g A numeric vector of instrumental variable.
#' @param x1 A numeric vector of the first measurement of X
#' @param x2 A numeric vector of the second measurement of X
#' @param y1 A numeric vector of the first measurement of Y
#' @param y2 A numeric vector of the second measurement of Y
#' @return A vector of causal effect estimate, bootstrap SE and 95% CI
#' @author Xinhui Liu

rcmiv <- function(g, x1, x2, y1, y2) {
  beta_gy1 <- summary(lm(y1 ~ g)) $ coef[2, 1]
  beta_gy2 <- summary(lm(y2 ~ g)) $ coef[2, 1]
  beta_x1x2 <- summary(lm(x2 ~ x1)) $ coef[2, 1]
  beta_y1y2 <- summary(lm(y2 ~ y1)) $ coef[2, 1]
  beta_gx1 <- summary(lm(x1 ~ g)) $ coef[2, 1]
  beta_gx2 <- summary(lm(x2 ~ g)) $ coef[2, 1]
  if (beta_gy1 / beta_gx1 > 0) {
    sign <- 1
  }
  if (beta_gy1 / beta_gx1 < 0) {
    sign <- (-1)
  }
  b <- sign * (sqrt((beta_gy1 * beta_x1x2 * beta_gy2) / (beta_gx1 *
                    beta_y1y2 * beta_gx2)))
  data <- data.frame(g, x1, x2, y1, y2)
  colnames(data) <- c("g", "x1", "x2", "y1", "y2")
  b_es <- rep(0, 2000)
  for (i in 1:2000) {
    nr <- nrow(data)
    select <- sample(c(1:nr), nr, replace = T)
    data_use <- data[select, ]
    g <- data_use $ g
    x1 <- data_use $ x1
    x2 <- data_use $ x2
    y1 <- data_use $ y1
    y2 <- data_use $ y2
    beta_gy1 <- summary(lm(y1 ~ g)) $ coef[2, 1]
    beta_gy2 <- summary(lm(y2 ~ g)) $ coef[2, 1]
    beta_x1x2 <- summary(lm(x2 ~ x1)) $ coef[2, 1]
    beta_y1y2 <- summary(lm(y2 ~ y1)) $ coef[2, 1]
    beta_gx1 <- summary(lm(x1 ~ g)) $ coef[2, 1]
    beta_gx2 <- summary(lm(x2 ~ g)) $ coef[2, 1]
    if (beta_gy1 / beta_gx1 > 0) {
      sign <- 1
    }
    if (beta_gy1 / beta_gx1 < 0) {
      sign <- (-1)
    }
    b_es[i] <- sign * (sqrt((beta_gy1 * beta_x1x2 * beta_gy2) / (beta_gx1 *
                  beta_y1y2 * beta_gx2)))
  }
  ### SE\CI###
  b_es <- na.omit(b_es)
  lengthbs <- length(b_es)
  for (j in 1 : lengthbs) {
    if (abs(b_es[j]) > 1) {
      b_es[j] <- NA
    }
  }
  b_es <- na.omit(b_es)
  se <- round(sd(b_es), 3)
  cin <- paste0("(", round((b - 1.96 * se), 3), ",",
                round((b + 1.96 * se), 3), ")")
  ciq <- paste0("(", round(quantile(b_es, probs = 0.025), 3), ",",
                round(quantile(b_es, probs = 0.975), 3), ")")
  result <- list(b, se, cin, ciq)
  names(result) <- c("Estimate", "SE", "95% CIn", "95% CIq")
  return(result)
}


#' IV approach when X&Y have dependent measurement error with continuous IV
#' @name rcmivcor
#' @param g A numeric vector of instrumental variable.
#' @param x1 A numeric vector of the first measurement of X
#' @param x2 A numeric vector of the second measurement of X
#' @param y1 A numeric vector of the first measurement of Y
#' @param y2 A numeric vector of the second measurement of Y
#' @return A vector of causal effect estimate, bootstrap SE and 95% CI
#' @author Xinhui Liu

rcmivcor <- function(g, x1, x2, y1, y2) {
  beta_gy1 <- summary(lm(y1 ~ g)) $ coef[2, 1]
  beta_gy2 <- summary(lm(y2 ~ g)) $ coef[2, 1]
  beta_x1x2 <- summary(lm(x2 ~ x1)) $ coef[2, 1]
  beta_y1y2 <- summary(lm(y2 ~ y1)) $ coef[2, 1]
  beta_gx1 <- summary(lm(x1 ~ g)) $ coef[2, 1]
  beta_gx2 <- summary(lm(x2 ~ g)) $ coef[2, 1]
  beta_x1y1 <- summary(lm(y1 ~ x1)) $ coef[2, 1]
  beta_x2y2 <- summary(lm(y2 ~ x2)) $ coef[2, 1]
  k2 <- (beta_x1y1 * beta_gx2 * beta_gy2 - beta_x2y2 * beta_gx1 * beta_gy1) /
    (beta_gx2 * beta_gy2 - beta_gx1 * beta_gy1)
  if (beta_gy1 / beta_gx1 > 0) {
    sign <- 1
  }
  if (beta_gy1 / beta_gx1 < 0) {
    sign <- (-1)
  }
  b <- sign * (sqrt((beta_gy1 * (beta_x1x2 - k2) * beta_gy2) /
      (beta_gx1 * (beta_y1y2 - k2) * beta_gx2)))
  data <- data.frame(g, x1, x2, y1, y2)
  colnames(data) <- c("g", "x1", "x2", "y1", "y2")
  b_es <- rep(0, 2000)
  for (i in 1:2000) {
    nr <- nrow(data)
    select <- sample(c(1 : nr), nr, replace = T)
    data_use <- data[select, ]
    g <- data_use $ g
    x1 <- data_use $ x1
    x2 <- data_use $ x2
    y1 <- data_use $ y1
    y2 <- data_use $ y2
    beta_gy1 <- summary(lm(y1 ~ g)) $ coef[2, 1]
    beta_gy2 <- summary(lm(y2 ~ g)) $ coef[2, 1]
    beta_x1x2 <- summary(lm(x2 ~ x1)) $ coef[2, 1]
    beta_y1y2 <- summary(lm(y2 ~ y1)) $ coef[2, 1]
    beta_gx1 <- summary(lm(x1 ~ g)) $ coef[2, 1]
    beta_gx2 <- summary(lm(x2 ~ g)) $ coef[2, 1]
    beta_x1y1 <- summary(lm(y1 ~ x1)) $ coef[2, 1]
    beta_x2y2 <- summary(lm(y2 ~ x2)) $ coef[2, 1]
    k2 <- (beta_x1y1 * beta_gx2 * beta_gy2 - beta_x2y2 * beta_gx1 * beta_gy1) /
      (beta_gx2 * beta_gy2 - beta_gx1 * beta_gy1)
    if (beta_gx1 / beta_gy1 > 0) {
      sign <- 1
    }
    if (beta_gx1 / beta_gy1 < 0) {
      sign <- (-1)
    }
    if ((beta_gy1 * (beta_x1x2 - k2) * beta_gy2) / (beta_gx1 *
          (beta_y1y2 - k2) * beta_gx2) > 0) {
      b_es[i] <- sign * (sqrt((beta_gy1 * (beta_x1x2 - k2) * beta_gy2) /
              (beta_gx1 * (beta_y1y2 - k2) * beta_gx2)))
    }
  }
  ### SE\CI###
  b_es <- na.omit(b_es)
  lengthbes <- length(b_es)
  for (j in 1 : lengthbes) {
    if (abs(b_es[j]) > 1) {
      b_es[j] <- NA
    }
  }
  b_es <- na.omit(b_es)
  se <- round(sd(b_es), 3)
  cin <- paste0("(", round((b - 1.96 * se), 3), ",",
                round((b + 1.96 * se), 3), ")")
  ciq <- quantile(b_es, probs = c(0.025, 0.975))
  ciq <- paste0("(", round(ciq[1], 3), ",", round(ciq[2], 3),
                ")")
  result <- list(b, se, cin, ciq)
  names(result) <- c("Estimate", "SE", "95% CIn", "95% CIq")
  return(result)
}


#' IV approach when X and Y subject to misclassification in binary IV model
#' @name rcmivbinary
#' @param G A numeric vector of binary instrumental variable.
#' @param X1 A numeric vector of the first measurement of binary exposure X
#' @param X2 A numeric vector of the second measurement of binary exposure X
#' @param Y1 A numeric vector of the first measurement of binary outcome Y
#' @param Y2 A numeric vector of the second measurement of binary outcome Y
#' @return A vector of causal effect estimate, bootstrap SE and 95% CI
#' @author Xinhui Liu

rcmivbinary <- function(gbi, x1, x2, y1, y2) {
  data2 <- data.frame(gbi, x1, x2, y1, y2)
  colnames(data2) <- c("gbi", "x1", "x2", "y1", "y2")
  nn <- nrow(data2)
  data2$a <- ifelse(data2 $ x1 == 1 & data2 $ x2 == 1
                    & data2 $ gbi == 1, 1, 0)
  data2$b <- ifelse(data2 $ x1 == 1 & data2 $ x2 == 1
                    & data2 $ gbi == 0, 1, 0)
  data2$c <- ifelse(data2 $ x1 + data2 $ x2 == 1
                    & data2 $ gbi == 1, 1, 0)
  data2$d <- ifelse(data2 $ x1 + data2 $ x2 == 1
                    & data2 $ gbi == 0, 1, 0)
  data2$e <- ifelse(data2 $ x1 == 0 & data2 $ x2 == 0
                    & data2 $ gbi == 1, 1, 0)
  data2$f <- ifelse(data2 $ x1 == 0 & data2 $ x2 == 0
                    & data2 $ gbi == 0, 1, 0)

  data2$g <- ifelse(data2 $ y1 == 1 & data2 $ y2 == 1
                    & data2 $ gbi == 1, 1, 0)
  data2$h <- ifelse(data2 $ y1 == 1 & data2 $ y2 == 1
                    & data2 $ gbi == 0, 1, 0)
  data2$i <- ifelse(data2 $ y1 + data2 $ y2 == 1 & data2
                    $ gbi == 1, 1, 0)
  data2$j <- ifelse(data2 $ y1 + data2 $ y2 == 1 & data2
                    $ gbi == 0, 1, 0)
  data2$k <- ifelse(data2 $ y1 == 0 & data2 $ y2 == 0
                    & data2 $ gbi == 1, 1, 0)
  data2$l <- ifelse(data2 $ y1 == 0 & data2 $ y2 == 0
                    & data2 $ gbi == 0, 1, 0)
  data2$sumx1x21 <- ifelse(data2 $ x1 + data2 $ x2 == 1, 1, 0)
  data2$sumy1y21 <- ifelse(data2 $ y1 + data2 $ y2 == 1, 1, 0)

  data2$g1 <- ifelse(data2 $ gbi == 1, 1, 0)
  data2$g0 <- ifelse(data2 $ gbi == 0, 1, 0)
  a <- sum(data2 $ a) / nn
  b <- sum(data2 $ b) / nn
  e <- sum(data2 $ e) / nn
  f <- sum(data2 $ f) / nn

  g <- sum(data2 $ g) / nn
  h <- sum(data2 $ h) / nn
  k <- sum(data2 $ k) / nn
  l <- sum(data2 $ l) / nn

  sumx1x21 <- sum(data2 $ sumx1x21) / nn
  sumy1y21 <- sum(data2 $ sumy1y21) / nn

  g1 <- sum(data2 $ g1) / nn
  g0 <- sum(data2 $ g0) / nn


  alpha <- 0.5 + 0.5 * sqrt(1 - 2 * sumx1x21)
  beta <- 0.5 + 0.5 * sqrt(1 - 2 * sumy1y21)

  c11 <- (1 / 2) * (1 + (((a / g1) - (e / g1)) /
                           (2 * alpha - 1)))
  d11 <- (1 / 2) * (1 + (((b / g0) - (f / g0)) /
                           (2 * alpha - 1)))
  a11 <- (1 / 2) * (1 + (((g / g1) - (k / g1)) /
                           (2 * beta - 1)))
  b11 <- (1 / 2) * (1 + (((h / g0) - (l / g0)) /
                           (2 * beta - 1)))

  bbbb <- (a11 - b11) / (c11 - d11)


  data_use <- data.frame(gbi, x1, x2, y1, y2)
  colnames(data_use) <- c("gbi", "x1", "x2", "y1", "y2")
  b_es <- rep(0, 2000)
  for (kk in 1:2000) {
    nr <- nrow(data_use)
    select <- sample(c(1 : nr), nr, replace = T)
    data2 <- data_use[select, ]
    nn <- nrow(data2)
    data2$a <- ifelse(data2 $ x1 == 1 & data2 $ x2 == 1
                      & data2 $ gbi == 1, 1, 0)
    data2$b <- ifelse(data2 $ x1 == 1 & data2 $ x2 == 1
                      & data2 $ gbi == 0, 1, 0)
    data2$c <- ifelse(data2 $ x1 + data2 $ x2 == 1
                      & data2 $ gbi == 1, 1, 0)
    data2$d <- ifelse(data2 $ x1 + data2 $ x2 == 1
                      & data2 $ gbi == 0, 1, 0)
    data2$e <- ifelse(data2 $ x1 == 0 & data2 $ x2 == 0
                      & data2 $ gbi == 1, 1, 0)
    data2$f <- ifelse(data2 $ x1 == 0 & data2 $ x2 == 0
                      & data2 $ gbi == 0, 1, 0)

    data2$g <- ifelse(data2 $ y1 == 1 & data2 $ y2 == 1
                      & data2 $ gbi == 1, 1, 0)
    data2$h <- ifelse(data2 $ y1 == 1 & data2 $ y2 == 1
                      & data2 $ gbi == 0, 1, 0)
    data2$i <- ifelse(data2 $ y1 + data2 $ y2 == 1 & data2
                      $ gbi == 1, 1, 0)
    data2$j <- ifelse(data2 $ y1 + data2 $ y2 == 1 & data2
                      $ gbi == 0, 1, 0)
    data2$k <- ifelse(data2 $ y1 == 0 & data2 $ y2 == 0
                      & data2 $ gbi == 1, 1, 0)
    data2$l <- ifelse(data2 $ y1 == 0 & data2 $ y2 == 0
                      & data2 $ gbi == 0, 1, 0)
    data2$sumx1x21 <- ifelse(data2 $ x1 + data2 $ x2 == 1, 1, 0)
    data2$sumy1y21 <- ifelse(data2 $ y1 + data2 $ y2 == 1, 1, 0)

    data2$g1 <- ifelse(data2 $ gbi == 1, 1, 0)
    data2$g0 <- ifelse(data2 $ gbi == 0, 1, 0)
    a <- sum(data2 $ a) / nn
    b <- sum(data2 $ b) / nn
    e <- sum(data2 $ e) / nn
    f <- sum(data2 $ f) / nn

    g <- sum(data2 $ g) / nn
    h <- sum(data2 $ h) / nn
    k <- sum(data2 $ k) / nn
    l <- sum(data2 $ l) / nn

    sumx1x21 <- sum(data2 $ sumx1x21) / nn
    sumy1y21 <- sum(data2 $ sumy1y21) / nn

    g1 <- sum(data2 $ g1) / nn
    g0 <- sum(data2 $ g0) / nn


    alpha <- 0.5 + 0.5 * sqrt(1 - 2 * sumx1x21)
    beta <- 0.5 + 0.5 * sqrt(1 - 2 * sumy1y21)

    c11 <- (1 / 2) * (1 + (((a / g1) - (e / g1)) /
                             (2 * alpha - 1)))
    d11 <- (1 / 2) * (1 + (((b / g0) - (f / g0)) /
                             (2 * alpha - 1)))
    a11 <- (1 / 2) * (1 + (((g / g1) - (k / g1)) /
                             (2 * beta - 1)))
    b11 <- (1 / 2) * (1 + (((h / g0) - (l / g0)) /
                             (2 * beta - 1)))

    b_es[kk] <- (a11 - b11) / (c11 - d11)
  }
  ### SE\CI###
  b_es <- na.omit(b_es)

  se <- round(sd(b_es), 3)
  cin <- paste0("(", round((bbbb - 1.96 * se), 3), ",",
                round((bbbb + 1.96 * se), 3), ")")
  ciq <- quantile(b_es, probs = c(0.025, 0.975))
  ciq <- paste0("(", round(ciq[1], 3), ",", round(ciq[2], 3), ")")
  result <- list(bbbb, se, cin, ciq)
  names(result) <- c("Estimate", "SE", "95% CIn", "95% CIq")
  return(result)
}

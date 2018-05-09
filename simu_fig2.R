CompPall2 <- function(p.vec,
                      N.vec) {

  n <- length(p.vec)
  r <- 2
  p.selected <- sort(p.vec, decreasing = T)[-n]

  n0 <- n -1
  p.fisher <- pchisq(-2 * sum(log(p.selected)), 
                    2 * n0,
                    lower.tail = F)
  p.simes <- min(n0/(n0:1) * p.selected)
  p.stou0 <- pnorm(sum(qnorm(1 - p.selected)) / sqrt(n0), lower.tail = F)


  N.uni <- unique(N.vec)
  gu.pu <- sapply(N.uni, function(N) {
                  idx.out <- which(p.vec == min(p.vec[N.vec == N]))
                  p.selected <- p.vec[-idx.out]
                  N.selected <- N.vec[-idx.out]
                  return(pnorm(sum(sqrt(N.selected) * 
                                   qnorm(1 - p.selected))/sqrt(sum(N.selected)), 
                               lower.tail = F))})
  p.stou <- max(gu.pu)
  return(c(p.fisher, p.simes, p.stou0, p.stou))
}

CompPower <- function(mu, sd, B = 20000, 
                      n = 8, r0 = 2,
                      alpha = 0.05,
                      N.vec = c(rep(100, 3), rep(500, 3), rep(1000, 2))) {
  if (sd > 150 * mu)
    return(rep(NA, 8))
  a <- (mu/sd)^2
  b <- mu / sd^2
  u.mat <- matrix(rgamma(r0 * B, a, b), nrow = B)
 # print(c(mu, sd))
 # print(c(mean(u.mat), sd(as.vector(u.mat))))
 # print("")
  u.mat <- cbind(matrix(0, B, n - r0), u.mat)
  u.mat <- t(sqrt(N.vec) * t(u.mat))
 # print(dim(u.mat))
  p.mat <- 2 * pnorm(abs(matrix(rnorm(n * B), nrow = B) + u.mat), lower.tail = F)
#  p.mat <- pnorm(matrix(rnorm(n * B), nrow = B) + u.mat, lower.tail = F)
  
  pc.pvec <- t(apply(p.mat, 1, CompPall2, N.vec))
  rej <- pc.pvec <= alpha
  return(c(colMeans(rej),
           apply(rej, 2, sd)/sqrt(B)))
}

library(reshape)

GenMap <- function(mu.range = c(0.03, 0.38), 
                   sd.range = c(0.001, 0.5), 
                   mm = 40,
                   r0 = 2) {
  mu.seq <- seq(mu.range[1], mu.range[2], length.out = mm)
  sd.seq <- seq(sd.range[1], sd.range[2], length.out = mm)


  idx <- melt(matrix(0, mm, mm))
  idx <- idx[, 1:2]

  require(parallel)
  cl <- makeCluster(3)

  clusterExport(cl, ls(.GlobalEnv))
  clusterExport(cl, ls(), environment())

  result <- parSapply(cl, 1:nrow(idx), function(i) {
                      require(shaRNG)
                    i %splitRngEval% CompPower(mu.seq[idx[i, 1]], 
                                               sd.seq[idx[i, 2]], r0 = r0)})

  result <- cbind(idx, t(result))
  return(list(result = result, 
              mu.seq = mu.seq,
              sd.seq = sd.seq))

}

PowerPlot <- function(result,
                      method.idx = 1,
                      zlim1 = zz1,
                      larger.margin = c(T, T, T, T),
                      ylab = "") {
  temp <- c(2, 1, 1, 2)
  temp[larger.margin] <- 4
  par(mar = temp)
  mu.seq <- result$mu.seq
  sd.seq <- result$sd.seq
  result <- result$result

  aa <- result[, 1:2]
  result <- result[, -(1:2)]
  if (method.idx != 1) {
    value <- result[, method.idx] - result[, 1]
    print(range(value))
    value[value > zlim1[2]] <- zlim1[2]
    value[value < zlim1[1]] <- zlim1[1]
    aa <- cbind(aa, value)
    zlim <- zlim1
    col <- my_palette
  } else {
    value <- result[, method.idx]
    aa <- cbind(aa, value)
    zlim <- c(0, 1)
    col <- my_palette1
  }

  zmat <- cast(aa, X1~X2)[, -1]
 # print(head(zmat))
#  zmat <- zmat[, -1]
#  print(zmat)
  main <- ""
  if (larger.margin[3]) {
    if (larger.margin[2])
      main <- "Power of Fisher's BHPC"
    else if (larger.margin[4])
      main <- "Power Gain of WS GBHPC"
    else
      main <- "Power Gain of Simes BHPC"
  }
  if (!larger.margin[2])
    yaxt <- "n" 
  else 
    yaxt <- "s"
  xlim <- range(mu.seq)
  xlim[2] <- xlim[2] + 0.005
  print(xlim)
  if (larger.margin[4])
    xlim[2] <- xlim[2] + (xlim[2] - xlim[1])/40
  xlab <- ""
  if (larger.margin[1])
    xlab <- expression(bold(mu[0]))
  image(mu.seq, sd.seq, as.matrix(zmat), 
        zlim = zlim,
        xlim = xlim,
        col = col, xlab = xlab, ylab = ylab,
        yaxt = yaxt,  
        main = main)


}
legend.col <- function(col, lev){

  opar <- par

  n <- length(col)

  bx <- par("usr")

  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n

  xx <- rep(box.cx, each = 2)

  par(xpd = TRUE)
  for(i in 1:n){

    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])

  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25, cex.axis = 1)
  par <- opar
}



#my_palette <- gray(0:200/200)
my_palette <- colorRampPalette(c("green", "white", "blue"))(n = 299)
#my_palette1 <- heat.colors(200)[200:1]
my_palette1 <- heat.colors(290)[c(290:101, seq(100, 1, by = -10))]
zz1 <- c(-0.1, 0.1)
zz2 <- c()

mu.range2 <- c(0.02, 0.3)
sd.range2 <- c(0.001, 0.3)

mu.range4 <- c(0.03, 0.15)
sd.range4 <- c(0.001, 0.15)

mu.range6 <- c(0.02, 0.11)
sd.range6 <- c(0.001, 0.11)


DoPlot <- function(result) {
  layout(matrix(1:9, nrow = 3, byrow = T), c(1.3, 0.9, 1), c(1.1, 0.9, 1))
  PowerPlot(result, 1, larger.margin = c(F, T, T, T))
  legend.col(my_palette1, c(0, 1))
  PowerPlot(result, 2, larger.margin = c(F, F, T, F))
  #PowerPlot(result, 3)
  PowerPlot(result, 4, larger.margin = c(F, F, T, T))
  legend.col(my_palette, zz1)
}

DoPlotWrapper <- function(result1, result2, result3){
  layout(matrix(1:9, nrow = 3, byrow = T), c(1.3, 0.9, 1), c(1.1, 0.9, 1))
  PowerPlot(result1, 1, larger.margin = c(F, T, T, T),
            ylab = expression(bold(paste(sigma[0], " (", r[0], " = 2)"))))
  legend.col(my_palette1, c(0, 1))
  PowerPlot(result1, 2, larger.margin = c(F, F, T, F))
  #PowerPlot(result, 3)
  PowerPlot(result1, 4, larger.margin = c(F, F, T, T))
  legend.col(my_palette, zz1)
  PowerPlot(result2, 1, larger.margin = c(F, T, F, T),
            ylab = expression(bold(paste(sigma[0], " (", r[0], " = 4)"))))
  legend.col(my_palette1, c(0, 1))
  PowerPlot(result2, 2, larger.margin = c(F, F, F, F))
  PowerPlot(result2, 4, larger.margin = c(F, F, F, T))
  legend.col(my_palette, zz1)
  PowerPlot(result3, 1, larger.margin = c(T, T, F, T),
            ylab = expression(bold(paste(sigma[0], " (", r[0], " = 6)"))))
  legend.col(my_palette1, c(0, 1))
  PowerPlot(result3, 2, larger.margin = c(T, F, F, F))
  PowerPlot(result3, 4, larger.margin = c(T, F, F, T))
  legend.col(my_palette, zz1)
}


###########################################################################################
## This can tak a long time
result1 <- GenMap(r0 = 2)
result2 <- GenMap(r0 = 4)
result3 <- GenMap(r0 = 6)

DoPlotWrapper(result1, result2, result3)


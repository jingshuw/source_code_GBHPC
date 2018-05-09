############ Code for Figure 1 ##########

###### left, upper panel ######
alpha = 0.1
n.req <- 10
dens <- 20
layout(matrix(1:4, nrow = 2, byrow= 2))
par(mar = c(2, 2.2, 3.5, 2))


plot((1:10)/10, (1:10)/10, type = "n", xlim = c(0, 1), ylim = c(0, 1),
	 xaxs = "i", xaxt = "n", yaxs = "i", yaxt = "n", 
   main = expression(bold(paste("Rejection region of (", 
										   p[1], ",",p[2], ")", sep = ""))),
	 col.main = "blue"
   );

for (i in 1:n.req) {
  if (i != 1 && i != n.req)
    lty <- 3
  else
    lty <- 1
  rect((i - 1) * alpha, (i-1) * alpha, 
       i * alpha, i * alpha, lty = lty, density = dens, col = "red")
}
axis(1, c(0, alpha, 0.5, 1 - alpha, 1), label = T,  cex.axis = 1.1)
axis(2, c(0, alpha, 0.5, 1 - alpha, 1), label = T,  cex.axis = 1.1,
     las = 2)


#### right, upper panel ####

limit <- 3
plot((-limit):limit, (-limit):limit, type = "n", xaxs = "i", yaxs = "i", bty = "n",
	 xaxt = "n", yaxt = "n", xlab = expression(Z[1]), 
	 ylab = expression(Z[2]), main = expression(bold(paste("Rejection region of (", 
										   Z[1], ",",Z[2], ")", sep = ""))),
	 col.main = "blue"
   );

pos1 <- qnorm(1 - alpha/2);
pos2 <- qnorm(0.5 + 0:(n.req - 1) * alpha/2)
pos2 <- c(pos2, limit)
segments(pos1, pos1, c(pos1, limit), c(limit, pos1), col = "RED");
segments(pos1, -pos1, c(pos1, limit), c(-limit, -pos1), col = "RED");
segments(-pos1, pos1, c(-pos1, -limit), c(limit, pos1), col = "RED");
segments(-pos1, -pos1, c(-pos1, -limit), c(-limit, -pos1), col = "RED")
for (i in 1:n.req) {
  if (i != 1 && i != n.req)
    lty <- 3
  else
    lty <- 1
  if (i == n.req)
    border <- F
  else
    border <- T
  pos.l <- pos2[i]
  pos.r <- pos2[i + 1]
  rect(pos.l, pos.l, pos.r, pos.r, density = dens, border = border, 
       lty = lty, col = "RED")
  rect(pos.l, -pos.l, pos.r, -pos.r, density = dens, border = border, 
       lty = lty, col = "RED")
  rect(-pos.l, pos.l, -pos.r, pos.r, density = dens, border = border, 
       lty = lty, col = "RED")
  rect(-pos.l, -pos.l, -pos.r, -pos.r, density = dens, border = border, 
       lty = lty, col = "RED")
}

axis(1, c(-2, -1, 1, 2), pos = 0, label = T,  cex.axis = 1.3)
axis(2, seq(-2, 2, 1), pos = 0, label = T,  cex.axis = 1.3, las = 2)
segments(-limit, 0, limit, 0)
segments(0, -limit, 0, limit)


##### bottom two panels ######

Power1 <- function(a1, a2, alpha = 0.2) {
  thres <- qnorm(1 - alpha/2)
  beta1 <- pnorm(thres, a1, lower.tail = F) + pnorm(-thres, a1) 
  beta2 <- pnorm(thres, a2, lower.tail = F) + pnorm(-thres, a2)
  return(beta1 * beta2)
}

Power2 <- function(a1, a2, alpha = 0.2) {
  thres.v <- qnorm(1 - (n.req:1 -1) * alpha/2)
  thres.v <- c(-thres.v[length(thres.v):1], 0, thres.v)
  n <- length(thres.v)
  pos.mat1 <- cbind(thres.v[-n], thres.v[-1], thres.v[-n], thres.v[-1])
  pos.mat2 <- cbind(thres.v[-n], thres.v[-1], thres.v[(n-1):1], thres.v[n:2])
  pos.mat <- rbind(pos.mat1, pos.mat2)
  pp <- apply(pos.mat, 1, function(v) {
                  beta1 <- pnorm(v[2], a1) - pnorm(v[1], a1)
                  beta2 <- pnorm(v[4], a2) - pnorm(v[3], a2)
                  return(beta1 * beta2)
        })
  return(sum(pp))

}

PowerPlot <- function(low, high, 
                      pow.fun,
                      cols,
                      alpha = 0.2, B = 100,
                      x.adjust= F,
                      idx = 1,
                      p.thres = c(0.15, 0.5),
                      tol = c(0.005, 0.01),
                      lty = 2:3) {
  p.seq <- seq(low, high, length.out = B)
  mat <- sapply(p.seq, function(p1) {
        return(sapply(p.seq, function(p2) pow.fun(p1, p2, alpha)))
       })
  require(reshape)
  m.p <- length(p.thres)
  pp.list <- list()
  for (i in 1:m.p) {
    tt <- abs(mat - p.thres[i]) < tol[i]
    temp <- melt(tt)
    pp <- temp[temp$value == 1, ]
    pp[, 1] <- p.seq[pp[, 1]]
    pp[, 2] <- p.seq[pp[, 2]]
    pp <- pp[, 1:2]
    pp.list[[i]] <- pp 
  }
  mat <- rbind(cbind(mat[B:1, B:1], mat[B:1, ]), cbind(mat[, B:1], mat))
  p.seq <- c(-p.seq[B:1], p.seq)
  xlim <- range(p.seq)
  if (x.adjust)
    xlim[2] <- xlim[2] + (xlim[2] - xlim[1])/30
  image(p.seq, p.seq, mat, 
        zlim = c(0, 1),
        xlim = xlim,
        col = cols)
  if (idx == 1)
    title(expression(bold(paste("power of ", varphi))))
  if (idx == 2)
  title(expression(bold(paste("power of ", tilde(varphi)))))

  sapply(1:m.p, function(i) {
         pp <- pp.list[[i]]
         m <- nrow(pp)
         points(pp, type = "l", lty = lty[i])
         points(t(c(-1, 1) * t(pp)), type = "l", lty = lty[i])
         points(t(c(-1, -1) * t(pp)), type = "l", lty = lty[i])
         points(t(c(1, -1) * t(pp)), type = "l", lty = lty[i])
        })
  


  return(list(mat = mat, pp = pp))
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


par(mar = c(2, 2, 3.5, 3))
my_palette <- heat.colors(100)[100:1]
a1 <- PowerPlot(0.01, 5, Power1, my_palette, alpha = alpha)
a2 <- PowerPlot(0.01, 5, Power2, my_palette, alpha = alpha,
          x.adjust = T, idx = 2)
legend.col(my_palette, c(0, 1))



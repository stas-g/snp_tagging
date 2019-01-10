foo1 <- function(zeta, r1, r2, r12) {
  cnd1 <- -2 * r1 * r2 * r12 + r1^2 + r2^2 + r12^2 <= 1
  cnd2 <- r1 + r2 < sqrt(2 * (1 + r12))
  ind <- cnd1 * cnd1 == 1
  zeta[ind] <- r1[ind] <- r2[ind] <- r12[ind] <- NA
  
  1 - pnorm(0.5 * abs(zeta) * sqrt( 2 * (1 + r12) - (r1 + r2)^2 ))
}






cnd1 <- function(r1, r2, r12) -2 * r1 * r2 * r12 + r1^2 + r2^2 + r12^2 <= 1
cnd2 <- function(r1, r2, r12) r1 + r2 < sqrt(2 * (1 + r12))

zeta <- c(0.01, 0.05, 0.1, 0.3, 0.5)

r <- seq(-1, 1, by = 0.01)
r.mat <- data.table(expand.grid(r1 = r, r2 = r, r12 = 0))
r.mat[, r.sum := r1 + r2]
r.mat[, c1 := cnd1(r1, r2, r12)]
r.mat[, c2 := cnd2(r1, r2, r12)]
r.mat <- r.mat[c1 & c2]


r.mat <- r.mat[!duplicated(r.mat[, c("r.sum", "zeta")], MARGIN = 1), ]

r.mat$c1


# SNP 3 is causal


x <- z <- seq(-1, 1, by = 0.01)

y1 <- x^2 + z^2
y2 <- x + z

plot(range(y1, y2), range(y1, y2), type = 'n', xlab = 'x2 + z2', ylab = 'x + z')
rect(-2, -2, 2, sqrt(2), col = alpha('firebrick1', 0.5), border = NA)
rect(-2, -2, 1, 2, col = alpha('dodgerblue', 0.5), border = NA)
abline(h = sqrt(2), lty = 2, col = 'firebrick1', lwd = 2)
abline(v = 1, lty = 2, col = 'dodgerblue', lwd = 2) 

points(r.mat$r1^2 + r.mat$r2^2, r.mat$r1 + r.mat$r2, col = "purple")
lines(y1, y2, col = "seagreen3", type = 'l', lwd = 3)


c1 <- r.mat$r1^2 + r.mat$r2^2 
c2 <- r.mat$r1 + r.mat$r2
















```{r}
foo1 <- function(zeta, r1, r2, r12) {
  cnd1 <- -2 * r1 * r2 * r12 + r1^2 + r2^2 + r12^2 <= 1
  cnd2 <- r1 + r2 < sqrt(2 * (1 + r12))
  ind <- cnd1 * cnd1 == 1
  zeta[ind] <- r1[ind] <- r2[ind] <- r12[ind] <- NA
  
  1 - pnorm(0.5 * abs(zeta) * sqrt( 2 * (1 + r12) - (r1 + r2)^2 ))
}

cnd1 <- function(r1, r2, r12) -2 * r1 * r2 * r12 + r1^2 + r2^2 + r12^2 <= 1
cnd2 <- function(r1, r2, r12) r1 + r2 < sqrt(2 * (1 + r12))

zeta <- c(0.01, 0.05, 0.1, 0.3, 0.5)
r <- seq(-1, 1, by = 0.1)
r.mat <- data.table(expand.grid(r1 = r, r2 = r, r12 = 0, zeta = zeta))
r.mat[, r.sum := r1 + r2]
r.mat[, c1 := cnd2(r1, r2, r12)]
r.mat[, c2 := cnd2(r1, r2, r12)]
r.mat <- r.mat[c1 & c2]

r.mat <- r.mat[!duplicated(r.mat[, c("r.sum", "zeta")], MARGIN = 1), ]

r.mat$c1

scatter3D(r.mat$r.sum, r.mat$zeta, )
```
# SNP 3 is causal


x <- z <- seq(-1, 1, by = 0.01)

y1 <- x^2 + z^2
y2 <- x + z

plot(range(y1, y2), range(y1, y2), type = 'n', xlab = 'x2 + z2', ylab = 'x + z')
rect(-2, -2, 2, sqrt(2), col = alpha('firebrick1', 0.5), border = NA)
rect(-2, -2, 1, 2, col = alpha('dodgerblue', 0.5), border = NA)
abline(h = sqrt(2), lty = 2, col = 'firebrick1', lwd = 2)
abline(v = 1, lty = 2, col = 'dodgerblue', lwd = 2) 

lines(y1, y2, col = "seagreen3", type = 'l', lwd = 3)
points(r.mat$r1^2 + r.mat$r2^2, r.mat$r1 + r.mat$r2, col = "purple")










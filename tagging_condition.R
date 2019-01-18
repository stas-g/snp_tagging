library(ggplot2)

plot(range(y1, y2), range(y1, y2), type = 'n', xlab = 'x2 + z2', ylab = 'x + z')
rect(-2, -2, 2, sqrt(2), col = alpha('firebrick1', 0.5), border = NA)
rect(-2, -2, 1, 2, col = alpha('dodgerblue', 0.5), border = NA)
abline(h = sqrt(2), lty = 2, col = 'firebrick1', lwd = 2)
abline(v = 1, lty = 2, col = 'dodgerblue', lwd = 2) 

points(r.mat$r1^2 + r.mat$r2^2, r.mat$r1 + r.mat$r2, col = "purple")
lines(y1, y2, col = "seagreen3", type = 'l', lwd = 3)


#-----------------------------------------------------------------------------------------
#SNPs 1 & 2 are causal
r12.vec <- c(-1, -0.5, 0, 0.5, 1)
zeta <- seq(0, 10, by = 1)

foo <- function(x, a = 0.01) {
  r <- seq(-1, 1, by = a)
  dat <- data.table(expand.grid(r, r))
  colnames(dat) <- c('r1', 'r2')
  cnd1 <- -2 * x * dat$r1 * dat$r2 + dat$r1^2 + dat$r2^2 + x^2 <= 1
  cnd2 <- round(dat$r1 + dat$r2, 2) < sqrt(2 * (1 + x))
  dat <- dat[cnd1 & cnd2, ]
  dat[, rsum := r1 + r2]
  dat
}

dat <- lapply(r12.vec, FUN = function(r12) {
  dat <- foo(r12)
  if(nrow(dat) == 0) return(NA)
  
  val.mat <- lapply(zeta, FUN = function(z) {
    v <- 0.5 * z * sqrt(2 * (1 + r12) - dat$rsum^2)
    1 - pnorm(v)
  }) %>% do.call(cbind, .)
  colnames(val.mat) <- paste0("zeta=", zeta)
  
  dat <- data.frame(dat, val.mat)
  dat <- melt(dat, id = c("r1", "r2", "rsum"))
  dat$r12 <- r12
  dat 
}) %>% do.call(rbind, .)

dat <- dat[!is.na(dat$r12), ]

ggplot(dat, aes(x = rsum^2, y = value)) + geom_line(aes(colour = variable)) + 
  facet_grid(r12 ~ .) +
  labs(x = expression((r[1] + r[2])^"2"), y = "Prob(choose SNP3, SNPs 1 & 2 are causal)") + theme_bw()


#-----------------------------------------------------------------------------------------
#SNP 3 is causal

r12.vec <- c(-1, -0.5, 0, 0.5, 1)
zeta <- seq(0, 10, by = 1)

foo <- function(x, a = 0.05) {
  r <- seq(-1, 1, by = a)
  dat <- data.table(expand.grid(r, r))
  colnames(dat) <- c('r1', 'r2')
  cnd1 <- -2 * x * dat$r1 * dat$r2 + dat$r1^2 + dat$r2^2 + x^2 <= 1
  cnd2 <- -2 * x * dat$r1 * dat$r2 + dat$r1^2 + dat$r2^2 <= 1
  dat <- dat[cnd1 & cnd2, ]
  dat
}

dat <- lapply(r12.vec, FUN = function(r12) {
  dat <- foo(r12)
  if(nrow(dat) == 0) return(NA)
  
  val.mat <- lapply(zeta, FUN = function(z) {
    prob <- 1 - pnorm(-0.5 * z * sqrt(round(1 - dat$r1^2 - dat$r2^2 + 2 * dat$r1 * dat$r2 * r12, 4)))
  }) %>% do.call(cbind, .)
  colnames(val.mat) <- paste0("zeta=", zeta)
  
  dat <- data.frame(dat, val.mat)
  dat <- melt(dat, id = c("r1", "r2"))
  dat$r12 <- r12
  dat 
}) %>% do.call(rbind, .)

dat <- dat[!is.na(dat$r12), ]

ggplot(dat, aes(x = r1, y = r2)) + geom_point(aes(colour = value)) + 
  facet_grid(r12 ~ variable) +
  labs(x = expression(r[1]), y = expression(r[2])) + 
  scale_colour_gradient("Prob(choose SNP3)", low = "steelblue", high = "darkorange1") +
  theme_bw()



#-----------------------------------------------------------------------------------------
#SNPs 1 & 2 are causal (BIC)

r12 <- c(-1, -0.5, 0, 0.5, 1)
zeta <- seq(0, 10, by = 1)

foo <- function(r12, a = 0.01) {
  r <- seq(-1, 1, by = a)
  dat <- data.table(expand.grid(r, r))
  colnames(dat) <- c('r1', 'r2')
  cnd1 <- -2 * r12 * dat$r1 * dat$r2 + dat$r1^2 + dat$r2^2 + r12^2 <= 1
  cnd2 <- round(dat$r1 + dat$r2, 2) < sqrt(2 * (1 + r12))
  dat <- dat[cnd1 & cnd2, ]
  dat[, rsum := r1 + r2]
  dat
}

dat <- lapply(r12, FUN = function(r) {
  dat <- foo(r)
  if(nrow(dat) == 0) return(NA)
  
  val.mat <- lapply(zeta, FUN = function(z) {
    v <- 0.5 * z * sqrt(2 * (1 + r) - dat$rsum^2)
    1 - pnorm(v)
  }) %>% do.call(cbind, .)
  colnames(val.mat) <- paste0("zeta=", zeta)
  
  dat <- data.frame(dat, val.mat)
  dat <- melt(dat, id = c("r1", "r2", "rsum"))
  dat$r12 <- r
  dat 
}) %>% do.call(rbind, .)

dat <- dat[!is.na(dat$r12), ]

ggplot(dat, aes(x = rsum^2, y = value)) + geom_line(aes(colour = variable)) + 
  facet_grid(r12 ~ .) +
  labs(x = expression((r[1] + r[2])^"2"), y = "Prob(choose SNP3, SNPs 1 & 2 are causal)") + theme_bw()































##


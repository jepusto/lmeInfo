set.seed(201)

beta <- c(0.25, 0.1, 2, 1)

s <- c(rep(2,5), rep(-2,5), rep(0,5))
p <- c(-1, -1, 3, 3, -2,
        2.5, 2.5, -2, -2, -0.5,
       -0.5, -0.5, -0.5, 0.6, 0.6)
x1 <- c(rep(2, 5), rep(1, 5), rep(2,5))
x2 <- rnorm(15)

y <- 0.2 + cbind(s, p, x1, x2) %*% beta + 0.1*rnorm(15)

dat <- data.frame(
  Source=c(rep("A",5),
           rep("B",5),
           rep("C",5)),
  Path=c("D", "D", "E", "E", "F",
         "E", "E", "F", "F", "G",
         "D", "D", "D", "G", "G"),
  Log10.Yield.=x1, Log10.Delta.=x2, y=y)

source("auxiliary/yield_im.r")

H=1
dat1 = list(H)
dat1[[1]] = dat
mod = list(H)

params = list()
params$iset = 1
params$dat = dat1
params$nmu = 2

require(nlme)
require(lmeInfo)
lm1 <- gls(y ~ Log10.Delta. + Log10.Yield.,
           data = dat,
           method = "REML")

mod[[1]] = lm1
params$ivc = "Residual"
params$vc = list(H)
params$vc[[1]] = mod[[1]]$sigma^2
params$mode = "REML"

im1 = yim(params)

im1.f = cbind(im1$I11[[1]], im1$I12[[1]])
im1.f = rbind(im1.f, cbind(t(im1$I12[[1]]), im1$I22[[1]]))
vcov1 = solve(im1.f)

print(vcov1)
print(lm1$varBeta)

im1.vc = Fisher_info(lm1, type="expected")

print(im1$I33[[1]])
print(im1.vc)

lm1a <- gls(y ~ Log10.Delta. + Log10.Yield.,
            data = dat,
            method = "ML")

mod[[1]] = lm1a
params$ivc = "Residual"
params$vc = list(H)
params$vc[[1]] = mod[[1]]$sigma^2
params$mode = "ML"

im1a = yim(params)

im1a.f = cbind(im1a$I11[[1]], im1a$I12[[1]])
im1a.f = rbind(im1a.f, cbind(t(im1a$I12[[1]]), im1a$I22[[1]]))
vcov1a = solve(im1a.f)

print(vcov1a)
print(lm1a$varBeta)

im1a.vc = Fisher_info(lm1a, type="expected")

print(im1a$I33[[1]])
print(im1a.vc)

lm2 <- lme(y ~ Log10.Delta. + Log10.Yield.,
           data = dat,
           random = ~1 | Source,
           method = "REML")

mod[[1]] = lm2
params$ivc = "Source"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Source[2]),
                   mod[[1]]$sigma)^2
params$mode = "REML"

im2 = yim(params)

im2.f = cbind(im2$I11[[1]], im2$I12[[1]])
im2.f = rbind(im2.f, cbind(t(im2$I12[[1]]), im2$I22[[1]]))
vcov2 = solve(im2.f)

print(vcov2)
print(lm2$varFix)

im2.vc = Fisher_info(lm2, type="expected")

print(im2$I33[[1]])
print(im2.vc)

lm2a <- lme(y ~ Log10.Delta. + Log10.Yield.,
            data = dat,
            random = ~1 | Source,
            method = "ML")

mod[[1]] = lm2a
params$ivc = "Source"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Source[2]),
                   mod[[1]]$sigma)^2
params$mode = "ML"

im2a = yim(params)

im2a.f = cbind(im2a$I11[[1]], im2a$I12[[1]])
im2a.f = rbind(im2a.f, cbind(t(im2a$I12[[1]]), im2a$I22[[1]]))
vcov2a = solve(im2a.f)

print(vcov2a)
print(lm2a$varFix)

im2a.vc = Fisher_info(lm2a, type="expected")

print(im2a$I33[[1]])
print(im2a.vc)

lm3 <- lme(y ~ Log10.Delta. + Log10.Yield.,
           data = dat,
           random = ~1 | Path,
           method = "REML")

mod[[1]] = lm3
params$ivc = "Path"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Path[2]),
                   mod[[1]]$sigma)^2
params$mode = "REML"

im3 = yim(params)

im3.f = cbind(im3$I11[[1]], im3$I12[[1]])
im3.f = rbind(im3.f, cbind(t(im3$I12[[1]]), im3$I22[[1]]))
vcov3 = solve(im3.f)

print(vcov3)
print(lm3$varFix)

im3.vc = Fisher_info(lm3, type="expected")

print(im3$I33[[1]])
print(im3.vc)

lm3a <- lme(y ~ Log10.Delta. + Log10.Yield.,
            data = dat,
            random = ~1 | Path,
            method = "ML")

mod[[1]] = lm3a
params$ivc = "Path"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Path[2]),
                   mod[[1]]$sigma)^2
params$mode = "ML"

im3a = yim(params)

im3a.f = cbind(im3a$I11[[1]], im3a$I12[[1]])
im3a.f = rbind(im3a.f, cbind(t(im3a$I12[[1]]), im3a$I22[[1]]))
vcov3a = solve(im3a.f)

print(vcov3a)
print(lm3a$varFix)

im3a.vc = Fisher_info(lm3a, type="expected")

print(im3a$I33[[1]])
print(im3a.vc)

lm4 <- lme(y ~ Log10.Delta. + Log10.Yield.,
           data = dat,
           random = ~1 | Source/Path,
           method = "REML")


mod[[1]] = lm4
params$ivc = "Full"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Source[2]),
                   as.numeric(intervals(mod[[1]])$reStruct$Path[2]),
                   mod[[1]]$sigma)^2
params$mode = "REML"

im4 = yim(params)

im4.f = cbind(im4$I11[[1]], im4$I12[[1]])
im4.f = rbind(im4.f, cbind(t(im4$I12[[1]]), im4$I22[[1]]))
vcov4 = solve(im4.f)

print(vcov4)
print(lm4$varFix)

im4.vc = Fisher_info(lm4, type="expected")


lm4_sep <- lme(y ~ Log10.Delta. + Log10.Yield.,
               data = dat,
               random = list(~1 | Source, ~ 1 | Path),
               method = "REML")

Fisher_info(lm4, type = "expected")
Fisher_info(lm4_sep, type = "expected")

print(im4$I33[[1]])
print(im4.vc)

lm4a <- lme(y ~ Log10.Delta. + Log10.Yield.,
            data = dat,
            random = ~1 | Source/Path,
            method = "ML")

mod[[1]] = lm4a
params$ivc = "Full"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Source[2]),
                   as.numeric(intervals(mod[[1]])$reStruct$Path[2]),
                   mod[[1]]$sigma)^2
params$mode = "ML"

im4a = yim(params)

im4a.f = cbind(im4a$I11[[1]], im4a$I12[[1]])
im4a.f = rbind(im4a.f, cbind(t(im4a$I12[[1]]), im4a$I22[[1]]))
vcov4a = solve(im4a.f)

print(vcov4a)
print(lm4a$varFix)

im4a.vc = Fisher_info(lm4a, type="expected")

print(im4a$I33[[1]])
print(im4a.vc)



lm5 <- lme(y ~ 1,
           data = dat,
           random = ~1 | Source/Path,
           method = "REML")


mod[[1]] = lm5
params$ivc = "Full"
params$vc = list(H)
params$vc[[1]] = c(as.numeric(intervals(mod[[1]])$reStruct$Source[2]),
                   as.numeric(intervals(mod[[1]])$reStruct$Path[2]),
                   mod[[1]]$sigma)^2
params$mode = "REML"

im5 = yim(params)

im5.vc = Fisher_info(lm5, type="expected")
print(im5$I33[[1]])
print(im5.vc)

# Parameter list:
parameters <- list(mu    = 5,
                   phi   = matrix(c(.8, .1), ncol = 1),
                   theta = matrix(c(-.42, -.2), ncol = 1),
                   omega = .9,
                   alpha = matrix(c(.1,.05,.2,.1), ncol = 2),
                   beta  = matrix(c(.6,.2), ncol = 1))
# Simulate data:
m <- c(25, 1)
burnin <- 300
n <- 2500
W <- create.neighbourhood.array(m = m, sp = 3)
set.seed(1234)
y <- simSTARMAGARCH(parameters, n = n, m = m, W=W, burnin = burnin, type= "queen",
               torus = TRUE)
# Create likelihood object:
f <- CreateLikelihood(y, W=W,
                 init = apply(y,1,var), parameters=parameters)

# Optimize the likelihood using:
fit <- fitSTARMAGARCH(f, data= y, print = FALSE)
summary(fit)

plot(fit)
plot_garch(fit)
# Checking residuals:
qqnorm(residuals(fit), pch = 19)
qqline(residuals(fit), col = "red")

# ---------------------------------
# ------ MUST BE FIXED ------------
# ---------------------------------
#Set parameters$beta[2,1] to zero:
#map <- parameterlist2maptemplate(parameters)
#map$beta[2,1]<-as.factor(NA)
#parameters$beta[2,1]<-0
#f <- CreateLikelihood(y, W=W,
#                      init = apply(y,1,var), parameters=parameters, map = map)
#
# Optimize the likelihood using:
#fit <- fitSTARMAGARCH(f, data= y, parameters, print = FALSE)

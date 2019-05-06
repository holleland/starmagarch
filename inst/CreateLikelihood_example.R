# Parameter list:
parameters <- list(mu    = 5,
                   phi   = matrix(c(.8, .1), ncol = 1),
                   theta = matrix(c(-.42, -.2), ncol = 1),
                   omega = .9,
                   alpha = matrix(c(.1,.05,.2,.1), ncol = 2),
                   beta  = matrix(c(.6,.2), ncol = 1))
m <- c(20, 1)
burnin <- 150
n <- 500

# Simulate data:
set.seed(1234)
W <- create.neighbourhood.array(m = m, sp = 3, type = "queen", torus = TRUE)
y <- simSTARMAGARCH(parameters, n = n, m = m, W = W, burnin = burnin)
# Create likelihood object:
f <- CreateLikelihood(y, W = W, init = apply(y, 1, var),
                      parameters = parameters, silent = FALSE)

# Optimize the likelihood using:
fit <- fitSTARMAGARCH(f, data = y, print = FALSE)
summary(fit)

plot(fit)
plot_garch(fit)
# Checking residuals:
qqnorm(residuals(fit), pch = 19)
qqline(residuals(fit), col = "red")

# ---------------------------------
# ------ MUST BE FIXED ------------
# ---------------------------------
# If you want to set a parameter to a fixed value,
# for instance the mean "mu" to zero, use the map argument:
parameters$mu = 0
f <- CreateLikelihood(y - 5, W = W,
                      init = apply(y, 1, var),
                      parameters = parameters,
                      map = list(mu = factor(NA)))

# Optimize the likelihood using:
fit <- fitSTARMAGARCH(f, data= y, parameters, print = FALSE)
summary(fit)

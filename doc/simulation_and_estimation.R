## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----parameters, eval=TRUE----------------------------------------------------
parameters <- list(
  mu    = 5,
  phi   = matrix(c(.8, .1), ncol = 1), #temporal order 1, spatial order 2
  theta = matrix(c(-.42, -.2), ncol = 1),  #temporal order 1, spatial order 2
  omega = .9,
  alpha = matrix(c(.09, .13,
                   .03, .07), ncol = 2), # temporal order 2, spatial order 2
  beta  = matrix(c(.5, .15), ncol = 1) #temporal order 1, spatial order 2
)
sum(parameters$alpha,parameters$beta)

## ----neighbourmatrix, eval = TRUE---------------------------------------------
# spatial dimension: 
m <- c(25, 1)
library(starmagarch)
W <- create.neighbourhood.array(m = m, sp = 3, type = "rook", torus =TRUE) 

## ----simulation, eval = TRUE, fig.width = 10----------------------------------
set.seed(1234)
y <- simSTARMAGARCH(parameters, n = 2500, m = m, W=W, burnin = 500)   # 1
#y <- simSTARMAGARCH(parameters, n = 2500, m = m, burnin = 500,        # 2
                    #type= "rook", torus = TRUE) 

## ----plotting_data------------------------------------------------------------
library(ggplot2)
library(reshape2)
ggplot(data=melt(y), aes(Var2,Var1))+ geom_raster(aes(fill = value))+
  scale_fill_gradient2(low="blue", mid="white", high = "red", midpoint = mean(y)) + xlab("Time") + ylab("Space") + 
  theme_bw()

## ----initial_parameters, eval = TRUE------------------------------------------
initial.parameters <- list(
  mu    = mean(y),
  phi   = matrix(c(.7, .01), ncol = 1), 
  theta = matrix(c(.01, .01), ncol = 1), 
  omega = 1,
  alpha = matrix(c(.01, .01,
                   .01, .01), ncol = 2), 
  beta  = matrix(c(.01, .01), ncol = 1)
)

## ----creatingLikelihood, eval = TRUE------------------------------------------
# Create likelihood object:
# f <- CreateLikelihood(y, W=W,
#                  init = apply(y,1,var), parameters=initial.parameters)
# To save time, we use the true parameters as initial values: 
f <- CreateLikelihood(y, W=W,
                 init = apply(y,1,var), parameters=parameters) 


## ----fitting_model, eval = TRUE-----------------------------------------------
fit <- fitSTARMAGARCH(f, data = y, print = FALSE)
summary(fit)

## ----plotting_fit,eval = TRUE, fig.show='hold', fig.height = 5----------------
plot(fit)
plot_garch(fit)

## ----comparing_results--------------------------------------------------------
round(rbind(unlist(parameters), round(coef(fit),4)), 4)

## ----STGARCH, eval = FALSE----------------------------------------------------
#  parameters <- list(
#    mu    = 0,
#    phi   = matrix(0, ncol = 1), #temporal order 1, spatial order 2
#    theta = matrix(0, ncol = 1),  #temporal order 1, spatial order 2
#    omega = .9,
#    alpha = matrix(c(.09, .13,
#                     .03, .07), ncol = 2), # temporal order 2, spatial order 2
#    beta  = matrix(c(.5, .15), ncol = 1) #temporal order 1, spatial order 2
#  )
#  y <- simSTARMAGARCH(parameters, n = 2500, m = c(25,1), W = W)
#  
#  initial.parameters <- list(
#    mu    = mean(y),
#    phi   = matrix(c(.7, .01), ncol = 1),
#    theta = matrix(c(.01, .01), ncol = 1),
#    omega = 1,
#    alpha = matrix(c(.01, .01,
#                     .01, .01), ncol = 2),
#    beta  = matrix(c(.01, .01), ncol = 1)
#  )
#  f <- CreateLikelihood(y, W=W,
#                   init = apply(y,1,var), parameters = initial.parameters)
#  fit <- fitSTARMAGARCH(f, data= y, print = FALSE)
#  summary(fit)
#  #>          Estimates          SD     Zscore       Pvalue
#  #> mu      0.01272265 0.009833892  1.2937556 9.787493e-02
#  #> phi1   -0.37973220 0.282238217 -1.3454315 8.924296e-02
#  #> phi2    0.37536087 0.783248661  0.4792359 3.158854e-01
#  #> theta1  0.37367411 0.283705318  1.3171206 9.389911e-02
#  #> theta2 -0.36493089 0.791497331 -0.4610639 3.223764e-01
#  #> omega   0.98304477 0.066042385 14.8850586 2.060627e-50
#  #> alpha1  0.08262752 0.006244003 13.2331007 2.825359e-40
#  #> alpha2  0.14105291 0.011764199 11.9900142 2.004181e-33
#  #> alpha3  0.03092677 0.007764663  3.9830156 3.402315e-05
#  #> alpha4  0.06686734 0.014844508  4.5045170 3.326203e-06
#  #> beta1   0.51746041 0.033736850 15.3381362 2.126149e-53
#  #> beta2   0.10477532 0.048024420  2.1817093 1.456550e-02
#  #>
#  #>
#  #> Standardized residual standard error:  1 .

## ----STGARCH_with_map, eval = FALSE-------------------------------------------
#  # Use "map"" to fix parameters to zero:
#  initial.parameters$mu <- 0
#  initial.parameters$phi <- matrix(0, ncol=1)
#  initial.parameters$theta <- matrix(0, ncol=1)
#  map <- list()
#  map$mu <- as.factor(NA)
#  map$phi <-  as.factor(NA)
#  map$theta <- as.factor(NA)
#  dim(map$theta) <- dim(map$phi) <- c(1,1)
#  map
#  #> $mu
#  #> [1] <NA>
#  #> Levels:
#  #>
#  #> $phi
#  #>      [,1]
#  #> [1,] <NA>
#  #> Levels:
#  #>
#  #> $theta
#  #>      [,1]
#  #> [1,] <NA>
#  #> Levels:
#  f <- CreateLikelihood(y, W=W,
#                   init = apply(y,1,var), parameters = initial.parameters, map = map)
#  
#  fit <- fitSTARMAGARCH(f, data= y, print = FALSE)
#  summary(fit)
#  #>         Estimates          SD    Zscore       Pvalue
#  #> omega  0.98281266 0.066025733 14.885297 2.053285e-50
#  #> alpha1 0.08264733 0.006242327 13.239826 2.583416e-40
#  #> alpha2 0.14108928 0.011764874 11.992417 1.946865e-33
#  #> alpha3 0.03103404 0.007765063  3.996624 3.212613e-05
#  #> alpha4 0.06653875 0.014844846  4.482279 3.692497e-06
#  #> beta1  0.51727367 0.033725151 15.337920 2.133229e-53
#  #> beta2  0.10519927 0.048014155  2.191005 1.422570e-02
#  #>
#  #>
#  #> Standardized residual standard error:  1 .
#  #>  AIC:  148332.7   BIC:  148396
#  
#  round(rbind(unlist(parameters)[-(1:3)], fit$coefficients),4)
#  #>       omega alpha1 alpha2 alpha3 alpha4  beta1  beta2
#  #> [1,] 0.9000 0.0900 0.1300  0.030 0.0700 0.5000 0.1500
#  #> [2,] 0.9828 0.0826 0.1411  0.031 0.0665 0.5173 0.1052

## ----STARMA_with_map, eval = FALSE--------------------------------------------
#  parameters <- list(
#    mu    = 5,
#    phi   = matrix(c(.8, .1), ncol = 1), #temporal order 1, spatial order 2
#    theta = matrix(c(-.42, -.2), ncol = 1),  #temporal order 1, spatial order 2
#    omega = 30,
#    alpha = matrix(0, ncol = 1), # temporal order 2, spatial order 2
#    beta  = matrix(0, ncol = 1) #temporal order 1, spatial order 2
#  )
#  y <- simSTARMAGARCH(parameters, n = 2500, m = c(25,1), W = W)
#  initial.parameters <- list(
#    mu    = mean(y),
#    phi   = matrix(c(.7, .01), ncol = 1),
#    theta = matrix(c(.01, .01), ncol = 1),
#    omega = 1,
#    alpha = matrix(c(.01, .01,
#                     .01, .01), ncol = 2),
#    beta  = matrix(c(.01, .01), ncol = 1)
#  )
#  f <- CreateLikelihood(y, W=W,
#                   init = apply(y,1,var), parameters = initial.parameters)
#  fit <- fitSTARMAGARCH(f, data= y, print = FALSE)
#  summary(fit)
#  #>            Estimates           SD        Zscore       Pvalue
#  #> mu      5.0714676601  0.076651590  6.616259e+01 0.000000e+00
#  #> phi1    0.7977517884  0.005226702  1.526301e+02 0.000000e+00
#  #> phi2    0.1119027906  0.010694371  1.046371e+01 6.339603e-26
#  #> theta1 -0.4178454614  0.007787193 -5.365803e+01 0.000000e+00
#  #> theta2 -0.2233686304  0.016775178 -1.331543e+01 9.414894e-41
#  #> omega  29.8237638244 15.355764525  1.942187e+00 2.605724e-02
#  #> alpha1  0.0000000100  0.002948183  3.391920e-06 4.999986e-01
#  #> alpha2  0.0005736162          NaN           NaN          NaN
#  #> alpha3  0.0001235004  0.002756660  4.480074e-02 4.821331e-01
#  #> alpha4  0.0001525837  0.009304228  1.639939e-02 4.934579e-01
#  #> beta1   0.0000000100          NaN           NaN          NaN
#  #> beta2   0.0000000100          NaN           NaN          NaN
#  #>
#  #>
#  #> Standardized residual standard error:  1 .
#  #>  AIC:  194676     BIC:  194784.5
#  
#  # Use "map"" to fix parameters to zero:
#  initial.parameters$alpha <- matrix(0, ncol=1)
#  initial.parameters$beta <- matrix(0, ncol=1)
#  map <- list()
#  map$alpha <-  as.factor(NA)
#  map$beta <- as.factor(NA)
#  dim(map$alpha) <- dim(map$beta) <- c(1,1)
#  map
#  #> $alpha
#  #>      [,1]
#  #> [1,] <NA>
#  #> Levels:
#  #>
#  #> $beta
#  #>      [,1]
#  #> [1,] <NA>
#  #> Levels:
#  
#  f <- CreateLikelihood(y, W=W,
#                   init = apply(y,1,var), parameters = initial.parameters, map = map)
#  fit <- fitSTARMAGARCH(f, data= y, print = FALSE)
#  summary(fit)
#  #>         Estimates          SD    Zscore       Pvalue
#  #> mu      5.0714315 0.076639614  66.17246 0.000000e+00
#  #> phi1    0.7977492 0.005225273 152.67130 0.000000e+00
#  #> phi2    0.1118951 0.010696428  10.46098 6.525081e-26
#  #> theta1 -0.4178449 0.007784903 -53.67374 0.000000e+00
#  #> theta2 -0.2233527 0.016775596 -13.31414 9.577871e-41
#  #> omega  29.8436922 0.168888961 176.70600 0.000000e+00
#  #>
#  #>
#  #> Standardized residual standard error:  1 .
#  #>  AIC:  194664     BIC:  194718.2
#  
#  round(rbind(unlist(parameters)[1:6],coef(fit)), 4)
#  #>          mu   phi1   phi2  theta1  theta2   omega
#  #> [1,] 5.0000 0.8000 0.1000 -0.4200 -0.2000 30.0000
#  #> [2,] 5.0714 0.7977 0.1119 -0.4178 -0.2234 29.8437


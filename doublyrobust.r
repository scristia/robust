#### Simulate covariates
i=1
mu.est <- matrix(rep(NA, 1000*4), ncol=4)
mu.est.rob <- matrix(rep(NA, 1000*4), ncol=4)
# mu.est <- t(replicate(1000,{
for(i in 1:1000) {
    set.seed(i)
    Z <- cbind(1, matrix(rnorm(1000*4), nrow=1000, ncol=4))

    x1 <- exp(Z[,2]/2)
    x2 <- Z[,3]/(1 + exp(Z[,2])) + 10
    x3 <- (Z[,2] * Z[,4] / 25 + 0.6)^3
    x4 <- (Z[,4] + Z[,5] + 20)^2

    X <- cbind(1, x1, x2, x3, x4)

    ### True outcome
    #     y <- Z %*% c(210, 24.4, 13.7, 13.7, 13.7) + rnorm(1000, 0, 1)
    y <- Z %*% c(210, 24.4, 13.7, 13.7, 13.7)
    y[1:100] <- y[1:100] + rnorm(100, 0, 7^2)
    y[101:1000] <- y[101:1000] + rnorm(900, 0, 1)

    ### True propensity score
    expit <- function(x) exp(x)/(1 + exp(x))
    pi0 <- expit(-Z[,2] + 0.5*Z[,3] - 0.25*Z[,4] - 0.1*Z[,5])

    R <- rbinom(1000, 1, pi0)
    #### True
    beta.h <- solve((t(Z) %*% Z)) %*% t(Z) %*% y

    ### complete cases (correct regression model)
    beta.or <- solve(t(Z[R==1,]) %*% Z[R==1,]) %*% t(Z[R==1,]) %*% y[R==1]

    ### complete cases (incorrect regression model)
    beta.or.wrong <- solve(t(X[R==1,]) %*% X[R==1,]) %*% t(X[R==1,]) %*% y[R==1]

    fit.robust <- rlm(y[R==1] ~ Z[R==1,-1])
    beta.robust <- fit.robust$coeff

    yfit.true <- y
    yfit.true[R==1] <- fit.robust$fitted.values

    fit.robust.wrong <- rlm(y[R==1] ~ X[R==1,-1])
    beta.robust.wrong<- fit.robust.wrong$coeff

    yfit.wrong <- y
    yfit.wrong[R==1] <- fit.robust.wrong$fitted.values

    ### estimate propensity score
    gamma.h <- coef(glm(R ~ Z[,-1], family="binomial"))
    ps <- expit(Z %*% gamma.h)
    mu.ipw <- mean(R * y / ps)

    ### incorect propensity score
    ps.wrong <- expit(X %*% coef(glm(R ~ X[,-1], family="binomial")))

    #### Doubly Robust estimator
    mu.dr <- mean(R*y/ps - (R - ps)/ps * Z %*% beta.or)
    mu.or.wrong <- mean(R*y/ps - (R - ps)/ps * X %*% beta.or.wrong)
    mu.ps.wrong <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.or)
    mu.both.wrong <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.or.wrong)

#     mu.dr.rob <- mean(R*yfit.true/ps - (R - ps)/ps * Z %*% beta.robust)
#     mu.or.wrong.rob <- mean(R*yfit.wrong/ps - (R - ps)/ps * X %*% beta.robust.wrong)
#     mu.ps.wrong.rob <- mean(R*yfit.true/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust)
#     mu.both.wrong.rob <- mean(R*yfit.wrong/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong)

    w.fit <- rep(1, length(y))
    w.fit[R==1] <- fit.robust$w
    w.fit.wrong <- rep(1, length(y))
    w.fit.wrong[R==1] <- fit.robust.wrong$w

    mu.dr.rob <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta.robust, w.fit)
    mu.or.wrong.rob <- weighted.mean(R*y/ps - (R - ps)/ps * X %*% beta.robust.wrong, w.fit.wrong)
    mu.ps.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust, w.fit)
    mu.both.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong, w.fit.wrong)

    mu.est[i,] <- c(mu.dr, mu.or.wrong, mu.ps.wrong, mu.both.wrong)
    mu.est.rob[i,] <- c(mu.dr.rob, mu.or.wrong.rob, mu.ps.wrong.rob, mu.both.wrong.rob)
    i <- i+1
}
# }))

colnames(mu.est) <- c("mu.dr", "mu.or.wrong", "mu.ps.wrong", "mu.both.wrong")
bias <- 210 - colMeans(mu.est)
mse <- bias^2 + apply(mu.est, 2, var)
bias
sqrt(mse)

colnames(mu.est.rob) <- c("mu.dr.rob", "mu.or.wrong.rob", "mu.ps.wrong.rob", "mu.both.wrong.rob")
bias <- 210 - colMeans(mu.est.rob)
mse <- bias^2 + apply(mu.est.rob, 2, var)
bias
sqrt(mse)



### seed=484
i=484

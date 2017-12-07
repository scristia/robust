library(MASS)
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
    #     y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7) + rnorm(1000, 0, 1)
    y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7)
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

#     fit.robust <- rlm(Z[R==1,], y[R==1], psi=psi.hampel)
    fit.robust <- rlm(Z[R==1,], y[R==1], psi=psi.hampel, maxit=50)
    beta.robust <- fit.robust$coeff

    yfit.true <- y
    yfit.true[R==1] <- fit.robust$fitted.values

#     fit.robust.wrong <- rlm(X[R==1,], y[R==1], psi=psi.hampel)
    fit.robust.wrong <- rlm(X[R==1,], y[R==1], psi=psi.hampel, maxit=50)
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

#     mu.dr <- mean(R*y/ps - (R - ps)/ps * Z %*% beta.robust)
#     mu.or.wrong <- mean(R*y/ps - (R - ps)/ps * X %*% beta.robust.wrong)
#     mu.ps.wrong <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust)
#     mu.both.wrong <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong)

    w.fit <- rep(1, length(y))
    w.fit[R==1] <- fit.robust$w
    w.fit.wrong <- rep(1, length(y))
    w.fit.wrong[R==1] <- fit.robust.wrong$w

#     mu.dr.rob <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta.or, w.fit)
#     mu.or.wrong.rob <- weighted.mean(R*y/ps - (R - ps)/ps * X %*% beta.or.wrong, w.fit.wrong)
#     mu.ps.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.or, w.fit)
#     mu.both.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.or.wrong, w.fit.wrong)

    mu.dr.rob <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta.robust, w.fit)
    mu.or.wrong.rob <- weighted.mean(R*y/ps - (R - ps)/ps * X %*% beta.robust.wrong, w.fit.wrong)
    mu.ps.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust, w.fit)
    mu.both.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong, w.fit.wrong)

#     mu.dr.rob <- mean(R*y/ps - (R - ps)/ps * Z %*% beta.robust)
#     mu.or.wrong.rob <- mean(R*y/ps - (R - ps)/ps * X %*% beta.robust.wrong)
#     mu.ps.wrong.rob <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust)
#     mu.both.wrong.rob <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong)

    mu.est[i,] <- c(mu.dr, mu.or.wrong, mu.ps.wrong, mu.both.wrong)
    mu.est.rob[i,] <- c(mu.dr.rob, mu.or.wrong.rob, mu.ps.wrong.rob, mu.both.wrong.rob)
    i <- i+1
}
# }))

colnames(mu.est) <- c("mu.dr", "mu.or.wrong", "mu.ps.wrong", "mu.both.wrong")
bias <- 210 - colMeans(mu.est[-484,])
mse <- bias^2 + apply(mu.est[-484,], 2, var)
bias
sqrt(mse)

colnames(mu.est.rob) <- c("mu.dr.rob", "mu.or.wrong.rob", "mu.ps.wrong.rob", "mu.both.wrong.rob")
bias <- 210 - colMeans(mu.est.rob[-484,])
mse <- bias^2 + apply(mu.est.rob[-484,], 2, var)
bias
sqrt(mse)

### seed=484

### No Outliers
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
    y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7) + rnorm(1000, 0, 1)

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

    fit.robust <- rlm(Z[R==1,], y[R==1], psi=psi.hampel, maxit=50)
    beta.robust <- fit.robust$coeff

    yfit.true <- y
    yfit.true[R==1] <- fit.robust$fitted.values

    fit.robust.wrong <- rlm(X[R==1,], y[R==1], psi=psi.hampel, maxit=50)
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
    #     y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7) + rnorm(1000, 0, 1)
    y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7)
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

#     fit.robust <- rlm(Z[R==1,], y[R==1], psi=psi.hampel)
    fit.robust <- rlm(Z, y, psi=psi.hampel, maxit=50)
    beta.robust <- fit.robust$coeff

    yfit.true <- y
    yfit.true <- fit.robust$fitted.values

#     fit.robust.wrong <- rlm(X[R==1,], y[R==1], psi=psi.hampel)
    fit.robust.wrong <- rlm(X, y, psi=psi.hampel, maxit=50)
    beta.robust.wrong<- fit.robust.wrong$coeff

    yfit.wrong <- y
    yfit.wrong <- fit.robust.wrong$fitted.values

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

#     mu.dr <- mean(R*y/ps - (R - ps)/ps * Z %*% beta.robust)
#     mu.or.wrong <- mean(R*y/ps - (R - ps)/ps * X %*% beta.robust.wrong)
#     mu.ps.wrong <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust)
#     mu.both.wrong <- mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong)

    w.fit <- rep(1, length(y))
    w.fit <- fit.robust$w
    w.fit.wrong <- rep(1, length(y))
    w.fit.wrong <- fit.robust.wrong$w

#     mu.dr.rob <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta.or, w.fit)
#     mu.or.wrong.rob <- weighted.mean(R*y/ps - (R - ps)/ps * X %*% beta.or.wrong, w.fit.wrong)
#     mu.ps.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.or, w.fit)
#     mu.both.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.or.wrong, w.fit.wrong)

    mu.dr.rob <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta.robust, w.fit)
    mu.or.wrong.rob <- weighted.mean(R*y/ps - (R - ps)/ps * X %*% beta.robust.wrong, w.fit.wrong)
    mu.ps.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * Z %*% beta.robust, w.fit)
    mu.both.wrong.rob <- weighted.mean(R*y/ps.wrong - (R - ps.wrong)/ps.wrong * X %*% beta.robust.wrong, w.fit.wrong)


    mu.est[i,] <- c(mu.dr, mu.or.wrong, mu.ps.wrong, mu.both.wrong)
    mu.est.rob[i,] <- c(mu.dr.rob, mu.or.wrong.rob, mu.ps.wrong.rob, mu.both.wrong.rob)
    i <- i+1
}

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




#####################
## data: Y (observations), R (missingness indicator), Z (auxilliary variables)
## parameters = mu, beta1 .. betap, gamma1 .. gammap
## transformations = u1, u2, u3, A
## weights w
## constants: a, p(dim(parameter space))

##### Startings values
w <- rep(1, length(y))
beta <- lm(y[R==1] ~ Z[R==1,-1], weights=w[R==1])$coeff
gamma <- coef(glm(R ~ Z[,-1], family="binomial"), weights=w)
ps <- expit(Z %*% gamma)
mu <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta, w)
a <- 1.48

beta.start <- beta
gamma.start <- gamma
ps.start <- ps
mu.start <- mu

eps = 0.0001
maxit=100

expit <- function(x) exp(x)/(1 + exp(x))
robust.dr <- function(y, R, Z, X, mu.start, beta.start, gamma.start, correct="both", maxit=50, eps=0.0001) {
    mu.update <- mu.start
    beta.update <- beta.start
    gamma.update <- gamma.start
    ps <- expit(Z %*% gamma)

    if(correct=="OR" | correct=="both") mu.or <- Z %*% beta
    if(correct=="PS" | correct == "none") {
        mu.or <- X %*% beta
    }

    for(i in 1:maxit) {
        ee.u1 <- function(y, R, mu.or, ps, w, mu) {
            w*(R*y/ps - (R - ps)/ps * mu.or - mu)
        }
        u1 <- ee.u1(y, R, mu.or, ps, w, mu.update)

        ee.u2 <- function(R, Z, ps, w) {
            apply(Z, 2, function(z) w * (R - ps)*z)
        }
        if(correct=="PS" | correct=="both") u2 <- ee.u2(R, Z, ps, w)

        if(correct=="OR" | correct == "none") {
            u2 <- ee.u2(R, X, ps, w)
        }

        ee.u3 <- function(y, R, Z, beta, w) {
            apply(Z, 2, function(z) w * R * (y - Z %*% beta)*z)
        }

        if(correct=="OR" | correct == "both") u3 <- ee.u3(y, R, Z, beta, w)

        if(correct=="PS" | correct == "none") {
            u3 <- ee.u3(y, R, X, beta, w)
        }

        ### x0 and x2 look to be redundant in u2/u3 leading to singular matrix. Why?
        ### testing with this
#         u2 <- u2[, -1]
#         u3 <- u3[, -1]

        A.list <- lapply(1:1000, function(x) {
                             U <- c(u1[x],  u2[x,], u3[x,])
                             U %*% t(U)
        })
        A <- Reduce("+", A.list) / length(A.list)

        ### A singular when OR and PS models both incorrect.
        U <- cbind(u1, u2, u3)
        findWeights <- function(U, A, a) {
#             norm <-  apply(U, 1, function(u) sqrt(t(u) %*% solve(A) %*% u))
            norm <-  apply(U, 1, function(u) sqrt(t(u) %*% chol2inv(chol(A)) %*% u))
            w0 <- a * sqrt(ncol(A)) / norm
            ifelse(w0 > 1, 1, w0)
            w <- psi.hampel(norm/(a*median(norm)))

            ### there needs to be some sort of quality assurance check here
            ### because the norm behaves funny for very extreme values.
            ### I just give these values weight 0 for now.

            ## extreme outliers get weight 0 on odd iterations and weight 1 on even.
            ## fix so if an observation gets 0 weight, it always keeps that weight.
            w[which(norm == 0)] <- 0
            w
        }
        w <- findWeights(U, A, a)


        theta.prev <- c(beta.update, gamma.update, mu.update)
        if(correct=="OR"|correct=="both") {
            beta.update <- lm(y[R==1] ~ Z[R==1,-1], weights=w[R==1])$coeff
            mu.or <- Z %*% beta.update
        }

        if(correct=="PS"|correct=="none") {
            beta.update <- lm(y[R==1] ~ X[R==1,-1], weights=w[R==1])$coeff
            mu.or <- X %*% beta.update
        }

        if(correct=="PS"|correct=="both") {
            gamma.update <- coef(glm(R ~ Z[,-1], family="binomial", weights=w))
            ps <- expit(Z %*% gamma.update)
        }
        if(correct=="OR"|correct=="none") {
            gamma.update <- coef(glm(R ~ X[,-1], family="binomial", weights=w))
            ps <- expit(X %*% gamma.update)
        }

        mu.update <- weighted.mean(R*y/ps - (R - ps)/ps * mu.or, w)
        theta <- c(beta.update, gamma.update, mu.update)

        converged <- sqrt(sum(((theta - theta.prev)^2)))/sqrt(sum(((theta)^2))) < eps
        if(converged) break
    }
    mu.update
}

mu.both.correct  <- rep(NA, 500)
mu.ps.correct <- rep(NA, 500)
mu.or.correct <-  rep(NA, 500)
mu.both.incorrect <-  rep(NA, 500)
mu <- rep(NA, 500)

for(i in 1:500) {
    set.seed(i)
    Z <- cbind(1, matrix(rnorm(1000*4), nrow=1000, ncol=4))

    x1 <- exp(Z[,2]/2)
    x2 <- Z[,3]/(1 + exp(Z[,2])) + 10
    x3 <- (Z[,2] * Z[,4] / 25 + 0.6)^3
    x4 <- (Z[,4] + Z[,5] + 20)^2

    X <- cbind(1, x1, x2, x3, x4)

    ### True outcome
    y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7) + rnorm(1000, 0, 1)
    ###   ^^ no outliers ^^ -- vv outliers vv
#     y <- Z %*% c(210, 27.4, 13.7, 13.7, 13.7)
#     y[1:100] <- y[1:100] + rnorm(100, 0, 7^2)
#     y[101:1000] <- y[101:1000] + rnorm(900, 0, 1)

    expit <- function(x) exp(x)/(1 + exp(x))
    pi0 <- expit(-Z[,2] + 0.5*Z[,3] - 0.25*Z[,4] - 0.1*Z[,5])

    R <- rbinom(1000, 1, pi0)
    ##### Startings values
    w <- rep(1, length(y))
    beta <- lm(y[R==1] ~ Z[R==1,-1], weights=w[R==1])$coeff
    gamma <- coef(glm(R ~ Z[,-1], family="binomial"), weights=w)
    ps <- expit(Z %*% gamma)
    mu <- weighted.mean(R*y/ps - (R - ps)/ps * Z %*% beta, w)
    a <- 1.48

    beta.start <- beta
    gamma.start <- gamma
    ps.start <- ps
    mu.start <- mu

    eps = 0.0001
    maxit=16
    ###################

    mu.both.correct[i] <- robust.dr(y, R, Z, X, mu.start, beta.start, gamma.start, correct="both", maxit=maxit, eps=0.001)
    mu.ps.correct[i] <- robust.dr(y, R, Z, X, mu.start, beta.start, gamma.start, correct="PS", maxit=maxit, eps=0.001)
    mu.or.correct[i] <- robust.dr(y, R, Z, X, mu.start, beta.start, gamma.start, correct="OR", maxit=maxit, eps=0.001)
    mu.both.incorrect[i] <- robust.dr(y, R, Z, X, mu.start, beta.start, gamma.start, correct="none", maxit=maxit, eps=0.001)
}

# bias <- 210 - mean(mu.or.correct[-484])
# mse <- bias^2 + var(mu.or.correct[-484])
bias.both.correct <- 210 - mean(mu.both.correct)
mse.both.correct <- bias.both.correct^2 + var(mu.both.correct)
bias.both.correct
sqrt(mse.both.correct)

bias.ps.correct <- 210 - mean(mu.ps.correct)
mse.ps.correct <- bias.ps.correct^2 + var(mu.ps.correct)
bias.ps.correct
sqrt(mse.ps.correct)

bias.or.correct <- 210 - mean(mu.or.correct)
mse.or.correct <- bias.or.correct^2 + var(mu.or.correct)
bias.or.correct
sqrt(mse.or.correct)

bias.both.incorrect <- 210 - mean(mu.both.incorrect)
mse.both.incorrect <- bias.both.incorrect^2 + var(mu.both.incorrect)
bias.both.incorrect
sqrt(mse.both.incorrect)

bias <- c(bias.both.correct, bias.ps.correct, bias.or.correct, bias.both.incorrect)
rmse <- sqrt(c(mse.both.correct, mse.ps.correct, mse.or.correct, mse.both.incorrect))

tab <- cbind(bias, rmse)
rownames(tab) <- c("Both Correct", "OR Wrong", "PS Wrong", "Both Wrong")

tab

#### This doesn't work. Debug.
#### A singular. Happens on first iteration when all weights are 1. Must be a bug,
#### since this is computationally equivalent to the standard (nonrobust) DR estimator.
mu.both.incorrect <- robust.dr(y, R, Z, X, mu.start, beta.start, gamma.start, correct="none", maxit=1, eps=0.001)

debug(robust.dr)
robust.dr(y, R, Z, X, mu.start, beta.start, gamma.start, correct="none", maxit=1, eps=0.001)

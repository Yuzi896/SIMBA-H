## Generate data
# source("Main.r")
K <- 2 ## number of clusters
p <- c(300, 500) ## number of features for each view
n <- 100
nbrcov <- 4
OutcType <- "continuous"
# OutcType <- "binary"
# OutcType <- "censored"
dat <- iGenDataHMM(
    OrderFeatured = TRUE, K = K, NbrPlat = 2, OutcType = OutcType,
    MeanY = c(rep(1,floor(K/2)),rep(-1,K-floor(K/2))), pc = nbrcov,
    betacov = 1.0, MeanOverExpre = c(1:K) / 2,
    MeanUnderExpre = -c(1:K) / 2, p = p, seed = 1, n = n,
    SigmaOE = rep(1, K), SigmaUE = rep(1, K)
)

#### Generate the test data
ntest <- 50
datTest <- iGenDataHMM(
    OrderFeatured = TRUE, K = K, NbrPlat = 2, OutcType = OutcType,
    MeanY =  c(rep(1,floor(K/2)),rep(-1,K-floor(K/2))), pc = nbrcov,
    betacov = 1.0, MeanOverExpre = c(1:K) / 2,
    MeanUnderExpre = -c(1:K) / 2, p = p, seed = 2, n = ntest,
    SigmaOE = rep(1, K), SigmaUE = rep(1, K)
)



#ntrain <- 1:(N / 2)
#dataListT <- lapply(dat$X, function(x) x[ntrain, ])
#dataListTest <- lapply(dat$X, function(x) x[-ntrain, ])
#CovT <- dat$Cov[ntrain, ]
#CovTest <- dat$Cov[-ntrain, ]

K1 <- K
system.time(
    result <- iBiclust(
        Method = 1, dataList = dat$X, Youtcome = dat$Y, OutcType = OutcType,
        Cov = dat$Cov, sample = 3000, burnin = 1000, K = K1, normalized = FALSE,
        seed = 1, hyperSigma = c(1, 1), alphaP = rep(1, K1),
        delt = rep(1, 3), TruncVal = 0.1, dataListTest = datTest$X,
        CovTest = datTest$Cov
    )
)

if (OutcType == "continuous") {
    MSE <- mean(apply(result$YpredTestSample, 2, mean) - mean(datTest$Y))
    MSE.LB <- mean(apply(result$YpredTestSample, 2, function(x) quantile(x, probs = 0.025)) - mean(datTest$Y))
    MSE.UB <- mean(apply(result$YpredTestSample, 2, function(x) quantile(x, probs = 0.975)) - mean(datTest$Y))
    print(MSE)
} else if (OutcType == "binary") {
    MCE <- 1 - sum(diag(table(datTest$Y, apply(result$YpredTestSample, 2, mean) > 0))) / ntest
    print(MCE)
} else if (OutcType == "censored") {
    CINDEX <- Cindex(Tpredict = apply(result$YpredTestSample, 2, mean), Ttest = datTest$Y[, 1], deltatest = datTest$Y[, 2])
    print(paste("Cindex=", CINDEX))
}

quantile(result$CoefySample[, 1], probs = c(0.025, 0.5, 0.975))
quantile(result$CoefySample[, 2], probs = c(0.025, 0.5, 0.975))
summary(result$Betasample[, 1])
quantile(result$Betasample[, 1], probs = c(0.025, 0.5, 0.975))
quantile(result$Betasample[, 2], probs = c(0.025, 0.5, 0.975))
result$MuMean[1, , ] ## Means within clusters for platform 1
result$MuMean[2, , ] ## Mean within clusters for platform 2
result$TransMatMean[1, , ] ## platform 1
result$TransMatMean[2, , ] ## platform 2
result$probZMean

result$KappaMean[[1]][1, , ] ## prob values for platform 1 in cluster 1
result$KappaMean[[1]][2, , ] ## prob values for platform 1 in cluster 2
plot(result$logpostsample,type = "l")
print(result$DIC)

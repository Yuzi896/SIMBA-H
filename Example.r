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

EvaluateTest(result,datTest$Y,OutcType,threshold = 0,return_all = TRUE)

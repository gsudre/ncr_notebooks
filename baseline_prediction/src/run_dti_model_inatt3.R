args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Arguments: seed[-1] cpuDiff tuneLength mymod root_fname", call.=FALSE)
} else {
  if (args[1] == '-1') {
    myseed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  } else {
    myseed = as.integer(args[1])
  }
  cpuDiff = as.integer(args[2])
  tuneLength = as.integer(args[3])
  mymod = args[4] # AdaBoost.M1, AdaBag, ada
  root_fname = args[5]
}

###########

fname = sprintf('%s_%04d.log', root_fname, myseed)
sink(fname, append=FALSE, split=TRUE)
source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')

tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
rm_me = (tract_data$fa_avg < .4 | tract_data$ad_avg < 1.18 | tract_data$rd_avg > .65 | tract_data$rd_avg < .5 |
           tract_data$norm.trans > .45 | tract_data$norm.rot > .008 | tract_data$goodSlices < 45 |
           tract_data$goodSlices > 70)
tract_data = tract_data[!rm_me, ]
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf = gf[gf$BASELINE=='BASELINE', ]
my_ids = intersect(gf$MRN, tract_data$MRN)
merged = mergeOnClosestDate(gf, tract_data, my_ids)
rm_me = abs(merged$dateX.minus.dateY.months) > 12
merged = merged[!rm_me, ]

phen_vars = c(which(grepl("^FA_", colnames(merged))),
              which(grepl("^AD_", colnames(merged))),
              which(grepl("^RD_", colnames(merged))),
              which(grepl("^MO_", colnames(merged)))
)
X = merged[, phen_vars]
y = merged$inatt3_named
y = factor(y, levels=c('low', 'medium', 'high'))

# save X and y if not already done so
fname = sprintf('%s_Xy.RData', root_fname)
if(!file.exists(fname)){
  save(X, y, file=fname, compress=T)
}

set.seed(myseed)
split <- createDataPartition(y, p = .8, list = FALSE)
Xtrain <- X[ split, ]
ytrain <- y[ split ]
Xtest  <- X[-split, ]
ytest = y[-split]

pp = preProcess(Xtrain, method=c('YeoJohnson', 'center', 'scale'))
filtXtrain = predict(pp, Xtrain)
nzv = nearZeroVar(filtXtrain)
nzv
if (length(nzv) > 0) {
    filtXtrain = filtXtrain[, -nzv]
}
correlations = cor(filtXtrain)

highCorr = findCorrelation(correlations, cutoff=.75)
length(highCorr)
noncorrXtrain = filtXtrain[, -highCorr]
noncorrXtest = predict(pp, Xtest)[, -highCorr]

library(pROC)

set.seed(myseed)
index <- createMultiFolds(ytrain, k = 5, times = 5)

library(doMC)
ncpus <- detectBatchCPUs()
njobs <- ncpus - cpuDiff
registerDoMC(njobs)
selFunc = 'best'  # maybe try oneSE and tolerance as well?

set.seed(myseed)
fullCtrl <- trainControl(method = "repeatedcv",
                         index = index,
                         search='grid',
                         summaryFunction = multiClassSummary,
                         classProbs = TRUE)

ptm <- proc.time()
m1 <- train(noncorrXtrain, ytrain,
            method = mymod,
            trControl = fullCtrl,
            tuneLength = tuneLength,
            metric = 'Mean_ROC')
print(proc.time() - ptm)
m1
getTrainPerf(m1)
pred = predict(m1, noncorrXtest)
postResample(pred, ytest)
multiclass.roc(as.numeric(ytest), as.numeric(pred))

fname = sprintf('%s_%04d.RData', root_fname, myseed)
save_list = c('m1', 'myseed', 'index', 'split')
save(list=save_list, file=fname)

sink()

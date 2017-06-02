args = commandArgs(trailingOnly=TRUE)
s=as.numeric(args[1])
ntimes = 50
model = args[2]
mysx = args[3]
myseed = 1234
tuneLength = 10
cpuDiff = 0
adhdOnly = F

out_fname = sprintf('/data/NCR_SBRB/loocv/%s_sx_slope_%s/s%03d.log', mysx, model, s)
# out_fname = sprintf('~/tmp/s%03d.log', s)
sink(out_fname, append=FALSE, split=TRUE)

get_SX_slope = function (df, ids) {
  inatt = c()
  hi = c()
  mrns = c()
  for (s in ids) {
    idx = df$MRN==s
    # proceed if we have more than one observation in the data
    if (sum(idx) >= 2) {
      inatt = c(inatt, lm(SX_inatt ~ age, data=df[idx,])$coefficients[[2]])
      hi = c(hi, lm(SX_HI ~ age, data=df[idx,])$coefficients[[2]])
      mrns = c(mrns, s)
    }
  }
  res = cbind(inatt, hi)
  colnames(res) = c('inatt', 'hi')
  rownames(res) = mrns
  return(res)
}
source('~/ncr_notebooks/baseline_prediction/src/load_voting_data.R')

library(caret)
library(doMC)
ncpus <- detectBatchCPUs()
njobs <- ncpus - cpuDiff
registerDoMC(njobs)
if (adhdOnly) {
  idx = which(merged$DX_BASELINE!='NV')
  ids = merged[idx,]$MRN
} else {
  idx = 1:nrow(merged)
  ids = merged$MRN
}
slopes = get_SX_slope(gf, ids)
y = slopes[, mysx]
X = cbind(prs, geospatial, neuropsych, struct_rois, dti_tracts)[idx, ]

Xtrain <- X[-s, ]
ytrain <- y[-s]
Xtest  <- X[s, ]
ytest  <- y[s]

# remove training examples that are all NaN. Note that this will never happen for testing!
rm_me = rowSums(is.na(Xtrain)) == ncol(Xtrain)
Xtrain = Xtrain[!rm_me,]
ytrain = ytrain[!rm_me]

set.seed(myseed)
fullCtrl <- trainControl(method = "boot",
                         number = ntimes,
                         savePredictions="final",
                         summaryFunction=defaultSummary)
set.seed(myseed)
mymod = train(Xtrain, ytrain,
               tuneLength=tuneLength,
               trControl=fullCtrl,
               method=model,
               preProcess = c('medianImpute'))

means = colMeans(slopes)
print(sprintf('NIR %s:', mysx))
RMSE(means[mysx], slopes[, mysx])

pred = predict(mymod, newdata=Xtest)

mymod
myseed
tuneLength
ntimes
adhdOnly
pred
sink()








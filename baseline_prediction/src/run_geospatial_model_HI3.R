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

geo_data = read.csv('~/data/baseline_prediction/stripped/geospatial.csv')
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf = gf[gf$BASELINE=='BASELINE', ]
merged = merge(gf, geo_data, by='MRN')
# some variables are being read as numeric...
merged$Home_Price = as.numeric(merged$Home_Price)
merged$Fam_Income = as.numeric(merged$Fam_Income)
merged$Crime_Rate = as.numeric(merged$Crime_Rate)

phen_vars = c('SES', 'Home_Type', 'Home_Price', 'Fam_Income', 'Pop_BPL', 'Fam_BPL', 'Pub_School',
              'Crime_Rate', 'Green_Space', 'Park_Access', 'Air_Quality', 'Obesity_Rate',
              'Food_Index', 'Exercise_Access', 'Excessive_Drinking')
keep_me = c()
for (v in phen_vars) {
  keep_me = c(keep_me, which(colnames(merged) == v))
}
X = merged[, keep_me]
y = merged$HI3_named
y = factor(y)

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

pp = preProcess(Xtrain, method=c('YeoJohnson', 'center', 'scale', 'knnImpute'))
filtXtrain = predict(pp, Xtrain)
nearZeroVar(filtXtrain)
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
            metric = 'ROC')
print(proc.time() - ptm)
m1
pred = predict(m1, noncorrXtest)
postResample(pred, ytest)
roc(as.numeric(ytest), as.numeric(pred))

fname = sprintf('%s_%04d.RData', root_fname, myseed)
save_list = c('m1', 'myseed', 'index', 'split')
save(list=save_list, file=fname)

sink()

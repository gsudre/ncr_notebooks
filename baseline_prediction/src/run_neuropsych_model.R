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

beery_data = read.csv('~/data/baseline_prediction/stripped/beeryVMI.csv')
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf = gf[gf$BASELINE=='BASELINE', ]
my_ids = intersect(gf$MRN, beery_data$Medical.Record...MRN)
mbeery = mergeOnClosestDate(gf, beery_data, my_ids, y.date='record.date.collected', y.id='Medical.Record...MRN')
rm_me = abs(mbeery$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mbeery), nrow(mbeery)-sum(rm_me)))
mbeery = mbeery[!rm_me, ]
mbeery$dateClinical.minus.dateBeery.months = mbeery$dateX.minus.dateY.months
mbeery$dateX.minus.dateY.months = NULL

cpt_data = read.csv('~/data/baseline_prediction/stripped/cpt.csv')
my_ids = intersect(gf$MRN, cpt_data$MRN)
mcpt = mergeOnClosestDate(gf, cpt_data, my_ids)
rm_me = abs(mcpt$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mcpt), nrow(mcpt)-sum(rm_me)))
mcpt = mcpt[!rm_me, ]
mcpt$dateClinical.minus.dateCPT.months = mcpt$dateX.minus.dateY.months
mcpt$dateX.minus.dateY.months = NULL

iq_data = read.csv('~/data/baseline_prediction/stripped/iq.csv')
my_ids = intersect(gf$MRN, iq_data$Medical.Record...MRN)
miq = mergeOnClosestDate(gf, iq_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(miq$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(miq), nrow(miq)-sum(rm_me)))
miq = miq[!rm_me, ]
miq$dateClinical.minus.dateIQ.months = miq$dateX.minus.dateY.months
miq$dateX.minus.dateY.months = NULL

wisc_data = read.csv('~/data/baseline_prediction/stripped/wisc.csv')
my_ids = intersect(gf$MRN, wisc_data$Medical.Record...MRN)
mwisc = mergeOnClosestDate(gf, wisc_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(mwisc$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mwisc), nrow(mwisc)-sum(rm_me)))
mwisc = mwisc[!rm_me, ]
mwisc$dateClinical.minus.dateWISC.months = mwisc$dateX.minus.dateY.months
mwisc$dateX.minus.dateY.months = NULL

wj_data = read.csv('~/data/baseline_prediction/stripped/wj.csv')
my_ids = intersect(gf$MRN, wj_data$Medical.Record...MRN)
mwj = mergeOnClosestDate(gf, wj_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(mwj$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mwj), nrow(mwj)-sum(rm_me)))
mwj = mwj[!rm_me, ]
mwj$dateClinical.minus.dateWJ.months = mwj$dateX.minus.dateY.months
mwj$dateX.minus.dateY.months = NULL

merged = merge(mwj, mbeery, by='MRN', all.x = T, all.y = T)
merged = merge(merged, miq, by='MRN', all.x = T, all.y = T)
merged = merge(merged, mwisc, by='MRN', all.x = T, all.y = T)
merged = merge(merged, mcpt, by='MRN', all.x = T, all.y = T)

phen_vars = c('FSIQ',
              # CPT
              'N_of_omissions', 'N_commissions', 'hit_RT', 'hit_RT_SE', 'variability_of_SE', 'N_perservations',
              'hit_RT_block_change', 'hit_RT_SE_block_change', 'hit_RT_ISI_change', 'hit_RT_SE_ISI_change',
              # WISC
              'Raw.score..DSF', 'Raw.score..DSB', 'Raw.score..SSF', 'Raw.score..SSB',
              # WJ
              'PS',
              # Beery
              'Standard.score')
keep_me = c()
for (v in phen_vars) {
  keep_me = c(keep_me, which(colnames(merged) == v))
}
X = merged[, keep_me]
y = merged$DX_BASELINE
y[y != 'NV'] = 'ADHD'
y = factor(y, levels = c('NV', 'ADHD'))

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
                         summaryFunction = twoClassSummary,
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

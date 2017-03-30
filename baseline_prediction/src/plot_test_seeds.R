load('~/data/baseline_prediction/results/bw/TMP_Xy.RData')
load('~/data/baseline_prediction/results/bw/TMP_1234.RData')
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
preds = predict(m1, newdata=noncorrXtest)
probs = predict(m1, newdata=noncorrXtest, type="prob")
ts = data.frame(obs=ytest, pred=preds, probs)
res = multiClassSummary(ts, lev=levels(ts$obs))
res

# par(mfrow=c(2, 3))
# boxplot(compile_metrics(mylist, 1, run_models), ylab='logLoss', las=2)
# boxplot(compile_metrics(mylist, 2, run_models), ylab='AUC', las=2)
# boxplot(compile_metrics(mylist, 3, run_models), ylab='Accuracy', las=2, ylim=c(.4, 1))
# library(Hmisc)
# res = binconf(x=max(table(groups)), n=length(groups), alpha=.05)
# abline(h=res[1], col='red')
# abline(h=res[2], col='red', lty=2)
# abline(h=res[3], col='red', lty=2)
# boxplot(compile_metrics(mylist, 4, run_models), ylab='Kappa', las=2)
# boxplot(compile_metrics(mylist, 5, run_models), ylab='Sensitivity', las=2)
# boxplot(compile_metrics(mylist, 6, run_models), ylab='Specificity', las=2)
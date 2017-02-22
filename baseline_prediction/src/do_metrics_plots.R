compile_metrics = function(res, col, hdr) {
  b = vector()
  for (d in res) {
    b = cbind(b, d[,col])
  }
  colnames(b) = hdr
  return(b)
}
mylist = list(rndForestAll, lrAll, lsvmAll, rsvmAll, xgbAll, gbmAll)
hdr = run_models
par(mfrow=c(2, 3))
boxplot(compile_metrics(mylist, 1, hdr), ylab='logLoss', las=2)
boxplot(compile_metrics(mylist, 2, hdr), ylab='AUC', las=2)
boxplot(compile_metrics(mylist, 3, hdr), ylab='Accuracy', las=2)
library(Hmisc)
res = binconf(x=max(table(groups)), n=length(groups), alpha=.05)
abline(h=res[1], col='red')
abline(h=res[2], col='red', lty=2)
abline(h=res[3], col='red', lty=2)
boxplot(compile_metrics(mylist, 4, hdr), ylab='Kappa', las=2)
boxplot(compile_metrics(mylist, 5, hdr), ylab='Sensitivity', las=2)
boxplot(compile_metrics(mylist, 6, hdr), ylab='Specificity', las=2)
# loading clinical data
gf = read.csv('~/data/prs/clinical_06192017.csv')
gf = gf[gf$ADHD_current_yes_no!='exclude',]
gf[gf$ADHD_current_yes_no=='no' & is.na(gf$SX_HI),]$SX_HI = 0
gf[gf$ADHD_current_yes_no=='no' & is.na(gf$SX_inatt),]$SX_inatt = 0

# loading PRS data
pgc = read.csv('~/data/prs/PRS2017_original_clump_default.csv')
df = merge(gf, pgc)
# remove duplicated MRNs
df = df[!duplicated(df$MRN),]

# loading brain structural
brain = read.csv('~/data/prs/neuropsych_07282017.csv')
mydata = merge(df, brain, by='MRN')

# choosing mediators
Ms = c(33:ncol(mydata))

out_fname = '~/data/prs/results/model4_p3_neuropsych_hi.csv'
X = mydata$PROFILES.0.3.profile
Y = mydata$SX_HI
nboot = 1000
ncpus = 12


run_model4 = function(X, M, Y, nboot=1000) {
  library(lavaan)
  idx = is.na(X) | is.na(Y) | is.na(M)
  Y = Y[!idx]
  if (!is.factor(Y)) {
    Y = scale(Y)
  }
  run_data = data.frame(X = scale(X[!idx]),
                        Y = Y,
                        M = scale(M[!idx]))
  
  hayes4 <- ' # direct effect
  Y ~ c*X
  direct := c
  
  # regressions
  M ~ a*X
  Y ~ b*M
  
  # indirect effect (a*b)
  indirect := a*b
  
  # total effect
  total := c + (a*b)'
  
  # fit model
  sem4 <- sem(model = hayes4,
              data = run_data,
              se = "bootstrap",
              bootstrap = nboot)
  
  res = parameterEstimates(sem4,
                           boot.ci.type = "bca.simple",
                           level = .95, ci = TRUE,
                           standardized = FALSE)
  res2 = res[7:9, 5:10]
  rownames(res2) = res[7:9,]$label
  # adding n
  res2 = cbind(rep(nrow(run_data), nrow(res2)), res2)
  colnames(res2)[1] = 'n'
  return(res2)
}

run_wrapper = function(m, run_model, mydata, nboot, X, Y) {
  cat('\t', sprintf('M=%s', colnames(mydata)[m]), '\n')
  tmp = run_model(X, mydata[, m], Y, nboot=nboot)
  res = c()
  for (i in 1:nrow(tmp)) {
    for (j in 1:ncol(tmp)) {
      res = c(res, tmp[i, j])
      names(res)[length(res)] = sprintf('%s_%s', rownames(tmp)[i], colnames(tmp)[j])
    }
  }
  return(res)
}

# no need to change anything below here. The functions remove NAs and zscore variables on their own
library(parallel)
cl <- makeCluster(ncpus)

m1_res = parLapply(cl, Ms, run_wrapper, run_model4, mydata, nboot, X, Y)
# m1_res = lapply(Ms, run_wrapper, run_model4, mydata, nboot, X, Y)
all_res = do.call(rbind, m1_res)
rownames(all_res) = colnames(mydata)[Ms]
  
stopCluster(cl)
write.csv(all_res, file=out_fname)

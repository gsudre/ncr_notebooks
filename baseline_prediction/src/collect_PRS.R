# Collects the polygenic risk scores from PROFILE files, created by PRSice
res_dir = '~/data/baseline_prediction/prs/results/'

# read in result files
files = dir(path = res_dir, pattern = 'profile$')
prs = c()
for (f in files) {
  a = read.table(sprintf('%s/%s', res_dir, f), header=1)
  prs = cbind(prs, a$SCORE)
}
colnames(prs) = files

# cleaning weird looking NSB strings
clean_nsb = function(nsb) {
  if (grepl('@', nsb)) {
    nsb = strsplit(nsb, '@')[[1]][1]
    nsb = strsplit(nsb, '-')[[1]][3]
  }
  return(nsb)
}
nsbs = as.numeric(sapply(as.character(a$IID), clean_nsb))

# transform to a nice looking dataframe
prs_df = as.data.frame(cbind(nsbs, prs))
colnames(prs_df)[1] = 'NSB_GENOTYPE_INDEX'

# some NSBs have been genotyped twice. I checked a few of them and their PRS correlation is 
# > .99 between genotype waves, so let's just pick the first one
dups = duplicated(prs_df$NSB_GENOTYPE_INDEX)
prs_df = prs_df[!dups, ]

# merge NSB to MRNs to clip to only baseline prediction project
mrns = read.csv('~/data/baseline_prediction/stripped/genotype_wave.csv')
mrns = mrns[!is.na(mrns$NSB_GENOTYPE_INDEX), ]

m = merge(mrns, prs_df)

write.csv(m, file='~/data/baseline_prediction/stripped/PRS.csv', row.names=F)
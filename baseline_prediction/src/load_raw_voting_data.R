source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf_base = gf[gf$BASELINE=='BASELINE' & gf$age <= 12, ]
my_ids = gf_base$MRN

print('Loading DTI tracts')
tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
rm_me = (tract_data$fa_avg < .4 | tract_data$ad_avg < 1.18 | tract_data$rd_avg > .65 | tract_data$rd_avg < .5 |
           tract_data$norm.trans > .45 | tract_data$norm.rot > .008 | tract_data$goodSlices < 45 | tract_data$goodSlices > 70)
tract_data = tract_data[!rm_me, ]
merged = mergeOnClosestDate(gf_base, tract_data, my_ids)
phen_vars = c(which(grepl("^FA_", colnames(merged))),
              which(grepl("^AD_", colnames(merged))),
              which(grepl("^RD_", colnames(merged))),
              which(grepl("^MO_", colnames(merged))),
              which(grepl("^age$", colnames(merged))),
              which(grepl("^SEX", colnames(merged)))
)
rm_me = abs(merged$dateX.minus.dateY.months) > 12
X = merged[, phen_vars]
X[which(rm_me), ] = NA
dti_tracts = as.data.frame(X)
if (sum(gf_base$MRN != merged$MRN) != 0) {
  print('ERROR merging!')
}

print('Loading structural ROIs')
struct_data = read.csv('~/data/baseline_prediction/stripped/structural.csv')
rm_me = (struct_data$mprage_score > 2)
struct_data = struct_data[!rm_me, ]
merged = mergeOnClosestDate(gf_base, struct_data, my_ids)
phen_vars = c(which(grepl("^lh_", colnames(merged))),
              which(grepl("^rh_", colnames(merged))),
              which(grepl("^age$", colnames(merged))),
              which(grepl("^SEX", colnames(merged)))
)
X = merged[, phen_vars]
rm_me = abs(merged$dateX.minus.dateY.months) > 12
X[which(rm_me), ] = NA
struct_rois = as.data.frame(X)
if (sum(gf_base$MRN != merged$MRN) != 0) {
  print('ERROR merging!')
}

print('Loading neuropsych')
beery_data = read.csv('~/data/baseline_prediction/stripped/beeryVMI.csv')
mbeery = mergeOnClosestDate(gf_base, beery_data, my_ids, y.date='record.date.collected', y.id='Medical.Record...MRN')
rm_me = abs(mbeery$dateX.minus.dateY.months) > 12
mbeery[which(rm_me),] = NA
X = mbeery$Standard.score
cpt_data = read.csv('~/data/baseline_prediction/stripped/cpt.csv')
mcpt = mergeOnClosestDate(gf_base, cpt_data, my_ids)
rm_me = abs(mcpt$dateX.minus.dateY.months) > 12
mcpt[which(rm_me),] = NA
phen_vars = c('N_of_omissions', 'N_commissions', 'hit_RT', 'hit_RT_SE', 'variability_of_SE', 'N_perservations',
              'hit_RT_block_change', 'hit_RT_SE_block_change', 'hit_RT_ISI_change', 'hit_RT_SE_ISI_change')
keep_me = sapply(phen_vars, function(d) which(colnames(mcpt) == d))
X = cbind(X, mcpt[, keep_me])
colnames(X)[1] = 'BeeryVM.StdScore'
iq_data = read.csv('~/data/baseline_prediction/stripped/iq.csv')
miq = mergeOnClosestDate(gf_base, iq_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(miq$dateX.minus.dateY.months) > 12
miq[which(rm_me),] = NA
X = cbind(X, miq$FSIQ)
colnames(X)[ncol(X)] = 'FSIQ'
wisc_data = read.csv('~/data/baseline_prediction/stripped/wisc.csv')
mwisc = mergeOnClosestDate(gf_base, wisc_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(mwisc$dateX.minus.dateY.months) > 12
mwisc[which(rm_me),] = NA
phen_vars = c('Raw.score..DSF', 'Raw.score..DSB', 'Raw.score..SSF', 'Raw.score..SSB')
keep_me = sapply(phen_vars, function(d) which(colnames(mwisc) == d))
X = cbind(X, mwisc[, keep_me])
wj_data = read.csv('~/data/baseline_prediction/stripped/wj.csv')
mwj = mergeOnClosestDate(gf_base, wj_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(mwj$dateX.minus.dateY.months) > 12
mwj[which(rm_me),] = NA
phen_vars = c('Raw.Score..VM', 'Raw.Score..DS', 'PS')
keep_me = sapply(phen_vars, function(d) which(colnames(mwj) == d))
X = cbind(X, mwj[, keep_me])
X = cbind(X, mwj[, c('age', 'SEX')])
neuropsych = as.data.frame(X)
if (sum(gf_base$MRN != merged$MRN, na.rm=T) != 0) {
  print('ERROR merging!')
}

print('Loading GeoSpatial')
geo_data = read.csv('~/data/baseline_prediction/stripped/geospatial.csv')
merged = merge(gf_base, geo_data, by='MRN', sort=F)
# some variables are not being read as numeric...
merged$Home_Price = as.numeric(merged$Home_Price)
merged$Fam_Income = as.numeric(merged$Fam_Income)
merged$Crime_Rate = as.numeric(merged$Crime_Rate)
phen_vars = c('SES', 'Home_Type', 'Home_Price', 'Fam_Income', 'Pop_BPL', 'Fam_BPL', 'Pub_School',
              'Crime_Rate', 'Green_Space', 'Park_Access', 'Air_Quality', 'Obesity_Rate',
              'Food_Index', 'Exercise_Access', 'Excessive_Drinking',
              'poverty_education', 'social_environment', 'physical_health_environment',
              'physical_envronment', 'age', 'SEX')
keep_me = sapply(phen_vars, function(d) which(colnames(merged) == d))
X = merged[, keep_me]
geospatial = as.data.frame(X)
if (sum(gf_base$MRN != merged$MRN) != 0) {
  print('ERROR merging!')
}

print('Loading PRS')
prs_data = read.csv('~/data/baseline_prediction/stripped/PRS_original_clump_default.csv')
# we don't need the extra BASELINE or GenotypeWave columns
prs_data = prs_data[, -c(3, 4)]
# remove people with more than one genotype
prs_data = prs_data[!duplicated(prs_data$MRN), ]
merged = merge(gf_base, prs_data, by='MRN', sort=F, all.x=T)
# somehow all.x=T is canceling sort=F, so we need to resort things here
X = c()
mycols = c(which(grepl("^SEX", colnames(merged))), which(grepl("age", colnames(merged))),
           which(grepl("^PROFILES", colnames(merged))))
for (s in gf_base$MRN) {
  X = rbind(X, merged[merged$MRN==s, mycols])
}
prs = as.data.frame(X)
if (sum(gf_base$age != X$age) != 0) {
  print('ERROR merging!')
}

##################
# # print('Loading DTI voxels')
# # tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
# # load('~/data/baseline_prediction/dti/ad_voxelwise.RData')
# # dti_vdata = cbind(tract_data$maskid, ad_data)
# # merged = mergeOnClosestDate(gf_base, tract_data, my_ids)
# # rm_me = abs(merged$dateX.minus.dateY.months) > 12
# # dti_base_vdata = merge(merged$maskid, dti_vdata, by.x=1, by.y=1, all.y=F, all.x=T)
# # X = dti_base_vdata[, 2:ncol(dti_base_vdata)]
# # X[which(rm_me), ] = NA
# # library(parallel)
# # cl <- makeCluster(ncpus)
# # X_resid = parSapply(cl, X, get_needed_residuals, 'y ~ df$age + I(df$age^2) + df$SEX', .1, merged)
# # stopCluster(cl)
# # brain_ad = as.data.frame(X_resid)
# # 
# # tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
# # load('~/data/baseline_prediction/dti/fa_voxelwise.RData')
# # dti_vdata = cbind(tract_data$maskid, fa_data)
# # merged = mergeOnClosestDate(gf_base, tract_data, my_ids)
# # rm_me = abs(merged$dateX.minus.dateY.months) > 12
# # dti_base_vdata = merge(merged$maskid, dti_vdata, by.x=1, by.y=1, all.y=F, all.x=T)
# # X = dti_base_vdata[, 2:ncol(dti_base_vdata)]
# # X[which(rm_me), ] = NA
# # library(parallel)
# # cl <- makeCluster(ncpus)
# # X_resid = parSapply(cl, X, get_needed_residuals, 'y ~ df$age + I(df$age^2) + df$SEX', .1, merged)
# # stopCluster(cl)
# # brain_fa = as.data.frame(X_resid)
# # 
# # tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
# # load('~/data/baseline_prediction/dti/rd_voxelwise.RData')
# # dti_vdata = cbind(tract_data$maskid, rd_data)
# # merged = mergeOnClosestDate(gf_base, tract_data, my_ids)
# # rm_me = abs(merged$dateX.minus.dateY.months) > 12
# # dti_base_vdata = merge(merged$maskid, dti_vdata, by.x=1, by.y=1, all.y=F, all.x=T)
# # X = dti_base_vdata[, 2:ncol(dti_base_vdata)]
# # X[which(rm_me), ] = NA
# # library(parallel)
# # cl <- makeCluster(ncpus)
# # X_resid = parSapply(cl, X, get_needed_residuals, 'y ~ df$age + I(df$age^2) + df$SEX', .1, merged)
# # stopCluster(cl)
# # brain_rd = as.data.frame(X_resid)
# 
# # print('Loading structural voxels')

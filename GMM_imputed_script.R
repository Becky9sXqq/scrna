## First step:
##after loading GMM model's all functions
##fitting for single cell matrix: row is genes column is cell
source('https://github.com/Becky9sXqq/scrna/raw/master/odgmm.R')
source('https://github.com/Becky9sXqq/scrna/raw/master/gmm_impute_core_functions')

GMM_impute_wrapper <- function(sce, seed = 42, comp = 2,ncores = 10 , rate = 0.8) {
  
library(parallel)
mu0 = 0
sigma0 = 0.1
set.seed(seed)
X_tab = as.matrix(assays(sce)[["counts"]])
comp = comp
  
## main step for fitting
mclapply(1:dim(X_tab)[1],function(i){
  X = X_tab[i,]#X_tab is the matrix need to be fitted by GMM
  comp = comp
  res = tryCatch(odgmm(X,comp), error=function(e) NA)
  alpha = tryCatch(res[[1]], error=function(e) NA)
  beta = tryCatch(res[[2]], error=function(e) NA)
  ws = tryCatch(res[[3]], error=function(e) NA)
  logml = tryCatch(res[[4]], error=function(e) NA)
  Z = tryCatch(res[[5]], error=function(e) NA)
  final_bic = tryCatch(res[[6]], error=function(e) NA)
  emus = tryCatch(alpha2mu(alpha), error=function(e) NA)
  esgs = tryCatch(beta2sigma(beta), error=function(e) NA)
  ews = tryCatch(ws, error=function(e) NA)
  x = tryCatch(seq(-2,max(X)+1,0.01), error=function(e) NA)
  ey = tryCatch(gmmpdf(x, emus, esgs, ews) + tail(ews,n=1)/max(X), error=function(e) NA)
  K = length(alpha)
  print(paste0('--the--',i,'--line--'))
  return(list(res = res , x = x , ey = ey , K = K ,
              emus = emus , ews = ews,
              esgs = esgs , components = length(ws)))
},mc.cores = ncores) -> GMM_res

##if GMM output contain 0 or NA, we need to filter them out 
X_tab = X_tab[which(unlist(lapply(GMM_res,function(x) length(x)))!=0),]
GMM_res[which(unlist(lapply(GMM_res,function(x) length(x)))!=0)] -> GMM_res

lapply(GMM_res,function(x){
  !unique(unlist(x$res)) %in% NA
}) -> GMM_res_succeed

unlist(lapply(GMM_res_succeed,function(x) unlist(unique(x)))) -> GMM_res_succeed

GMM_res_withoutNA = GMM_res[which(GMM_res_succeed == TRUE)]
X_tab = X_tab[which(GMM_res_succeed == TRUE),]

## Second step:
## GMM output into imputing

library(e1071)
n_gene = dim(X_tab)[1]
n_cell = dim(X_tab)[2]

X = matrix(0, nrow = n_cell, ncol = n_gene)
xvals <- vector("list", length = n_gene)

proc_res = function(res){
  Z = res[[5]]
  mus = res[[1]]  
  L = dim(Z)[2]
  y = apply(Z,1,which.max)
  y = y-1
  new_mus = c(mus[-1],mean(mus[-1]))
  return(list(label=y, mus=new_mus))
}      
              
for(i in seq(n_gene)){
  cat("proc gene",i,"\n")
  tmpres = proc_res(GMM_res_withoutNA[[i]]$res)
  X[,i] = tmpres[[1]]
  xvals[[i]] = tmpres[[2]]
}

# filter genes

keep_rate = apply(X>0,2,sum)/n_cell
gene_inds = which(keep_rate > rate)

# perform imputation

XX = X[,gene_inds]
XX_val = xvals[gene_inds]
res_after_impute = gmm_impute(XX,XX_vals)
X_impute = res_after_impute$impute_label_mat
X_val_impute = res_after_impute$impute_val_mat # NA means values that does not need to be imputed
              
orig_data = t(X_tab[gene_inds,])
valid_inds = is.na(X_val_impute)
X_val_impute[valid_inds] = orig_data[valid_inds]
GMM_imputed_X_tab = X_tab
t(X_val_impute) ->  GMM_imputed_X_tab[gene_inds,]         
return(GMM_imputed_X_tab)
}


res <- GMM_impute_wrapper(sce,42,2,10,0.8) 
              
assays(sce)[["GMM_impute"]] <- res
              
#finally all the results are saved in sce object    

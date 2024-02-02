library(magrittr)
library(purrr)
library(abind)
load_as_list <- function(fname){
  e <- new.env()
  load(fname, e)
  as.list(e)
}

for(prefix in c("flex_CIs","S_CIs","para_CIs","sflex_CIs","sS_CIs")){
filenames <- 1:100 %>% sprintf(paste(prefix,"_%03d.rda",sep=''), .)
comb_list <- lapply(filenames, load_as_list)
varnames <- names(comb_list[[1]])
comb_env <- new.env()
for(var in varnames){
  assign(var, 
  reduce(lapply(1:100, function(x) get(var,comb_list[[x]])), abind, along=1),
  comb_env
  )
}
save(list=ls(comb_env), file=paste(prefix,".rda",sep=''),envir=comb_env)
}
rm(list=ls());gc()
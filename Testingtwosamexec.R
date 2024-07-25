#.libPaths("/moto/stats/users/nd2560/rpackages") --- Location of packages
#source("/moto/home/nd2560/Testingtwosamprep.R") --- Location of source files
statdist=mainf() #--- Function to invoke  
for (e in commandArgs(TRUE)) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    temp = as.numeric(temp)
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
    print(statdist)
    print(temp)
  }
}
# Code for Table 1.
# clusterCall(cl, function(pkgs) {
#   for (req in pkgs) {
#     require(req, character.only=TRUE)
#   }
# }, c("Rcplex","Matrix"))
#rm(list=setdiff(ls(), "seed"))

#library(cobs)
# library(parallel)
# library(foreach)
# library(doParallel)
#library(np)

#print("#----------------------------------------------------------------------")
#cat("seed ", seed, "\n")
#print("#----------------------------------------------------------------------")

file.name.basic <- paste0("data", "seed_", seed,"methods.RData")
save(statdist, file =file.name.basic )
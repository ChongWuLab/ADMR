slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

# I-squared function
Isq = function(y,s){
    k          = length(y)
    w          = 1/s^2; sum.w  = sum(w)
    mu.hat     = sum(y*w)/sum.w
    Q          = sum(w*(y-mu.hat)^2)
    Isq        = (Q - (k-1))/Q
    Isq        = max(0,Isq)
    return(Isq)
}

# pr
library(devtools)
library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/3.6", install_github("MRCIEU/TwoSampleMR"))

library(TwoSampleMR)
library(data.table)

library(BEDMatrix)
job.id = job.id - 1
#job.id = 16

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
job <- as.numeric(job)

cat(job.id)

exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure_noclump/"

out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/rev_IGAP_proxy/"
cutoff = 5e-8

tmp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/tmp/"

#exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure/"


files = c("ieu-a-1239","ukb-b-11615","ukb-b-17557","ukb-b-17409","ukb-b-2002","ebi-a-GCST006250","ieu-a-1013","ukb-b-17729","ukb-d-II_NEOPLASM","ukb-b-16446","ukb-d-30130_irnt","ebi-a-GCST005367","ieu-a-89","ieu-a-1001","ukb-b-4461","ukb-b-16099","ukb-b-19925","ukb-b-10787","ukb-b-4710","ukb-b-8764","ukb-b-9685","ukb-b-17271")


#wgt_list[i] = "ubm-a-125.RData"
load(paste(out.dir,files[job.id],".RData",sep=""))
        
cat("Harmizoning data\n")
        

presso.res = NULL
tryCatch({
        presso.res = run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
        presso.res
}, error = function(e) {
           cat("ERROR msg :", conditionMessage(e), "\n")
})
       
out.name = paste(out.dir,files[job.id],".RData",sep="")
out.name2 = paste(out.dir,files[job.id],".RData",sep="")

res_loo <- mr_leaveoneout(dat)

save(mF,Isq,ind.pre.SNP,proxy.inf,instruments,exp_dat,out_dat,dat,mr.res,ivw.res,egger.res,presso.res,exp.name,res_loo,file = out.name)
save(mF,Isq,ind.pre.SNP,proxy.inf,instruments,exp_dat,out_dat,dat,mr.res,ivw.res,egger.res,presso.res,exp.name,res_loo,file = out.name2)

        #
cat("Finish Analysis\n")
        


files = c("ieu-a-1239","ukb-b-11615","ukb-b-17557","ukb-b-17409","ukb-b-2002","ebi-a-GCST006250","ieu-a-1013","ukb-b-17729","ukb-d-II_NEOPLASM","ukb-b-16446","ukb-d-30130_irnt","ebi-a-GCST005367","ieu-a-89","ieu-a-1001","ukb-b-4461","ukb-b-16099","ukb-b-19925","ukb-b-10787","ukb-b-4710","ukb-b-8764","ukb-b-9685","ukb-b-17271")
for(id in files) {
    load(paste(id,".RData",sep=""))
    cat(id,"  ",presso.res[[1]][2]$`MR-PRESSO results`$`Global Test`$Pvalue,"\n")
}

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

# pr
library(devtools)
library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/3.6", install_github("MRCIEU/TwoSampleMR"))

job.id = 2
job.id = job.id - 1


library(TwoSampleMR)

ao <- available_outcomes()

ao = as.data.frame(ao)
ao = ao[ao[,"access"] == "public",]
ao = ao[ao[,"population"]=="European",]

ao = ao[grepl("ebi-a",ao[,1]) |grepl("ieu-a",ao[,1]) | grepl("ubm-a",ao[,1])| grepl("ukb-b",ao[,1])| grepl("ukb-d",ao[,1]),]
exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure_noclump/"
files = list.files(path = exp.dir)
files = gsub(".RData","",files)
length(files)

ao = ao[!ao[,1] %in% files,]


wgt_list = ao
#bmi_gwas <- subset(gwas_catalog, Phenotype == "Body mass index")
for(i in 1:dim(wgt_list)[1]) {#dim(ao)[1]
    #gene.indx = 3
    tryCatch({
        exp_dat <- extract_instruments(outcomes = wgt_list[i,1],p1 = 5e-06,clump = FALSE,p2 = 5e-06,r2 = 0.01,kb = 1000)
        
        save(exp_dat,file = paste("/gpfs/research/chongwu/Chong/Application/MR_AD/exposure_noclump/",wgt_list[i,1],".RData",sep=""))
            cat(i,"\n")
        
        
    }, error = function(e) {
        cat("Error ",i,"\n")
        cat("ERROR :", conditionMessage(e), "\n")
    })
}


exp_dat <- extract_instruments(outcomes = "ukb-b-17988",p1 = 1,clump = FALSE,p2 = 5e-06,r2 = 0.01,kb = 1000)

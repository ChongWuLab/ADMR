slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
chr.id <- as.numeric(slurm_arrayid)

# pr
library(devtools)
library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/3.6", install_github("MRCIEU/TwoSampleMR"))

library(TwoSampleMR)
library(data.table)

dat = fread("/gpfs/research/chongwu/shared/summary_statistics/AD/IGAP_stage_1.txt")

dat = as.data.frame(dat)


#chr.id = 22
freq = fread(paste("/gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.freq.CHR",chr.id,".afreq",sep=""))
freq = as.data.frame(freq)

dat = dat[dat[,"Chromosome"]==chr.id,]
dat = dat[!duplicated(dat[,"MarkerName"]),]

freq = freq[!duplicated(freq[,2]),]
dat$eaf = NA
sum(dat[,"MarkerName"] %in% freq[,"ID"])

rownames(dat) = dat$MarkerName
rownames(freq) = freq$ID


freq2 = freq[dat$MarkerName,]

indx = dat[,"Effect_allele"] ==freq2[,"ALT"] & dat[,"Non_Effect_allele"] ==freq2[,"REF"]
indx[is.na(indx)]=FALSE
dat[indx,"eaf"] = freq2[indx,5]
      
indx = dat[,"Effect_allele"] ==freq2[,"REF"] & dat[,"Non_Effect_allele"] ==freq2[,"ALT"]
indx[is.na(indx)]=FALSE
dat[indx,"eaf"]= 1 - freq2[indx,5]
   
saveRDS(dat,paste("/gpfs/research/chongwu/Chong/Application/MR_AD/IGAP/IGAP_GWAS_with_freq_chr",chr.id,".rds",sep=""))

IGAP = NULL
for(chr.id in 1:22) {
    tmp = readRDS(paste("/gpfs/research/chongwu/Chong/Application/MR_AD/IGAP/IGAP_GWAS_with_freq_chr",chr.id,".rds",sep=""))
    IGAP = rbind(IGAP,tmp)
}

saveRDS(IGAP,"/gpfs/research/chongwu/Chong/Application/MR_AD/IGAP/IGAP_GWAS_with_freq.rds")

################################
### For using AD as exposure ###
################################

exp_dat2 <- readRDS("/gpfs/research/chongwu/Chong/Application/MR_AD/IGAP_MR.rds")

freq2 = freq[freq[,"ID"] %in% exp_dat2[,"SNP"],]

rownames(exp_dat2) = exp_dat2[,"SNP"]
rownames(freq2) = freq2[,"ID"]
freq2 = freq2[rownames(exp_dat2),]

indx = exp_dat2[,"effect_allele.exposure"]==freq2[,"ALT"]
indx[is.na(indx)]=FALSE
exp_dat2[indx,"eaf.exposure"] = freq2[indx,5]
      
indx = exp_dat2[,"effect_allele.exposure"]==freq2[,"REF"]
indx[is.na(indx)]=FALSE

exp_dat2[indx,"eaf.exposure"] = 1 - freq2[indx,5]
   
#saveRDS(exp_dat2,"/gpfs/research/chongwu/Chong/Application/MR_AD/IGAP_MR_with_freq.rds")


#
#export PATH=$PATH:/gpfs/research/chongwu/shared/software

#for i in {1..22}; do
#plink2 --bfile 1000G.EUR.ALLSNP.QC.CHR$i --freq --out 1000G.EUR.ALLSNP.QC.freq.CHR$i
#done


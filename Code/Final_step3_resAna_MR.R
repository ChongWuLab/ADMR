slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

# pr
library(devtools)
library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/3.6", install_github("MRCIEU/TwoSampleMR"))

library(TwoSampleMR)
library(data.table)

job.id = job.id - 1
#job.id = 7

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
job <- as.numeric(job)

cat(job.id)

exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure_noclump/"

#exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure_r2_01/"

#out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/res_nopalindromic/" #IGAP

#out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/res_IGAP_nooutlier/"
out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/rev_IGAP_proxy/"
#out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/res_IGAP/"

tmp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/tmp/"


wgt_list = list.files(path = out.dir)
length(wgt_list)

final.out = as.data.frame(matrix(NA,length(wgt_list),23))

colnames(final.out) = c("exposure","nsnp","IVW_b","IVW_se","IVW_pval","Wght_med_b","Wght_med_se","Wght_med_pval","Wght_mod_b","Wght_mod_se","Wght_mod_pval","egger_intercept","egger_intercept_se","egger_intercept_pval","PRESSO-gloabl_p","PRESSO-beta","PRESSO-sd","PRESSO-corrected_p","F statistics", "No. of Significant SNPs","No. SNPs included in the analysis","Proxies used","Excluded")


wgt_list = list.files(out.dir)
for(i in 1:length(wgt_list)) {
    tryCatch({
        #wgt_list[i]="ukb-b-12209.RData"
        #load(paste(exp.dir,wgt_list[i],sep=""))

        #wgt_list[i]="ieu-a-965.RData"
        load(paste(out.dir,wgt_list[i],sep=""))
        dat
        if(is.null(exp.name)) {
            next
        }
        final.out[i,19] = mF
        
        final.out[i,20] = length(ind.pre.SNP)
        final.out[i,21] = length(instruments)
        final.out[i,22] = length(ind.pre.SNP) - sum(instruments %in% ind.pre.SNP)
        final.out[i,23] = length(ind.pre.SNP) - length(instruments)

        
        if(!is.null(presso.res)) {
            final.out[i,15] = presso.res[[1]][[2]][[1]][[2]] # global P value
            final.out[i,16:18] = presso.res[[1]][[1]][2,c(3,4,6)]
            
        }
        
        final.out[i,1] = gsub(".RData","",wgt_list[i])
        final.out[i,2:5] = mr.res[mr.res$method=="Inverse variance weighted",c("nsnp","b","se","pval")]
        
        final.out[i,6:8] = mr.res[mr.res$method=="Weighted median",c("b","se","pval")]
        final.out[i,9:11] = mr.res[mr.res$method=="Weighted mode",c("b","se","pval")]
        
        final.out[i,12:14] = c(egger.res$b_i,egger.res$se_i,egger.res$pval_i)
        
    }, error = function(e) {
        cat("Error: ",i,"\n")
        cat("ERROR :", conditionMessage(e), "\n")
    })
}

final.out = final.out[final.out$nsnp>2,]

final.out = final.out[!is.na(final.out[,2]),]
dim(final.out)
rownames(final.out) = final.out[,1]
ao <- available_outcomes()
ao <- as.data.frame(ao)

ao = ao[ao[,1] %in% final.out[,1],]
rownames(ao) = ao[,1]
ao = ao[rownames(final.out),]

final.out$trait = ao[,2]
#final.out$id = ao[,1]

#final.out = final.out[!grepl("Alzheimer",final.out$trait),]

#final.out = final.out[final.out[,1]!="ukb-b-17557",]

final.out = final.out[!grepl("ubm-a",final.out[,1]),]

saveRDS(final.out,"/gpfs/research/chongwu/Chong/Application/MR_AD/res_IGAP1_summary.rds")
#files = ao1[,1]
#files = files[grepl("ebi-a",files) |grepl("ieu-a",files) |  grepl("ukb-b",files)| grepl("ukb-d",files)]
ao1 = ao[ao[,1] %in% final.out[,1],]

tmp2 = ao1[grepl("medication",ao1[,2]),]
tmp3 = ao1[grepl("ubm-a",ao1[,1]),]
tmp1 = ao1[!ao1[,1] %in% c(tmp2[,1],tmp3[,1]),]

exp.inf = rbind(tmp1,tmp2)

exp.inf = exp.inf[,c("id","trait","sample_size","nsnp","pmid","consortium","note")]

tmp = final.out[,19:23]
tmp = tmp[rownames(exp.inf),]
exp.inf = cbind(exp.inf,tmp)

#exp.inf = rbind(exp.inf,tmp3)
write.csv(exp.inf,"exposure_inf.csv") #supl Tables 1 & 2


# Supl Table 3
final.out = final.out[final.out[,1] %in% exp.inf[,1],]


out = final.out[,c(1,24,19,2:8,14)]

write.csv(out,"suplTable34.csv")

tmp = final.out[order(final.out[,"IVW_pval"]),]
out1 = out[!grepl("medication",out$trait),]
write.csv(out1,"suplTale3.csv")

out1 = out[grepl("medication",out$trait),]
write.csv(out1,"suplTale4.csv")
#False discovery rate 0.05


library(sgof)
BH(final.out[,"IVW_pval"], alpha = 0.05)


#final.out[final.out[,1]=="ukb-b-12209",]

indx = final.out[,"IVW_pval"] < 0.05/100#dim(final.out)[1]
tmp.p = sort(final.out[,"IVW_pval"])
BH(tmp.p[5:length(tmp.p)], alpha = 0.05)

res = final.out[order(final.out[,"IVW_pval"]),]
rownames(res) = res[,1]
res = res[res[,"IVW_pval"] < 0.05/100,]

#
tmp = res[res[,"Wght_med_pval"]<0.05,]
#rownames(tmp) = tmp$trait
tmp$IVW_OR = round(exp(tmp$IVW_b),2)

tmp$IVW_con = paste(round(exp(tmp$IVW_b-1.96 * tmp$IVW_se),2),"-",round(exp(tmp$IVW_b+1.96 * tmp$IVW_se),2),sep="")

tmp$med_OR = round(exp(tmp$Wght_med_b),2)
tmp$med_con = paste(round(exp(tmp$Wght_med_b-1.96 * tmp$Wght_med_se),2),"-",round(exp(tmp$Wght_med_b+1.96 * tmp$Wght_med_se),2),sep="")

tmp = tmp[,c("trait","F statistics","IVW_OR","IVW_con","IVW_pval","med_OR","med_con","Wght_med_pval","egger_intercept_pval","PRESSO-gloabl_p")]

write.csv(tmp,"table1.csv") #Table 1



################################
indx = final.out[,"IVW_pval"] < 0.005#dim(final.out)[1]
tmp.p = sort(final.out[,"IVW_pval"])
BH(tmp.p[5:length(tmp.p)], alpha = 0.05)

res = final.out[indx,]
res = res[order(res[,"IVW_pval"]),]
rownames(res) = res[,1]
#res = res[8:16,]

#
tmp = res[res[,"Wght_med_pval"]<0.05,] #table 2
tmp = tmp[9:dim(tmp)[1],]


tmp$IVW_OR = round(exp(tmp$IVW_b),2)
tmp$IVW_con = paste(round(exp(tmp$IVW_b-1.96 * tmp$IVW_se),2),"-",round(exp(tmp$IVW_b+1.96 * tmp$IVW_se),2),sep="")

tmp$med_OR = round(exp(tmp$Wght_med_b),2)
tmp$med_con = paste(round(exp(tmp$Wght_med_b-1.96 * tmp$Wght_med_se),2),"-",round(exp(tmp$Wght_med_b+1.96 * tmp$Wght_med_se),2),sep="")

tmp = tmp[,c("trait","F statistics","IVW_OR","IVW_con","IVW_pval","med_OR","med_con","Wght_med_pval","egger_intercept_pval","PRESSO-gloabl_p")]

write.csv(tmp,"table2.csv") #Table 2





tmp = final.out[grepl("medication",final.out[,"trait"]),]
#rownames(tmp) = tmp$trait
tmp$OR = exp(tmp$IVW_b)
tmp$lower = exp(tmp$IVW_b-1.96 * tmp$IVW_se)
tmp$upper = exp(tmp$IVW_b+1.96 * tmp$IVW_se)

tmp["ukb-b-15445",]


res2 = res[res[,"Wght_med_pval"]<0.05,]
write.table(res2,"/gpfs/research/chongwu/Chong/Application/MR_AD/significant_res.txt")
res[c("ieu-a-1239","ukb-b-16489","ebi-a-GCST006250","ukb-b-11615"),]

res[c("ebi-a-GCST005920","ebi-a-GCST005921","ebi-a-GCST005923"),]

system("export PATH=$PATH:/gpfs/research/chongwu/shared/software")



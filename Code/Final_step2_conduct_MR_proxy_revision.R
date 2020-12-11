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

APOEexclude <- as.character(args[[2]])

#APOEexclude = "exclude"
#job= 2
#job.id = 1
#cat(job.id)

exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure_noclump/"

if(APOEexclude == "exclude") {
    out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/rev_IGAP_proxy_no_APOE/"
    dir.create(out.dir)
    tmp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/tmp/"
    dir.create(tmp.dir)
    
} else {
    out.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/rev_IGAP_proxy/"
    dir.create(out.dir)
    tmp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/tmp2/"
    dir.create(tmp.dir)
}


cutoff = 5e-8
#exp.dir = "/gpfs/research/chongwu/Chong/Application/MR_AD/exposure/"


files = list.files(path = exp.dir)

files = files[grepl("ebi-a",files) |grepl("ieu-a",files) | grepl("ubm-a",files)| grepl("ukb-b",files)| grepl("ukb-d",files)]

finished = list.files(path = out.dir)
files = files[!files%in%finished ]

each = ceiling(length(files) / 290)
last.job.num = floor(length(files) / each)

if (job.id == last.job.num) wgt_list = files[(job.id * each + 1):length(files)] else wgt_list = files[(job.id * each + 1):(job.id * each + each)]


final.out = as.data.frame(matrix(NA,length(wgt_list),18))

colnames(final.out) = c("exposure","nsnp","IVW_b","IVW_se","IVW_pval","Wght_med_b","Wght_med_se","Wght_med_pval","Wght_mod_b","Wght_mod_se","Wght_mod_pval","egger_intercept","egger_intercept_se","egger_intercept_pval","PRESSO-gloabl_p","PRESSO-beta","PRESSO-sd","PRESSO-corrected_p")

# add information for
IGAP = readRDS("/gpfs/research/chongwu/Chong/Application/MR_AD/IGAP/IGAP_GWAS_with_freq.rds")
IGAP = IGAP[!is.na(IGAP$eaf),]

IGAP$Pvalue = as.numeric(IGAP$Pvalue)
# remove genome-wide significant SNPs
IGAP = IGAP[IGAP$Pvalue>5e-8,]

# remove genome-wide significant SNPs and APOE region

# Excluding APOE gene region (44,409,039–46,412,650 bp according to GRch37/Feb 2009).; https://doi.org/10.1038/s41380-018-0030-8
if(APOEexclude == "exclude") {
    
    indx = IGAP[,1]==19
    indx2 = IGAP[,2]<46412650 & IGAP[,2]>44409039
    indx = indx & indx2
    indx = !indx
    IGAP = IGAP[indx,]
}


for(i in 1:length(wgt_list)) {
    tryCatch({
        
        #wgt_list[i] = "ubm-a-125.RData"
        load(paste(exp.dir,wgt_list[i],sep=""))
        
        
        exp_dat0 = exp_dat
        exp_dat = exp_dat[exp_dat[,"pval.exposure"]<cutoff ,]
        
 
        if(is.null(exp_dat)) {
            
            out.name = paste(out.dir,wgt_list[i],sep="")
            exp.name = NULL
            save(exp.name,file = out.name)
            
            next
        }
        
        
        if(dim(exp_dat)[1]<3) {
            out.name = paste(out.dir,wgt_list[i],sep="")
            exp.name = NULL
            save(exp.name,file = out.name)
            
            next
        }
        
        dim(exp_dat)
        #############################################
        # Prepare SNP-exposure and SNP-outcome data
        #############################################
        chr.list = unique(as.numeric(exp_dat$chr.exposure))
        chr.list = sort(chr.list)
        chr.list = chr.list[!is.na(chr.list)]
        snp_id = NULL
        for (chr.id in chr.list) {
            tmp = fread(paste("/gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.CHR", chr.id, ".bim",sep=""))
            tmp2 = tmp[duplicated(tmp[,2]),]
            tmp = as.data.frame(tmp)
            tmp2 = as.data.frame(tmp2)
            tmp = tmp[!tmp[,2] %in% tmp2[,2],]
            
            tmp2 = exp_dat$SNP
            tmp2 = tmp2[tmp2 %in% tmp[,2]]
            snp_id = c(snp_id,tmp2)
        }
        exp_dat = exp_dat[exp_dat$SNP %in% snp_id,]
        
        dim(exp_dat)
        snp.list = NULL
        
        
        chr.list = unique(as.numeric(exp_dat$chr.exposure))
        chr.list = sort(chr.list)
        chr.list = chr.list[!is.na(chr.list)]
        
        
        for(chr.id in chr.list) {
            exp.name = gsub(".RData","",wgt_list[i])
            tmp.snpfile = paste(tmp.dir, exp.name,"_",chr.id,".txt", sep = "")
            exp_dat_tmp = exp_dat[exp_dat$chr.exposure==chr.id,,drop=FALSE]
            write.table(exp_dat_tmp$SNP, file = tmp.snpfile, quote = FALSE, col.names = FALSE, row.names = FALSE)
            
            tmp.ssfile = paste(tmp.dir, exp.name,"_",chr.id, "_expsoure.txt", sep = "")
            write.table(exp_dat_tmp, file = tmp.ssfile, quote = FALSE, col.names = TRUE, row.names = FALSE)
            
            system(paste0("plink --bfile ", "/gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.CHR", chr.id, " --memory 4000 --silent --chr ", chr.id, "  --extract ", tmp.snpfile, " --make-bed --out ", tmp.dir, exp.name,"_CHR",chr.id))
            
            command <- paste("plink --bfile ",tmp.dir, exp.name,"_CHR",chr.id, "  --memory 4000 --silent --clump ", tmp.ssfile, "  --clump-p1 ",cutoff," --clump-kb 10000 --clump-r2 0.001  --clump-snp-field SNP --clump-field pval.exposure --out ", tmp.dir, exp.name,"_CHR",chr.id, sep = "")
            system(command)
            
            tmp <- fread(paste(tmp.dir, exp.name,"_CHR",chr.id, ".clumped", sep = ""))
            tmp = as.data.frame(tmp)
            
            snp.list = c(snp.list,tmp[,"SNP"])
        }
        
        ##################################################################
        ind.pre.SNP = snp.list
        
        #exp_dat$mr_keep.exposure = TRUE
        avali.snp = exp_dat0[,"SNP"]
        avali.snp = avali.snp[avali.snp%in% IGAP[,"MarkerName"]]
        
        instruments = snp.list[snp.list %in% avali.snp]
        
        proxy.snp = snp.list[!snp.list %in% avali.snp]
        exp_tmp = exp_dat0[exp_dat0[,"SNP"] %in% proxy.snp,]
        proxy.snp = exp_tmp[,"SNP"]
        
        tmp.snpfile2 = paste(tmp.dir, exp.name,"_2.txt", sep = "")
        write.table(exp_dat0$SNP, file = tmp.snpfile2, quote = FALSE, col.names = FALSE, row.names = FALSE)
        
        proxy.inf = NULL
        if(length(proxy.snp)>0) {
            proxy.inf = as.data.frame(matrix(NA, length(proxy.snp),3))
            colnames(proxy.inf) = c("oringal","proxy","r2")
            proxy.inf[,1] = proxy.snp
            
            for(j in 1:length(proxy.snp)) {
                chr.id = exp_tmp[j,"chr.exposure"]
                tmp.snp = exp_tmp[j,"SNP"]
                tmp.matrixfile = paste(tmp.dir, exp.name,j,"_",tmp.snp, sep = "")
                P0 = exp_tmp[j,"pos.exposure"] - 1000000
                P1 = exp_tmp[j,"pos.exposure"] + 1000000
                
                if(P0 <0) {
                    P0  = 0
                }
                
                error.ind = FALSE
                tryCatch({
                    system(paste0("plink --bfile ", "/gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.CHR", chr.id, " --memory 4000 --silent --chr ", chr.id, " --from-kb ",P0/1e3," --to-kb ",P1/1e3, " --extract ", tmp.snpfile2, " --make-bed --out ",tmp.matrixfile))
                }, error = function(e) {
                    error.ind  = TRUE
                })
                
                if(error.ind) {
                    next
                }
                # read the data  --silent
                genos = BEDMatrix(tmp.matrixfile)
                geno <- as.matrix(genos)
                
                bim <- fread(paste0(tmp.matrixfile,".bim"))
                bim <- as.data.frame(bim)
                
                bim[bim[,2]==tmp.snp,]
                colnames(geno) = gsub("_.*","",colnames(geno))
                
                cor = cor(geno[,tmp.snp],geno,use = "pairwise.complete.obs")
                cor = cor[,cor^2>0.8,drop=FALSE]
                cor = cor[,!colnames(cor) %in%tmp.snp,drop=FALSE]
                cor = cor[,colnames(cor) %in% avali.snp,drop=FALSE]
                cor = cor[,order(cor[1,],decreasing=TRUE),drop=FALSE]
                
                if(dim(cor)[2]>=1) {
                    proxy.inf[j,2] = colnames(cor)[1]
                    proxy.inf[j,3] = cor[1,1]
                }
            }
            instruments = c(instruments,proxy.inf[,2])
            instruments = instruments[!is.na(instruments)]
        }
        
        cat("Finish instruments identification: len ",length(instruments),"\n")
        system(paste0("rm ",tmp.dir, exp.name,"*"))
        
        
        if(length(instruments)<2) {
            out.name = paste(out.dir,wgt_list[i],sep="")
            exp.name = NULL
            save(exp.name,file = out.name)
            next
        }
        
        exp_dat = exp_dat0[exp_dat0[,"SNP"] %in% instruments,]
        
        out_dat = IGAP[IGAP[,"MarkerName"] %in% exp_dat[,"SNP"],]
        out_dat$samplesize = 54162
        out_dat$outcome = "Alzheimer's disease (late onset)"
        
        out_dat = out_dat[,c(3,1,2,6,7,10,8,9,4,5,11)]
        colnames(out_dat) = c("SNP","chr","pos","beta.outcome","se.outcome","samplesize.outcome","pval.outcome","eaf.outcome","effect_allele.outcome","other_allele.outcome","outcome")
        out_dat$id.outcome = "ebi-a-GCST002245"
        out_dat$mr_keep.outcome = TRUE
        
        cat("IGAP data preparation\n")
        
        if(job==1) {
            dat <- harmonise_data(
            exposure_dat = exp_dat,
            outcome_dat = out_dat,
            action = 2
            )
        } else if(job==2) {
            dat <- harmonise_data(
            exposure_dat = exp_dat,
            outcome_dat = out_dat,
            action = 3
            )
        }
        
        
        if(dim(dat)[1]<2) {
            out.name = paste(out.dir,wgt_list[i],sep="")
            exp.name = NULL
            save(exp.name,file = out.name)
            
            next
        }
        if(sum(dat$mr_keep)<2) {
            out.name = paste(out.dir,wgt_list[i],sep="")
            exp.name = NULL
            save(exp.name,file = out.name)
            
            next
        }
        
        cat("Harmizoning data\n")
        
        dat = dat[dat$mr_keep,]
        b_exp = dat$beta.exposure
        b_out = dat$beta.outcome
        se_exp = dat$se.exposure
        se_out = dat$se.outcome
        
        BetaXG   = b_exp
        seBetaXG = se_exp
        
        Fall   = BetaXG^2/seBetaXG^2
        mF  = mean(Fall)
        
        #unweightedIsq
        Isq = Isq(BetaXG,seBetaXG) #unweighted
        
        ivw.res = mr_ivw_fe(b_exp,b_out,se_exp,se_out)
        
        egger.res = mr_egger_regression(b_exp,b_out,se_exp,se_out)
        
        # intercept result:
        #coef(egger.res[[11]])[1,]
        
        mr.res = mr(dat,method_list = c("mr_ivw","mr_ivw_fe","mr_weighted_median","mr_weighted_mode"))
        presso.res = NULL
        
        #tryCatch({
        #presso.res = run_mr_presso(dat, NbDistribution = 5000, SignifThreshold = 0.05)
        #}, error = function(e) {
        #   cat("ERROR msg :", conditionMessage(e), "\n")
        #})
        exp.name = gsub(".RData","",wgt_list[i])
        
        out.name = paste(out.dir,wgt_list[i],sep="")
        save(mF,Isq,ind.pre.SNP,proxy.inf,instruments,exp_dat,out_dat,dat,mr.res,ivw.res,egger.res,presso.res,exp.name,file = out.name)
        
        # output data
        write.csv(dat,file= paste(out.dir,exp.name,".csv",sep=""),col.names=TRUE,row.names=FALSE, quote=FALSE)
        
        #res_loo <- mr_leaveoneout(dat)
        cat("Finish Analysis\n")
        
        #directionality_test(dat)
        cat(wgt_list[i],"indx :",i,"\n")
    }, error = function(e) {
        id = wgt_list[i]
        save(job.id,id,file=paste("/gpfs/research/chongwu/Chong/Application/MR_AD/error_msg/msg_",job.id,"_",i,".RData",sep=""))
        cat(wgt_list[i]," wrong indx :",i,"\n")
        
        cat("ERROR2 :", conditionMessage(e), "\n")
    })
}


#dat2 = dat[dat$palindromic==FALSE,]
#mr.res2 = mr(dat2,method_list = mr_method_list()$obj)

system("export PATH=$PATH:/gpfs/research/chongwu/shared/software")


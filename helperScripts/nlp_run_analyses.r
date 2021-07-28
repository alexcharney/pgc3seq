#!/usr/bin/Rscript
 
##Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages(library(data.table))
suppressMessages(library(stats))
suppressMessages(library(rms))
suppressMessages(library(metafor))
args <- commandArgs(trailingOnly=TRUE)
hlp <- fread(args[[1]], header=F)
out <- args[[2]]
colnames(hlp) <- c("ancestryCol", "ancestryValue", "covarType", "pheType", "pheDx", "pheCol", "genScoreType", "genScoreCol")
dt <- fread("/sc/arion/projects/psychgen2/psychosis_nlp/files/BioMe_PsychosisNLP_MASTER_TABLE.tsv", na=c("", "NA"))
szprs.col <- grep("sczprs", colnames(dt), value=T)
cgprs.col <- grep("cogprs", colnames(dt), value=T)
edprs.col <- grep("eduprs", colnames(dt), value=T)
cnvct.col <- grep("cnvbur", colnames(dt), value=T)
rvprs.col <- grep("rvPRS", colnames(dt), value=T)
score.col <- list("prs"=c(szprs.col, cgprs.col, edprs.col), "rvprs"=rvprs.col, "cnv"=cnvct.col)
phers.col <- list( "psychosis"=list("raw"="psychosis.phers", "weighted"="psychosis.phers.weighted"),
                 "mania"=list("raw"="mania.phers", "weighted"="mania.phers.weighted"))
myicd.col <- list( "scz"=list("count"="scz.icd.count", "binary"="scz.icd.binary"),
                 "bip"=list("count"="bip.icd.count", "binary"="bip.icd.binary"),
                 "cog"=list("count"="cog.icd.count", "binary"="cog.icd.binary"),
                 "psychosis"=list("count"="psychosis.icd.count", "binary"="psychosis.icd.binary"))
myphe.col <- list("phers"=phers.col, "icd"=myicd.col)
genpc.col <- grep("^PC", colnames(dt), value=T)
admix.col <- c("ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI")
covar.col <- list("pca"=genpc.col, "admix"=admix.col)
ancat.col <- c("AncestryClass", "gill.ContinentalGrouping","gill.IBDcommunity")

##Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

res <- c()
for( i in 1:nrow(hlp)){
    cat(i, "of", nrow(hlp), "\n")
    cur <- hlp[i]
    phe <- cur$pheCol
    gen <- cur$genScoreCol
    cov <- covar.col[[cur$covarType]]
    anc <- cur$ancestryCol
    anv <- cur$ancestryValue
    if ( is.na(anc) ){
        cov <- covar.col[[cur$covarType]]
        mydata <- dt[!is.na(get(phe)) & !is.na(get(gen)) & !is.na(get(cov[1]))]
    } else {
        cov <- covar.col[[cur$covarType]][1:5]
        mydata <- dt[!is.na(get(phe)) & !is.na(get(gen)) & !is.na(get(cov[1])) & get(anc)==anv]
    }
    if (cur$genScoreType == "cnv") {
        cov <- c(cov, "chip")
    }
    formula1 <- as.formula(paste(phe, " ~ ", paste(c(gen, cov), collapse=" + "), sep=""))    
    formula1x <- as.formula(paste(phe, " ~ ", paste(c(paste0("I(", gen, "/1000)"), cov), collapse=" + "), sep=""))    
    formula2 <- as.formula(paste(phe, " ~ ", paste(cov, collapse=" + "), sep=""))
    if ( is.na(anc) ){    
        x <- summary(glm(formula1, data=mydata))
        m1 <- lrm(formula1x, data=mydata)
        m2 <- lrm(formula2, data=mydata)
        cur$p <- x$coef[gen,"Pr(>|t|)"]
        cur$t <- x$coef[gen,"t value"]
        cur$beta <- x$coef[gen,"Estimate"]
        cur$betase <- x$coef[gen,"Std. Error"]
        cur$r2 <- m1$stats["R2"] - m2$stats["R2"]
        cur$n <- nrow(mydata)
    } else {
        nind <- nrow(mydata)
        nphe <- uniqueN(mydata[[phe]])
        ngen <- uniqueN(mydata[[gen]])
        check1 <- nind > 10
        check2 <- nphe > 1
        check3 <- ngen > 1
        check4 <- sum(table(mydata[[phe]]) <= 1) == 0
        if (check1 & check2 & check3 & check4){
            x <- summary(glm(formula1, data=mydata))
            m1 <- lrm(formula1x, data=mydata)
            m2 <- lrm(formula2, data=mydata)
            cur$p <- x$coef[gen,"Pr(>|t|)"]
            cur$t <- x$coef[gen,"t value"]
            cur$beta <- x$coef[gen,"Estimate"]
            cur$betase <- x$coef[gen,"Std. Error"]
            cur$r2 <- m1$stats["R2"] - m2$stats["R2"]           
            cur$n <- nind
        } else {
            cur$p <- NA
            cur$t <- NA
            cur$beta <- NA
            cur$betase <- NA
            cur$r2 <- NA            
            cur$n <- NA
        }
    }
    res <- rbind(res,cur)
}

## WRITE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fwrite (res, na="NA", sep='\t', row=F, quo=F, file=out)
cat("DONE\n")



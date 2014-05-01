# process fpkm data.----------------

fpkm_expres <- fpkm[c("gene_id", "gene_short_name", grep("FPKM", colnames(fpkm), value=T))]
library(reshape2)
fpkm_expres <- melt(fpkm_expres, id.vars=c("gene_id", "gene_short_name"))
fpkm_ci <- fpkm[c("gene_id", grep("conf_hi", colnames(fpkm), value=T))]
fpkm_ci <- melt(fpkm_ci, id.vars=c("gene_id"))
dat <- cbind(fpkm_expres, fpkm_ci)
dat <- dat[c(1,2,3,4,7)]
colnames(dat) <- c("gene_id", "gene_short_name", "comparison", "fpkm", "ci_hi")
dat$error <- dat$ci_hi - dat$fpkm
dat$sdev <- dat$error / 3.94
dat$sem <- dat$sdev/sqrt(2)
dat$hi_error <- dat$fpkm + dat$sem
dat$lo_erro <- dat$fpkm - dat$sem

dat <- dat[c("gene_id", "gene_short_name", "comparison", "fpkm", "hi_error", "lo_erro")]
colnames(dat) <- c("ensembl_id", "gene", "sample", "fpkm", "hi_error", "lo_error")

dat$sample_good_name <- NULL
d1 <- dat
d1[(d1$sample == "embryonic_MN_FPKM"), "sample_good_name"] <- "Embryonic MNs"
d1[(d1$sample == "iMN_FPKM"), "sample_good_name"] <- "iMNs"
d1[(d1$sample == "mES_MN_FPKM"), "sample_good_name"] <- "ESC MNs"
d1[(d1$sample == "miPS_MN_FPKM"), "sample_good_name"] <- "iPSC MNs"
d1[(d1$sample == "mES_FPKM"), "sample_good_name"] <- "ESC"
d1[(d1$sample == "miPS_FPKM"), "sample_good_name"] <- "iPSC"
d1[(d1$sample == "Fib_d15_N3_FPKM"), "sample_good_name"] <- "MEF1"
d1[(d1$sample == "Fib_d15_FBS_FPKM"), "sample_good_name"] <- "MEF2"
d1[(d1$sample == "Fib_d15_FBS_1p_FPKM"), "sample_good_name"] <- "MEF3"
d1[(d1$sample == "Fib_d2_FBS_FPKM"), "sample_good_name"] <- "MEF4"

saveRDS(object=d1, file="iMN/data/genesFpkm.rds")

  

# process diff data------------------
diff <- readRDS("iMN/data/iMN_data.rds")

### swap diff names... 
d1 <- diff
d1$sample_1 <- as.character(d1$sample_1)
d1[(d1$sample_1 == "embryonic_MN"), "sample_1"] <- "Embryonic MNs"
d1[(d1$sample_1 == "iMN"), "sample_1"] <- "iMNs"
d1[(d1$sample_1 == "mES_MN"), "sample_1"] <- "ESC MNs"
d1[(d1$sample_1 == "miPS_MN"), "sample_1"] <- "iPSC MNs"
d1[(d1$sample_1 == "mES"), "sample_1"] <- "ESC"
d1[(d1$sample_1 == "miPS"), "sample_1"] <- "iPSC"
d1[(d1$sample_1 == "Fib_d15_N3"), "sample_1"] <- "MEF1"
d1[(d1$sample_1 == "Fib_d15_FBS"), "sample_1"] <- "MEF2"
d1[(d1$sample_1 == "Fib_d15_FBS_1p"), "sample_1"] <- "MEF3"
d1[(d1$sample_1 == "Fib_d2_FBS"), "sample_1"] <- "MEF4"

d1$sample_2 <- as.character(d1$sample_2)
d1[(d1$sample_2 == "embryonic_MN"), "sample_2"] <- "Embryonic MNs"
d1[(d1$sample_2 == "iMN"), "sample_2"] <- "iMNs"
d1[(d1$sample_2 == "mES_MN"), "sample_2"] <- "ESC MNs"
d1[(d1$sample_2 == "miPS_MN"), "sample_2"] <- "iPSC MNs"
d1[(d1$sample_2 == "mES"), "sample_2"] <- "ESC"
d1[(d1$sample_2 == "miPS"), "sample_2"] <- "iPSC"
d1[(d1$sample_2 == "Fib_d15_N3"), "sample_2"] <- "MEF1"
d1[(d1$sample_2 == "Fib_d15_FBS"), "sample_2"] <- "MEF2"
d1[(d1$sample_2 == "Fib_d15_FBS_1p"), "sample_2"] <- "MEF3"
d1[(d1$sample_2 == "Fib_d2_FBS"), "sample_2"] <- "MEF4"

d1$sample_1 <- factor(d1$sample_1, levels=c(
  "Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", "MEF1", "MEF2", "MEF3", "MEF4","ESC", "iPSC"), 
  ordered=T)

d1$sample_2 <- factor(d1$sample_2, levels=c(
  "Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", "MEF1", "MEF2", "MEF3", "MEF4","ESC", "iPSC"), 
  ordered=T)

saveRDS(d1, "iMN/data/iMN_diff.rds")

# process isoform data. -------------------
isoforms <- read.table("~/host/isoforms.fpkm_tracking", sep="\t", header=T)
# keep only fpkm columns, fix sample names and write as RDS.
isoforms <- isoforms[c("tracking_id", "gene_id", "gene_short_name", "length", grep("FPKM", colnames(isoforms), value=T))]
isoforms.m <- melt(isoforms, id.vars=c("tracking_id", "gene_id", "gene_short_name", "length"))

isoforms.m$sample_good_name <- NULL

isoforms.m[(isoforms.m$variable == "embryonic_MN_FPKM"), "sample_good_name"] <- "Embryonic MNs"
isoforms.m[(isoforms.m$variable == "iMN_FPKM"), "sample_good_name"] <- "iMNs"
isoforms.m[(isoforms.m$variable == "mES_MN_FPKM"), "sample_good_name"] <- "ESC MNs"
isoforms.m[(isoforms.m$variable == "miPS_MN_FPKM"), "sample_good_name"] <- "iPSC MNs"
isoforms.m[(isoforms.m$variable == "mES_FPKM"), "sample_good_name"] <- "ESC"
isoforms.m[(isoforms.m$variable == "miPS_FPKM"), "sample_good_name"] <- "iPSC"
isoforms.m[(isoforms.m$variable == "Fib_d15_N3_FPKM"), "sample_good_name"] <- "MEF1"
isoforms.m[(isoforms.m$variable == "Fib_d15_FBS_FPKM"), "sample_good_name"] <- "MEF2"
isoforms.m[(isoforms.m$variable == "Fib_d15_FBS_1p_FPKM"), "sample_good_name"] <- "MEF3"
isoforms.m[(isoforms.m$variable == "Fib_d2_FBS_FPKM"), "sample_good_name"] <- "MEF4"

head(isoforms.m)
isoforms.m <- isoformFPKM
isoforms.m$sample_good_name <- factor(isoforms.m$sample_good_name, levels=c(
  "Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", "MEF1", "MEF2", "MEF3", "MEF4","ESC", "iPSC"), 
  ordered=T)
                                      
saveRDS(isoforms.m, "iMN/data/isoform_fpkm.rds")

# process methylation data --------------------

meth <- read.table("~/host/all.promoters_4_1txt.txt", sep="\t", header=T)
meth$position <- paste0(meth$chrom,":", meth$chromStart, "-", meth$chromEnd)

meth$Embryonic_MN_meth <- rowMeans(meth[c(7,8,20)], na.rm=T)
meth$iMN_meth <- rowMeans(meth[c(15,21)], na.rm=T)
meth$es_MN_meth <- rowMeans(meth[c(11,12)], na.rm=T)
meth$ips_MN_meth <- rowMeans(meth[c(13,14)], na.rm=T)
meth$mef1_meth <- rowMeans(meth[c(27,28)], na.rm=T)
meth$ips_meth <- rowMeans(meth[c(52,53)], na.rm=T)
meth$es_meth <- rowMeans(meth[c(50,51)], na.rm=T)
meth$mef2_meth <- rowMeans(meth[c(32,33)], na.rm=T)
meth$mef3_meth <- rowMeans(meth[c(30,31)], na.rm=T)
meth$mef4_meth <- rowMeans(meth[c(34,35)], na.rm=T)
meth$mef5_meth <- rowMeans(meth[c(48,49)], na.rm=T)

meth2 <- select(meth, position, name, geneNames, ensemblIds, Embryonic_MN_meth, iMN_meth, 
                es_MN_meth, ips_MN_meth, mef1_meth, ips_meth, es_meth, mef2_meth, 
                mef3_meth, mef4_meth, mef5_meth)

# for now, just include average methylation levels, stats test are 1 sided - this presents some challenge 
# to display this data. 

meth2 <- melt(meth2, id.vars=c("position", "geneNames", "ensemblIds", "name"))

# fix sample names. 

meth2$variable <- as.character(meth2$variable)
meth2[(meth2$variable == "Embryonic_MN_meth"), "variable"] <- "Embryonic MNs"
meth2[(meth2$variable == "iMN_meth"), "variable"] <- "iMNs"
meth2[(meth2$variable == "es_MN_meth"), "variable"] <- "ESC MNs"
meth2[(meth2$variable == "ips_MN_meth"), "variable"] <- "iPSC MNs"
meth2[(meth2$variable == "mef1_meth"), "variable"] <- "MEF1"
meth2[(meth2$variable == "ips_meth"), "variable"] <- "iPSC"
meth2[(meth2$variable == "es_meth"), "variable"] <- "ESC"
meth2[(meth2$variable == "mef2_meth"), "variable"] <- "MEF2"
meth2[(meth2$variable == "mef3_meth"), "variable"] <- "MEF3"
meth2[(meth2$variable == "mef4_meth"), "variable"] <- "MEF4"
meth2[(meth2$variable == "mef5_meth"), "variable"] <- "MEF5"

# order samples. 
meth2$variable <- factor(meth2$variable, levels=c(
  "Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", "MEF1", "MEF2", "MEF3", "MEF4", "MEF5", "ESC", "iPSC"), 
  ordered=T)

saveRDS(meth2, "iMN/data/methylation.rds")


meth <- readRDS("iMN/data/methylation.rds")
  GOI <- "Sox17"
methDat <- meth
methPlot <- function(GOI, methDat){
  res <- methDat[(methDat$geneNames %in% GOI),]
  methP <- ggplot(res, aes(x=variable, y=value, fill=variable))+geom_bar(stat="identity")
  methP + scale_y_continuous(expand=c(0,0), limits=c(0,1.05))+theme_bw()+
    scale_fill_manual(values=c("#666699", "#0DB14B", "#008ccf", "#118697", 
                            "#744c28", "#8a5d3b", "#9a8478", "#594a41", "tan",
                            "#d8403f", "#f7931d"))+xlab("")+ylab("Methylation Ratio")+
    geom_hline(y=.2, color="#bdbdbd", linetype="dashed")+
    geom_hline(y=.6, color="#bdbdbd", linetype="dashed")+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12), legend.position="none")
}
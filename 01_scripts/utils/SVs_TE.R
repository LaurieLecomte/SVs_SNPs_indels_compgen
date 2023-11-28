#setwd("E:/coregonus/genotyping_SV/08_fasta_SV/")

library(tidyr)
library(dplyr)
library(tibble)


# 1. Access files in command line, import and format ----------------------
argv <- commandArgs(T)
#RM_TABLE <- argv[1]
RM_TABLE <- 'TE/SVs_CHR_POS_ID_REF_ALT.fasta.out.table'
  SV_INFOS <- argv[2]

SV_INFOS <- 'TE/SVs_CHR_POS_ID_REF_ALT.info'
#repeat_tab_full<-read.table("ALLDP1_MISS50_2all.variants.TEreformatted", header=TRUE, stringsAsFactors = FALSE, fill=T)
repeat_tab_full <- read.table(RM_TABLE, header=FALSE, stringsAsFactors = FALSE, fill=T,
                              col.names = c('SW_score', 'perc_div', 'perc_del', 'perc_ins', 'query', 
                                            'start_in_query', 'end_in_query', 'left_in_query', '+', 
                                            'matching_repeat', 'repeat_class_fam', 
                                            'start_in_repeat', 'end_in_repeat', 'left_in_repeat', 
                                            'ID', 'star'))

#info SV
SVinfo<-read.table(SV_INFOS, header = FALSE)
colnames(SVinfo) <-c("CHROM","POS", "ID", "REF", "ALT","ref_len", "alt_len", "LEN", "TYPE", "start_bed", "stop_bed", "chr_pos", "max_len")


head(SVinfo)
SVinfo<-mutate(SVinfo, SV_true_len=pmax(ref_len,alt_len, na.rm=T))




head(repeat_tab_full)
#stars in the last column indicate that there is a better match that overlap. After looking at some exemples I think it is worth removing the star lines

repeat_tab<-repeat_tab_full[-which(repeat_tab_full$star=="*"),]
repeat_tab$query_length<-repeat_tab$query_end-repeat_tab$query_begin

head(repeat_tab)
#repeat elements diversity
elements<-levels(repeat_tab$class_family)
cbind(table(repeat_tab$class_family))
#create a new factor more simple
repeat_tab$rep_class<-repeat_tab$class_family
repeat_tab$rep_class[grep("LINE", repeat_tab$class_family)]<-"line"
repeat_tab$rep_class[grep("DNA", repeat_tab$class_family)]<-"dnaTE"
repeat_tab$rep_class[grep("RNA", repeat_tab$class_family)]<-"rna"
repeat_tab$rep_class[grep("LTR", repeat_tab$class_family)]<-"ltr"
repeat_tab$rep_class[grep("SINE", repeat_tab$class_family)]<-"sine"
cbind(table(repeat_tab$rep_class))
head(repeat_tab)

#split Sv_name
repeat_tab<-separate(repeat_tab, CHR, c("SV_id","allele"), sep="-")
head(repeat_tab)

#View((repeat_tab))

#loop over annotated SV to simplify the annotation
SVlist<-levels(as.factor(repeat_tab$SV_id))
nSV<-length(SVlist)

SVTE<-matrix(ncol=9)
colnames(SVTE)<-c("Sv_id","allele","TE_class","TE_len", "TE_div","TE_maj","TE_maj_len","TE_maj_div","TE_len_tot")

for (i in(1:nSV))
{
  repeat_tab_i<-repeat_tab[repeat_tab$SV_id==SVlist[i],]
  
  #case in which only the REf or the ALT is annotated
    if (length(levels(as.factor(repeat_tab_i$allele)))==1)
    {
      len_by_class<-rownames_to_column(as.data.frame(cbind(by(repeat_tab_i$query_length, repeat_tab_i$class_family, sum),by(repeat_tab_i$perc_div, repeat_tab_i$class_family, mean))))#problem is output with rownames is annoying
      #len_by_class<-cbind(summarise(group_by(repeat_tab_i,class_family),len=sum(query_length), .groups="rowwise"))#problem is that it does not work when just one level of factor
      TE_maj_len<-max(len_by_class$V1)
      TE_maj<-len_by_class$rowname[len_by_class$V1==TE_maj_len]
      div_TE_maj<-len_by_class$V2[len_by_class$V1==TE_maj_len]
      TE_len_sum<-sum (len_by_class$V1)
      n_class<-dim(len_by_class)[1]
      svte_i<-cbind(rep(SVlist[i], n_class), rep(repeat_tab_i$allele[1], n_class), len_by_class, 
                    rep(TE_maj, n_class),rep(TE_maj_len, n_class),rep(div_TE_maj, n_class), rep(TE_len_sum, n_class))
      
      colnames(svte_i)<-c("Sv_id","allele","TE_class","TE_len","TE_div","TE_maj","TE_maj_len","TE_maj_div","TE_len_tot")
      SVTE<-rbind(SVTE, svte_i)
      
    }
    
    
    #case in which both the REf or the ALT is annotated
    if (length(levels(as.factor(repeat_tab_i$allele)))==2)
    {
      print (i)
      #print(repeat_tab_i)
      len_by_class<-rownames_to_column(as.data.frame(cbind(by(repeat_tab_i$query_length, list(repeat_tab_i$class_family,repeat_tab_i$allele), sum),by(repeat_tab_i$perc_div, list(repeat_tab_i$class_family,repeat_tab_i$allele), mean))))#problem is output with rownames is annoying
      
      #len_by_class<-rownames_to_column(as.data.frame(cbind(by(repeat_tab_i$query_length, list(repeat_tab_i$class_family,repeat_tab_i$allele), sum))))#problem is output with rownames is annoying
      #keep the longest
      len_by_class<-mutate(len_by_class, len_max=pmax(ALT,REF, na.rm=T))
      TE_maj_len<-max(len_by_class$len_max, na.rm=T)
      TE_maj<-len_by_class$rowname[len_by_class$len_max==TE_maj_len]
      div_TE_maj<-(len_by_class[len_by_class$len_max==TE_maj_len,4]+len_by_class[len_by_class$len_max==TE_maj_len,5])/2
      TE_len_sum<-sum (len_by_class$len_max, na.rm=T)
      n_class<-dim(len_by_class)[1]
      svte_i<-cbind(rep(SVlist[i], n_class), rep("BOTH", n_class), len_by_class[,c(1,6,5)], 
                    rep(TE_maj, n_class),rep(TE_maj_len, n_class),rep(div_TE_maj, n_class),rep(TE_len_sum, n_class))
      
      colnames(svte_i)<-c("Sv_id","allele","TE_class","TE_len","TE_div","TE_maj","TE_maj_len","TE_maj_div","TE_len_tot")
      SVTE<-rbind(SVTE, svte_i)
    }

    
}

SVTE2<-SVTE[2:dim(SVTE)[1],]
  
  



SVinfo_TE<-left_join(SVinfo, SVTE2)
SVinfo_TE$TE_prop<-round(SVinfo_TE$TE_maj_len/SVinfo_TE$SV_true_len,2)
head(SVinfo_TE)
SVinfo_TE$TE_prop[which(is.na(SVinfo_TE$TE_prop))]<-0
SVinfo_TE$TE_prop[SVinfo_TE$TE_prop>=1]<-1

#some check
hist(SVinfo_TE$TE_prop[SVinfo_TE$TE_prop<=1], breaks=80)
hist(SVinfo_TE$SV_true_len[SVinfo_TE$SV_true_len<=4000 & SVinfo_TE$SV_true_len>=200], breaks = 60)
SVinfo_TE[SVinfo_TE$TE_prop>=10,]

length(levels(as.factor(as.character(SVinfo_TE$Sv_id[which(SVinfo_TE$TE_prop>=0.10 & SVinfo_TE$TE_prop<=0.5)]))))
length(levels(as.factor(as.character(SVinfo_TE$Sv_id[which(SVinfo_TE$TE_prop>0 & SVinfo_TE$TE_prop<=0.10)]))))
length(levels(as.factor(as.character(SVinfo_TE$Sv_id[which(SVinfo_TE$TE_prop==0 )]))))
length(levels(as.factor(as.character(SVinfo_TE$Sv_id[which(SVinfo_TE$TE_prop<=0.50 )]))))
length(levels(as.factor(as.character(SVinfo_TE$Sv_id[which(SVinfo_TE$TE_prop>=0.50 )]))))
length(levels(as.factor(as.character(SVinfo_TE$Sv_id[which(SVinfo_TE$TE_prop>=0.75 )]))))

#removing if cover <50% of the length
SVinfo_TE$TE_class2<-SVinfo_TE$TE_class
SVinfo_TE$TE_class2[SVinfo_TE$TE_prop<=0.50]<-"no_TE"
SVinfo_TE$TE_maj2<-SVinfo_TE$TE_maj
SVinfo_TE$TE_maj2[SVinfo_TE$TE_prop<=0.50]<-"no_TE"

head(SVinfo_TE,10)

#create a new factor more simple
SVinfo_TE$TE_maj2_class<-SVinfo_TE$TE_maj2
SVinfo_TE$TE_maj2_class[grep("LINE", SVinfo_TE$TE_maj2)]<-"line"
SVinfo_TE$TE_maj2_class[grep("DNA", SVinfo_TE$TE_maj2)]<-"dnaTE"
SVinfo_TE$TE_maj2_class[grep("RNA", SVinfo_TE$TE_maj2)]<-"rna"
SVinfo_TE$TE_maj2_class[grep("LTR", SVinfo_TE$TE_maj2)]<-"ltr"
SVinfo_TE$TE_maj2_class[grep("SINE", SVinfo_TE$TE_maj2)]<-"sine"



SVinfo_TE_unique<-SVinfo_TE %>% distinct(Sv_id, .keep_all = T)
cbind(table(SVinfo_TE_unique$TE_maj2_class))

SVinfo_TE_unique_summary<-as.data.frame(summarise(group_by(SVinfo_TE_unique,TE_maj2_class), n_SV=n()))
SVinfo_TE_unique_summary_order<-SVinfo_TE_unique_summary[order(SVinfo_TE_unique_summary$n_SV,decreasing=T),]

pie(SVinfo_TE_unique_summary_order$n_SV, SVinfo_TE_unique_summary_order$TE_maj2_class, col=c("grey","black","orange","red","chartreuse","cyan","blue","purple","yellow","grey50","pink"))


write.table(SVinfo_TE, "ALLDP1_MISS50_2all.variants.FULLTE", row.names=F, quote=F, sep="\t")
write.table(SVinfo_TE_unique, "ALLDP1_MISS50_2all.variants.annotated", row.names=F, quote=F, sep="\t")

#write.table(SVinfo_TE, "LR_SR_GA.FULLTE", row.names=F, quote=F, sep="\t")

#write.table(SVinfo_TE_unique, "LR_SR_GA.annotated", row.names=F, quote=F, sep="\t")


SVinfo_TE_unique<-read.table("LR_SR_GA.annotated", header=T)



write.table(SVinfo_TE, "all_samples_DP1_MISS50_MAF5_2all.variants.FULLTE", row.names=F, quote=F, sep="\t")

write.table(SVinfo_TE_unique, "all_samples_DP1_MISS50_MAF5_2all.variants.annotated", row.names=F, quote=F, sep="\t")
#!/usr/bin/env Rscript

# Load libraries:
library("tidyverse")
library("optparse")

# Setup arguments
option_list = list(
  make_option(c("-b", "--bsc"), type="character", default=NULL, 
              help="BSC SV file", metavar="character"),
  make_option(c("-c", "--charite"), type="character", default=NULL, 
              help="charite SV file", metavar="character"),
  make_option(c("-d", "--curie"), type="character", default=NULL, 
              help="curie SV file", metavar="character"),
  make_option(c("-H", "--hartwig"), type="character", default=NULL, 
              help="hartwig SV file", metavar="character"),
  make_option(c("-O", "--oicr"), type="character", default=NULL, 
              help="oicr SV file", metavar="character"),
  make_option(c("-t", "--truth"), type="character", default=NULL, 
              help="truth file", metavar="character"),
  make_option(c("-o", "--outputDir"), type="character", default=NULL, 
              help="oicr SV file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bsc) && is.null(opt$charite) && is.null(opt$curie) && is.null(opt$hartwig) && is.null(opt$oicr)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# testing:
#bscTable=read.table("/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/TCGA/results/DO32237_curie/SV_benchmark_results.csv", header=TRUE, sep =",",stringsAsFactors = FALSE)
#bscTable$Center="BSC"
#svTable=bscTable[-c(4,8,12,16,20),]

#truthTsv=read.table("/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/TCGA/results/DO32237_curie/truth_sv_dataframe.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")

#colo=read.table("/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/TRUTH/COLO829_truth.tsv",header=TRUE,stringsAsFactors = FALSE,sep="\t")

# Check if arg is present and add to svTable
svTable=data.frame()

if (!is.null(opt$bsc)){
  bscTable = read.table(opt$bsc, header=TRUE, sep =",",stringsAsFactors = FALSE)
  bscTable$Center="BSC"
  svTable=rbind(svTable,bscTable[-c(4,8,12,16,20),])
}

if (!is.null(opt$charite)){
  chariteTable = read.table(opt$charite, header=TRUE, sep =",",stringsAsFactors = FALSE)
  chariteTable$Center="Charite"
  svTable=rbind(svTable,chariteTable[-c(4,8,12,16,20),])
}

if (!is.null(opt$curie)){
  curieTable = read.table(opt$curie, header=TRUE, sep =",",stringsAsFactors = FALSE)
  curieTable$Center="Curie"
  svTable=rbind(svTable,curieTable[-c(4,8,12,16,20),])
}

if (!is.null(opt$hartwig)){
  hartwigTable = read.table(opt$hartwig, header=TRUE, sep =",",stringsAsFactors = FALSE)
  hartwigTable$Center="Hartwig"
  svTable=rbind(svTable,hartwigTable[-c(4,8,12,16,20),])
}

if (!is.null(opt$oicr)){
  oicrTable = read.table(opt$oicr, header=TRUE, sep =",",stringsAsFactors = FALSE)
  oicrTable$Center="OICR"
  svTable=rbind(svTable,oicrTable[-c(4,8,12,16,20),])
}

# Load truth table
truthTsv=read.table(opt$truth,header=TRUE,stringsAsFactors = FALSE,sep=",")

# Transform svTable: anonymize, change bin name
svTable=svTable %>% 
  mutate(Center = case_when(Center == "BSC" ~ "Node 1",
                            Center == "Curie" ~ "Node 2",
                            Center == "Charite" ~ "Node 3",
                            Center == "Hartwig" ~ "Node 4",
                            Center == "OICR" ~ "Node 5")) %>% 
  mutate(Center=factor(Center,levels=c("Node 1", "Node 2", "Node 3", "Node 4","Node 5"))) %>% 
  mutate(Bin.0.50.bp = case_when(Bin.0.50.bp == "All results" ~ "All",
                                 Bin.0.50.bp == "Bin 0-50 bp" ~ "0-50",
                                 Bin.0.50.bp == "Bin 50-200 bp" ~ "50-200",
                                 Bin.0.50.bp == "Bin 200-1000 bp" ~ "200-1000",
                                 Bin.0.50.bp == "Bin 1000-100000000000000000000000000000 bp" ~ "> 1000",
                                 Bin.0.50.bp == "Bin NaN bp" ~ "NaN")) %>% 
  dplyr::rename(Bin = Bin.0.50.bp)

write.table(svTable,file = paste0(opt$outputDir ,"svTable.csv"),row.names = FALSE)

print("svTable:")
print(svTable,row.names = FALSE)

# Transform svTable and extract only useful info
tidySV=svTable %>% 
  filter(Bin != "All" & TIER == "tier3") %>% 
  select(Bin,Center,starts_with("TP")& -"TP") %>% 
  mutate(Bin = case_when(Bin == "NaN" ~ "NaN",
                         Bin == "0-50" ~ "0-50",
                         Bin == "50-200" ~ "50-200",
                         Bin == "200-1000" ~ "200-1000",
                         Bin == "> 1000" ~ "> 1000")) %>% 
  mutate(TP_DEL=as.numeric(TP_DEL),TP_INS=as.numeric(TP_INS),TP_DUP=as.numeric(TP_DUP),TP_INV=as.numeric(TP_INV),TP_BND=as.numeric(TP_BND))

colnames(tidySV)=c("Bin","Center","DEL","INS","DUP","INV","BND")

print("truth")
print(head(truthTsv))

# Transform Truth Table
truthTsv= truthTsv %>% 
  mutate(ID=paste(start_chrom,start,sep = "_")) %>% 
  select(ID,type,length) %>% 
  mutate(bin = case_when(length == 0 | length == 1 | length == NA | type == "BND" ~ "NaN",
                         length > 1 & length <= 50 ~ "0-50",
                         length > 50 & length <= 200 ~ "50-200",
                         length > 200 & length <= 1000 ~ "200-1000",
                         length > 1000 ~ "> 1000"))

truthTsv$bin=factor(truthTsv$bin,levels=c("NaN","0-50","50-200","200-1000","> 1000"))
truthTsv$type=factor(truthTsv$type,levels=c("BND","DEL","DUP","INS","INV"))

print("truth")
print(truthTsv)

# Count the number of TP of each category and put everything together

SV_count=t(table(truthTsv$type,truthTsv$bin))
SV_count = as.data.frame.matrix(SV_count) %>% 
  rownames_to_column() %>% 
  dplyr::rename(Bin = rowname) %>%
  mutate(Center="Truth") %>% 
  select(Bin,Center,DEL,INS,DUP,INV,BND) %>% 
  bind_rows(tidySV)

SV_count$Bin=fct_relevel(as.factor(SV_count$Bin), "0-50", "50-200","200-1000","> 1000","NaN")

print("General matrix")
print(SV_count)

# Plot general metrics:
## Plot TP FP FN raw 

tidySvTP=svTable %>% 
  filter(Bin == "All" & TIER == "tier3") %>% 
  select(-c("FP_original", "FP_tier", "Recall", "Precision", "F1.score", "TP_DEL", "TP_INS", "TP_DUP", "TP_INV", "TP_BND")) %>% 
  pivot_longer(cols=c(TP,FP,FN), names_to="metric",values_to="count")

ggplot(tidySvTP,aes(x= Center,y=as.numeric(count),fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Set1") +
  theme_bw()+
  labs(x = "Centers",y="SV counts") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplot_SV_general_tier3_all.png"),width=30,height=20,units='cm')

## Plot TP FP FN by SV type

tidySvTP_all=svTable %>% 
  filter(Bin == "All" & TIER == "tier3") %>% 
  select(Bin,Center,starts_with("TP")) %>% 
  pivot_longer(cols=starts_with("TP"),names_to = "metric",values_to = "count")

ggplot(tidySvTP_all,aes(x= Center,y=as.numeric(count),fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  labs(x = "Centers",y="SV counts") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplot_SV_TP_tier3_all.png"),width=30,height=20,units='cm')

## Plot Precision Recall & F1 score
tidySvF1=svTable %>% 
  filter(Bin == "All" & TIER == "tier3") %>% 
  select(-c(TIER,FP,FP_original,FP_tier,FN,TP_DEL,TP_INS,TP_DUP,TP_INV,TP_BND )) %>% 
  pivot_longer(cols = c(Recall,Precision,F1.score), names_to = "metric", values_to = "count") %>% mutate(count=as.numeric(count),TP=as.numeric(TP))


tidySvF1$metric=fct_relevel(as.factor(tidySvF1$metric), "Recall","Precision","F1.score")

ggplot(tidySvF1,aes(x= Center,y=count*100,fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  labs(x = "Centers",y="Score (%)") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplot_SV_F1_tier3_all.png"),width=30,height=20,units='cm')

## Compare to truth:

tidySV_count=SV_count %>% 
  pivot_longer(cols=c(DEL,INS,DUP,INV,BND),names_to = "type",values_to="count")

ggplot(tidySV_count,aes(x=Bin,y=count,fill=Center)) +
  geom_bar(stat="identity",position="dodge") +
  facet_wrap(~type) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplot_SV_truth_types.png"),width=30,height=20,units='cm')

ggplot(mapping = aes(x=Bin,y=count,fill=Center)) +
  geom_bar(data=tidySV_count[which(tidySV_count$Center != "Truth"),],stat="identity",position="dodge") +
  geom_bar(data=tidySV_count[which(tidySV_count$Center == "Truth"),],stat="identity",position="dodge",linetype="dashed", colour="black",fill=NA) +
  facet_wrap(~type,scales = "free") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  labs(x = "SV length (bp)",y="SV counts") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        strip.text.x = element_text(size = 14),
        legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplot_SV_truth_types.png"),width=30,height=20,units='cm')

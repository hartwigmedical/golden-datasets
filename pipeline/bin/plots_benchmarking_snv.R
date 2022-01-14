#!/usr/bin/env Rscript

# Load libraries:
# library("tidyverse")
library("DT")
library("optparse")
library("ggplot2")
library("dplyr")
library("tidyr")
library("forcats")
library("tibble")

# Setup arguments
option_list = list(
    make_option(c("-b", "--bsc"), type="character", default=NULL,
                help="BSC snv file", metavar="character"),
    make_option(c("-c", "--charite"), type="character", default=NULL,
                help="charite snv file", metavar="character"),
    make_option(c("-d", "--curie"), type="character", default=NULL,
                help="curie snv file", metavar="character"),
    make_option(c("-H", "--hartwig"), type="character", default=NULL,
                help="hartwig snv file", metavar="character"),
    make_option(c("-O", "--oicr"), type="character", default=NULL,
                help="oicr snv file", metavar="character"),
    make_option(c("-o", "--outputDir"), type="character", default=NULL,
                help="oicr snv file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bsc) && is.null(opt$charite) && is.null(opt$curie) && is.null(opt$hartwig) && is.null(opt$oicr)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Check if arg is present and add to snvTable
snvTable=data.frame()

if (!is.null(opt$bsc)){
    bscTable = read.table(opt$bsc, header=TRUE, sep =",")
    bscTable$Node="BSC"
    bscTable=bscTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]
    snvTable=rbind(snvTable,bscTable)
}

if (!is.null(opt$charite)){
    chariteTable = read.table(opt$charite, header=TRUE, sep =",")
    chariteTable$Node="Charite"
    chariteTable=chariteTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]
    snvTable=rbind(snvTable,chariteTable)
}
if (!is.null(opt$curie)){
    curieTable = read.table(opt$curie, header=TRUE, sep =",")
    curieTable$Node="Curie"
    curieTable=curieTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]
    snvTable=rbind(snvTable,curieTable)
}

if (!is.null(opt$hartwig)){
    hartwigTable = read.table(opt$hartwig, header=TRUE, sep =",")
    hartwigTable$Node="Hartwig"
    hartwigTable=hartwigTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]
    snvTable=rbind(snvTable,hartwigTable)
}

if (!is.null(opt$oicr)){
    oicrTable = read.table(opt$oicr, header=TRUE, sep =",")
    oicrTable$Node="OICR"
    oicrTable=oicrTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]
    snvTable=rbind(snvTable,oicrTable)
}

#setwd("/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/results/CURIE/COLO829_curie_211130_nextflow")

#snvTable=read.table("SNV/COLO829_curie_211130_nextflow.stats.csv",header=TRUE,sep=",")

#snvTable$Node="Curie"
#snvTable=snvTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]

print(snvTable)

# Compute F1 + order factors
snvTable=snvTable %>%
    mutate(recall = round(recall,3)) %>%
    mutate(precision = round(precision,3)) %>%
    mutate(F1 = round(2*(precision*recall)/(precision+recall),3)) %>%
    mutate(type=factor(type,levels=c("SNVs","indels","records"))) %>%
    mutate(Node = case_when(Node == "BSC" ~ "Node 1",
                            Node == "Curie" ~ "Node 2",
                            Node == "Charite" ~ "Node 3",
                            Node == "Hartwig" ~ "Node 4",
                            Node == "OICR" ~ "Node 5")) %>%
    mutate(Node=factor(Node,levels=c("Node 1", "Node 2", "Node 3", "Node 4", "Node 5")))

datatable(snvTable)

write.table(snvTable,file = paste0(opt$outputDir ,"_snvTable.csv"))

# Transform table + order factors
tidySnv=snvTable %>%
    select(-c(total.truth,total.query,precision,recall,F1)) %>%
    pivot_longer(cols = c(tp,fp,fn), names_to = "metric", values_to = "count") %>%
    mutate(metric=factor(metric,levels=c("tp","fp","fn")))

# Plot TP FP FN
ggplot(tidySnv,aes(x= Node,y=count,fill=metric)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~type,scales="free") +
    scale_fill_brewer(palette="Set1") +
    theme_bw()+
    labs(x = "Centers",y="SNV counts") +
    theme(legend.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          axis.text=element_text(size=11),
          axis.title=element_text(size=16,face="bold"),
          plot.title = element_text(size=22),
          strip.text.x = element_text(size = 14),
          legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplotTPFPFN.png"),width=30,height=20,units='cm')

# Transform table for Precision recall F1 + order factors
tidySnv=snvTable %>%
    pivot_longer(cols = c(recall,precision,F1), names_to = "metric", values_to = "count") %>%
    mutate(metric=factor(metric,levels=c("recall","precision","F1")))

# Plot Precision Recall F1
ggplot(tidySnv,aes(x= Node,y=count*100,fill=metric)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~type,scales="free") +
    scale_fill_brewer(palette="Set1") +
    theme_bw() +
    labs(x = "Centers",y="Score (%)") +
    theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 14),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave(paste0(opt$outputDir ,"barplotPrecisionRecallF1.png"),width=30,height=20,units='cm')

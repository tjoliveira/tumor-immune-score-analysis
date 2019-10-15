#!/usr/bin/env Rscript

library(ggplot2)

#setwd('/Users/mesm/Desktop/TIDEanalysis/datasets/')

args <- commandArgs(trailingOnly = TRUE)

results <- read.csv(args[1])
#results <- read.csv('TIDE_outputs/results.csv')

CTL <- results[,c('CTL')]

correlation_exclusion_CTL <- cor.test(results[,c("MDSC_Exclusion")],CTL,
                                      method='pearson')

sp2<-ggplot(results, aes(x=results[,c("MDSC_Exclusion")],
                         y=CTL)) +geom_point(colour='darkgreen') +ggtitle("MDSCs")+theme_gray() + theme(plot.title = element_text(hjust = 0.5))


pearson <- paste('italic(R) ==',toString(correlation_exclusion_CTL$estimate))

pvalue <-paste('italic(p) ==',toString(correlation_exclusion_CTL$p.value))

sp2 <- sp2 +labs(x = "Cor(tumor,immune suppressive cell signature)")+annotate("text", x = 0.2, y = 3, label = pearson,parse = TRUE, hjust = 1)

sp2 <- sp2 +annotate("text", x = 0.2, y = 2.75, label = pvalue,parse = TRUE, hjust = 1)

sp2

filename <- args[2]

ggsave(filename, plot = last_plot(),width =11.2, height = 7.4)

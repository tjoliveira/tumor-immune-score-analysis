#!/usr/bin/env Rscript

library(ggplot2)
library(gdata)

setwd('/Users/mesm/Desktop/TIDEanalysis_26/')

args <- commandArgs(trailingOnly = TRUE)

gene_expression <- read.csv(args[1])
gene_expression <- read.csv('CRC.expression_centered_29.csv')

results<- read.csv(args[2])
results <- read.csv('29_CRC_results.csv')

#gene_names <- unlist(strsplit(args[3],','))
gene_names <- unlist(strsplit('TGFB1,CD274,SOX10',','))

gene_expression <- gene_expression[gene_expression[,c(1)] %in% gene_names,]

gene_expression_transposed <- t(gene_expression[,-1])

gene_list <- gdata::drop.levels(as.vector(gene_expression[,1]))

colnames(gene_expression_transposed) <- gene_list

gene_expression_means <- colMeans(gene_expression_transposed)

df <- as.data.frame(gene_expression_transposed)

results <- results[order(results$Patient),]

df['CTL_Classification'] <- results[,c('CTL_Classification')]
df['CTL'] <- results[,c('CTL')]

df['Dysfunction_standardized'] <- results[,c('Dysfunction_standardized')]

df['Exclusion_standardized'] <- results[,c('Exclusion_standardized')]

df['Patient'] <- rownames(df)

df <- df[order(df$CTL_Classification),]

write.table(df,file=args[4],sep=',',row.names = FALSE)

write.table(df,file='29_focused_expression.csv',sep=',',row.names = FALSE)

for(i in gene_names){
  
  
  p <- ggplot(data=df, aes_string(x='Patient',y=i,fill='CTL_Classification'))
  
  p<- p+ theme(axis.text.x = element_text(angle = 90))+geom_bar(stat = "identity")
  
  p <- p + labs(y = paste(i,' (log2(FPKM+1) - avg)'))
  
  filename <- paste(args[5],i,'expression_CTL_plot.png')
    
  ggsave(filename, plot = p,width =11.2, height = 7.4)
  
  correlation_expression_dysfunction <- cor.test(df[,c(i)],results[,c('Dysfunction_standardized')],
                                                method='pearson')
  
  p2 <-ggplot(df, aes(x=df[,c(i)],
                      y=df[,c('Dysfunction_standardized')],
                      colour = CTL))+geom_point(size=2)+geom_text(aes(label=Patient),colour='black',size=2,hjust=0, vjust=1)+scale_colour_gradient2(low = "blue",high = "red",mid = "white")
  
  pearson <- paste('italic(R) ==',toString(correlation_expression_dysfunction$estimate))
  
  pvalue <-paste('italic(p) ==',toString(correlation_expression_dysfunction$p.value))
  
  axis1 <- paste(i, '(log2(FPKM+1)-avg)')
  
  p2 <- p2 +labs(x = axis1, y= "Dysfunction")+annotate("text", x = 2, y = 2, label = pearson,parse = TRUE, hjust = 1)
  
  p2 <- p2 +annotate("text", x = 2, y = 1.75, label = pvalue,parse = TRUE,, hjust = 1)
  
  filename2 <- paste(args[5],i,'expression_dysfunction_plot.png')
  
  ggsave(filename2, plot = p2,width =11.2, height = 7.4)
  
  correlation_expression_exclusion <- cor.test(df[,c(i)],results[,c('Exclusion_standardized')],                                               method='pearson')
  
  p3 <-ggplot(df, aes(x=df[,c(i)],
                      y=df[,c('Exclusion_standardized')],
                      colour = CTL))+geom_point(size=2)+geom_text(aes(label=Patient),colour='black',size=2,hjust=0, vjust=1)+scale_colour_gradient2(low = "blue",high = "red",mid = "white")
  
  pearson <- paste('italic(R) ==',toString(correlation_expression_exclusion$estimate))
  
  pvalue <-paste('italic(p) ==',toString(correlation_expression_exclusion$p.value))
  
  p3 <- p3 +labs(x = axis1, y= "Exclusion")+annotate("text", x = 2, y = 2, label = pearson,parse = TRUE, hjust = 1)
  
  p3 <- p3 +annotate("text", x = 2, y = 1.75, label = pvalue,parse = TRUE,, hjust = 1)
  
  filename3 <- paste(args[5],i,'expression_exclusion_plot.png')
  
  ggsave(filename3, plot = p3,width =11.2, height = 7.4)
  
}


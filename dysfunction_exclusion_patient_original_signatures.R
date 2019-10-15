#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(ggpubr)
library(compare)

#setwd('/Users/mesm/Desktop/TIDE_calculation')

args <- commandArgs(trailingOnly = TRUE)

readmat = function(mat) as.matrix(read.table(mat, sep='\t', header=T, check.names=F, quote=NULL))

writemat = function(result, output)
{
  write.table(result, file=output, sep='\t', quote=F)
}

gene_dysfunction_signatures <- read.table(args[1], dec='.',sep='\t', header=T, check.names=T,quote=NULL,colClasses=c('character',"numeric"))
#gene_dysfunction_signatures <- read.table('signature', dec='.',sep='\t', header=T, check.names=T, quote=NULL,colClasses=c('character',"numeric"))

check <- gene_dysfunction_signatures[,1]
gene_dysfunction_signatures <- as.matrix(gene_dysfunction_signatures[,c('Dysfunction')])
colnames(gene_dysfunction_signatures) <- c('Dysfunction')
rownames(gene_dysfunction_signatures) <- check

gene_exclusion_signatures <- read.table(args[2], dec='.',sep='\t', header=T, check.names=F, quote=NULL,colClasses=c('character',"numeric","numeric","numeric","numeric"))
#gene_exclusion_signatures <- read.table('exclusion.signature', dec='.',sep='\t', header=T, check.names=F, quote=NULL,colClasses=c('character',"numeric","numeric","numeric","numeric"))

MDSC_exclusion_signatures <- data.frame(gene_exclusion_signatures[,c('MDSC')])
colnames(MDSC_exclusion_signatures) <- c('MDSC')
rownames(MDSC_exclusion_signatures) <- rownames(gene_exclusion_signatures)
MDSC_exclusion_signatures <- na.omit(MDSC_exclusion_signatures)

CAF_exclusion_signatures <- data.frame(gene_exclusion_signatures[,c('CAF')])
colnames(CAF_exclusion_signatures) <- c('CAF')
rownames(CAF_exclusion_signatures) <- rownames(gene_exclusion_signatures)
CAF_exclusion_signatures <- na.omit(CAF_exclusion_signatures)

TAM_exclusion_signatures <- data.frame(gene_exclusion_signatures[,c('M2')])
colnames(TAM_exclusion_signatures) <- c('TAM')
rownames(TAM_exclusion_signatures) <- rownames(gene_exclusion_signatures)
TAM_exclusion_signatures <- na.omit(TAM_exclusion_signatures)

exclusion <-read.table(args[3], dec='.',sep='\t', header=T, check.names=T, quote=NULL,colClasses=c('character',"numeric"))
#exclusion <-read.table('average.signature.exclusion', dec='.',sep='\t', header=T, check.names=T, quote=NULL,colClasses=c('character',"numeric"))
average_exclusion_signatures <- data.frame(exclusion[,c('Exclusion')])
colnames(average_exclusion_signatures) <- c('average')
rownames(average_exclusion_signatures) <- exclusion$X
#average_exclusion_signatures <- na.omit(average_exclusion_signatures)

gene_expression_profiles <- read.table(args[4], sep='\t', header=T, check.names=F, quote=NULL)
#gene_expression_profiles <- read.table('colon.expression', sep='\t', header=T, check.names=F, quote=NULL)

gene_expression_profiles_4dysfunction <- gene_expression_profiles[rownames(gene_expression_profiles) %in% row.names(gene_dysfunction_signatures),]
gene_dysfunction_signatures <- gene_dysfunction_signatures[rownames(gene_expression_profiles_4dysfunction),,drop=FALSE]

print('Number of genes in gene expression profile for T cell dysfunction score:')
paste(length(gene_expression_profiles_4dysfunction[,1]),length(gene_dysfunction_signatures[,1]))

gene_expression_profiles_4exMDSC <- gene_expression_profiles[row.names(MDSC_exclusion_signatures),]
gene_expression_profiles_4exMDSC <- na.omit(gene_expression_profiles_4exMDSC)
MDSC_exclusion_signatures <- data.frame(MDSC_exclusion_signatures[rownames(gene_expression_profiles_4exMDSC),])

print('Number of genes in gene expression profile for T cell MDSC exclusion score:')
paste(length(gene_expression_profiles_4exMDSC[,1]),length(MDSC_exclusion_signatures[,1]))

gene_expression_profiles_4exCAF <- gene_expression_profiles[row.names(CAF_exclusion_signatures),]
gene_expression_profiles_4exCAF <- na.omit(gene_expression_profiles_4exCAF)
CAF_exclusion_signatures <- data.frame(CAF_exclusion_signatures[rownames(gene_expression_profiles_4exCAF),])

print('Number of genes in gene expression profile for T cell CAF exclusion score:')
paste(length(gene_expression_profiles_4exCAF[,1]),length(CAF_exclusion_signatures[,1]))

gene_expression_profiles_4exTAM <- gene_expression_profiles[row.names(TAM_exclusion_signatures),]
gene_expression_profiles_4exTAM <- na.omit(gene_expression_profiles_4exTAM)
TAM_exclusion_signatures <- data.frame(TAM_exclusion_signatures[rownames(gene_expression_profiles_4exTAM),])

print('Number of genes in gene expression profile for T cell TAM exclusion score:')
paste(length(gene_expression_profiles_4exTAM[,1]),length(TAM_exclusion_signatures[,1]))

gene_expression_profiles_4exaverage <- gene_expression_profiles[row.names(average_exclusion_signatures),]
gene_expression_profiles_4exaverage <- na.omit(gene_expression_profiles_4exaverage)
average_exclusion_signatures <- data.frame(average_exclusion_signatures[rownames(gene_expression_profiles_4exaverage),])

print('Number of genes in gene expression profile for T cell average exclusion score:')
paste(length(gene_expression_profiles_4exaverage[,1]),length(average_exclusion_signatures[,1]))

patients <- colnames(gene_expression_profiles)

results <- data.frame(matrix(nrow=length(patients),ncol=13))

results_high_ctl <- data.frame(matrix(ncol=13))

results_low_ctl <- data.frame(matrix(ncol=13))

colnames(results_high_ctl) <- c('Patient','Dysfunction','p-value','MDSC_Exclusion','p-value',
                                'CAF_Exclusion','p-value','TAM_Exclusion','p-value','Average_Exclusion',
                                'p-value','CTL','CTL_Classification')

colnames(results_low_ctl) <- c('Patient','Dysfunction','p-value','MDSC_Exclusion','p-value',
                               'CAF_Exclusion','p-value','TAM_Exclusion','p-value','Average_Exclusion',
                               'p-value','CTL','CTL_Classification')


gene_expression_ctl <- gene_expression_profiles[c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1'),]

CTL_per_patient <- colMeans(gene_expression_ctl)

mean_CTL_samples <- mean(CTL_per_patient)

index_high <- 1

index_low <- 1

for(i in patients){
  
  dysfunction <- cor.test(gene_dysfunction_signatures[,1],gene_expression_profiles_4dysfunction[,c(i)],method='pearson')
  MDSC_exclusion <- cor.test(MDSC_exclusion_signatures[,1],gene_expression_profiles_4exMDSC[,c(i)],method='pearson')
  CAF_exclusion <- cor.test(CAF_exclusion_signatures[,1],gene_expression_profiles_4exCAF[,c(i)],method='pearson')
  TAM_exclusion <- cor.test(TAM_exclusion_signatures[,1],gene_expression_profiles_4exTAM[,c(i)],method='pearson')
  Average_exclusion <- cor.test(average_exclusion_signatures[,1],gene_expression_profiles_4exaverage[,c(i)],method='pearson')
  
  if (CTL_per_patient[c(i)] > mean_CTL_samples & all(gene_expression_ctl[c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1'),c(i)]>0)){
    
    
    results_high_ctl[index_high,] <- list(i,dysfunction$estimate,dysfunction$p.value,
                                          MDSC_exclusion$estimate,MDSC_exclusion$p.value,
                                          CAF_exclusion$estimate,CAF_exclusion$p.value,
                                          TAM_exclusion$estimate,TAM_exclusion$p.value,
                                          Average_exclusion$estimate,MDSC_exclusion$p.value,
                                          CTL_per_patient[c(i)],'High CTL')
    
    index_high <- index_high +1
    
  } else {
    
    results_low_ctl[index_low,] <- list(i,dysfunction$estimate,dysfunction$p.value,
                                        MDSC_exclusion$estimate,MDSC_exclusion$p.value,
                                        CAF_exclusion$estimate,CAF_exclusion$p.value,
                                        TAM_exclusion$estimate,TAM_exclusion$p.value,
                                        Average_exclusion$estimate,MDSC_exclusion$p.value,
                                        CTL_per_patient[c(i)],'Low CTL')
    
    index_low <- index_low +1
  }
  
  
}

results <- rbind(results_high_ctl,results_low_ctl)

#dysfunction_standarized <- results[,c('Dysfunction'),drop=FALSE]/0.143325264
#colnames(dysfunction_standarized ) <- c('Dysfunction_standardized')
#results <- cbind(results,dysfunction_standarized)

dysfunction_standarized <- results[,c('Dysfunction'),drop=FALSE]/sd(results[,c('Average_Exclusion')])
colnames(dysfunction_standarized ) <- c('Dysfunction_standardized')
results <- cbind(results,dysfunction_standarized)

#exclusion_standardized <-results[,c('Average_Exclusion'),drop=FALSE]/0.077219244
#colnames(exclusion_standardized ) <- c('Exclusion_standardized')
#results <- cbind(results,exclusion_standardized)

exclusion_standardized <-results[,c('Average_Exclusion'),drop=FALSE]/sd(results[,c('Average_Exclusion')])
colnames(exclusion_standardized ) <- c('Exclusion_standardized')
results <- cbind(results,exclusion_standardized)

results['CTL_standardized'] <- scale(results[,c('CTL')],center=TRUE,scale=TRUE)

results['TIDE_score'] <- results['CTL_Classification']

for (i in 1:length(rownames(results['TIDE_score']))){
  if (results[c(i),c('TIDE_score')]=='High CTL'){
    results[c(i),c('TIDE_score')] <- results[c(i),'Dysfunction_standardized']
  } else{
    results[c(i),c('TIDE_score')] <- results[c(i),'Exclusion_standardized']
  }
}


filename1 <- args[5]

#write.table(results,'colon_results.csv',row.names=FALSE,col.names=TRUE, sep=',')
write.table(results,filename1,row.names=FALSE,col.names=TRUE, sep=',')



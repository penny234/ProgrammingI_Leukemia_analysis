#install.packages('matrixStats')
#source("http://bioconductor.org/biocLite.R")
#biocLite()
biocLite("affy")
biocLite("limma")
biocLite("hgu133plus2.db")
biocLite("annotate")
biocLite("genefilter")
install.packages("purrr")

library(data.table)
library(affy)
library(hgu133plus2.db)
library(annotate)
library(purrr)
library(genefilter)
library(limma)
library(matrixStats)

setwd('/Users/Penny/Desktop/dataset/ALL_P1')
#this will read in all the *.CEL files in the working directory 
affydata <- ReadAffy()
# Normalization
gene_expr_rma <- rma(affydata)

### Filter out probes with features exhibiting little variation, or a consistently low signal, across samples
gene_expr_filtered <-nsFilter(gene_expr_rma, require.extrez = F, remove.dupEntrez = F)


# Convert to data frame gene expresion rma data
dt_gene <- data.frame(exprs(gene_expr_filtered$eset))

# Look for gene expresion names
probe_names <- data.table(probe_names = row.names(dt_gene), A = 1:nrow(dt_gene))

gene_names <- data.table(gene_name = getSYMBOL(probe_names[, probe_names], "hgu133plus2"), 
                         A = 1:nrow(dt_gene))

probe_gene_names <- gene_names[probe_names, on = 'A']

# Number of duplicated genes
temp <- probe_gene_names[, .N, by = gene_name] %>%
  setorder(., -N)
temp

# Merge names with gene expression data
setDT(dt_gene)
dt_gene[, A := 1:.N]
dt_gene_names <- probe_gene_names[dt_gene, on = 'A' ]

# Melt DT
long_dt_gene_name<- melt(dt_gene_names, id.vars = c('gene_name', 'A', 'probe_names'),
                          variable.name = 'patient',
                          value.name = 'gene_expression')

# Calculate min, p25, median, p75, max
long_dt_gene_name[, `:=`(MIN = quantile(gene_expression)[[1]],
                         p25 = quantile(gene_expression)[[2]],
                         MEDIAN = quantile(gene_expression)[[3]],
                         p75 = quantile(gene_expression)[[4]],
                         MAX = quantile(gene_expression)[[5]],
                         Avg = mean(gene_expression) ),
                  by = gene_name]

# 12,132 unique genes
quantiles_gene<- unique(long_dt_gene_name[, .(gene_name,
                                               MIN,
                                               p25,
                                               MEDIAN,
                                               p75,
                                               MAX)])

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(long_dt_gene_name$MEDIAN, horizontal=TRUE, ylim=c(2.5,14), xaxt="n" ,col=rgb(0.8,0.8,0,0.5) , frame=F, outcex=0.5)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(long_dt_gene_name$MEDIAN,xaxt='n',xlab='median expression level of genes',ylab='frequency',main='histogram of median gene expression level')
axis(side=1, at=seq(2.5,13.5, 0.5))
abline(v=median(long_dt_gene_name$MEDIAN),col='red')
text(4.77, y=0, '4.77',pos=4, cex=0.8) 

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(long_dt_gene_name$Avg, horizontal=TRUE, ylim=c(2.5,14), xaxt="n" ,col=rgb(0.8,0.8,0,0.5) , frame=F, outcex=0.5)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(long_dt_gene_name$Avg,xaxt='n',xlab='average expression level of genes',ylab='frequency',main='histogram of average gene expression level')
axis(side=1, at=seq(2.5,13.5, 0.5))
abline(v=median(long_dt_gene_name$Avg),col='red')
text(4.93, y=0, '4.93',pos=4, cex=0.8) 



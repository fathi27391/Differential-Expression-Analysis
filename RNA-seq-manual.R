# Installation of all necessary libraries
install.packages("BiocManager")
install.packages("matrixTests")
# Installation of genefilter from Bioconductor Manager
BiocManager::install("genefilter")

# activate of all necessary libraries
library(BiocManager)
library(matrixTests)
library(genefilter)

# loading of the data by using read.csv method and load the data as matrix
norm = as.matrix(read.csv("lusc-rsem-fpkm-tcga_paired.csv", row.names = 1)) # normal gene exp. data
## row.names argument to make the first column as our row index
tumor = as.matrix(read.csv("lusc-rsem-fpkm-tcga-t_paired.csv", row.names = 1)) # tumor gene exp. data

# exploration of the data by using dim function that return (rows number and column numbers)
dim(norm)
dim(tumor)
## Both have the same number of rows(genes) and columns (samples)

# combine both matrix (normal and tumor) into single dataset, as both have the same row number
combined = cbind(norm, tumor)
dim(combined) # exploration of new data dimensions

# explore if there is any missing data
sum(is.null(combined))

# A simple Visulatization of the combined data to see its dirtribution using a histogram
hist(combined, col="orange", main="Histogram")
# for better visulaization: we better use log2 data
# adding one (+1) to the combined data to avoid getting infinity if we log Zero Value
hist(log2(combined+1), col='cyan', main='Histogram')

# filteration of data(genes) that has low mean count
combined = combined[rowMeans(combined) > 1,]

### Calculation of fold change ###

# fold change for normal(first 50 columns) and tumor from column 51 till the end of the combined data
norm.mean = apply(log2(combined+1)[,1:50], 1 ,mean)
tumor.mean = apply(log2(combined+1)[,51:dim(combined)[2]], 1 ,mean)

# the difference between logged mean equals to the fold change in tumor sample
fold = tumor.mean - norm.mean

# visualization of fold change
hist(fold, col='cyan') # most of our data is around Zero

### Differential Expression Statistical Analysis ###

# Create a phenotype table: each row represents the phenotype [normal or tumor]
# that is corresponding to columns in the count data .. both with the same order
phenotype = as.data.frame(factor(rep(c("norm","tumor"), c(50,50))))
colnames(phenotype) = "group"

# make hypothesis testing using t-test for each row(gene)
# we use rowttest from gene filter library
t_table = rowttests(combined, phenotype$group)

# use p.adjust method for fixing of p-value that we achieved from t-table
# use the FDR method[False Discovery Rate] in case of gene expression data
p.adj = p.adjust(t_table$p.value, 'fdr')

# make our results in a new dataset
result = as.data.frame(cbind(fold, p.adj))

# filter out DEG based on p.value < 0.05 significant fold changes
# and also based on fold change > 2 either in postive(upregulated) or negative(downregulated) direction
result.deg = result[result$p.adj < 0.05 & abs(result$fold) > 2,]

## saving our data for further analysis
write.csv(as.matrix(result.deg), file = 'result_deg.csv', quote=F, row.names = TRUE)

#### using non-parametric test if the data isn't normally distributed ####
w_table = row_wilcoxon_twosample(combined[,1:50], combined[,51:dim(combined)[2]])
# use fdr method for calculating adjusted p-value
p.adj_2 = p.adjust(w_table$pvalue, 'fdr')
# combine fold change with p.adj_2
result_2 = as.data.frame(cbind(fold, p.adj_2))
# DEG from wilcoxon statistical test
deg_2 = result_2[abs(result_2$fold) > 2 & result_2$p.adj_2 < 0.05,]

hist(deg_2$fold)
hist(result.deg$fold)

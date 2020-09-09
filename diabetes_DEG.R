# installation of DESeq2 from bio-conductor
BiocManager::install("DESeq2")
# activate DESeq2 Library
library(DESeq2)

# loading of the data as matrix with gene name as row index
data = as.matrix(read.csv('Diabetes_infection_count.csv', row.names = 1))
# loading phenotype table , each row represent the sample in the count matrix data
phenotype = read.csv('Diabetes_infection_pheno.csv', row.names = 1)

table(phenotype$CL4) # similar to value_counts in python: count each unique item in the column

# overview of the data
dim(data)
sum(is.null(data))

# Visulaization of data distribution
hist(data, col='cyan', main='Histogram')

# use log2 transformation for better visulization 
# add 1 to data to avoid log2 of zero values
logged = log2(data +1)
hist(logged, col='cyan', main='Histogram') # there is a peak of data at zero

# to use DESeq2: it is a must for both the column names in our data count and the
# row name in the phenotype table are of the same order
phenotype = phenotype[colnames(data),]

# second thing before running DESeq 2 is that it requires the count data to be integer
# 1) save gene names in a variable as row index will be deleted after converting data type
genes = row.names(data)
# 2) convert data into integer
data = apply(data, 2, as.integer)
# 3) restore original row names of the data
row.names(data) = genes

##### Use DESeq2 to determine DEG [Differentially Expressed Gene] #### 

# A) Specify how many condition you want to compare according to phenotype data
table(phenotype$CL4)
cond1 = "Healthy"
cond2 = "Infection"

# B) Create a DESeq Dataset object
dds = DESeqDataSetFromMatrix(countData = data, colData = phenotype, design = ~ CL4)
  # use DESeq function of DESeq Dataset object to perform DESeq workflow
dds_2 = DESeq(dds)

# C) get the DESeq result by using DESeq result method
res = results(dds_2, contrast = c("CL4", cond1, cond2))
  # convert result DESeq object into data frame & remove null values using complete.cases to filter out the data
res = as.data.frame(res[complete.cases(res),])

# D) from results data frame filter out the DEG using subsets only p-value <0.05 and fold change
deseq.deg = res[res$padj < 0.05 & abs(res$log2FoldChange) > log2(2),]

# E) save the DEG into csv file
write.csv(as.matrix(deseq.deg),file ='DEG-diabetes.csv', quote = F, row.names = T)

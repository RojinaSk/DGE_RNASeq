
setwd("your_path")

#install.packages("magrittr")
library(magrittr)
#install.packages("statmod") 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)



#----------------------------
# 1. Load and organize data.
#----------------------------

phenotype_info_file <- "targets.txt"
raw_counts_file <- "GSE60450_Lactation-GenewiseCounts.txt.gz"

#12 samples and 6 categories
targets <- read.delim(phenotype_info_file, stringsAsFactors=FALSE) #delimited files, defaulting to the TAB character for the delimiter
  
View(targets)

group <- paste(targets$CellType, targets$Status, sep = ".")
group <- as.factor(group)

#group <- targets %$% paste(CellType, Status, sep = ".") %>% factor()

group

#Length of a gene is the total number of bases in exons and UTRs for that gene.
GenewiseCounts <- read.delim(raw_counts_file , row.names="EntrezGeneID")
head(GenewiseCounts,1)

colnames(GenewiseCounts)=substring(colnames(GenewiseCounts), 1, 7)

#colnames(GenewiseCounts) %<>% substring(.,1,7)

head(GenewiseCounts,1)



#------------------------
# 2. MA plots
#------------------------

two_samples <- GenewiseCounts[, c(2, 3)] + 1 #2, 3: Replicate samples 12
two_samples <- log2(two_samples)
head(two_samples)

#two_samples <- GenewiseCounts[, c(2, 3)] %>% 
#  add(., 1) %>%
#  log2()

plotData <- data.frame(M = two_samples[, 1] - two_samples[, 2],
                       A = (two_samples[, 1] + two_samples[, 2])/2)

#log fold-change versus mean expression between two treatments
ggplot(plotData, aes(x = A, y = M)) +
  geom_point() +
  geom_smooth()



#------------------------
# 3. Create DGElist object and retrieve gene symbols.
#------------------------

y <- DGEList(counts = GenewiseCounts[,-1], 
             group = group,
             genes = GenewiseCounts[,1,drop=FALSE]) #drop=FALSE keeps the same dimension

y$genes$Symbol <- mapIds(org.Mm.eg.db,
                         keys = rownames(y),
                         keytype="ENTREZID", 
                         column="SYMBOL")

View(y)



#------------------------
# 4. Independent filtering.
#------------------------

#Filter genes whose symbols are not found.
keep <- !is.na(y$genes$Symbol) 

head(keep)
table(keep)

y <- y[keep, ]



# It is not possible to make reliable inference for genes if gene count is too small

# https://f1000research.com/articles/5-1438 :
# As a rule of thumb, we require that a gene have a count of at least 10~15 in at least some libraries
# filter on count-per-million (CPM) values so as to avoid favoring genes that are expressed in larger 
# libraries over those expressed in smaller libraries


minimum_counts_reqd <- 10
cutoff <- median(y$samples$lib.size)
cutoff = cutoff/10^6
cutoff = minimum_counts_reqd/cutoff
cutoff = round(cutoff, 1)
cutoff

#cutoff <- y$samples$lib.size %>% 
#  median() %>% 
#  divide_by(10^6) %>% 
#  divide_by(minimum_counts_reqd, .) %>% #check
#  round(1)

#... or simply set cutoff to 0.5

keep <- cpm(y, prior.count=1)
keep = keep > cutoff
keep = rowSums(keep)
keep = keep >= 2

head(keep)

#keep <- cpm(y) %>% #counts per million
#  is_greater_than( cutoff) %>%
#  rowSums() %>%
#  is_weakly_greater_than( 2)

y <- y[keep, , keep.lib.sizes=FALSE] # keep.lib.sizes=FALSE: the library sizes to be recomputed after the filtering



#-------------------------
# 5. Normalizing the counts.
#-------------------------

# normalize for RNA composition by a set of scaling factors 
# that minimize the log-fold changes between the samples for most genes
head(y$samples)

y <- calcNormFactors(y)
head(y$samples)

# Intuitively, we should expect similar adjustments for similar samples.
# Plot the normalization factors by sample.
ggplot(cbind(y$samples, replicate = factor(1:2)), 
       aes(x = group, y = norm.factors, fill = replicate)) +
  geom_col(position = position_dodge())

#"A normalization factor below one indicates that a small number of high count genes 
#...are monopolizing the sequencing, causing the counts for other genes to be lower 
#...than would be usual given the library size." Chen et al., 2016
#Looks like L.lactating samples contain a number of very highly upregulated genes.





#------------------------
# 6. Exploratory visualizations
#------------------------

#MDS plot
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group], gene.selection = "common")
legend("top", legend=levels(group) %>% substr(1, 3), 
       pch=pch, col=colors, ncol=2, cex = 0.5)
# leading log2FC: the root-mean-square average of the top largest log2-fold-changes between those two samples


#PCA plot
cpm <- cpm(y, log = TRUE)
rv <- apply(cpm,1,var) 

#Select genes with highest variance.
keep <- order(rv, decreasing = TRUE)[1:500]
selected <- cpm[keep, ] %>% t()

#Transpose is needed to ensure that each row is a vector.
pca <- prcomp(selected, scale=T, center = T)

stddev <- pca$sdev
pc1_var <- round(100*stddev[1]^2/sum(stddev^2))
pc2_var <- round(100*stddev[2]^2/sum(stddev^2))

PlotData <- data.frame(cbind(PC1 = pca$x[,1], PC2 = pca$x[,2]))
PlotData <- targets[, c("CellType", "Status")] %>%
  cbind(PlotData, .)

head(PlotData)

ggplot(PlotData, aes(x=PC1, y=PC2, color=CellType, shape=Status)) + 
  geom_point(size=4.5) +  
  xlab(paste("PC1:", pc1_var, "% variance")) + 
  ylab(paste("PC2:", pc2_var, "% variance"))



#--------------------
# 7. Fitting the model.
#--------------------

design <- model.matrix(~ 0 + group) %>%   
  set_colnames(levels(group))

design

y <- estimateDisp(y, design, robust=TRUE) 

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE) 
head(fit$coefficients)

plotQLDisp(fit) 



#--------------------
# 8. Hypothesis testing
#--------------------

#Example 1: Hypothesis: The mean expression of B.lactating is different from the mean expression of B.pregnant.
B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
res <- glmQLFTest(fit, contrast=B.LvsP)
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

#Example 2: Hypothesis: The mean expression of B.lactating is different from the mean expression of B.pregnant and 
#  the log2-fold-changes are significantly greater than 1.5
B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
res <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

#Example 3: Hypothesis: The change in expression between lactating and pregnant mice is the same for basal cells as it is for luminal cells
con <- makeContrasts(
  (L.lactating-L.pregnant)-(B.lactating-B.pregnant),
  levels=design)
res <- glmQLFTest(fit, contrast=con)
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

#Example 4: Hypothesis : Comparing virgin, pregnant, and lactating mouse in the luminal population
con <- makeContrasts(
  L.PvsL = L.pregnant - L.lactating,
  L.VvsL = L.virgin - L.lactating,
  L.VvsP = L.virgin - L.pregnant, levels=design)
res <- glmQLFTest(fit, contrast=con)
topTags(res)
is.de <- decideTestsDGE(res)
summary(is.de)


#other contrast examples: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/

#----------------------------
# 9. Save results in a table
#----------------------------
result <- as.data.frame( topTags(res, n = nrow(y$counts)) )
write.table(result, 
            file = "DE_result.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)





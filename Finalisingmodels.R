
# Adonis test
library(vegan)


### test out adonis 2 
library(devtools)

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
adonis.test2 <- adonis2(bray_ps_nemadults ~  Season + Sex + Id,
                            by="margin", 
                            data = metadf.bxadults)  

# reg adonis 
adonis.testadultsmulti<- adonis(bray_ps_nemadults ~ Id + Season + Sex, data = metadf.bxadults)  


R2 <- adonis.test.soay$aov.tab$R2 
Terms <- adonis.test.soay$aov %>% row.names()

# adonois2 

R2 <- adonis.test2$R2 
Terms <- adonis.test2%>% row.names()

adonisDFadults <- data.frame(R2, Terms) %>% 
  filter(Terms!="Total") 

adonisDFadults %>% 
  ggplot(aes(x=1, y=R2, fill=Terms)) + 
  geom_bar(stat="identity") + 
  labs(x='Variance Explained') + 
  scale_fill_brewer(palette = "Spectral") + 
  theme_bw(base_size=14) + 
  theme(axis.text.x = element_blank())


adonis.test.soay.pair <- pairwise.adonis2(bray_ps.bx.soay ~ Season, 
                                          by="margin", 
                                          data = metadf.bx.soay)

## this is not working but should if the data can be formatted - deseq2 good option as alternate 
coef <- coefficients(adonis.testadultssex)["Sex",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]     #error in evaluating the argument 'x' in selecting a method for function 'rev': non-numeric argument to mathematical function
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

##deseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

ds2 <- phyloseq_to_deseq2(ps_nem_filt_repeatsrem, ~ Sex) 
diagdds = DESeq(ds2, test="Wald", fitType="parametric") #alternative to Run DESeq analysis, comes up with same error. Sourced from https://joey711.github.io/phyloseq-extensions/DESeq2.html
#calculate geometric means prior to estimate size factors
#following is sourced here https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds2), 1, gm_mean)
ds2 = estimateSizeFactors(ds2, geoMeans = geoMeans) #comes up with Error in .local(object, ..., value) : all(!is.na(value)) is not TRUE
ds2 = DESeq(ds2, fitType="local")
#comes up with error, in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

 
#another alternative found here https://support.bioconductor.org/p/62246/#62250
ds2 <- ds2[ rowSums(counts(ds2)) > 5, ]
cts <- counts(ds2)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
ds2 <- estimateSizeFactors(ds2, geoMeans=geoMeans) #Error in .local(object, ..., value) : all(!is.na(value)) is not TRUE


# Run DESeq2 analysis (all taxa at once!)
dds <- DESeq(ds2) # comes up with error Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

# Investigate results
deseq.results <- as.data.frame(results(dds))
deseq.results$taxon <- rownames(results(dds))

# Sort (arrange) by pvalue and effect size
library(knitr)
deseq.results <- deseq.results %>%
  arrange(pvalue, log2FoldChange)

# Print the result table
# Let us only show significant hits
knitr::kable(deseq.results %>%
               filter(pvalue < 0.05 & log2FoldChange > 1.5),
             digits = 5)

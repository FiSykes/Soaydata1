rm(list=ls()) #clearing environment

#libraries needed

library(tidyverse) #important 
library(phyloseq) #important 
library(microbiome)  # important 
library(knitr) #for nice tables
library(ggpubr) #optional for plots
library(reshape2)
library(RColorBrewer)
library(microbiomeutilities) # important 
library(viridis)
library(tibble)
library(metagMisc) # important 
library(ggforce) 
library(pals)



#taking a look

view(otu_table(phyloseqNemabiomeSoaySpecies))
tax_table(phyloseqNemabiomeSoaySpecies)
view(sample_data(phyloseqNemabiomeSoaySpecies))

#prevalence 

phyloseq_prevalence_plot(phyloseqNemabiomeSoaySpecies, taxcolor = "SpeciesFull", facet=FALSE, point_alpha = 0.5, 
                         prev.trh = 0.05) +
  theme_bw(base_size = 12) + 
  geom_label(aes(label=SpeciesFull,))  #just for viz purposes the labels a bit clunky 

#Stuff post meeting with Amy

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame %>% pull(SpeciesFull) -> SpeciesNames # pulls the actual species names


taxa_names(phyloseqNemabiomeSoaySpecies) <- SpeciesNames  ## saves the object with new names

taxa_sums(phyloseqNemabiomeSoaySpecies)


ps_nem_filt_noncomp <- filter_taxa(phyloseqNemabiomeSoaySpecies,
                                function(x) {sum(x>1000) > 5}, prune=TRUE) #filtering equation, try as compositional too

taxa_sums(ps_nem_filt_noncomp)

#conversion to compositional pre filter

ps_nem_comp <- microbiome::transform(phyloseqNemabiomeSoaySpecies, "compositional") #back to sum to one 
otu_table(ps_nem_comp) %>% rowSums()

saveRDS(ps_nem_comp, "phyloseqSpeciesProportional.rds")
view(otu_table(ps_nem_comp))  #IT WORKED!!!! YAY!!

#Removing individual problem samples - samples predominated by species not on St Kilda
#1 is extreme 50% O.venulosum sample - IS A REPEAT SAMPLE 
ps_nem_samplesremoved <- prune_samples(!(sample_names(ps_nem_comp) =="10349_070319" ), ps_nem_comp) 
view(otu_table(ps_nem_samplesremoved))

#2 is problematic H.contortus sample - lamb, NON REPEAT SAMPLE 11932_180719
ps_nem_samplesremoved2 <- prune_samples(!(sample_names(ps_nem_samplesremoved) =="11932_180719" ), ps_nem_samplesremoved) 
view(otu_table(ps_nem_samplesremoved2))

#3 is  cooperia oncophora sample at over 16000 reads (most predominant species in this sample). Lamb, NOT A REPEAT 12035_270719
ps_nem_samprem<- prune_samples(!(sample_names(ps_nem_samplesremoved2) =="12035_270719" ), ps_nem_samplesremoved2) 
view(otu_table(ps_nem_samprem))


#now to filter! 
ps_nem_filt <- filter_taxa(ps_nem_samprem,
                                   function(x) {sum(x>0.02) > 3}, prune=TRUE) #filtering equation. Remove taxa not seen above 2% in a sample at least 3 times
#this filtering value is the lowest to only include the five most abundant species.
taxa_sums(ps_nem_filt)
view(otu_table(ps_nem_filt))
saveRDS(ps_nem_filt, "phyloseqSpeciesProportionalFiltered.rds")



#Further analysis taken from TraditionalAnalysis tutorial
#Diversity Measures  -----------------------------------------------------
  
  
  ### ALPHA 
  #replace 'pseq' with my dataset
tab <- alpha(phyloseqNemabiomeSoaySpecies, index = "all")  # all alpha diversity metrics, doesn't work with compositional data, can't get to work at all?
tab %>% head()

tab <- alpha(ps_nem_filt, index = "diversity_inverse_simpson")  
tab %>% head()

diet_even <- evenness(ps_nem_filt_noncomp, index = "all")  # evenness metrics, this works completely fine
diet_even %>% head() 

metaps_nem <-meta(ps_nem_filt_noncomp)
metaps_nem$simpson_inverse <- tab$diversity_inverse_simpson # or gini_simpson, or shannon 
metaps_nem$simpson_gini <- tab$diversity_gini_simpson #these last 2 won't work unless I do the tab thing and get it to show all
metaps_nem$shannon <- tab$diversity_shannon

metaps_nem <- as.data.frame(metaps_nem)

## download 
library(ggforce)
library(ggpubr)
metaps_nem %>% 
  ggplot(aes(x=sex, y=shannon, colour=sex)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  theme_bw(base_size = 16) -> DiversityBySex 

metaps_nem %>% 
  ggplot(aes(x=AgeClass, y=simpson_inverse, colour=AgeClass)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova") #p= 0.00014. Sig difference in simpson_inverse diversity between AgeClasses. Adults more diverse

metaps_nem %>% 
  ggplot(aes(x=Season, y=simpson_inverse, colour=Season)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova") # p = 3e -10. Sig diff in index between seasons. March is most diverse with this dataset

metaps_nem$Sex<-as.factor(metaps_nem$Sex) 
levels(metaps_nem$Sex)

factor_readlabel

metaps_nem %>% 
  ggplot(aes(x=as.factor(Sex), y=simpson_inverse, colour=as.factor(Sex))) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova")

## can probably ignore this !!! 
# conditioning on variable type 
# create a list of pairwise comaprisons
smtypeDiet <- levels(as.factor(metaDiet$sex)) # get the variables for chosen variable sex 

smtypeDiet2 <- levels(as.factor(metaDiet$bmi_group)) # get the variables for chosen variable bmi

# make a pairwise list that we want to compare.
smtype.pairs <- combn(seq_along(smtypeDiet), 2, simplify = FALSE, FUN = function(i)smtypeDiet[i])
smtype.pairs2 <- combn(seq_along(smtypeDiet2), 2, simplify = FALSE, FUN = function(i)smtypeDiet2[i])

print(smtype.pairs) # check your comparisons 

DiversityPlot <- ggviolin(metaDiet, x = "sex", y = "simpson_inverse",
                          add = "boxplot", fill = "sex", 
                          palette = c("#FF6EB4", "#63B8FF"),
                          legend = "right") +stat_compare_means(method = "t.test")

DiversityPlot2 <- ggviolin(metaDiet, x = "bmi_group", y = "simpson_inverse",
                           add = "boxplot", fill = "bmi_group", 
                           palette = c("#FF6EB4", "#63B8FF", "purple"),
                           legend = "right") +stat_compare_means(method = "anova")



### BETA 

## Beta Diversity 
#Need to use compostional data set

bx.ord_pcoa_bray_ps_nem <- ordinate(ps_nem_filt, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nem) + theme_bw()

# Axis 1 and 2 are of interest - pcoa 
beta_ps_nem <- plot_ordination(ps_nem_filt, 
                             bx.ord_pcoa_bray_ps_nem, 
                             color="AgeClass",  #factor of interest 
                             #label = "subject",
) + 
  geom_point(aes(shape = AgeClass), size= 4) +  #factor of interest 
  #geom_label(aes(label=sampleId), size=1) + 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)

plot(beta_ps_nem)
# Axis 1 and 2 are of interest.
beta.ps_nem <- plot_ordination(ps_nem_filt, 
                             bx.ord_pcoa_bray_ps_nem, 
                             color="Season", #factor of interest 
                             #label = "sample_id",
) + 
  geom_point(aes(shape = Season), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)


beta.ps_nem_sex <- plot_ordination(ps_nem_filt, #need to change sex from numbers to labels to get this to work
                               bx.ord_pcoa_bray_ps_nem, 
                               color= as.factor("Sex"), #factor of interest 
                               #label = "sample_id",
) + 
  geom_point(aes(shape =as.factor(Sex), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse(), #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)
  
plot(beta.ps_nem_sex),

x=as.factor(Sex), y=simpson_inverse, colour=as.factor(Sex))

## PERMANOVA
metadf.bx <- data.frame(sample_data(pseq_1))
bray_ps.bxnDiet <- phyloseq::distance(physeq = diet.rel, method = "bray")

set.seed(995)
# Adonis test
library(vegan)
adonis.testDiet <- adonis(bray_ps.bxnDiet ~ bmi_group, data = metadf.bx)

adonis.testDiet  ## matches mcmcglmm output very well again

# Note the assumption of similar multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be potentially explained by that.
dist <- vegdist(t(abundances(diet.rel)))
anova(betadisper(dist, metadf.bx$bmi_group)) # signif different spreads among groups 


## dissimilarity stuff #can't currently get this to work
b.Lean <- as.data.frame(divergence(subset_samples(pseq_1, bmi_group == "lean"))) %>% 
  rename(div=1) %>%  mutate(bmi="lean")
b.Over <- as.data.frame(divergence(subset_samples(pseq_1, bmi_group == "overweight"))) %>% 
  rename(div=1) %>% mutate(bmi="overewight")
b.Obese <-as.data.frame(divergence(subset_samples(pseq_1, bmi_group == "obese"))) %>% 
  rename(div=1) %>% mutate(bmi="obese")


div_df <- bind_rows(b.Lean, b.Over, b.Obese)

ggpubr::ggboxplot(div_df, "bmi", "div", 
                  ylab = "Divergence", 
                  xlab = "Season", 
                  add = "jitter",
                  fill = "bmi",
                  palette = c("#a6cee3", "#b2df8a", "purple")) + stat_compare_means(method = "anova")




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
                                   function(x) {sum(x>0.02) > 3}, prune=TRUE) #filtering equation. Remove taxa not seen above 2% in a sample at least 3 times across the dataset.
#this filtering value is the lowest to only include the five most abundant species.
taxa_sums(ps_nem_filt)
view(otu_table(ps_nem_filt))
saveRDS(ps_nem_filt, "phyloseqSpeciesProportionalFiltered.rds")

#making a repeats only data frame 
reps <- sample_data(ps_nem_filt) %>% as_tibble %>% filter(Repeat =="Y") %>% pull(SampleId)
ps_nem_filt_repeats <- prune_samples((sample_names(ps_nem_filt)%in% reps), ps_nem_filt)
view(sample_data(ps_nem_filt_repeats))
saveRDS(ps_nem_filt_repeats, "RepeatsonlyDF.rds")

#filtering out repeat samples for main datset
view(sample_data(ps_nem_filt))

firstSamples <-
  meta(ps_nem_filt) %>% group_by(IdSeason) %>% mutate(Jandate=as.Date(DateCollect) %>% format(.,"%j")) %>% #converts to numeric value 1-365
  mutate(SampleRank = order(order(Jandate))) %>% #this makes a rank of each sample within season
  filter(SampleRank==1) %>% pull(SampleId)

ps_nem_filt_repeatsrem <- prune_samples((sample_names(ps_nem_filt) %in% firstSamples), ps_nem_filt)
saveRDS(ps_nem_filt_repeatsrem, "FilteredwithRepeatsRem.rds")
view(sample_data(ps_nem_filt_repeatsrem))
view(sample_data(ps_nem_perf))

ps_nem_perf <- sample_data(ps_nem_filt_repeatsrem)$SexCorr <- ifelse(sample_data(ps_nem_filt_repeatsrem)$Sex ==1, "Female", "Male")
view(sample_data(ps_nem_perf))
#Further analysis taken from TraditionalAnalysis tutorial
#Diversity Measures  -----------------------------------------------------
  
  
  ### ALPHA 
  #replace 'pseq' with my dataset
tab <- alpha(phyloseqNemabiomeSoaySpecies, index = "all")  # all alpha diversity metrics, doesn't work with compositional data, can't get to work at all?
tab %>% head()

tab <- alpha(ps_nem_filt_lambs, index = c("diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon"))  
tab %>% head()

diet_even <- evenness(ps_nem_perf, index = "all")  # evenness metrics, this works completely fine
diet_even %>% head() 

metaps_nem <-meta(ps_nem_perf)
metaps_nem$simpson_inverse <- tab$diversity_inverse_simpson # or gini_simpson, or shannon 
metaps_nem$simpson_gini <- tab$diversity_gini_simpson #these last 2 won't work unless I do the tab thing and get it to show all
metaps_nem$shannon <- tab$diversity_shannon

metaps_nem <- as.data.frame(metaps_nem)
metaps_nem %>% mutate(Sex=case_when(Sex== 1 ~ "Female" , Sex== 2 ~ "Male")) #comes up with error, 'non numeric argument to binary operator?'
#did this just work :0

metaps_nem_adults <- meta(ps_nem_filt_adults)
metaps_nem_adults <- as.data.frame(metaps_nem_adults)
metaps_nem_adults$simpson_inverse <- tab$diversity_inverse_simpson
metaps_nem_adults$shannon <- tab$diversity_shannon

metaps_nem_lambs <- meta(ps_nem_filt_lambs)
metaps_nem_lambs <- as.data.frame(metaps_nem_lambs)
metaps_nem_lambs$simpson_inverse <- tab$diversity_inverse_simpson
metaps_nem_lambs$shannon <- tab$diversity_shannon
## download 
library(ggforce)
library(ggpubr)

#shannon plots

metaps_nem_lambs%>% 
  ggplot(aes(x=SexCorr, y=shannon, colour=SexCorr)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  labs(x="Sex", y="Shannon diversity index") +
  theme_bw(base_size = 16)+ stat_compare_means(method = "anova") -> DiversityBySexLambs 

metaps_nem_adults%>% 
  ggplot(aes(x=SexCorr, y=shannon, colour=SexCorr)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  labs(x="Sex", y="Shannon diversity index") +
  theme_bw(base_size = 16)+ stat_compare_means(method = "anova") -> DiversityBySexAdults

metaps_nem_adults%>% 
  ggplot(aes(x=Season, y=shannon, colour=Season)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) + 
  labs(x="Season", y="Shannon diversity index") +
  theme_bw(base_size = 16)  

metaps_nem_sexnonbin%>% 
  ggplot(aes(x=AgeClass, y=shannon, colour=AgeClass)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  scale_color_discrete ("Age Class") +
  labs(x="Age Class", y="Shannon diversity index") +
  theme_bw(base_size = 16) +stat_compare_means(method = "anova") #p = 0.00061

metaps_nem %>% 
  ggplot(aes(x=as.factor(Age), y=shannon, colour=as.factor(Age))) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  labs(x = "Age (years)", y= "Shannon'Diversity Index") +
  scale_color_discrete ("Age (years)") +
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova")

#inverse_simpson plots

metaps_nem %>% 
  ggplot(aes(x=as.factor(Age), y=simpson_inverse, colour=as.factor(Age))) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  labs(x = "Age (years)", y= "Simpson's Inverse Diveristy Index") +
  scale_color_discrete ("Age (years)") +
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova") #p= 0.0047. Doesn't seem to be much of a change after 3 years old. Even at 10 doesn't seem to be much difference. 
#Highest diversity seems to be at 2! Though lower number of data points. Try this out as continunous plot
plot_age <- ggplot(metaps_nem_sexnonbin, aes(x = Age, y = simpson_inverse, colour = Age)) +
  geom_point(size = 3, shape = 1) +
  labs(x = "Age (years)", y= "Simpson's Inverse Diveristy Index") +
  geom_smooth(method = "lm", linetype = 2, se = TRUE) #sugests diversity increases with Age. No sample data of sheep at 1 year old.

metaps_nem %>% 
  ggplot(aes(x=AgeClass, y=simpson_inverse, colour=AgeClass)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  labs(x = "AgeClass", y= "Simpson's Inverse Diveristy Index")+
  scale_color_discrete ("Age Class") +
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova") #p= 3e-04. Sig difference in simpson_inverse diversity between AgeClasses. Adults more diverse

metaps_nem%>%mutate(Season = fct_relevel(Season, "Mar-19", "May-19", "Jul-19"))

metaps_nem_adults %>% 
  ggplot(aes(x=Season, y=simpson_inverse, colour=Season)) + 
  geom_sina() + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +  #Needed to change to Mar-19 instead of March
  geom_boxplot(fill="transparent")  + 
  labs(x = "Season (2019)", y= "Simpson's Inverse Diveristy Index") +
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova")
      # p = 9.5e-08. Sig diff in index between seasons. March is most diverse with this dataset. Post changing season order cannot generate anova?


metaps_nem$Sex<-as.factor(metaps_nem$Sex) 
levels(metaps_nem$Sex)
metaps_nem_sexnonbin <- metaps_nem%>%mutate(Sex=case_when(Sex==1 ~ 'Female', Sex == 2 ~ 'Male'))

factor_readlabel

metaps_nem_lambs %>% 
  ggplot(aes(x=SexCorr, y=simpson_inverse, colour=SexCorr)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  theme_bw(base_size = 16)  +stat_compare_means(method = "anova") #Females are more diverse than males p= 2.1e-0.5

metaps_nem_adults%>% 
  ggplot(aes(x=SexCorr, y=simpson_inverse, colour=SexCorr)) + 
  geom_sina() + 
  geom_boxplot(fill="transparent")  + 
  labs(x="Sex", y="Simpson's Diversity Index") +
  theme_bw(base_size = 16)+ stat_compare_means(method = "anova") 

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

bx.ord_pcoa_bray_ps_nem <- ordinate(ps_nem_filt_repeatsrem, "PCoA", "bray")
bx.ord_pcoa_bray_ps_nemlambs <- ordinate(ps_nem_filt_lambs, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nem) + theme_bw()
plot_sree(bx.ord_pcoa_bray_ps_nemlambs) + theme_bw()


# Axis 1 and 2 are of interest - pcoa 
beta_ps_nemlambssex <- plot_ordination(ps_nem_filt_lambs, 
                               bx.ord_pcoa_bray_ps_nemlambs, 
                               color="SexCorr",  #factor of interest 
                               #label = "subject",
) + 
  geom_point(aes(shape = SexCorr), size= 4) +  #factor of interest 
  #geom_label(aes(label=sampleId), size=1) + 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)


beta_ps_nem <- plot_ordination(ps_nem_filt_repeatsrem, 
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
beta.ps_nem_season <- plot_ordination(ps_nem_filt_repeatsrem, 
                             bx.ord_pcoa_bray_ps_nem, 
                             color="Season", #factor of interest 
                             #label = "sample_id",
) + 
  geom_point(aes(shape = Season), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)

ps_nem_sexnonbin <- as.phyloseq(metaps_nem_sexnonbin)

ps_nem_sexnonbin <- ps_nem_filt_repeatsrem%>%mutate(Sex=case_when(Sex==1 ~ 'Female', Sex == 2 ~ 'Male')) #can't be applied to phyloseq object


beta.ps_nem_sexadults <- plot_ordination(ps_nem_filt_adults, #need to change sex from numbers to labels within a phyloseq object. Won't work with metadata.
                               bx.ord_pcoa_bray_ps_nem, 
                               color= "SexCorr", #factor of interest 
                               #label = "sample_id",
) + 
  geom_point(aes(shape =SexCorr), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120) 


plot(beta.ps_nem_sexadults)


## PERMANOVA. Used to detect whether different groups (e.g male vs female) have a serious effect on overall worm composition
metadf.bx <- data.frame(sample_data(ps_nem_filt_repeatsrem))
bray_ps_nem <- phyloseq::distance(physeq = ps_nem_filt_repeatsrem, method = "bray")

metadf.bxadults <- data.frame(sample_data(ps_nem_filt_adults))
bray_ps_nemadults <- phyloseq::distance(physeq = ps_nem_filt_adults, method = "bray")

metadf.bxlambs <- data.frame(sample_data(ps_nem_filt_lambs))
bray_ps_nemlambs <- phyloseq::distance(physeq = ps_nem_filt_lambs, method = "bray")

set.seed(995)
# Adonis test
library(vegan)
adonis.testlambssex <- adonis(bray_ps_nemlambs ~ SexCorr, data = metadf.bxlambs)

adonis.testlambssex



adonis.testadultssex <- adonis(bray_ps_nemadults ~ SexCorr, data = metadf.bxadults)

adonis.testadultssex


adonis.test <- adonis(bray_ps_nem ~ Season, data = metadf.bx)

adonis.test  ## so for season, there is not much of a correlation between groups r2= 0.15593? Though a sig diff between groups p = 0.001?
#does this result suggest a difference between just 2 of the groups, or that every group is significantly different from each other?
#can try for AgeClass
adonis.test.Age <- adonis(bray_ps_nem ~ AgeClass, data = metadf.bx)

adonis.test.Age  #for ageClass, slightly greater correlation between groups r2= 0.29584? Also a sig difference between the groups p = 0.001

#can try for sex when things sorted out?
adonis.test.Sex <- adonis(bray_ps_nem ~ SexCorr, data = metadf.bx)

adonis.test.Sex #very little correlation between groups, r2 = 0.01666, p value = 0.002

#the following looks to check that the varianc of homogeneity assumptions hold, to ensure relability of the results of the permanova
dist <- vegdist(t(abundances(ps_nem_filt_repeatsrem)))
anova(betadisper(dist, metadf.bx$Season)) # signif different distance spreads among season groups, explained by lambs being present? p = 1.508e-08
anova(betadisper(dist, metadf.bx$SexCorr))  # gives p value of 0.009281. signif different distance spreads between sexes
anova(betadisper(dist, metadf.bx$AgeClass)) #  AgeClass p = 0.0001214

distadults <- vegdist(t(abundances(ps_nem_filt_adults)))
anova(betadisper(distadults, metadf.bxadults$SexCorr)) 

distlambs <- vegdist(t(abundances(ps_nem_filt_lambs)))
anova(betadisper(distlambs, metadf.bxlambs$SexCorr)) 

#microbiome site shows determining the top factors after this step (showing coefficeitn for the top taxa separating groups
#is that needed here?
coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

## dissimilarity stuff #can't currently get this to work
b.Male <- as.data.frame(divergence(subset_samples(ps_nem_filt_repeatsrem, Season == "Mar-19"))) %>% 
  rename(div=1) %>%  mutate(Season ="Mar-19")
b.May <- as.data.frame(divergence(subset_samples(ps_nem_filt_repeatsrem, Season == "May-19"))) %>% 
  rename(div=1) %>% mutate(Season ="May-19")
b.July <-as.data.frame(divergence(subset_samples(ps_nem_filt_repeatsrem, bmi_group == "Jul-19"))) %>% 
  rename(div=1) %>% mutate(Season = "Jul-19")


div_df <- bind_rows(b.March, b.May, b.July)

ggpubr::ggboxplot(div_df, "Season", "div", 
                  ylab = "Divergence", 
                  xlab = "Season", 
                  add = "jitter",
                  fill = "Season",
                  palette = c("#a6cee3", "#b2df8a", "purple")) + stat_compare_means(method = "anova")










#Some abundance plots!
#taken from phyloseq site
title = "Abundance plot"
plot_bar(ps_nem_filt_lambs, "Sex", "Abundance", "SpeciesFull", title=title)


#microbiome site one, not sure how to work yet
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
p <- plot_composition(ps_nem_filt_repeatsrem,
                      taxonomic.level = "SpeciesFull",
                      sample.sort = "SampleId",
                      x.label = "SampleId") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  #facet_wrap(~AgeClass) +
  labs(x = "SampleId", y = "Relative abundance (%)",
       title = "Relative abundance data",
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  #theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        #strip.background = element_rect(fill="white")) +
  theme_ipsum(grid="Y")
print(p + theme(axis.text.x = element_text(angle = 90)))

#another one, doesn't look great
fullplot <- plot_composition(ps_nem_filt_repeatsrem, sample.sort = NULL, 
                          otu.sort = NULL,
                          x.label = "SampleId", # sample type
                          plot.type = "barplot", 
                          verbose = FALSE) + 
  theme_bw() + scale_fill_brewer("SpeciesFull", palette = "Paired")
# we can rotate x axis labels 
print(fullplot + theme(axis.text.x = element_text(angle = 90)))

#another one, comes up with 'object of type 'closure' is not subsettable 
p <- ps_nem_filt_lambs %>%
  plot_composition(sample.sort = "Sex", otu.sort = "abundance") +
  # Set custom colors
  scale_fill_manual(values = default_colors("SpeciesFull")[taxa(ps_nem_filt_lambs)]) +
  scale_y_continuous(label = scales::percent)

print(p)
#another microbiome site
p <- plot_composition(ps_nem_filt_adults,
                      average_by = "Sex", transform = "compositional") #will only work for something that is not classed as character.
print(p)

#The one in Amy's phyloseq  script
psdat.species <- tax_glom(ps_nem_filt_repeatsrem, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.species <- psmelt(psdat.species) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

ps_nem_filt_repeatsrem%>%mutate(Season = fct_relevel(Season, "Mar-19", "May-19", "Jul-19")) ##doesn't work for a phyloseq object

RelSpecies<-ggplot(ps.melt.species, 
                   aes(x = Sample, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpecies) 



psdat.specieslambs <- tax_glom(ps_nem_filt_lambs, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.specieslambs <- psmelt(psdat.specieslambs) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesLambs<-ggplot(ps.melt.specieslambs, 
                   aes(x = Sample, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  #facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesLambs) 



psdat.speciesadults <- tax_glom(ps_nem_filt_adults, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.speciesadults <- psmelt(psdat.speciesadults) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesAdults<-ggplot(ps.melt.speciesadults, 
                        aes(x = Sample, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  #facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesAdults) 


#splitting the datasets into adult only and lamb only 

adults <- sample_data(ps_nem_filt_repeatsrem) %>% as_tibble %>% filter(AgeClass =="Adult") %>% pull(SampleId)
ps_nem_filt_adults <- prune_samples((sample_names(ps_nem_filt_repeatsrem)%in% adults), ps_nem_filt_repeatsrem)
view(sample_data(ps_nem_filt_adults))
saveRDS(ps_nem_filt_adults, "Adultsonly.rds")

lambs <- sample_data(ps_nem_filt_repeatsrem) %>% as_tibble %>% filter(AgeClass =="Lamb") %>% pull(SampleId)
ps_nem_filt_lambs <- prune_samples((sample_names(ps_nem_filt_repeatsrem)%in% lambs), ps_nem_filt_repeatsrem)
view(sample_data(ps_nem_filt_lambs))
saveRDS(ps_nem_filt_lambs, "Lambsonly.rds")
view(otu_table(ps_nem_filt_lambs))



#splitting adults into seasonal datasets
march <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Season =="Mar-19") %>% pull(SampleId)
ps_nem_filt_march <- prune_samples((sample_names(ps_nem_filt_adults)%in% march), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_march))


may <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Season =="May-19") %>% pull(SampleId)
ps_nem_filt_may <- prune_samples((sample_names(ps_nem_filt_adults)%in% may), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_may))

july <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Season =="Jul-19") %>% pull(SampleId)
ps_nem_filt_july <- prune_samples((sample_names(ps_nem_filt_adults)%in% july), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_july))

#seasonal abundance plot
psdat.speciesmarch <- tax_glom(ps_nem_filt_march, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.speciesmarch <- psmelt(psdat.speciesmarch) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesMarch<-ggplot(ps.melt.speciesmarch, 
                         aes(x = Sample, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  #facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesMarch) 


psdat.speciesmay <- tax_glom(ps_nem_filt_may, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.speciesmay <- psmelt(psdat.speciesmay) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesMay<-ggplot(ps.melt.speciesmay, 
                        aes(x = Sample, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  #facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesMay) 



psdat.speciesjuly <- tax_glom(ps_nem_filt_july, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.speciesjuly <- psmelt(psdat.speciesjuly) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesJuly<-ggplot(ps.melt.speciesjuly, 
                      aes(x = Sample, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  #facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19", "Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesJuly) 

ps.melt.speciesseason%>%mutate(Season = fct_relevel(Season, "Mar-19", "May-19", "Jul-19"))

psdat.speciesseason <- tax_glom(ps_nem_filt_adults, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.speciesseason <- psmelt(psdat.speciesseason) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesSeason<-ggplot(ps.melt.speciesseason, 
                          aes(x = Season, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ Season , scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  scale_x_discrete(limits=c("Mar-19", "May-19","Jul-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesSeason) #need to figure out how to get labels

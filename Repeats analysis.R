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

view(otu_table(ps_nem_filt_repeats))
tax_table(ps_nem_filt_repeats)
view(sample_data(ps_nem_filt_repeats))
#all samples are adults during March and may only 
ps_nem_sex <- sample_data(ps_nem_filt_repeats)$SexCorr <- ifelse(sample_data(ps_nem_filt_repeats)$Sex ==1, "Female", "Male")

#abundance plot all repeats

psdat.speciesrepeats <- tax_glom(ps_nem_filt_repeats, taxrank = "SpeciesFull")  #ask about exactly what this does again
ps.melt.speciesrepeats <- psmelt(psdat.speciesrepeats) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))

RelSpeciesrepeats<-ggplot(ps.melt.speciesrepeats, 
                       aes(x = Id, y = Abundance, fill = SpeciesFull)) + 
  geom_bar(stat = "identity", position = "fill") + 
  #facet_wrap(~ Age ) + #, scales = "free_y") +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  #scale_x_discrete(limits=c("Mar-19", "May-19")) +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpeciesrepeats) 

p <- plot_composition(ps_nem_filt_repeats,
                      average_by = "Age", transform = "compositional")


## Beta Diversity 
#Need to use compostional data set

bx.ord_pcoa_bray_ps_nemrepeats <- ordinate(ps_nem_filt_repeats, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nemrepeats) + theme_bw()


# Axis 1 and 2 are of interest - pcoa 
beta_ps_nemrepeatsSeason <- plot_ordination(ps_nem_filt_repeats, 
                               bx.ord_pcoa_bray_ps_nemrepeats, 
                               color="Season",  #factor of interest 
                               #label = "subject",
) + 
  geom_point(aes(shape = Season), size= 4) +  #factor of interest 
  #geom_label(aes(label=sampleId), size=1) + 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)

plot(beta_ps_nemrepeatsSeason)

beta.ps_nem_repeatssex <- plot_ordination(ps_nem_filt_repeats, #need to change sex from numbers to labels within a phyloseq object. Won't work with metadata.
                                   bx.ord_pcoa_bray_ps_nemrepeats, 
                                   color= "SexCorr", #factor of interest 
                                   #label = "sample_id",
) + 
  geom_point(aes(shape =SexCorr), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120) 


plot(beta.ps_nem_repeatssex)


## PERMANOVA. Used to detect whether different groups (e.g male vs female) have a serious effect on overall worm composition
metadf.bxrepeats <- data.frame(sample_data(ps_nem_filt_repeats))
bray_ps_nemrepeats <- phyloseq::distance(physeq = ps_nem_filt_repeats, method = "bray")

set.seed(995)
# Adonis test
library(vegan)
adonis.testSeason <- adonis(bray_ps_nemrepeats ~ Season, data = metadf.bxrepeats)

adonis.testSeason  ## so for season, there is not much of a correlation between groups r2= 0.15593? Though a sig diff between groups p = 0.001?
#does this result suggest a difference between just 2 of the groups, or that every group is significantly different from each other?
#can try for AgeClass
adonis.test.Age <- adonis(bray_ps_nem ~ AgeClass, data = metadf.bx)

adonis.test.Age   


adonis.test.Sex <- adonis(bray_ps_nem ~ SexCorr, data = metadf.bx)

adonis.test.Sex 

#the following looks to check that the variance of homogeneity dispersal assumption holds, to ensure relability of the results of the permanova/adonis
#a non signif p value indicates that there are no differences within groups and that therefore they are homogenous
dist <- vegdist(t(abundances(ps_nem_filt_repeats)))
anova(betadisper(dist, metadf.bxrepeats$Season))# p = 1.033e-05, therefore dispersions are heterogenous, violating the assumption of adonis test
anova(betadisper(dist, metadf.bxrepeats$SexCorr))  
anova(betadisper(dist, metadf.bx$AgeClass)) 


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



#Investigating within season repeatability----------------------------------------------------------------------

#Need to run ordination and permanova against Season Id

bx.ord_pcoa_bray_ps_nemrepeats <- ordinate(ps_nem_filt_repeats, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nemrepeats) + theme_bw()


# Axis 1 and 2 are of interest - pcoa 
beta_ps_nemrepeatsSeason <- plot_ordination(ps_nem_filt_repeats, 
                                            bx.ord_pcoa_bray_ps_nemrepeats, 
                                            color="Tag",  #factor of interest 
                                            #label = "subject",
) + 
  #geom_point(aes(shape = Tag), size= 4) +  #factor of interest 
  #geom_label(aes(label=Tag), size=2) +
  #geom_point(aes(color=Tag), size=4) +
  scale_colour_brewer(palette="PuOr")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)

plot(beta_ps_nemrepeatsSeason)



## PERMANOVA. Used to detect whether different groups (e.g male vs female) have a serious effect on overall worm composition
metadf.bxrepeats <- data.frame(sample_data(ps_nem_filt_repeats))
bray_ps_nemrepeats <- phyloseq::distance(physeq = ps_nem_filt_repeats, method = "bray")

set.seed(995)
# Adonis test
library(vegan)
#Against Season Id
adonis.testSeason <- adonis(bray_ps_nemrepeats ~ IdSeason, data = metadf.bxrepeats)

adonis.testSeason

distSeasonId <- vegdist(t(abundances(ps_nem_filt_repeats)))
anova(betadisper(distSeasonId, metadf.bxrepeats$Season))

#Against Individual Tag
adonis.testTag <- adonis(bray_ps_nemrepeats ~ Tag, data = metadf.bxrepeats)

adonis.testTag

distTag <- vegdist(t(abundances(ps_nem_filt_repeats)))
anova(betadisper(distTag, metadf.bxrepeats$Tag))

#Investigating across season repeatability--------------------------------------------------------------------------

#Need to filter main dataset to sheep Ids that appear 3 times to represent individuals that were sampled at each season
seasonrepeat <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Tag >2) %>% pull(SampleId)
ps_nem_filt_seasonrepeat <- prune_samples((sample_names(ps_nem_filt_adults)%in% seasonrepeat), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_seasonrepeat))


ps_nem_filt_noncomp <- filter_taxa(phyloseqNemabiomeSoaySpecies,
                                   function(x) {sum(x>1000) > 5}, prune=TRUE)
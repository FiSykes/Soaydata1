library(tidyverse) #important 
library(phyloseq) #important 
library(microbiome)  # important 
library(knitr)
library(ggpubr)
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
view(sample_data(phyloseqSpeciesProportionalCorrected))

#histograms

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Bunostomum_trigonocephalum") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Trichostrongylus_axei") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Teladorsagia_circumcincta") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Trichostrongylus_vitrinus") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Chabertia_ovina") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Trichostrongylus_colubriformis") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Haemonchus_contortus") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Oesophagostomum_venulosum") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Cooperia_oncophora") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Cooperia_fuelleborni") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Cooperia_curticei") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Nematodirus_battus") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Trichostrongylus_rugatus") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Muellerius_capillaris") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Dictyocaulus_filaria") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much 

tax_table(phyloseqNemabiomeSoaySpecies) %>% as.data.frame() %>% filter(SpeciesFull=="Oesophagostomum_dentatum") %>% rownames -> N_bat
otu_table(phyloseqNemabiomeSoaySpecies)[,N_bat] %>% rowSums() %>% hist() # really not very much

#Looking at sequence depths

SeqDepth = rowSums(otu_table(phyloseqNemabiomeSoaySpecies))
hist(SeqDepth)
sample_data(phyloseqNemabiomeSoaySpecies)$SeqDepth = SeqDepth

min(SeqDepth) #minimum sequence depth is 4007
max(SeqDepth) #maximum sequence depth is 169871
view(sort(SeqDepth))


ggbarplot(meta(phyloseqNemabiomeSoaySpecies), "SampleId", "SeqDepth") +     #creates a fun plot
  theme_bw(base_size=9) 

meta(phyloseqNemabiomeSoaySpecies) %>% 
  ggplot(aes(x=AgeClass, y=SeqDepth, colour=SampleType)) + 
  geom_sina() + ggsave("SeqDepthSina.png", units="mm", height=150, width=150)     #creates a fun plot

#Used the following to remove species at 0, this object is named ps_nem_soay

any(taxa_sums(phyloseqNemabiomeSoaySpecies) ==0)   #comes back as TRUE
ps_nem_soay <- prune_taxa(taxa_sums(phyloseqNemabiomeSoaySpecies) > 0, phyloseqNemabiomeSoaySpecies) #this has taken out species with 0 reads and named the object as ps_nem_soay
view(otu_table(ps_nem_soay))

#summary

summarize_phyloseq(phyloseqNemabiomeSoaySpecies)
sample(phyloseqNemabiomeSoaySpecies)



#Making plots!!

colnames(ps_nem_soay_comptest@tax_table)[8] <- "SpeciesFull"


psdat.species <- tax_glom(ps_nem_soay_comptest, taxrank = "Species")  #ask about exactly what this does again
ps.melt.species <- psmelt(psdat.species) %>% 
  mutate(SampleId=fct_reorder(Sample, AgeClass))


RelSpecies<-ggplot(ps.melt.species, 
                   aes(x = Sample, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ Age) +
  #scale_fill_viridis_d(option="plasma")+
  theme_bw(base_size = 8) + # labs(subtitle = "A") + 
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=2), 
        strip.background = element_rect(fill="white"))
plot(RelSpecies)  #doesn't distinguish between different Trichostrongylus species, lumps them together???
#What is a better alternative to having all the sample Ids squished together? Tried splitting into ageclass using facet function but looks weird


phyloseq_prevalence_plot(phyloseqNemabiomeSoaySpecies, taxcolor = "SpeciesFull", facet=TRUE, point_alpha = 0.5, prev.trh = 0.05) +
  theme_bw(base_size = 12)+theme(legend.position = "none")  #makes a prevalence plot sorta like the one Amy made

phyloseq_prevalence_plot(phyloseqNemabiomeSoaySpecies, taxcolor = "SpeciesFull", facet=FALSE, point_alpha = 0.5, 
                         prev.trh = 0.05) +
  theme_bw(base_size = 12) + 
  geom_label(aes(label=SpeciesFull))  #just for viz purposes the labels a bit clunky 


#Found this on the microbiome website, not currently sure how to get to work. Should produce something relatively similar to first plot
p <- plot_composition(ps_nem_soay_comptest,
                      taxonomic.level = "Genus",
                      sample.sort = "nationality",
                      x.label = "nationality") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data",
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  theme_ipsum(grid="Y")
print(p)  

title = "Abundance plot"
plot_bar(ps_nem_soay_comptest, "Season", "Abundance", "Genus", title=title) #phyloseq uses this function to plot abundance plots?
#makes some weird looking plots??? 

#converting data to a compostional dataset

ps_nem_soay_comptest <- microbiome::transform(Thresholdtest, "compositional") #back to sum to one 
otu_table(ps_nem_soay_comptest) %>% rowSums()

saveRDS(ps_nem_soay_comptest, "phyloseqSpeciesProportionalTest.rds")
view(otu_table(ps_nem_soay_comptest))  #IT WORKED!!!! YAY!!


#Trying out some thresholds/filtering related stuff

any(taxa_sums(ps_nem_soay_comptest) <0.05)   #comes back as TRUE
ps_nem_prunetest <- prune_taxa(taxa_sums(ps_nem_soay_comptest) > 0.05, ps_nem_soay_comptest)
view(otu_table(ps_nem_prunetest)) #worked. Got rid of everything up to N.battus. Does <0.05 mean less than that across all samples,yes.hmmm 


GP = filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE) # found this on the phyloseq site. Means Remove taxa not seen more than 3 times in at least 20% of the samples.    

Thresholdtest <- filter_taxa(phyloseqNemabiomeSoaySpecies, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
view(otu_table(Thresholdtest)) #this actually leads to the removal of all the problematic species :o Even if reduced to in 5% of samples


Thresholdtest3 <- filter_taxa(phyloseqNemabiomeSoaySpecies, function(x) sum(x > 0) > (0.05*length(x)), TRUE)
view(otu_table(Thresholdtest3)) # this, removing taxa with less than 5% prevalence also removes all problematic species

Thresholdtest4 <- filter_taxa(phyloseqNemabiomeSoaySpecies, function(x) {sum(x > 1000) > 5}, prune = TRUE)
view(otu_table(Thresholdtest4)) #I thiiiiink this translates as remove samples that don't appear at least 5 times at more than 1000 reads.
#Also removed all problematic species


#Try subsetting repeat samples

ps_nem_repeats <- prune_samples(!(sample_names(ps.melt.species) %in% RepeatTech), ps.melt.species)  
#not sure how to get this to work

ps_nem_soay <- prune_taxa(taxa_sums(phyloseqNemabiomeSoaySpecies) > 0, phyloseqNemabiomeSoaySpecies ) #abundance filtering, change 0 to a specific number of reads

ps_nem_soay <- filter_taxa(ps_nem_soay, function(x){sum(x > 0) > 5 }, prune = TRUE) #prevalence filtering. In this form means that taxa appearing less than 5 times will be removed. 


#Trying the thresholds out
Filteringtest <- prune_taxa(taxa_sums(phyloseqNemabiomeSoaySpecies) > 5000, phyloseqNemabiomeSoaySpecies ) #abundance filtering
view(otu_table(Filteringtest))
Thresholdtest2 <- filter_taxa(Filteringtest, function(x){sum(x > 0) > 5 }, prune = TRUE) #prevalence filtering
view(otu_table(Thresholdtest2))
#sadly because some samples for Haemonchus and O. venulosum are at 15000 it means these do not get filtered out. Is there a way to combine the two functions so they work at the same time?

#Thinking about alpha diversity (how many species present)
# Absolute abundances for the single most abundant taxa in each sample # gives the summary in a table not sure how to interpret
tab <- dominance(phyloseqNemabiomeSoaySpecies, index = "all")
kable(head(tab))
#rarity and low abundance taxa, again not sure how to interpret
tab <- rarity(phyloseqNemabiomeSoaySpecies, index = "all")
kable(head(tab))


plot_richness(ps.melt.species, x="SampleType", measures=c("Chao1", "Shannon")) #

plot_richness(ps.melt.species)





#Analysis questions to ask that are not already in script
#Minimum sequence depth is 4007, so take it samples with less than 4000 reads taken out this time?
#is it possible to merge the metadata dataset yolanda sent with the main dataset to incorporate Larval Numbers?
#confused by the plot alpha diversity tutorial on the phyloseq website

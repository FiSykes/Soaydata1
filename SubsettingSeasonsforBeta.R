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

#March and JUly beta diversity-------------------------------------------------------------------------------- 

marchandjuly <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Season =="May-19") %>% pull(SampleId)
ps_nem_filt_marchjuly <- prune_samples(!(sample_names(ps_nem_filt_adults)%in% marchandjuly), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_marchjuly))

bx.ord_pcoa_bray_ps_nemmarchjuly <- ordinate(ps_nem_filt_marchjuly, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nemmarchjuly) + theme_bw()


beta.ps_nem_marchjuly <- plot_ordination(ps_nem_filt_marchjuly, 
                                      bx.ord_pcoa_bray_ps_nemmarchjuly, 
                                      color="Season", #factor of interest 
                                      #label = "sample_id",
) + 
  geom_point(aes(shape = Season), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)



metadf.bxmarchjuly <- data.frame(sample_data(ps_nem_filt_marchjuly))
bray_ps_nemmarchjuly <- phyloseq::distance(physeq = ps_nem_filt_marchjuly, method = "bray")
set.seed(995)


#adonis and betadisper code
adonis.test.marchjuly <- adonis(bray_ps_nemmarchjuly ~ Season, data = metadf.bxmarchjuly)

adonis.test.marchjuly 

#the following looks to check that the varianc of homogeneity assumptions hold, to ensure relability of the results of the permanova
dist <- vegdist(t(abundances(ps_nem_filt_marchjuly)))
anova(betadisper(dist, metadf.bxmarchjuly$Season))






#March and May beta diversity-------------------------------------------------------------------

marchandmay <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Season =="Jul-19") %>% pull(SampleId)
ps_nem_filt_marchmay <- prune_samples(!(sample_names(ps_nem_filt_adults)%in% marchandmay), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_marchmay))

bx.ord_pcoa_bray_ps_nemmarchmay <- ordinate(ps_nem_filt_marchmay, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nemmarchmay) + theme_bw()


beta.ps_nem_marchmay <- plot_ordination(ps_nem_filt_marchmay, 
                                         bx.ord_pcoa_bray_ps_nemmarchmay, 
                                         color="Season", #factor of interest 
                                         #label = "sample_id",
) + 
  geom_point(aes(shape = Season), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)



metadf.bxmarchmay <- data.frame(sample_data(ps_nem_filt_marchmay))
bray_ps_nemmarchmay <- phyloseq::distance(physeq = ps_nem_filt_marchmay, method = "bray")
set.seed(995)


#adonis and betadisper code
adonis.test.marchmay <- adonis(bray_ps_nemmarchmay ~ Season, data = metadf.bxmarchmay)

adonis.test.marchmay 

#the following looks to check that the varianc of homogeneity assumptions hold, to ensure relability of the results of the permanova
dist <- vegdist(t(abundances(ps_nem_filt_marchmay)))
anova(betadisper(dist, metadf.bxmarchmay$Season))






#July and May beta diversity-------------------------------------------------------------------

julyandmay <- sample_data(ps_nem_filt_adults) %>% as_tibble %>% filter(Season =="Mar-19") %>% pull(SampleId)
ps_nem_filt_julymay <- prune_samples(!(sample_names(ps_nem_filt_adults)%in% julyandmay), ps_nem_filt_adults)
view(sample_data(ps_nem_filt_julymay))

bx.ord_pcoa_bray_ps_nemjulymay <- ordinate(ps_nem_filt_julymay, "PCoA", "bray")
#bx.ord_nmds_brayDiet <- ordinate(ps_nem_filt, "NMDS", "bray") #not using MMDS analysis for this I don't think

#Scree plot
plot_scree(bx.ord_pcoa_bray_ps_nemjulymay) + theme_bw()


beta.ps_nem_julymay <- plot_ordination(ps_nem_filt_julymay, 
                                        bx.ord_pcoa_bray_ps_nemjulymay, 
                                        color="Season", #factor of interest 
                                        #label = "sample_id",
) + 
  geom_point(aes(shape = Season), size= 4) + # factor of interest 
  scale_colour_brewer(palette="Paired")+ theme_bw() +
  theme(plot.title = element_text(hjust = 0, size = 12))+
  stat_ellipse() #+ggsave("SeasonOrdination.tiff", units="mm", width=150, height=120)



metadf.bxjulymay <- data.frame(sample_data(ps_nem_filt_julymay))
bray_ps_nemjulymay <- phyloseq::distance(physeq = ps_nem_filt_julymay, method = "bray")
set.seed(995)


#adonis and betadisper code
adonis.test.julymay <- adonis(bray_ps_nemjulymay ~ Season, data = metadf.bxjulymay)

adonis.test.julymay 

#the following looks to check that the varianc of homogeneity assumptions hold, to ensure relability of the results of the permanova
dist <- vegdist(t(abundances(ps_nem_filt_julymay)))
anova(betadisper(dist, metadf.bxjulymay$Season))

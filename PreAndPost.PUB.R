#Mason S Ward
#email: msw5688@psu.edu
#Co authors: Jon Sweetman (PhD), Sara Hermann (PhD), Alice Belskis (MS)
#Penn State University
#Department of Ecosystem Science & Management
#Started analysis in April 2024
#Research question: Does pre-drying & post-refilling of vernal ponds 
#have an effect on temporal changes in macroinvertebrate diversity?


#Load all applicable libraries - will probably not use the majority of these:
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(vegan)
library(tidyverse)
library(adespatial)
library(betapart)
library(lme4)
library(cowplot)
library(viridis)
library(sf)
library(ggspatial)
library(reshape2)

#Load in dataset from Github - 
##Username: masward -> Repository -> macro_dry
macro.pre <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/macro_predry.csv")
macro.post <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/macro_postdry.csv")
#Check the structure of each dataframe - looks good

#Vernal pond mapping----
##Need to have applicable shape files for mapping downloaded. 
fig.theme <- theme(text = element_text(family = "Times New Roman", color = "black"),
                   axis.text = element_text(color = "black", size = 10),
                   axis.title = element_text(size = 12),
                   panel.grid = element_blank(),
                   legend.text = element_text(size = 11),
                   legend.title = element_text(size = 11),
                   plot.caption.position = "plot", plot.caption = element_text(hjust =0, size = 11))

vernal.map_df <- data.frame(state_forest  = c("Bald_Eagle", "Bald_Eagle","Moshannon","Moshannon", "Sproul"),
                            SITE_ID  = c("BDEG_2", "BDEG_3", "MSHN_1", "MSHN_2","SPRL_3"),
                            lat = c("40.80081", "40.80368", "40.84595", "40.93583", "41.15978"),
                            long = c("-77.506195", "-77.499264", "-78.097878", "-78.032284", "-77.899456"))

##Use shapefile to create map!
pond_points <- st_as_sf(vernal.map_df, coords = c("long", "lat"), crs = 4269)  #4269 is the code for NAD83 obtained from here https://epsg.io/4269 


shp_forest <- st_read("/Users/masonward/Downloads/DCNR_BOF_StateForests202306/DCNR_BOF_StateForests202306.shp") #shape file of state forests
shp_forest #display information of the shapefile
st_crs(shp_forest)
st_bbox(shp_forest)
shp_boundaries <- st_read("/Users/masonward/Downloads/DCNR_BOF_Bndry_SFM201703/DCNR_BOF_Bndry_SFM201703.shp") #shape file of the boundaries of state forests
shp_boundaries #display information of the shapefile
st_crs(shp_boundaries)
st_bbox(shp_boundaries)

#Subset state forests by study area (both forests and boundaries)
shp_forest_subset <- subset(shp_forest, SF_Name %in% c('Moshannon', 'Bald Eagle','Sproul')) 
shp_boundaries_subset <- subset(shp_boundaries, DistrictNa %in% c('Moshannon', 'Bald Eagle', 'Sproul'))

#Plot and check how everything is looking
(ForestMap <- ggplot() +geom_sf(data = shp_forest, aes(fill = SF_Name)) +
    geom_sf(data = shp_boundaries) +
    ggtitle("Map of Pennsylvania State Forests") +
    ylab(expression("Latitude ("*degree*")" )) + 
    xlab(expression("Longitude ("*degree*")" )) +
    guides(fill=guide_legend(title="State Forests")) +
    fig.theme)
ggsave(filename = "forestmap.png", plot = ForestMap)

ggplot() +
  geom_sf(data = pond_points)

#Overlay study sites with the subset of state forests/boundaries
map <- (VernalMap <- ggplot() + geom_sf(data = shp_forest_subset, aes(fill = SF_Name)) + 
    geom_sf(data = shp_boundaries_subset, aes(fill = DistrictNa)) + 
    geom_sf(data = pond_points) + 
    scale_fill_manual(values = c("#E5F5F9", "#D4B9DA", "#D9D9D9")) + 
    ggtitle("Map of vernal pond study sites") +
    theme_bw() + 
    ylab(expression("Latitude ("*degree*")" )) + 
    xlab(expression("Longitude ("*degree*")" )) +
    guides(fill=guide_legend(title="State Forests")) +
    coord_sf() +
    annotation_scale()) +
  fig.theme
ggsave(filename = "vernalmap.png", plot = map)
#Edit map in Powerpoint to add on identification tags

#Macroinvertebrate Diversity - Alpha ----
##Calculate richness and Shannon's Div using the vegan package
richness_pre <- specnumber(macro.pre[3:37]) #specify the columns to calculate
richness_post <- specnumber(macro.post[3:37]) #specify the columns to calculate
#Check for normality
shapiro.test(richness_pre) #not normal
shapiro.test(richness_post) #not normal
#Use wilcoxen test - non-parametric t test
wilcox.test(richness_pre, richness_post, exact = FALSE)
#W = 20, p-value = 0.1425; not significant


shannon_pre <- diversity(macro.pre[3:37]) #specify the columns to calculate
shannon_post <- diversity(macro.post[3:37]) #specify the columns to calculate
#Check for normality
shapiro.test(shannon_pre) #not normal
shapiro.test(shannon_post) #not normal
#Use wilcoxen test - non-parametric t test
wilcox.test(shannon_pre, shannon_post, exact = FALSE)
#W = 19, p-value = 0.2101; not significant


#Create new dataframe with values computed from vegan package -> did this in excel
##uploaded to Github as "macro_rich_shdiv"
macro_rich_shannon <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/macro_rich_shdiv.csv")

new_x_label <- as_labeller(c("pre_dry" = "pre_dry",
                             "post_dry" = "pst_dry"))

#Create boxplots of richness and shannon diversity 
richness_plot <- ggplot(macro_rich_shannon, aes(x = fct_reorder(week_col, richness, .desc = TRUE), y=richness)) + #reordering so pre- goes before post-
  geom_boxplot(aes(fill = fct_reorder(week_col, richness, .desc = TRUE)), position = position_dodge(width = .75)) + # same as before - reordering so pre- goes before post-
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw(base_size = 14) +
  xlab(" ") +
  ylab("Richness") +
  theme_cowplot(12) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800"),
                    guide = guide_legend(title = "Pre- and Post dry"),
                    labels = c("pre_dry", "pst_dry")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

shannon_plot <- ggplot(macro_rich_shannon, aes(x = fct_reorder(week_col, shannon, .desc = TRUE), y=shannon)) + #reordering so pre- goes before post-
 geom_boxplot(aes(fill = fct_reorder(week_col, shannon, .desc = TRUE)), position = position_dodge(width = .75)) + #same as before - reordering so pre- goes before post-
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw(base_size = 14) +
  xlab("") +
  scale_x_discrete(labels = new_x_label) +
  ylab("Shannon") +
  theme_cowplot(12) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800"),
                    guide = guide_legend(title = "Pre- and Post dry"),
                    labels = c("pre_dry", "pst_dry")) 


alpha_plot <- ggarrange(richness_plot, shannon_plot,
                        labels = c("A", "B"), hjust = -0.5,
                        ncol = 1, nrow = 2,
                       legend = 'none',
                      heights = c(4,4)) +
  theme(plot.margin = NULL)
alpha_plot
ggsave(filename = "fig2updated.jpg", bg = "white", width=6, height=6, plot = alpha_plot)

ggsave(filename = "fig2updated.jpg", bg = "white", plot = alpha_plot)


#Macroinvertebrate Diversity - Temporal Beta ----
pa.macro.pre <- macro.pre %>%
  mutate_if(is.numeric, ~1* (. !=0))
pa.macro.post <- macro.post %>%
  mutate_if(is.numeric, ~1* (. !=0))

#Use 'adespatial' package with function TBI - Legendre 2019
tbi.sp <- TBI(pa.macro.pre[,-1:-2], pa.macro.post[,-1:-2], method ="jaccard", nperm = 9999, test.t.perm = TRUE)
tbi.sp$BCD.summary
#General trend of losses among sites;
#Not significant based among t test
#For all sites, there is losses but also not sig -> 77% L v 23% G of dissimilarity 
#tpaired.krandtest(pa.macro.pre[,-1:-2], pa.macro.post[,-1:-2], nperm = 9999, list.all = FALSE)

#create table with results
sites <- as.data.frame(pa.macro.post[,2])
TBI <- as.data.frame(tbi.sp[[1]])
p.TBI <- as.data.frame(tbi.sp[[2]])
BCD <- as.data.frame(tbi.sp[[4]])
p.adjusted <- as.data.frame(tbi.sp[[3]])
#Bind to site names for TBI results
TBI.results <- cbind(sites, TBI, p.TBI, p.adjusted, BCD)
colnames(TBI.results) <- c("sites", "TBI", "p.TBI", "p.adjusted", "losses", "gains", 
                      "D=(B+C)/(2A+B+B", "Change")

#Plotting TBI results
TBI_results <- ggplot(TBI.results, aes(x = losses, y = gains, fill = sites, shape = Change, size = TBI, colour = sites)) +
  geom_point(alpha = 0.65, aes(colour = sites), position = position_dodge(width = 0.9)) +
  scale_size_continuous(range = c(4,10)) +
  scale_shape_manual(values = c(15, 16),
                     labels = c("Losses", "Gains")) +
#  scale_color_manual(guide = "Sites") +
  theme_cowplot(12) +
  theme(aspect.ratio = 1) +
  xlab("\nSpecies Losses") + ylab("Species Gains\n") +
  xlim(0, 1) + ylim(0, 1) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray") +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray") +
  annotate("text", x=.25, y=.75, label= "Gains > Losses", fontface = "bold") +
  annotate("text", x=.75, y=.25, label= "Losses > Gains", fontface = "bold") +
  annotate('text', x = .90, y = .865, label = '1-1 line', size = 4, angle='45')
TBI_results


ggsave(filename = "fig3updated.jpg", bg = "white", width=8, height=10, plot = TBI_results)


tbi.sp$p.TBI[tbi.sp$p.TBI<0.05] #0.0075
tbi.sp$TBI[]
hist(tbi.sp$TBI)


#Use beta.temp for temporal beta diversity
beta.temp <- beta.temp(pa.macro.pre[,-1:-2], pa.macro.post[,-1:-2], index.family="jaccard")
mean(beta.temp[,1], na.rm=T) #turnover, 0.4031746
mean(beta.temp[,2], na.rm=T) #nested, 0.3377666
mean(beta.temp[,3]) # total, 0.7409412

betapre <- beta.multi(pa.macro.pre[,-1:-2], index.family="jaccard")
betapost <- beta.multi(pa.macro.post[,-1:-2], index.family="jaccard")

betapre <- betapart::beta.multi(pa.macro.pre[,-1:-2], index.family="jaccard")
betapost <- betapart::beta.multi(pa.macro.post[,-1:-2], index.family="jaccard")
beta.prepost <- cbind(betapre,betapost)

#Spatial beta diversity 
jac.pre <- beta.pair(pa.macro.pre[,-1:-2], index.family="jac")
plot(hclust(jac.pre$beta.jac, method = "average"), hang = -1, main = '', sub = '')
jac.post <- beta.pair(pa.macro.post[,-1:-2], index.family="jac")
plot(hclust(jac.post$beta.jac, method = "average"), hang = -1, main = '', sub = '')


# make matrix where each column has all pairs (remove diagonal)
jac.pre.turn <- as.matrix(jac.pre$beta.jtu)
jac.pre.nest <- as.matrix(jac.pre$beta.jne)
jac.pre.tot <- as.matrix(jac.pre$beta.jac)
diag(jac.pre.turn) <- NA
diag(jac.pre.nest) <- NA
diag(jac.pre.tot) <- NA

jac.post.turn <- as.matrix(jac.post$beta.jtu)
jac.post.nest <- as.matrix(jac.post$beta.jne)
jac.post.tot <- as.matrix(jac.post$beta.jac)
diag(jac.post.turn) <- NA
diag(jac.post.nest) <- NA
diag(jac.post.tot) <- NA

#Turnover - Pre- and Post- 
#Testing if pre- spatial turnover is statistically different from temporal turnover
dis.site.turn.pre =vector()
for (i in 1:NCOL(jac.pre.turn)) {
  dis.site.turn.pre[i] <- mean(jac.pre.turn[,i], na.rm=TRUE)
}
mean(dis.site.turn.pre);sd(dis.site.turn.pre)
mean(beta.temp[,1]);sd(beta.temp[,1])

wilcox.test(beta.temp[,1], dis.site.turn.pre, exact = FALSE)
#ns; W = 5, p-value = 0.1425


#Testing if post- spatial turnover is statistically different from temporal turnover
dis.site.turn.post =vector()
for (i in 1:NCOL(jac.post.turn)) {
  dis.site.turn.post[i] <- mean(jac.post.turn[,i], na.rm=TRUE)
}
mean(dis.site.turn.post);sd(dis.site.turn.post)
mean(beta.temp[,1]);sd(beta.temp[,1])

wilcox.test(beta.temp[,1], dis.site.turn.post, exact = FALSE)
#ns; W = 5, p-value = 0.1425


#Nestedness - Pre- and Post- 
#Testing if pre- spatial nestedness is statistically different from temporal nestedness
dis.site.nest.pre =vector()
for (i in 1:NCOL(jac.pre.nest)) {
  dis.site.nest.pre[i] <- mean(jac.pre.nest[,i], na.rm=TRUE)
}
mean(dis.site.nest.pre);sd(dis.site.nest.pre)
mean(beta.temp[,2]);sd(beta.temp[,2])

wilcox.test(beta.temp[,2], dis.site.nest.pre)
#ns; W = 10, p-value = 0.6905


#Testing if post- spatial nestedness is statistically different from temporal nestedness
dis.site.nest.post =vector()
for (i in 1:NCOL(jac.post.nest)) {
  dis.site.nest.post[i] <- mean(jac.post.nest[,i], na.rm=TRUE)
}
mean(dis.site.nest.post);sd(dis.site.nest.post)
mean(beta.temp[,2]);sd(beta.temp[,2])

wilcox.test(beta.temp[,2], dis.site.nest.post)
#ns; W = 11, p-value = 0.8413


#Total - Pre- and Post- 
#Testing if pre- spatial total is statistically different from temporal total
dis.site.tot.pre =vector()
for (i in 1:NCOL(jac.pre.tot)) {
  dis.site.tot.pre[i] <- mean(jac.pre.tot[,i], na.rm=TRUE)
}
mean(dis.site.tot.pre);sd(dis.site.tot.pre)
mean(beta.temp[,3]);sd(beta.temp[,3])

wilcox.test(beta.temp[,3], dis.site.tot.pre)
#ns; W = 16, p-value = 0.5476


#Testing if post- spatial total is statistically different from temporal total
dis.site.tot.post =vector()
for (i in 1:NCOL(jac.post.tot)) {
  dis.site.tot.post[i] <- mean(jac.post.tot[,i], na.rm=TRUE)
}
mean(dis.site.tot.post);sd(dis.site.tot.post)
mean(beta.temp[,3]);sd(beta.temp[,3])

wilcox.test(beta.temp[,3], dis.site.tot.post)
#ns; W = 5, p-value = 0.1508


#Combine datasets for figure -> Do in excel 
betadiv <- read_csv("https://raw.githubusercontent.com/masward/macro_dry/main/betadiv.csv")
#Separate nest, turnover, and total into 6 plots -> easier to do in Excel
tmp_pre_turn <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/tmp_pre_turn.csv")
tmp_pre_nest <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/temp_pre_nest.csv")
tmp_pre_total <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/tmp_pre_total.csv")
tmp_pst_turn <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/tmp_pst_turn.csv")
tmp_pst_nest <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/temp_pst_nest.csv")
tmp_pst_total <- read.csv("https://raw.githubusercontent.com/masward/macro_dry/main/tmp_pst_total.csv")

#Remove temporal data from pre- since pre- and post contain the same
#tmp_pre_turn <- tmp_pre_turn[-c(1:5), ]
#tmp_pre_nest <- tmp_pre_nest[-c(1:5), ]
#tmp_pre_total <- tmp_pre_total[-c(1:5), ]

#Combine pre- and pst- datasets for each variable - remove temporal duplicates
#Turnover
pre_pst_turn <- rbind(tmp_pre_turn, tmp_pst_turn)
pre_pst_turn <- pre_pst_turn[-c(11:15),]
#Nestedness
pre_pst_nest <- rbind(tmp_pre_nest, tmp_pst_nest)
pre_pst_nest <- pre_pst_nest[-c(11:15),]
#Total
pre_pst_total <- rbind(tmp_pre_total, tmp_pst_total)
pre_pst_total <- pre_pst_total[-c(11:15),]


pre_pst_turn_plot <- ggplot(pre_pst_turn, aes(x = component, y = index, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  ylab("Dissimilarity") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.0))

pre_pst_nest_plot <- ggplot(pre_pst_nest, aes(x = component, y = index, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  ylab("Dissimilarity") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.0))

pre_pst_total_plot <- ggplot(pre_pst_total, aes(x = component, y = index, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  ylab("Dissimilarity") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.0))
ggsave(filename = "turnoverupdated.jpg", bg = "white", width=12, height=6, plot = pre_pst_turn_plot)
ggsave(filename = "nestednessupdated.jpg", bg = "white", width=12, height=6, plot = pre_pst_nest_plot)
ggsave(filename = "totalupdated.jpg", bg = "white", width=12, height=6, plot = pre_pst_total_plot)




#Figure 4—6 graphing ----
prenest_plot <- ggplot(tmp_pre_nest, aes(y = index, x = component, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  ylab("Dissimilarity") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.0))


pretotal_plot <- ggplot(tmp_pre_total, aes(y = index, x = component, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  ylab("Dissimilarity") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, NA))


preturn_plot <- ggplot(tmp_pre_turn, aes(y = index, x = component, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  labs(y = "Dissimilarity") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, NA))


pstnest_plot <- ggplot(tmp_pst_nest, aes(y = index, x = component, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  labs(y = " ") +
  theme_cowplot(12) +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.0))


psttotal_plot <- ggplot(tmp_pst_total, aes(y = index, x = component, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  labs(y = " ") +
  theme_cowplot(12) +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, NA))

pstturn_plot <- ggplot(tmp_pst_turn, aes(y = index, x = component, fill = type)) +
  geom_boxplot() +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL),
         shape = "none") +
  theme_cowplot(12) +
  labs(y = " ") +
  theme(legend.title = element_blank(), 
        text =  element_text(size = 18),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, NA))


pre.turn.nest.total.fig <- ggarrange(preturn_plot, prenest_plot, pretotal_plot,
                                     labels = c("A", "B", "C"), hjust = -0.5,
                                     nrow = 1, legend = 'none', 
                                     widths = 1, heights = 1)
pre.turn.nest.total.fig #Maybe we should split into 3 separate graphs?

pst.turn.nest.total.fig <- ggarrange(pstturn_plot, pstnest_plot, psttotal_plot,
                                     labels = c("A", "B", "C"), hjust = -0.5,
                                     nrow = 1, legend = 'none', 
                                     widths = 1, heights = 1)
pst.turn.nest.total.fig #Maybe we should split into 3 separate graphs?



turnoverfig <- ggarrange(preturn_plot, pstturn_plot,
                         labels = c("A", "B"), hjust = -0.5,
                         nrow = 1, legend = 'none', 
                         widths = 1, heights = 1)
turnoverfig

nestedfig <- ggarrange(prenest_plot, pstnest_plot,
                         labels = c("A", "B"), hjust = -0.5,
                         nrow = 1, legend = 'none', 
                         widths = 1, heights = 1)
nestedfig

totalfig <- ggarrange(pretotal_plot, psttotal_plot,
                         labels = c("A", "B"), hjust = -0.5,
                         nrow = 1, legend = 'none', 
                         widths = 1, heights = 1)
totalfig

ggsave(filename = "fig4updated.jpg", bg = "white", width=12, height=6, plot = turnoverfig)
ggsave(filename = "fig5updated.jpg", bg = "white", width=12, height=6, plot = nestedfig)
ggsave(filename = "fig6updated.jpg", bg = "white", width=12, height=6, plot = totalfig)



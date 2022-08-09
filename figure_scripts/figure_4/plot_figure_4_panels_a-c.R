# library(ape)
# library(ggtree)
# library(phytools)
# library(tidytree)
# library(dplyr)
# library(stringr)
# library(vegan)
# library(ggplot2)
# library(tidyr)
# library(ggridges)
# library(ggbeeswarm)  



library(micro.gen.extra)
library(UpSetR)


#### 1 - get tree ####
tree_cols <- c("0" = "black",
               "1" = "#bab0ac",
               "2" = "#59a14f",
               "3" = "#e15759",
               "4" = "#76b7b2",
               "5" = "#4e79a7",
               "6" = "#f28e2b",
               "7" = "#9c755f",
               "8" = "#edc948",
               "9" = "#b07aa1",
               "10" = "#D7B5A6",
               "11" = "#ff9da7",
               "12" = "#a0a0a0",
               "13" = "#B3E2CD",
               "14" = "black",
               "0" = "black")

serratia_tree <- ggtree::read.tree("../../figshare_data/tree/snp_sites_aln_panaroo.aln.treefile")
rooted_serratia_tree <- phytools::midpoint.root(serratia_tree)
t1 <- ggtree::ggtree(rooted_serratia_tree, ladderize = T, right = T)


## fastani clusters:

fastani_data <- read.csv("../../figshare_data/clusters/95_to_99_ani_clusters.csv",
                         stringsAsFactors = F)

fastani_95 <- fastani_data %>% 
  select(Name,X95) %>% 
  rename(strain = Name, phylogroup = X95) %>% #phylogroups set by 95% ANI
  mutate(strain = gsub("#","_",strain))

##### dplyr example to get the groups from this dataframe to put into groupOTU:

dplyr_genus_groups <- fastani_data %>%
  select(Name, X95) %>% # only take the name of the strains and the clusters when cut off at 95% ANI
  mutate(Name = gsub("#","_",Name)) %>% # remove hashes from WSI lane id's
  group_by(X95) %>%  
  mutate(strains = list(Name)) %>% # get list of strains in each 95% ANI cluster
  select(X95,strains) %>%
  unique() %>% # remove the strain name/ lane id column and remove duplicate rows
  pull(strains, name = X95) # take the list of the lists of strain names/lane ids and name each item of the list by the fastANI cluster



genus_groups_tree <- ggtree::groupOTU(rooted_serratia_tree, dplyr_genus_groups) #, overlap='abandon',connect = T)
#get tibble data:
tree_tibble <- as_tibble(genus_groups_tree)


## plot tree second time with colours:

t2 <-  ggtree::ggtree(genus_groups_tree, ladderize = T, right = T, size=0.5, aes(color=group)) +
  scale_color_manual(values = c(tree_cols)) +
  ggtree::geom_treescale(x=0, y=0, offset = 10, width = 0.1, linesize = 0.5) +
  theme(legend.position = "none")



### add in fastbaps data: 
fastbaps <- read.csv("../../figshare_data/clusters/fastbaps_clusters.csv",
                     row.names = 1,
                     stringsAsFactors = T)


## use level three
fastbaps_l3 <- fastbaps %>%
  select(Level.3) %>%
  rename(cluster = Level.3) %>%
  tibble::rownames_to_column(var = "strain") %>%
  mutate(cluster = as.factor(cluster))

##get the groups of the fastbaps:

fastbaps_groups <- fastbaps_l3 %>%
  group_by(cluster) %>%
  summarise(strains = list(strain)) %>%
  pull(strains, name = cluster)

#get the clade names:
clade_labels <- micro.gen.extra::get_clade_nodes(genus_groups_tree, fastbaps_l3)


#draw the tree:
new_tree <- micro.gen.extra::add_clades_to_tree(t2,clade_labels)


#### 2 - read in GC content data: ####

intergenic <- read.csv("../../figshare_data/gc_content/all_intergenic_gc_values.csv",
                       header = F,
                       stringsAsFactors = F,col.names = c("strain","intergenic_gc"))

coding <- read.csv("../../figshare_data/gc_content/all_coding_gc_values.csv",
                   header = F,
                   stringsAsFactors = F,col.names = c("strain","coding_gc")) %>%
  mutate(strain = gsub("_coding","",strain))


##merge:
gc_data_w_phylogroup <- intergenic %>%
  left_join(coding, by = "strain") %>%
  mutate(strain = gsub("#","_",strain)) %>% # remove hashes from lane id's
  left_join(fastani_95)
  

#write out for source data

write.csv(file = "figure_4b.csv",
          x = gc_data_w_phylogroup,
          quote = F,
          row.names = F)

#plot coding vs. non coding - panel b

p2 <- ggplot(gc_data_w_phylogroup, aes(intergenic_gc, coding_gc)) + 
  geom_smooth(method=lm, colour = "black", alpha = 0.5) +
  ggpubr::stat_cor(method = "pearson",#p.accuracy = 0.0000000000000000000000000000001, r.accuracy = 0.01,#  +
                   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))
  ) +
  geom_point(aes(colour = as.factor(phylogroup))) +
  scale_color_manual(values = c(tree_cols)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("GC content: Non-coding (%)") +
  ylab("GC content: Coding (%)")


### get gc content of the whole genome :
serratia_metadata <- read.csv("../../figshare_data/metadata/serratia_metadata.csv", header = T,
                              quote = "",
                              stringsAsFactors = F,
                              comment.char = "")

#add in the phlogroup
whole_genome_gc <- subset.data.frame(serratia_metadata, select = c("File_prefix","GC...."))  %>%
  rename(gc_perc = GC...., strain = File_prefix) %>%
  merge(fastani_95)

#write out for source data

write.csv(file = "figure_4a.csv",
          x = whole_genome_gc,
          quote = F,
          row.names = F)


## plot whole genome gc histogram - panel a:
p1 <- ggplot(whole_genome_gc, aes(gc_perc, fill = as.factor(phylogroup))) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  scale_fill_manual(values = c(tree_cols)) +
  theme(legend.position = "none") +
  xlab("Whole genome average GC content (%) [binwidth = 0.1]") +
  ylab("Number of genomes")


# get average GC per species for factoring the order of panel c:
species_average_gc <- whole_genome_gc %>%
  group_by(phylogroup) %>%
  summarise(mean_gc = mean(gc_perc), n = n())


##genome and gc:
genome_and_gc <- subset.data.frame(serratia_metadata,select=c("File_prefix","GC...."))
colnames(genome_and_gc) <- c("strain","whole_genome_gc")

#merge in fastbaps data into the whole_genome_gc dataframe:
whole_genome_gc <- whole_genome_gc %>%
  left_join(fastbaps_l3)

## get the average gc values and arrange the dataframe to this :
gc_fastbaps_cluster <- whole_genome_gc %>%
  group_by(cluster) %>%
  summarise(whole_genome_gc = mean(gc_perc), n = n(), sd = sd(gc_perc)) %>%
  arrange(-as.numeric(whole_genome_gc))



#### 3 - read in the gc details for each gene in the panaroo dataset: ####

#first untar:
untar(tarfile = "../../figshare_data/gc_content/panaroo_gc_details.tar.gz",
      exdir = "../../figshare_data/gc_content/")

gc_details <- read.csv("../../figshare_data/gc_content/panaroo_gc_details.csv", stringsAsFactors = F, quote = "", comment.char = "") %>%
  mutate(strain = gsub("#","_",strain))

#now remove file

file.remove("../../figshare_data/gc_content/panaroo_gc_details.csv")

## add in the fastani data (phylogroup) and the fastbaps data (lineage):
gc_details_w_metadata <- gc_details %>%
  left_join(fastani_95 , by = "strain") %>%
  left_join(fastbaps_l3 %>% rename(lineage = cluster), by = "strain")

### read in the core genes names:
core_genes <- read.csv("../../figshare_data/pangenome/genus_core_gene_names.csv",stringsAsFactors = F) %>%
  mutate(type = "core") %>%
  rename(group = Gene)

##merge this in - keep all!:
gc_details_with_type <- gc_details_w_metadata %>%
  full_join(core_genes, by = "group")   %>%
  mutate(type = ifelse(is.na(type), "non-core",type)) # replace all NA in "type" column with "non-core"


### group by both strain and type to get summarised mean values to plot:
GC_per_strain_un_normalised <- gc_details_with_type %>%
  group_by(strain,phylogroup,lineage,type) %>%
  summarise(GC_gene = mean(GC_gene), GC1 = mean(GC1), GC2 = mean(GC2), GC3 = mean(GC3), n = n())

#chceck numbers of genes:
# ggplot(subset.data.frame(GC_per_strain_un_normalised, type =="core"), aes(x = n, fill = as.factor(phylogroup))) + 
#   geom_histogram(binwidth=1) + 
#   xlim(1650,1700) +
#   scale_fill_manual(values = c(tree_cols))

#inflate:
GC_un_normalised_inflated <- GC_per_strain_un_normalised %>%
  tidyr::gather(measurement, gc_perc, GC1:GC3 )

#inflate for entire dataset - not just the mean!
GC_un_normalised_inflated_whole <- gc_details_with_type %>%
  tidyr::gather(measurement, gc_perc, GC1:GC3 )

### plot the GC as un-normalised data:
#first factor the data so that it is in ascending order of average GC content:


#now factor:
GC_un_normalised_inflated_whole$lineage <- factor(GC_un_normalised_inflated_whole$lineage, levels = c(gc_fastbaps_cluster$cluster))


#ridges_plot <- 
p3 <- ggplot(GC_un_normalised_inflated_whole, aes(gc_perc, y = as.factor(lineage), fill = as.factor(phylogroup))) +
  #geom_histogram(binwidth = 0.1) +
  ggridges::geom_density_ridges(alpha = 0.8, aes(weight=length,height=..density..),
                      stat = "density") +
  scale_fill_manual(values =c(tree_cols)) +
  theme(legend.position = "none") +
  xlab("GC (%)") +
  ylab("Lineage") +
  xlim(25,100) +
  facet_grid(cols = vars(GC_un_normalised_inflated_whole$measurement),
             rows = vars(GC_un_normalised_inflated_whole$type),
             scales = "free") +
  theme_bw() +
  theme(legend.position = "none") 


##plot p1, p2 and p3:
library(gridExtra)
library(grid)
library(ggpubr)


library(cowplot)

svg("plot_3.svg", width = 9.843, height = 13.7795)
ggdraw() +
  draw_plot(p1, x = 0, y = .8, width = .5, height = .2) +
  draw_plot(p2,x = 0.5, y = 0.8, width = 0.5, height = .2) +
  draw_plot(p3,x = 0, y = 0.4, width = 1, height = .4) 
  #draw_plot(ridges_plot,x = 0, y = 0, width = 1, height = .5) +
dev.off()






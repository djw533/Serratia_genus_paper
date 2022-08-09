library(micro.gen.extra)


#### 1 - Plot tree #####
##first get tree:

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

### read in tree

serratia_tree <- ape::read.tree("../../figshare_data/tree/snp_sites_aln_panaroo.aln.treefile")

rooted_serratia_tree <- phytools::midpoint.root(serratia_tree)


#ape::write.tree(rooted_serratia_tree, file = "~/Documents/Serratia_genus_paper/figures_part_2/1_tree/rooted_serratia_tree.newick")

t1 <- ggtree::ggtree(rooted_serratia_tree, ladderize = T, right = T)



## fastani clusters:

fastani_data <- read.csv("../../figshare_data/clusters/95_to_99_ani_clusters.csv",
                         stringsAsFactors = F)

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

#for each level of fastbaps clustering:

trees <- list()

for (i in 1:length(colnames(fastbaps))) {

  level = colnames(fastbaps)[i]
  
  lvl_fastbaps <- fastbaps %>%
    select(!!level) %>%
    rename(cluster = !!level) %>%
    tibble::rownames_to_column(var = "strain") %>%
    mutate(cluster = as.factor(cluster))
  
  ##get the groups of the fastbaps:
  fastbaps_groups <- lvl_fastbaps %>%
    group_by(cluster) %>%
    summarise(strains = list(strain)) %>%
    pull(strains, name = cluster)
  
  #get the clade names:
  clade_labels <- micro.gen.extra::get_clade_nodes(genus_groups_tree, lvl_fastbaps)
  
  #draw the tree:
  highlighted_tree <- micro.gen.extra::add_clades_to_tree(t2,clade_labels, highlight_clades = TRUE)
  
  trees[[i]] <- highlighted_tree

}

library(ggpubr)
library(gridExtra)    


arranged_trees <- ggarrange(plotlist = trees,
          labels = c("a","b","c","d"),
          ncol = 2, nrow = 2)


svg("arranged_trees.svg", width = 5.77*1.75, height = 5.77*1.75)
arranged_trees
dev.off()





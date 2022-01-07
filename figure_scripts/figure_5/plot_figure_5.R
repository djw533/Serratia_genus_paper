
#### A. Libraries ####

# library(ggplot2)
# library(ggtree)
# library(phytools)
# library(tidyr)
# library(dplyr)
# library(tibble)
# library(castor)
# library(gggenes)
# library(purrr)

library(micro.gen.extra)
library(ggtree)


#### B. Functions ####

# #discard these
# test.overlap = function(start,end,all_starts,all_ends) {
#   rowSums(mapply(function(a,b) start > a & end < b, 
#                  all_starts, all_ends)) > 0
# }
# 
# single_test.overlap = function(start,end,start_2,end_2) {
#   return(start > start_2 & end < end_2)
# }



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


### just get the marcescens tree:
marc_tree <- castor::get_subtree_at_node(genus_groups_tree, 888 - length(genus_groups_tree$tip.label))$subtree

#get the marcescens only details
#first only want the baps clusters for the tips that are in the marcescens tree
marc_fastbap_df <- data.frame(Cluster = fastbaps[marc_tree$tip.label,], row.names = marc_tree$tip.label) %>%
  select(Cluster.Level.3) %>%
  tibble::rownames_to_column( var = "strain") %>%
  rename(cluster=  Cluster.Level.3)

marc_clade_labels <- micro.gen.extra::get_clade_nodes(marc_tree,marc_fastbap_df)

#### create t3 plot:
t3 <- ggtree::ggtree(marc_tree, ladderize = T, right = T) 

#add clades for the fastbaps clusters to the marcescens tree
t3_w_clades <- micro.gen.extra::add_clades_to_tree(t3,marc_clade_labels)

#get groups:
marc_groups_tree <- groupOTU(marc_tree, genus_groups) #, overlap='abandon',connect = T)



#### 2. Take copies of all complete loci between tRNA-Pro(ggg) and tRNA-Ser(tgt) ####

# first untar the directory:

untar(tarfile = "../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs.tar.gz",
      exdir = "../../figshare_data/tetrathionate_gentisate/")

cluster_names <- list.files("../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs/", pattern = "gff")

#### 3. Load in gggenes of the data ####

gggenes_df <- gggenes_df_from_gff_dir("../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs/")

#### 4. Load in panaroo data ####

genes_and_groups <- read.csv("../../figshare_data/pangenome/gene_presence_absence_w_header.csv",
                             header = T,
                             stringsAsFactors = F,
                             check.names = F,
                             comment.char = "")

pivotted <- subset.data.frame(tidyr::pivot_longer(genes_and_groups, names_to = "strain", values_to = "gene", cols = `26968_7_100`:SJC1071),
                              gene != "")

#drop the strain:
#pivotted <- pivotted[ , !(names(pivotted) %in% c("strain"))]


#### 5. Then the tRNA data ####

tRNAs_and_panaroo <- read.csv("../../figshare_data/tRNAs/tRNAs.csv",
                  header = F) %>%
  select(V10,V13,V1) %>%
  rename(gene = V10, Gene = V13,strain = V1) %>%
  mutate(gene = gsub("ID=","",gene), Gene = gsub("product=","",Gene), strain = gsub("#","_",strain)) %>%
  bind_rows(pivotted)


#### 6. Merge tRNA/panaroo data with the gggenes data ####

gggenes_w_tRNAs_panaroo <- gggenes_df %>%
  left_join(tRNAs_and_panaroo, by = "gene") 

#### 7. Then phaster data ####

#first just summary

phaster_data <- read.csv("../../figshare_data/tetrathionate_gentisate/phaster_summary.csv")


phaster_data <- phaster_data %>% # add in 
  mutate(common_phage = unlist(purrr::map(phaster_data$most_common_phage, ~ unlist(strsplit(.x,"[(]")[1])[1])))

#read in cluster name and contig name:
phaster_data_w_strain <- read.csv("../../figshare_data/tetrathionate_gentisate/cluster_and_contig_names.csv", 
                                  header = F, 
                                  comment.char = "") %>% 
  setNames(c("strain","contig")) %>%
  mutate(strain = gsub("#","_",gsub(".fasta","",(gsub("_cluster_1","",gsub("../1_fasta/","",strain))))))  %>% # an awful lot of find and replace.....
  left_join(phaster_data)  %>% # join the phaster_data in now 
  left_join(gggenes_w_tRNAs_panaroo, by = "strain") %>%
  mutate(phage_gene = ifelse(phage_start < start & phage_end > end, T, F))


genes_w_phage_data <- phaster_data_w_strain %>%
  select(strain,gene,start,end,direction,type,Gene,phage_gene,common_phage) %>%
  distinct() # remove any duplicates 


#then the details
phaster_data_details <- read.csv("../../figshare_data/tetrathionate_gentisate/parsed_phaster_details.csv")

phaster_data_details_w_strain <- read.csv("../../figshare_data/tetrathionate_gentisate/cluster_and_contig_names.csv", 
                                          header = F, 
                                          comment.char = "") %>% 
  setNames(c("strain","contig")) %>%
  mutate(strain = gsub("#","_",gsub(".fasta","",(gsub("_cluster_1","",gsub("../1_fasta/","",strain))))))  %>% # an awful lot of find and replace.....
  right_join(phaster_data_details, by = "contig") %>%
  right_join(gggenes_w_tRNAs_panaroo, by = "strain") %>% # merge in the gggenes data
  mutate(matched_gene = ifelse(direction.x == direction.y & start.x == start.y & end.x == end.y, TRUE,FALSE)) %>% # do the gene coordinates match?
  filter(matched_gene == TRUE) %>% # subset to only get matched genes
  select(gene,homologue,details) %>%  # and clean up by only keeping columns of interest
  distinct()  # make sure no duplicates

gggenes_w_tRNAs_panaroo_and_phaster  <- genes_w_phage_data %>%
  left_join(phaster_data_details_w_strain, by = "gene") %>% #merge into genes_w_phage_data
  select(gene,homologue,details,phage_gene,common_phage) %>%
  distinct() %>%
  right_join(gggenes_w_tRNAs_panaroo)

#### 8. Read in manually annotated pathway data ####

manual_annotation_file <- read.csv("../../figshare_data/tetrathionate_gentisate/manually_annotated_tet_gentisate_pathways.csv") %>%
  #select(Gene,simple_name,colour_to_plot,description_to_plot) %>%
  distinct() # remove duplicates

gggenes_w_tRNAs_panaroo_phaster_and_manual_annotation <- gggenes_w_tRNAs_panaroo_and_phaster %>%
  left_join(manual_annotation_file) 

#### now can work with this to create different dfs to plot the genes with:

genes_to_plot <- gggenes_w_tRNAs_panaroo_phaster_and_manual_annotation %>%
  select(strain,start,end,direction,type,Gene,colour_to_plot,phage_gene,common_phage,gene,simple_name) %>%
  distinct() %>% #remove duplicates 
  mutate(type = ifelse(is.na(phage_gene), type, 
                       ifelse(phage_gene == FALSE, type, "Prophage")))  %>% # set whether it's a prophage or not
  mutate(synteny_colour = ifelse(type == "Prophage" , common_phage, colour_to_plot))
  
  

#### 9. Plot all clusters against the tree for the supplementary plots ####

# get the tips in the tree that we don't have full data for

tips_w_no_data <- data.frame(table (append(unique(genes_to_plot$strain),marc_tree$tip.label))) %>%
  rename(strain = Var1, num = Freq) %>%
  filter(num == 1 ) %>%
  select(strain) %>%
  mutate(label2 = "#")
  

t3_highlighted <- t3 %<+% tips_w_no_data + geom_tiplab(aes(label = label2), color = "red", size = 16)

t3_highlighted_w_clades <- add_clades_to_tree(t3_highlighted,marc_clade_labels,T,T)


#make some colours:
gene_colours <- micro.gen.extra::random_n_colours(length(unique(genes_to_plot$colour_to_plot)),T)
names(gene_colours) <- unique(genes_to_plot$colour_to_plot)


#set the feature colours:
feature_colours <- c("CDS" = "#000000",
                     "ncRNA" =  "#00ff00",
                     "tRNA" = "#ff0000",
                     "Prophage" = "#0000ff")


p1 <- facet_plot(t3_highlighted_w_clades, panel='clusters',
                 mapping = aes(xmin = start, xmax = end, fill = colour_to_plot, color = type, x=start, forward = direction, y=y),
                 data=genes_to_plot, geom=gggenes::geom_gene_arrow, arrow_body_height = grid::unit(1.8, "mm"), arrowhead_height = grid::unit(2.5, "mm"), arrowhead_width = grid::unit(2.5,"mm"),size = 0.2) +
  #theme(panel.grid.minor = element_line(colour="grey", size=0.1)) +
  scale_colour_manual(values = c(feature_colours)) +
  scale_fill_manual(values = c(gene_colours)) + 
  #scale_y_continuous(minor_breaks = seq(1, 550, 2)) +
  #geom_treescale(width = 2, offset = -1.2) +
  theme(legend.position = "bottom")

ggsave(plot = p1,
       filename = "clusters_w_marc_tree_common_genes.pdf", dpi=300, width=30, height=50, units=c("in"), limitsize = F)



#### 10. Remove tips lacking data from the tree ####


marc_tree_dropped_tips <- micro.gen.extra::remove_tips_from_tree(marc_tree, tips_w_no_data$strain)

t4 <- ggtree(marc_tree_dropped_tips, right = TRUE, ladderize = T)

# now add the fastbaps data:
#first only want the baps clusters for the tips that are in this pruned tree
selected_marc_fastbap_df <- data.frame(Cluster = fastbaps[marc_tree_dropped_tips$tip.label,], row.names = marc_tree_dropped_tips$tip.label) %>%
  select(Cluster.Level.3) %>%
  tibble::rownames_to_column( var = "strain") %>%
  rename(cluster=  Cluster.Level.3)

#get the clade labels and the nodes to annotate these at
selected_marc_clade_labels <- micro.gen.extra::get_clade_nodes(marc_tree_dropped_tips,selected_marc_fastbap_df)

#add the clade labels
marc_tree_dropped_tips_w_clades <- micro.gen.extra::add_clades_to_tree(t4,selected_marc_clade_labels)

# draw with the clusters:
p2 <- facet_plot(marc_tree_dropped_tips_w_clades, panel='clusters',
                 mapping = aes(xmin = start, xmax = end, fill = colour_to_plot, color = type, x=start, forward = direction, y=y),
                 data=genes_to_plot, gggenes::geom_gene_arrow, arrow_body_height = grid::unit(1.8, "mm"), arrowhead_height = grid::unit(2.5, "mm"), arrowhead_width = grid::unit(2.5,"mm"),size = 0.2) +
  #theme(panel.grid.minor = element_line(colour="grey", size=0.1)) +
  scale_colour_manual(values = c(feature_colours)) +
  scale_fill_manual(values = c(gene_colours)) + 
  #scale_y_continuous(minor_breaks = seq(1, 550, 2)) +
  #geom_treescale(width = 2, offset = -1.2) +
  theme(legend.position = "bottom")

ggsave(plot = p2,
       filename = "clusters_w_selected_marc_tree_common_genes.pdf", dpi=300, width=30, height=50, units=c("in"), limitsize = F)


#### 11. Cluster using sets of genes present ####

#first fix some of the names of the genes that have incorrectly been termed as phage genes:

genes_to_plot_corrected  <- genes_to_plot %>%
  mutate(synteny_colour = ifelse(colour_to_plot == "cAMP phosphodiesterase" & ! is.na(colour_to_plot), colour_to_plot, synteny_colour)) %>%
  mutate(synteny_colour = ifelse(simple_name == "dinI" & ! is.na(simple_name), colour_to_plot, synteny_colour)) %>%
  mutate(synteny_colour = ifelse(colour_to_plot == "tetrathionate reduction I (to thiosulfate)" & ! is.na(colour_to_plot), colour_to_plot, synteny_colour)) %>%
  mutate(more_names = ifelse(synteny_colour == "tetrathionate reduction I (to thiosulfate)", simple_name, synteny_colour)) %>%
  mutate(more_names = ifelse(synteny_colour == "gentisate degradation I", simple_name, more_names)) %>%
  #mutate(more_names = ifelse(colour_to_plot == "cAMP phosphodiesterase", synteny_colour, more_names)) %>%
  mutate(more_names = ifelse(synteny_colour == "cAMP phosphodiesterase", simple_name, more_names)) 
  
  
  
#then filter this so can be used as input into the plot_synteny call
more_complex_plotting <- genes_to_plot_corrected %>%
  select(gene,more_names) %>%
  distinct() %>%
  rename(Gene = more_names)


#create a list of the "unique" patterns of blocks of genes in a 5' -> 3' direction .....
genes_between_tRNAs_clusters <- genes_to_plot_corrected %>%
  select(strain,start,synteny_colour) %>% 
  distinct() %>%
  filter(! is.na(strain)) %>%
  filter(! is.na(synteny_colour)) %>%
  group_by(strain) %>%
  arrange(start) %>%
  select(strain, synteny_colour) %>%
  distinct() %>%
  summarise(gene_groups = paste(synteny_colour, collapse = ",")) %>%
  arrange(gene_groups)


### .... then use this to create a heatmap of presence/absence of the order of possession of each of these different "unique" patterns of gene blocks

gene_sets_presence_absence <- genes_to_plot %>%
  select(strain,start,synteny_colour) %>% 
  distinct() %>%
  filter(! is.na(strain)) %>%
  filter(! is.na(synteny_colour)) %>%
  group_by(strain) %>%
  arrange(start) %>%
  select(strain, synteny_colour) %>%
  distinct() %>%
  arrange(synteny_colour) %>%
  mutate(presence = "Y") %>%
  tidyr::pivot_wider(names_from = synteny_colour, values_from = presence,values_fill = NA) %>%
  tibble::column_to_rownames(var = "strain")


#change "Y" to the colname, so that a coloured heatmap can be produced, where presence/absence correlates with the colours of the genes in the figures
temp_df <- gene_sets_presence_absence

for (i in 1:ncol(temp_df)) {
  text_string <- colnames(temp_df)[i]
  temp_df[,i] <- gsub("Y",text_string,temp_df[,i])
}


###plot against the tree:

##order the colnames into desired order, and in multiple dataframes:

#first main metabolic groups:
metabolic_groups <- temp_df %>%
  select("tetrathionate reduction I (to thiosulfate)","gentisate degradation I","3-chlorobenzoate degradation III (via gentisate)")

#then phages
phages <-  temp_df %>%
  select("PHAGE_Acinet_YMC11/11/R3177_NC_041866","PHAGE_Entero_mEp460_NC_019716","PHAGE_Erwini_phiEt88_NC_015295","PHAGE_Klebsi_phiKO2_NC_005857","PHAGE_Salmon_118970_sal3_NC_031940")

#then extras
extras <- temp_df %>%
  select("Acyltransferase 3","cAMP phosphodiesterase","Colicin V secretion")


#set the colours:
##read in the colours for these:

more_complex_plotting_gene_colours <- read.csv("colours.csv",
                                          header = F) %>%
  rename(gene_name = V1,colour = V2) %>%
  pull(colour, name = gene_name)


t4 %>% 
  #metabolic groups
  gheatmap(metabolic_groups, color = NULL,
                colnames_position = "top",
           width = 0.1 * ncol(metabolic_groups),
                colnames_angle = 60,
                hjust = 0, font.size = 3) %>%
  #extras
  gheatmap(extras, color = NULL,
           colnames_position = "top",
           colnames_angle = 60,
           width = 0.1 * ncol(extras),
           offset = 0.05,
           hjust = 0, font.size = 3) %>%
  #phages
  gheatmap(phages, color = NULL,
           width = 0.1 * ncol(phages),
           colnames_position = "top",
           colnames_angle = 60,
           offset = 0.12,
           hjust = 0, font.size = 3) + 
  scale_fill_manual(values = c(more_complex_plotting_gene_colours)) + 
  ylim(0,500)
  


##if want it in order??

left_phages <- temp_df %>%
  select("PHAGE_Acinet_YMC11/11/R3177_NC_041866","PHAGE_Entero_mEp460_NC_019716","PHAGE_Erwini_phiEt88_NC_015295","PHAGE_Salmon_118970_sal3_NC_031940")

position_one <- temp_df %>%
  select("gentisate degradation I", "Colicin V secretion") 

position_two <- temp_df %>% 
  select("tetrathionate reduction I (to thiosulfate)") 

position_three <- temp_df %>%
  select("cAMP phosphodiesterase","Acyltransferase 3")
  
right_phages <- phages <-  temp_df %>%
  select("PHAGE_Klebsi_phiKO2_NC_005857")

dfs_to_plot <- list(left_phages,position_one,position_two,position_three,right_phages)
  
temp_plot <- t4
tree_width <- max(t4$data$x)
offset_num = 0
col_width = 0.1 # is as proportion of the width of the tree!

for (df_name in dfs_to_plot) {
  
  width_num <- col_width * ncol(df_name)
  
  temp_plot <- temp_plot %>%
    gheatmap(df_name, color = NULL,
             width = width_num,
             colnames_position = "top",
             colnames_angle = 60,
             offset = offset_num,
             hjust = 0, font.size = 3)
  
  #set offset:
  offset_num = offset_num + (width_num * tree_width) + (col_width * tree_width)

}

p3 <- temp_plot + 
  scale_fill_manual(values = c(more_complex_plotting_gene_colours)) + 
  ylim(0,500)

ggsave(plot = p3,
       filename = "marc_tree_plasticity_zone_presence_absence.png", dpi = 600, height = 8, width = 16, limitsize = F)


#now create a table of these:
genes_clusters_table <- data.frame(table(genes_between_tRNAs_clusters$gene_groups)) %>%
  arrange(-Freq)

genes_clusters <- data.frame(gene_groups = genes_clusters_table$Var1,
                             cluster_num = rownames(genes_clusters_table)) %>%
  inner_join(genes_between_tRNAs_clusters) %>%
  mutate(cluster_num = as.character(cluster_num)) %>%
  tibble::column_to_rownames(var = "strain") %>%
  select(cluster_num)

#now plot on the marcescens tree:

#t4 %>% gheatmap(genes_clusters, color = NULL)

### and heatmap with various different columns:

genes_clusters_heatmap <- genes_clusters %>%
  tibble::rownames_to_column(var = "strain") %>% 
  mutate(Presence = "Y") %>% 
  tidyr::pivot_wider(names_from = cluster_num, values_from = Presence,values_fill = NA) %>%
  tibble::column_to_rownames(var = "strain")

#heatmap:
# 
# t4 %>% gheatmap(genes_clusters_heatmap,color = NULL,colnames_position = "top", colnames_angle = 90, hjust = 0 )  + 
#   ylim(0,450) 


### now collapse the tree:

clusters_and_strains <- genes_clusters %>%
  rename(Cluster = cluster_num) %>%
  tibble::rownames_to_column(var = "strain") %>%
  group_by(Cluster) %>%
  summarise(list_strains = list(strain),n = n()) 


#then use these to get a named vector set
clusters_and_strains_list <- clusters_and_strains$list_strains
names(clusters_and_strains_list) <- clusters_and_strains$Cluster


#work out which internal nodes should be collapsed (nodes should collapse if all descending nodes have the same gene organisation
# i.e. they are the same "cluster" in terms of synteny in the plasticity zone - defined above)
internal_nodes_and_subtree_tips <- subset.data.frame(ggtree::ggtree(marc_tree_dropped_tips)$data, isTip != TRUE, select = "node") %>%
  mutate(node = node - length(marc_tree_dropped_tips$tip.label)) %>%
  mutate(tips = purrr::map(node, ~ castor::get_subtree_at_node(marc_tree_dropped_tips, .x )$subtree$tip.label)) %>%
  # compare the tips in the subtrees with the strain names in each of the blast clusters
  mutate(cluster = purrr::map(tips, function(x) purrr::map(seq_along(clusters_and_strains_list), function(y)  ifelse(all(unlist(x) %in% unlist(clusters_and_strains_list[y])), names(clusters_and_strains_list)[y], NA)))) %>%
  tidyr::unnest(cols = c(cluster)) %>% # expand out as map returns a list
  filter(! is.na(cluster)) %>% # remove all NA from the above mapped mutate (2 lines above)
  # now need to try and get rid of nodes that are descendents of others that also are in the same cluster - so do similar to before
  mutate(is_subtree =   purrr::map(seq_along(tips), function(x) purrr::map(seq_along(tips), function(y)  ifelse(all(unlist(tips[x]) %in% unlist(tips[y])) & x != y , TRUE, FALSE)))) %>%
  mutate(is_subtree_2 = purrr::map(is_subtree, ~ ifelse(TRUE %in% .x, TRUE, FALSE))) %>%
  filter(is_subtree_2 == FALSE) # remove all FALSE from the above mapped mutate (2 lines above) %>%


#try collapsing:

collapsed_t4 <- collapse_tree(t4,
                              internal_nodes_and_subtree_tips$node + length(marc_tree_dropped_tips$tip.label),
                              "mixed",
                              FALSE)

##scaling....
temp_scaled <- t4

for ( node_num in (internal_nodes_and_subtree_tips$node + length(marc_tree_dropped_tips$tip.label))) {
  temp_scaled <- scaleClade(temp_scaled, node = node_num, scale = .001)
}

#set some random colours:
random_clusters_colours <- micro.gen.extra::random_n_colours(length(unique(clusters_and_strains$Cluster)),TRUE)
names(random_clusters_colours) <- unlist(unique(clusters_and_strains$Cluster))


collapsed_t4_w_colours <- collapse_tree(temp_scaled,
                                        internal_nodes_and_subtree_tips$node + length(marc_tree_dropped_tips$tip.label),
                                        "max",
                                        random_clusters_colours[unlist(internal_nodes_and_subtree_tips$cluster)])



# collapsed_t4_w_colours <- collapse_tree(t4,
#                                         internal_nodes_and_subtree_tips$node + length(marc_tree_dropped_tips$tip.label),
#                                         "max",
#                                         random_clusters_colours[unlist(internal_nodes_and_subtree_tips$cluster)])

##now colour the tips as points?

###add cluster number fo the df within the ggtree object created above:

collapsed_t4_w_colours$data  <- collapsed_t4_w_colours$data %>%
  left_join(genes_clusters %>% tibble::rownames_to_column(var = "label"), by = "label")

#and also get the node number for those nodes that were collapsed   
collapsed_t4_w_colours$data <- internal_nodes_and_subtree_tips %>%
  mutate(correct_node_num = node + length(marc_tree_dropped_tips$tip.label),
         cluster = as.character(cluster),
         n = as.integer(purrr::map(tips, ~ length(.x)))) %>%
  select(correct_node_num, cluster,n) %>%
  rename(node = correct_node_num) %>%
  right_join(collapsed_t4_w_colours$data, by = "node") 



collapsed_t4_w_colours_df <- genes_clusters %>%
  tibble::rownames_to_column(var = "label")


tips_coloured_collapsed_t4_w_colours <- collapsed_t4_w_colours + ggtree::geom_tippoint(aes(color = as.factor(cluster_num))) +
  ggtree::geom_tiplab(aes(label = as.factor(cluster_num))) +
  ggtree::geom_point2(aes(subset=(! is.na(cluster)), size = n), shape = 23, fill = "black", alpha = 0.5)  + 
  scale_color_manual(values = random_clusters_colours)


#get named sets of the clades to label: internal_nodes_and_subtree_tips$node + length(marc_tree_dropped_tips$tip.label

clade_labels_for_collapsed_tree <- internal_nodes_and_subtree_tips$node + length(marc_tree_dropped_tips$tip.label)
names(clade_labels_for_collapsed_tree) <- unlist(internal_nodes_and_subtree_tips$cluster)
  
  
with_labels <- add_clades_to_tree(drawn_tree = tips_coloured_collapsed_t4_w_colours, clade_labels = clade_labels_for_collapsed_tree,FALSE,F)

ggsave(plot = with_labels,
       filename = "collapsed_tree_condensed.svg", height = 20, width = 12 )

### 12 - plot synteny of these defined clusters ####

plot_order <- c(1:length(unique(genes_clusters$cluster_num)))

example_plots <- list()

for (plot in plot_order) {
  temp_df <- subset.data.frame(genes_clusters, cluster_num == plot) %>%
    tibble::rownames_to_column(var = "strain")
  example_plots <- append(example_plots, paste0(gsub("26968_7_","26968_7#",sample(temp_df$strain,1)),"_cluster_1_cluster_1_cluster_1.gff",sep=""))
}

names(example_plots) <- plot_order

##save this out

example_plots_df <- data.frame(cluster = names(unlist(example_plots)),
                               name = unlist(example_plots)) %>%
  write.csv(file = "example_plots_sets_of_genes.csv",
            row.names = F,
            quote = F)


#get gggenes df :
gggenes_input <- micro.gen.extra::gggenes_df_from_gff_list(paste("../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs/",example_plots, sep = "")) %>%
  left_join(more_complex_plotting %>% rename(colour_variable = Gene))

#create dnaseqs:

genoplot_data <- list()

genoplot_data$dna_seqs <- micro.gen.extra::create_dnaseqs(gggenes_df = gggenes_input,
                                                          set_gene_fill_colours = T,
                                                          gene_colours = more_complex_plotting_gene_colours)

##change the names of these:
names(genoplot_data$dna_seqs) <- names(example_plots)

#get fasta sequences:
fasta_sequences <- micro.gen.extra::fasta_from_gff_list(gff_list = paste("../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs/",example_plots, sep = ""),
                                     gff_names = names(example_plots),
                                     fasta_dir = "output_fasta_sequences",
                                     clean_up = FALSE)



#set the blast comparisons:
genoplot_data$comparisons <- micro.gen.extra::create_blast_comparisons(blast_order = gsub(".fasta","",list.files("output_fasta_sequences/")),
                                                                       fasta_dir = "output_fasta_sequences")


#plot the comparisons:
svg("synteny_of_clusters.svg", width = 12, height = 10,pointsize = 8)
genoPlotR::plot_gene_map(genoplot_data$dna_seqs, genoplot_data$comparisons,
                         scale=T,
                         legend = T,
                         arrow_head_len = 400,
                         annotation_height = 0.5)
dev.off()



#now plot a key of the genes (as arrows):

colours_to_plot <- data.frame(Gene = unique(more_complex_plotting$Gene),
                              value = sample.int(100, length(unique(more_complex_plotting$Gene))))

p4 <- ggplot(colours_to_plot) +
  geom_col(aes(x = Gene, y = value, fill = Gene)) +
  scale_fill_manual(values = c(more_complex_plotting_gene_colours))

l1 <- cowplot::ggdraw(cowplot::get_legend(p4))

ggsave(filename = "legend.svg",
       plot = l1,
       dpi = 600)



#### 13 - clean up ####
for (filename in list.files(path = "../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs",
                            full.names = TRUE)) {
  file.remove(filename)
}

#remove the directory:
unlink("../../figshare_data/tetrathionate_gentisate/plasticity_zone_gffs/",
       recursive = TRUE)



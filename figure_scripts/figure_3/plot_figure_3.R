

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



#### 2 now look at pathways: ####

#first untar:
untar(tarfile = "../../figshare_data/pathways/metabolic_pathways_rxns_with_genes_all_serratia.tar.gz",
      exdir = "../../figshare_data/pathways/")


pathways <- read.table("../../figshare_data/pathways/metabolic_pathways_rxns_with_genes_all_serratia.tsv",
                       header = T,
                       sep="\t",
                       quote = "",
                       comment.char = "",
                       stringsAsFactors = F) %>%
  filter(Reason.to.Keep != "PWY-HAS-NON-COMPUTATIONAL-EVIDENCE")

#now remove it:
file.remove("../../figshare_data/pathways/metabolic_pathways_rxns_with_genes_all_serratia.tsv")

##correct ones which were changed:
pathways$genome <- gsub("GCA_000465615_2","GCA_000465615_2_GS_De_Novo_Assembly_genomic",pathways$genome)
pathways$genome <- gsub("GCA_000468075_1","GCA_000468075_1_Serratia_fanticola_AU-P3_3__genomic",pathways$genome)
pathways$genome <- gsub("GCA_000477615_1","GCA_000477615_1_Serratia_fanticola_genomic",pathways$genome)
pathways$genome <- gsub("GCA_000805875_1","GCA_000805875_1_SmRM66262_v1_0_genomic",pathways$genome)
pathways$genome <- gsub("GCA_900175135_1","GCA_900175135_1_MFPA44A14-05_10_genomic",pathways$genome)


### now correct all the names back to what they should be:
correct_pathways <- data.frame(correct = genus_groups_tree$tip.label,
                         stringsAsFactors = F) %>%
  mutate(wrong = gsub("[.]","_",correct)) %>%
  left_join(data.frame(pathway_genome = unique(pathways$genome),stringsAsFactors = F), by = c("wrong" = "pathway_genome")) %>%
  right_join(pathways, by = c("wrong" = "genome"), keep = T)


# filter this to only keep the "correct" genome names (i.e formatted correctly)
pathways_and_genomes <- unique(subset.data.frame(correct_pathways, select = c("Pathway.Frame.id","correct"))) %>%
  rename("genome" = "correct")

#factor:
pathways_and_genomes <- within(pathways_and_genomes, Pathway.Frame.id <- factor(Pathway.Frame.id, levels=names(sort(table(Pathway.Frame.id),  decreasing=TRUE))))


### create a bar chart of the number of pathways into a line graph with continuous variables to go underneath the gheatmap tree:
pathway_numbers <- data.frame(table(unlist(pathways_and_genomes$Pathway.Frame.id))) %>%
  tibble::rownames_to_column(var = "number") %>%
  mutate(number = as.numeric(number)) %>%
  rename("pwy" = "Var1", "num" = "number", "num_pwys" = "Freq") %>%
  mutate(proportion = num / max(num) )


##plot the number of pathways / genome as line graph:
number_of_pathways_plot <- ggplot(pathway_numbers, aes(num_pwys, num)) + geom_line() +
  theme_classic() +
  xlab("Number of pathways") +
  ylab("Number of genomes") +
  scale_y_continuous(breaks = c(0,100,200,300,400,500,600,700)) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(num,0)),alpha = 0.5)

#save this to put underneath the gheatmap
ggsave(plot = number_of_pathways_plot,
       filename = "number_of_pathways_plot.svg",
       width = 12,
       height = 6)



#### now do presence/absence matrix of the different pathways to create a gheatmap:
df2 <- pathways_and_genomes %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = Pathway.Frame.id, values_from = present) %>%
  tibble::column_to_rownames(var = "genome")


### now sort colnames so that they are in descending order of prevalence:
orders <- subset.data.frame(as.data.frame(table(unlist(pathways_and_genomes$Pathway.Frame.id))), Freq > 0)

orders_to_sort <- as.character(unlist(orders$Var1))

df2<- df2[, orders_to_sort]

##plot pathways against the tree

pathways_gheatmap_plot <- new_tree %>% ggtree::gheatmap(df2, color = NULL,
                colnames_position = "top",
                colnames_offset_y = 100,
                width = 2,
                colnames_angle = 90,
                offset = 0.1,
                hjust = 0) +
  scale_fill_gradient(low = "#7EACE6", high = "#7EACE6", na.value = NA)+
  ylim(0,max(new_tree$data$y)+20) +
  theme(legend.position = "none")

ggsave(plot = pathways_gheatmap_plot,
       filename = "pathways_against_serratia_tree.png", 
       dpi = 600, 
       units = c("cm"), 
       width = 14*2, 
       height = 8.5*2)


#### 3 - create an upset for the supplementary plot ####

position_number <- function(x, levels, labels) { # function to count within the upsetR df generation below (taken from micro.gen.extra::create_upset_dfs)
  #and for x:
  if ( ! is.numeric(x)) {
    stop("x must be numeric")
  }
  #function
  for (i in 1:(length(levels)-1)) {
    min = levels[i]
    max = levels[i+1]
    if (x == 0) {
      label = "Absent"
      return(label)
    }
    if (x >= min & x < max) {
      label = labels[i]
      return(label)
    }
    if ( x == levels[length(levels)]) { # if x is at the top limits of the levels - then but in the top level ( fix for above x < max check)
      label = labels[length(labels)]
      return(label)
    }
  }
}


#create upset dataframes
levels = c(0,0.15,0.95,1)
labels = c("Cloud","Shell","Core")

pathway_occurence_df <- pathways %>%
  select(genome,Pathway.Frame.id) %>%
  rename(strain = genome, pathway = Pathway.Frame.id) %>%
  distinct() %>%
  mutate(presence = 1) %>% # add a nominal value for each present pathway in each genome - i.e 1 if present
  tidyr::pivot_wider(names_from = "pathway", values_from = "presence", values_fill = 0) %>%
  tidyr::pivot_longer(cols = !strain,names_to = "pathway", values_to = "presence") %>%
  left_join(fastbaps_l3 %>% mutate(strain = gsub("[.]","_",strain)), by = "strain") %>% # add groups in (using lineages)
  group_by(cluster,pathway) %>%
  summarise(perc_in_group = sum(presence) / length(strain)) %>% # get incidence of the gene within each lineage
  #filter(perc_in_group > 0) %>% # remove genes not seen in individual groups of strains
  mutate(occurence_label = purrr::map(perc_in_group, 
                                      ~ position_number(.x, 
                                                        levels, 
                                                        labels)
                                      )
         ) ## map over incidence and assign group for the gene

#create the three upset dfs from this:
upset_dfs <- vector(mode = "list", length = length(labels))
names(upset_dfs) <- labels

for (i in 1:length(labels)) {
  
  l = labels[i]
  
  upset_df <- pathway_occurence_df %>%
    mutate(perc_in_group = ifelse(occurence_label == l, 1, 0 )) %>%
    select(pathway,cluster,perc_in_group) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = perc_in_group)
  
  #change the colnames back:
  colnames(upset_df) <- gsub(glue::glue("{l}[.]"),"",colnames(upset_df))
  
  upset_dfs[i] <-list(as.data.frame(upset_df))
  
}



#get lineage order for plot:
lineage_order <- t2$data %>%
  filter(isTip== TRUE) %>%
  rename(strain = label) %>%
  left_join(fastbaps_l3, by = "strain") %>%
  group_by(cluster) %>%
  summarise(max_y = max(y)) %>%
  arrange(max_y) %>%
  mutate(cluster = as.character(cluster)) %>%
  pull(cluster)

###get metadata for colouring:
### get the species for each fastbaps cluster (lineage):
metadata_for_upset <- fastani_data %>%
  rename(strain = Name, phylogroup = X95) %>%
  mutate(strain = gsub("#","_",strain)) %>%
  left_join(fastbaps_l3, by = "strain") %>%
  select(cluster,phylogroup) %>%
  unique() %>%
  rename(sets = cluster, species = phylogroup) %>%
  mutate(species = gsub("11","S.rubidaea", species),
         species = gsub("12","S.rubidaea_like", species),
         species = gsub("13","S.marcescens_like", species),
         species = gsub("10","S.odorifera", species),
         species = gsub("4","S.liquefaciens", species),
         species = gsub("6","S.ficaria", species),
         species = gsub("3","S.fonticola", species),
         species = gsub("5","S.entomophila", species),
         species = gsub("1","S.proteamaculans_subs.quinovorans", species),
         species = gsub("9","S.proteamaculans", species),
         species = gsub("8","S.plymuthica", species),
         species = gsub("7","S.grimesii", species),
         species = gsub("2","S.marcescens", species))


#plot UpSet
svg("Core_and_softcore_pathways.svg", width = 10, height = 8)
UpSetR::upset(upset_dfs$Core, sets = c(lineage_order), order.by = "freq", keep.order = T,
              set.metadata = list(data = metadata_for_upset, plots = list(
                list(type = "matrix_rows",
                     column = "species",
                     colors = c(S.ficaria = "#f28e2b",
                                S.fonticola = "#e15759",
                                S.liquefaciens = "#76b7b2",
                                S.marcescens = "#59a14f",
                                S.plymuthica = "#edc948",
                                S.proteamaculans = "#b07aa1",
                                S.rubidaea = "#ff9da7",
                                S.grimesii = "#9c755f",
                                S.entomophila = "#4e79a7",
                                S.proteamaculans_subs.quinovorans = "#bab0ac",
                                S.odorifera = "#76b7b2",
                                S.marcescens_like = "#B3E2CD",
                                S.rubidaea_like = "#a0a0a0"),
                     alpha = 0.65))))
dev.off()
          

#### 4 - plot the selected pathways in marcescens for panel b ####


#get marcescnes tree first:
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


five_pathways <- subset.data.frame(df2, select = c("QUINATEDEG-PWY","PWY-2503","PWY0-1465","PWY-6223","PWY-5358"))
##re-set colnames:
colnames(five_pathways) <- c("Quinate","Benzoate","D-malate","Gentisate","Tetrathionate")
#

##convert values in five_pathways to factors:
five_pathways[five_pathways == 1] <- "Y"


#read in metadata for extra columns in the plot

serratia_metadata <- read.csv("../../figshare_data/metadata/serratia_metadata.csv", header = T,
                              quote = "",
                              stringsAsFactors = F,
                              comment.char = "")


biotype_df <- serratia_metadata %>%
  select(File_prefix,new_biotype) %>%
  rename(Biotype = new_biotype) %>%
  mutate(Biotype = ifelse(Biotype == "", NA, Biotype)) %>%
  tibble::column_to_rownames(var = "File_prefix")

source_df <- serratia_metadata %>%
  select(File_prefix,Generic_host) %>%
  rename(Source = Generic_host) %>%
  mutate(Source = ifelse(Source == "", NA, Source)) %>%
  tibble::column_to_rownames(var = "File_prefix")

prodigiosin_df <- read.csv("../../figshare_data/prodigiosin/pig_cluster_hamburger_cluster_statistics.csv") %>%
  select(gene_cluster) %>%
  distinct() %>%
  mutate(gene_cluster = gsub("#","_",gsub("_cluster_1","",gene_cluster))) %>%
  mutate(Prodigiosin = "R") %>%
  tibble::column_to_rownames(var = "gene_cluster")


#colours
cols <- c("Serratia entomophila" = "#4e79a7",
          "Serratia ficaria" = "#f28e2b",
          "Serratia fonticola" = "#e15759",
          "Serratia liquefaciens" = "#76b7b2",
          "Serratia marcescens" = "#59a14f",
          "Serratia plymuthica" = "#edc948",
          "Serratia proteamaculans" = "#b07aa1",
          "Serratia rubidaea" = "#ff9da7",
          "Serratia grimesii" = "#9c755f",
          "Serratia proteamaculans subs. quinovorans" = "#bab0ac",
          "Serratia proteamaculans\nsubs. quinovorans" = "#bab0ac",
          "Serratia sp." = "#d3d3d3",
          "Unknown" = "grey",
          "uncultured Serratia" = "#a0a0a0",
          "Rhodotorula mucilaginosa/Serratia"="#BF5B17",
          "Rhodotorula mucilaginosa/\nSerratia"="#BF5B17",
          "Serratia aquatilis"="#FF7F00",
          "Serratia myotis"="#E5C494",
          "Serratia nematodiphila"="#FFED6F",
          "Enterobacter liquefaciens"="#F781BF",
          "Serratia quinivorans"="#D95F02",
          "Serratia Serratia"="#CAB2D6",
          "Serratia sp.XT-30"="#E6F5C9",
          "Serratia symbiont"="#386CB0",
          "Serratia oryzae" = "#a0a0a0",
          "Serratia odorifera"="#F781BF",
          "Serratia symbiotica"="#B3E2CD",
          "Serratia ureilytica"="#A6CEE3",
          "Serratia vespertilionis"="#FFFFCC",
          "1" = "#59a14f",
          "2" = "#76b7b2",
          "3" = "#f28e2b",
          "4" = "#e15759",
          "5" = "#4e79a7",
          "6" = "#bab0ac",
          "7" = "#b07aa1",
          "8" = "#edc948",
          "9" = "#9c755f",
          "10" = "#ff9da7",
          "11" = "#a0a0a0",
          "12" = "#B3E2CD",
          "13" = "#76b7b2",
          "14" = "black",
          "0" = "black",
          "Grimont collection (Institut Pasteur)" = "#59a14f", 
          "Study w/ < 10 genomes" = "#76b7b2",
          "Chan et al; 2016" = "#f28e2b",
          "Roach et al; 2015" = "#e15759",
          "Moradigaravand et al; 2015 (BSAC)" = "#4e79a7",
          "Hemarajata et al; 2018" = "#bab0ac",
          "Anderson et al; 2017" = "#b07aa1",
          "Zhang et al; 2018" = "#edc948",
          "Stoesser et al; 2017" = "#9c755f",
          "Chen et al; 2017" = "#ff9da7",
          "Raymann et al; 2018" = "#a0a0a0",
          "Weingarten et al; 2018" = "#B3E2CD",
          "Martineau et al; 2018" = "#76b7b2",
          "Parks et al; 2015" = "#A6CEE3",
          "NCTC 3000" = "#FFFFCC",
          "UK hospitals" = "#F781BF",
          "Plant" = "#59a14f",
          "N/A" = "white",
          "Costelytra zealandica" = "#edc948",
          "Animal" = "#76b7b2",
          "Insect" = "#f28e2b",
          "Water" = "#e15759",
          "Dairy" = "#4e79a7",
          "Human" = "#9c755f",
          "Environmental" = "#B3E2CD",
          "Rhizosphere" = "#76b7b2",
          "Fish" = "#A6CEE3",
          "Crustacean" = "#FFFFCC",
          "Food" = "#e866b2",
          "Mammal" = "red",
          "Bird" = "purple",
          "Hospital" = "pink",
          "Air" = "brown",
          "AA" = "grey",
          "Y" = "#C9AB00",
          "A1a?" = "#e6194b",
          "A2" = "#ffe119",
          "A1b" = "#3cb44b",
          "A2a" = "#4363d8",
          "A4a" = "#f58231",
          "A5" = "#911eb4",
          "A6a" = "#46f0f0",
          "TCT" = "#808080",
          "R" = "#e15759")



selected_pathways_plot <- t3_w_clades %>% 
  ggtree::gheatmap(biotype_df, 
           color = NULL, 
           colnames_position = "top",
           colnames_offset_y = 0,
           colnames_angle = 90,
           hjust = 0,
           width = 0.25) %>%
  ggtree::gheatmap(source_df, 
           color = NULL, 
           colnames_position = "top",
           colnames_offset_y = 0,
           colnames_angle = 90,
           hjust = 0,
           offset = 0.01 + (max(t3_w_clades$data$x) * 0.25), 
           width = 0.25) %>% 
  ggtree::gheatmap(five_pathways, color = NULL,
           colnames_position = "top",
           colnames_offset_y = 0,
           width = 1.5,
           colnames_angle = 90,
           offset = 0.02 + (max(t3_w_clades$data$x) * 0.5),
           hjust = 0) %>%
  ggtree::gheatmap(prodigiosin_df, color = NULL,
           colnames_position = "top",
           colnames_offset_y = 0,
           width = 0.25,
           colnames_angle = 90,
           offset = 0.03 + (max(t3_w_clades$data$x) * 2) ,
           hjust = 0) +
  ylim(0,max(t3_w_clades$data$y) + 50) + 
  scale_fill_manual(values = c(cols),na.value = NA) +
  theme(legend.position = "none")

ggsave("selected_pathways.svg", height = 3, width = 8)







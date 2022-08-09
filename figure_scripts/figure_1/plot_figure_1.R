library(micro.gen.extra)

#also requires phytools

## colours:
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
           S.odorifera = "#D7B5A6",
           S.marcescens_like = "#B3E2CD",
           S.rubidaea_like = "#a0a0a0")



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
          
          "United Kingdom"="#cfd9ff",
          "USA"="#ffd8b1",
          "Spain"="#ffe119",
          "Canada"="#4363d8",
          "South Africa"="#f58231",
          "Italy"="#911eb4",
          "Greece"="#46f0f0",
          "Singapore"="#f032e6",
          "Denmark"="#bcf60c",
          "USSR"="#fabebe",
          "Portugal"="#008080",
          "Netherlands"="#e6beff",
          "Malaysia"="#9a6324",
          "Taiwan"="#fffac8",
          "Sweden"="#800000",
          "India"="#aaffc3",
          "Germany"="#808000",
          "Russia"="#3cb44b",
          "China"="#000075",
          "France"="#808080",
          "Argentina"="#1b4d29",
          "Australia"="#665900",
          "South Korea"="#ffa15e",
          "Mexico"="#e6194b",
          "Switzerland"="#594955",
          "Brazil"="#806568",
          "Venezuela"="#3c0d4a",
          "Japan"="#262626",
          "Romania"="#d8e3bc",
          "New Zealand" = "#e15759",
          "Sicily" = "#f28e2b",
          "Tunisia" = "#b07aa1")



serratia_tree <- ggtree::read.tree("../../figshare_data/tree/snp_sites_aln_panaroo.aln.treefile")
rooted_serratia_tree <- phytools::midpoint.root(serratia_tree)
ape::write.tree(rooted_serratia_tree, file = "rooted_serratia_tree.newick")

t1 <- ggtree::ggtree(rooted_serratia_tree, ladderize = T, right = T)


#add bootstraps:
t1$data <- t1$data %>%
  mutate(bootstrap = as.numeric(ifelse(isTip == FALSE & label != "Root", label, NA))) %>%
  mutate(bootstrap_class = cut(bootstrap, breaks = c(0,80,90,95,99,100),
                               labels = c("<= 80","<= 90","<= 95","<= 99","100"))) %>%
  mutate(num_descending_tips = as.numeric(purrr::map2(node, bootstrap, ~ ifelse(!is.na(.y), length(castor::get_subtree_at_node(rooted_serratia_tree, 
                                                                              as.numeric(.x) - length(rooted_serratia_tree$tip.label))$subtree$tip.label), NA))))




view <- t1$data

t1 + ggtree::geom_point2(aes(subset=!isTip & ! is.na(bootstrap) & bootstrap >= 80 & num_descending_tips > 20, color = bootstrap_class), na.rm = T)

t1 + ggtree::geom_point2(aes(subset=!isTip & ! is.na(bootstrap), color = bootstrap), na.rm = T)


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

#get df of fastabap and fastani:
phylogroups_and_lineages <- fastbaps_l3 %>%
  left_join(fastani_data %>%
              select(Name, X95) %>%
              rename(strain = Name, phylogroup = X95) %>%
              mutate(strain = gsub("#","_",strain)),
            by = "strain")

##get the groups of the fastbaps:

fastbaps_groups <- fastbaps_l3 %>%
  group_by(cluster) %>%
  summarise(strains = list(strain)) %>%
  pull(strains, name = cluster)

#get the clade names:
clade_labels <- micro.gen.extra::get_clade_nodes(genus_groups_tree, fastbaps_l3)

#draw the tree:
new_tree <- micro.gen.extra::add_clades_to_tree(t2,clade_labels)

###add bootstraps to this

new_tree$data <- new_tree$data %>%
  #add nodes corresponding to each lineage:
  left_join(data.frame(node = unlist(clade_labels), ancestral_to_lineage = names(clade_labels)), by = "node") %>%
  #add lineage and fastani data to the df
  left_join(phylogroups_and_lineages %>% rename(label = strain, lineage = cluster), by = "label") %>%
  #sort the bootstrap values
  mutate(bootstrap = as.numeric(ifelse(isTip == FALSE & label != "Root", label, NA))) %>%
  mutate(bootstrap_class = cut(bootstrap, breaks = c(0,80,90,95,99,100),
                               labels = c("<= 80","<= 90","<= 95","<= 99","100"))) %>%
  mutate(num_descending_tips = as.numeric(purrr::map2(node, bootstrap, ~ 
                                                        ifelse(!is.na(.y), 
                                                               length(castor::get_subtree_at_node(rooted_serratia_tree,
                                                                                                  as.numeric(.x) - length(rooted_serratia_tree$tip.label))$subtree$tip.label), NA))))
#draw new_tree with bootstraps:
library(ggnewscale)

#first all bootstraps as a color gradient
new_tree + 
  ggnewscale::new_scale_color() +
  ggtree::geom_point2(aes(subset=!isTip & ! is.na(bootstrap), color = bootstrap), na.rm = T, alpha = 0.75) + 
  theme(legend.position = "right")


#then just all bootstraps >=80 with at least 20 descending tips
new_tree <- new_tree +   
  ggnewscale::new_scale_color() +
  ggtree::geom_point2(aes(subset=!isTip & ! is.na(bootstrap) & bootstrap >= 80 & num_descending_tips > 10, color = bootstrap_class), na.rm = T) +
  scale_color_grey(start = 0.7, end = 0.3) +
  theme(legend.position = "right")


#then only the bootstraps
new_tree +
  ggnewscale::new_scale_color() +
  ggtree::geom_point2(aes(subset=!isTip & ! is.na(ancestral_to_lineage), color = bootstrap_class)) +
  scale_color_grey(start = 0.7, end = 0.3) +
  theme(legend.position = "right")


### start adding metadata columns - selecting species, source, dataset.
serratia_metadata <- read.csv("../../figshare_data/metadata/serratia_metadata.csv", header = T,
                              quote = "",
                              stringsAsFactors = F,
                              comment.char = "")


paper_df <- serratia_metadata %>%
  select(File_prefix,Short.species,Paper,Generic_host,Country) %>%
  mutate(Short.species = ifelse(Short.species == "Unknown",NA,Short.species)) %>%
  mutate(Country = ifelse(Country == "None",NA,Country)) %>%
  mutate(Country = ifelse(Country == "Unknown",NA,Country)) %>%
  rename(`Labelled species` = Short.species, Dataset = Paper, Source = Generic_host) %>%
  tibble::column_to_rownames(var = "File_prefix")
  
  

##plot a gheatmap:



new_tree %>% ggtree::gheatmap(paper_df, color = NULL,
                      colnames_offset_y = 5,
                      colnames_angle = 90,
                      hjust = 0,
                      colnames_position = "top",
                      width = 0.5,
                      offset = ( max(new_tree$data$x) * 0.05 )) + 
  scale_fill_manual(values = c(cols)) +
  theme(legend.position = "right") +
  ggplot2::ylim(0, max(new_tree$data$y) + 50)
  

ggsave("tree.svg",units = c("in"), width = 5.88*2, height = 7.00*2)

## do legend/key separately:

# library(cowplot)
# library(glue)
# library(stringr)

legends_list <- list()
i = 0

for (colname in colnames(paper_df)) {
  temp_gheatmap_df <- paper_df %>% select(colname)
  
  ### split the text if value longer than 20 spaces:
  # temp_gheatmap_df <- temp_gheatmap_df %>% 
  #   mutate(test_col = map(.[[1]]), ~ ifelse(length(.x) > 15, gsub(" ","\n",)))
  
  ##get legend
  temp_legend <- cowplot::get_legend(new_tree %>% ggtree::gheatmap(temp_gheatmap_df, color = NULL) + 
                                       scale_fill_manual(values = c(cols),
                                                          labels = function(x) str_wrap(x, width = 25)) +
                                       theme(legend.position = "bottom",
                                             legend.spacing.y = unit(0.1,"cm")) +
                                       ggplot2::guides(color = F,
                                                       fill = guide_legend(title.position = "left",
                                                                           ncol = 5,
                                                                           byrow = TRUE)) +
                                       ggplot2::labs(fill = colname))
  #now plot the legend
  cowplot::ggdraw(temp_legend)

  ggsave(filename = glue::glue("legend_{colname}.svg"), 
         units = c("in"), 
         width = 20, 
         height = 10)
  
  
}


### get maximal root-tip distances:

fastani_95 <- fastani_data %>%
  select(Name,X95) %>%
  rename(strain = Name, phylogroup = X95) %>%
  mutate(strain = gsub("#","_",strain))
#get root nodes for each phylogroup:

ani_clade_labels <- micro.gen.extra::get_clade_nodes(genus_groups_tree, 
                                                 fastani_95 %>%
                                                   rename(cluster=phylogroup))

ani_clade_labels_df <- data.frame(node = unlist(ani_clade_labels), phylogroup_ancestral_node = names(ani_clade_labels))


root_tip_distances <- new_tree$data %>%
  #add in phylogroup_nodes 
  left_join(ani_clade_labels_df, by = "node") %>%
  filter(isTip == TRUE | ! is.na(phylogroup_ancestral_node)) %>%
  mutate(phylogroup = ifelse(is.na(lineage), phylogroup_ancestral_node, phylogroup)) %>%
  #group
  group_by(phylogroup) %>%
  #get max to mix in x:
  mutate(max_dist = max(x) - min(x)) %>%
  select(phylogroup,max_dist) %>%
  unique()
  
#phylogroup names, numbers and colours
names_and_numbers <- data.frame(colour = tree_cols, numbers = names(tree_cols)) %>%
  left_join(data.frame(colour = colors, species = names(colors))) %>%
  left_join(root_tip_distances %>% rename(numbers = phylogroup)) %>%
  #join in phylogenetic diversity
  left_join( data.frame(numbers = names(dplyr_genus_groups)) %>%
               mutate(tips = purrr::map(numbers, ~ dplyr_genus_groups[[.x]])) %>%
               #now get the phylogenetic diversity - using caper
               mutate(pd = unlist(purrr::map(tips, ~ caper::pd.calc(genus_groups_tree, tip.subset = .x)))))


#read in the core and accessory numbers of genes:
core_and_accessory_size <- read.csv("/home/djwilliams/Documents/Serratia_genus_paper/core_and_accessory_gene_count.csv") %>%
  rename(species = Species) %>%
  left_join(names_and_numbers)

ggplot(core_and_accessory_size) + 
  geom_smooth(aes(max_dist, Accessory), method = "lm") +
  geom_point(aes(max_dist,Accessory, size = Num..genomes, color = numbers)) + 
  scale_color_manual(values = c(tree_cols)) + 
  theme_bw()

ggplot(core_and_accessory_size) + 
  geom_smooth(aes(pd, Accessory), method = "lm") +
  geom_point(aes(pd, Accessory, size = Num..genomes, color = numbers)) + 
  scale_color_manual(values = c(tree_cols)) + 
  theme_bw()






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


#ape::write.tree(rooted_serratia_tree, file = "../../figshare_data/tree/rooted_serratia_tree.newick")

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


#### 2 - read in the pangenome ####

csv_pangenome_df <- read.table("../../figshare_data/pangenome/gene_presence_absence_w_header.csv",
                           header = T,
                           sep=",",
                           quote = "",
                           comment.char = "",
                           stringsAsFactors = F,
                           check.names = F)



pangenome_df <- read.table("../../figshare_data/pangenome/gene_presence_absence.Rtab",
                           header = T,
                           sep="\t",
                           quote = "",
                           comment.char = "",
                           stringsAsFactors = F,
                           check.names = F)

pangenome_df_transposed <- pangenome_df %>%
  tidyr::pivot_longer(!Gene,names_to="strain",values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Gene, values_from = presence) %>%
  tibble::column_to_rownames(var = "strain")


## read in eggnog data:
filtered_annotation <- read.csv("../../figshare_data/pangenome/eggnog_annotations.csv",
                                header = T,
                                quote = "",
                                stringsAsFactors = F,
                                comment.char = "",
                                check.names = F) %>%
  select(query_name,description, COG_letter) %>%
  rename(Gene = query_name, desc = description)

#set an extra naming column:


upset_dfs <- create_upset_dfs(gene_presence_absence_file = "../../figshare_data/pangenome/gene_presence_absence.Rtab",
                              groups = fastbaps_l3 %>% rename(group = cluster))


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
### get the species for eachc poppunk cluster:
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



#plot
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



#### 3 - read in twilight data ####

twilight_analysis <- read.delim("../../figshare_data/pangenome/classification.tab",
                                sep = "\t")

#add this to the upsetR data:

core_upset_with_twilight <- upset_dfs$Core %>% left_join(twilight_analysis, by = c("Gene"="gene_name"))


## plot the upsetR plot for the core - but include gals colours for the different classes for the genes:


library(UpSetR)


core_genes_required <- core_upset_with_twilight %>% 
  #tibble::column_to_rownames(var = "Gene") %>%
  #select(Gene,`1`:`23`,specific_class)  %>% 
  rowwise() %>% 
  mutate(total_number = sum(c_across(`1`:`23`))) %>%
  filter(total_number > 0)

svg("twilight_colours_upset_include_singletons.svg", width = 13, height = 8)
upset(core_upset_with_twilight, sets = c(lineage_order), order.by = "freq", keep.order = T,
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
                        S.odorifera = "#D7B5A6",
                        S.marcescens_like = "#B3E2CD",
                        S.rubidaea_like = "#a0a0a0"),
             alpha = 0.65))),
      queries = list(
        list(query = elements,
             params = list("specific_class", c("Absent in large lineages","Core and intermediate","Core and rare","Core, intermediate and rare","Intermediate and rare","Collection rare","Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#d3d3d3", active = T),
        ###information storage and processing
        list(query = elements,
             params = list("specific_class", c("Core and intermediate","Core and rare","Core, intermediate and rare","Intermediate and rare","Collection rare","Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#238b45", active = T),
        ### metabolism
        list(query = elements,
             params = list("specific_class", c("Core and rare","Core, intermediate and rare","Intermediate and rare","Collection rare","Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#74c476", active = T),
        list(query = elements,
             params = list("specific_class", c("Core, intermediate and rare","Intermediate and rare","Collection rare","Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#bae4b3", active = T),
        list(query = elements,
             params = list("specific_class", c("Intermediate and rare","Collection rare","Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#edf8e9", active = T),
        ###cellular processing
        list(query = elements,
             params = list("specific_class", c("Collection rare","Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#b1300b", active = T),
        list(query = elements,
             params = list("specific_class", c("Multi-lineage rare","Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#d95f0e", active = T),
        list(query = elements,
             params = list("specific_class", c("Lineage specific rare","Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#fec44f", active = T),
        list(query = elements,
             params = list("specific_class", c("Collection intermediate","Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#7d1158", active = T),
        list(query = elements,
             params = list("specific_class", c("Multi-lineage intermediate","Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#c51b8a", active = T),
        list(query = elements,
             params = list("specific_class", c("Lineage specific intermediate","Collection core","Multi-lineage core","Lineage specific core")), color = "#fa9fb5", active = T),
        list(query = elements,
             params = list("specific_class", c("Collection core","Multi-lineage core","Lineage specific core")), color = "#08519c", active = T),
        list(query = elements,
             params = list("specific_class", c("Multi-lineage core","Lineage specific core")), color = "#8c96c6", active = T),#,
        list(query = elements,
             params = list("specific_class", c("Lineage specific core")), color = "#542788", active = T)))
dev.off()



##plot a legend:


colours_df = data.frame(
  Class = c( "Lineage specific core","Multi-lineage core", "Collection core","Lineage specific intermediate",
             "Multi-lineage intermediate","Collection intermediate","Lineage specific rare", "Multi-lineage rare" ,
             "Collection rare", "Intermediate and rare","Core, intermediate and rare","Core and rare", "Core and intermediate",
             "Absent in large lineages"),
  Colour = c("#542788","#8c96c6","#08519c","#fa9fb5","#c51b8a","#7d1158",
             "#fec44f","#d95f0e","#b1300b","#edf8e9","#bae4b3","#74c476","#238b45", "#d3d3d3"), stringsAsFactors = F
)

colours_df$Class <- factor(colours_df$Class, level = c(colours_df$Class))

g1 <- ggplot(colours_df, aes(Class, fill = Class)) + geom_bar() +
  scale_fill_manual(values = c(colours_df %>% pull(Colour, name = Class))) +
  guides(fill = guide_legend(ncol = 2))


library(cowplot)
l1 <- get_legend(g1)

ggdraw(l1)
ggsave("legend.svg")

#### 4 - now extract out the intersections of interest #### 

# sup fig 5 - entomophila - L17,18,19
# sup fig 6 - liq. complex - L2, 3, 5, 6, 7, 8
# sup fig 7 - marc L11,10,14
# sup fig 8 - marc L9

clusters_as_characters <- as.character(sort(unique(fastbaps_l3$cluster)))

#entomophila
entomophila_heatmap <- plot_intersection_on_tree(ggtree_object = new_tree,
                          intersection_members = c("17","18","19"),
                          upset_df = upset_dfs$Core,
                          eggnog_annotation = filtered_annotation,
                          pangenome_df_transposed = pangenome_df_transposed,
                          change_names = TRUE)

# 11, 10, 14
lineages_11_10_14 <- plot_intersection_on_tree(ggtree_object = new_tree,
                                                intersection_members = c("11","10","14"),
                                                upset_df = upset_dfs$Core,
                                                eggnog_annotation = filtered_annotation,
                                                pangenome_df_transposed = pangenome_df_transposed,
                                               change_names = TRUE)


# 9 
lineage_9 <- plot_intersection_on_tree(ggtree_object = new_tree,
                                       intersection_members = c("9"),
                                       upset_df = upset_dfs$Core,
                                       eggnog_annotation = filtered_annotation,
                                       pangenome_df_transposed = pangenome_df_transposed,
                                       change_names = TRUE)


# prot + grim + liq
prot_grim_liq <- plot_intersection_on_tree(ggtree_object = new_tree,
                                               intersection_members = c("3","2","5","6","7","8"),
                                               upset_df = upset_dfs$Core,
                                               eggnog_annotation = filtered_annotation,
                                               pangenome_df_transposed = pangenome_df_transposed,
                                           change_names = TRUE)

#selected
selected <- list(entomophila_heatmap,
                      lineages_11_10_14,
                      lineage_9,
                      prot_grim_liq)

names(selected) <- c("entomophila_heatmap",
                          "lineages_11_10_14",
                          "lineage_9",
                          "prot_grim_liq")

#save out colnames and read them back in
for (i in 1:length(selected)) {

  heatmap <- selected[i]
  name <- names(selected)[i]

  temp_df <- data.frame(name = colnames(heatmap[[1]]))

  write.csv(temp_df,
            file = glue::glue("colnames_for_heatmaps/{name}_colnames.csv"),
            quote = F,
            row.names = F)
}

##apply new names after shortening them:
for (i in 1:length(selected)) {
  
  heatmap <- selected[i]
  name <- names(selected)[i]
  
  new_colnames <- read.csv(glue::glue("colnames_for_heatmaps/new_{name}_colnames.csv"))$name
  
  colnames(selected[[i]]) <- new_colnames
  
}

# save original names against the changed names in a dataframe so the reader can check back to find genes of interest

lineage_numbers <- list(c("17","18","19"),
                        c("11","10","14"),
                        c("9"),
                        c("3","2","5","6","7","8"))

names(lineage_numbers) <- c("entomophila_heatmap",
                     "lineages_11_10_14",
                     "lineage_9",
                     "prot_grim_liq")


for (i in 1:length(lineage_numbers)) {
  name <- names(lineage_numbers)[[i]]
  nums <- lineage_numbers[[i]]
  
  temp_heatmap <- plot_intersection_on_tree(ggtree_object = new_tree,
                                            intersection_members = nums,
                                            upset_df = upset_dfs$Core,
                                            eggnog_annotation = filtered_annotation,
                                            pangenome_df_transposed = pangenome_df_transposed)
  
  names_comparison <- data.frame(original_names = colnames(temp_heatmap),
                                 new_names = colnames(selected[[i]]))
  
  #write this out
  
  write.csv(x = names_comparison,
            file = glue::glue("{name}_colnames.csv"),
            row.names = F,
            quote = F)
  

}


### now plot these supplementary plots

counter = 0

for (i in selected) {

  counter <- counter + 1

  new_tree %>% ggtree::gheatmap(i, color = NULL,
                                colnames_position = "top",
                                colnames_angle = 45,
                                colnames_offset_y = 5,
                                hjust = 0,
                                width = 3,
                                offset = 0.1,
                                font.size = 0.8) +
    scale_fill_gradient(low = "white", high = "blue") +
    theme(legend.position = "none",
          text = element_text(size = 0.5),
          rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent")) +
    ylim(0,max(t2$data$y)+300) +
    xlim(0,(max(t2$data$x) + (max(t2$data$x)*4)))

  name = names(selected)[counter]


  ggsave(glue::glue("{name}_plot.png", dpi = 600, height = 15*2, width = 20*2, units = c("cm")))

  ggsave(glue::glue("{name}_plot.svg", height = (15*2)/2.54, width = (20*2)/2.54, units = c("cm")))
  

}



#### 5 - read in gc details ####

whole_genome_gc <- read.csv("../../figshare_data/metadata/serratia_metadata.csv", header = T,
                              quote = "",
                              stringsAsFactors = F,
                              comment.char = "") %>%
  select(File_prefix, GC....) %>%
  rename(strain = File_prefix, wg_GC = GC....) %>%
  left_join(fastani_data %>% rename(strain = Name, phylogroup = X95) %>% select(strain,phylogroup) %>%
              mutate(strain = gsub("#","_",strain)), by = "strain") %>%
  group_by(phylogroup) %>%
  summarise(wg_GC = mean(wg_GC), n = n(), sd = sd(wg_GC)) %>%
  arrange(-as.numeric(wg_GC))
  
jkl <- pangenome_df %>% filter(Gene == "group_16143") %>%
  tidyr::pivot_longer(cols = !Gene,
                      names_to = "strain",
                      values_to = "presence")

csv_jkl <- csv_pangenome_df %>% filter(Gene == "group_16143") %>%
  tidyr::pivot_longer(cols = !Gene,
                      names_to = "strain",
                      values_to = "presence")


#read gc data

#first untar:
untar(tarfile = "../../figshare_data/gc_content/panaroo_gc_details.tar.gz",
      exdir = "../../figshare_data/gc_content/")

gc_details <- read.csv("../../figshare_data/gc_content/panaroo_gc_details.csv", stringsAsFactors = F, quote = "", comment.char = "") %>%
  mutate(strain = gsub("#","_",strain))

#now remove file

file.remove("../../figshare_data/gc_content/panaroo_gc_details.csv")


#replace hashes:
gc_details$strain <- gsub("#","_",gc_details$strain) 

# try to look at the GC3 values for the entomophila genes that are also seen across the gneus:


#re-pull out entomophila genes:
#entomophila
entomophila_heatmap <- plot_intersection_on_tree(ggtree_object = new_tree,
                                                 intersection_members = c("17","18","19"),
                                                 upset_df = upset_dfs$Core,
                                                 eggnog_annotation = filtered_annotation,
                                                 pangenome_df_transposed = pangenome_df_transposed)



names_comparison <- data.frame(original_names = colnames(entomophila_heatmap),
                               new_names = colnames(selected$entomophila_heatmap))

genes_and_strains_in_entomophila_heatmap <- entomophila_heatmap %>%
  tibble::rownames_to_column(var = "strain") %>%
  tidyr::pivot_longer(cols = !strain, names_to = "group", values_to = "presence") %>%
  filter(presence > 0 ) %>%
  #now join in the GC data:
  left_join(gc_details, by = c("strain","group")) %>%
  #add gals classifications:
  left_join(twilight_analysis %>% rename(group = gene_name) %>% select(group, specific_class), by = "group") %>%
  #merge in fastbaps cluster (lineage)
  left_join(fastbaps_l3 %>% rename(lineage = cluster), by = "strain") %>%
  #add phylogroup
  left_join(fastani_data %>% rename(strain = Name, phylogroup = X95) %>% select(strain,phylogroup) %>%
              mutate(strain = gsub("#","_",strain)), by = "strain") %>%
  #collapse the gc values into a single column:
  tidyr::pivot_longer(cols = c(GC1,GC2,GC3), names_to = "measurement", values_to = "gc_perc") 

#genes_and_strains_in_entomophila_heatmap$phylogroup <- factor(genes_and_strains_in_entomophila_heatmap$phylogroup, levels = c(whole_genome_gc$phylogroup))


# calculate average for all entomophila
average_entomophila <- genes_and_strains_in_entomophila_heatmap %>%
  select(strain,group,gc_perc,measurement,lineage) %>%
  mutate(gc_perc = as.numeric(gc_perc)) %>%
  filter(measurement == "GC3") %>%
  filter(!is.na(gc_perc)) %>%
  filter(lineage %in% c(17,18,19)) %>%
  group_by(group) %>%
  summarise(gc_perc = mean(gc_perc))


# n_over_2 <- c("26968_7_107","26968_7_178","26968_7_219","26968_7_249","26968_7_80","26968_7_83",
#               "26968_7_89","GCA_900478125.1_28193_B02_genomic","GCA_900635625.1_31436_B01_genomic")

#now create heatmap of distance away for GC3 from the 
GC3_entomophila_genes_data <- genes_and_strains_in_entomophila_heatmap %>%
  #select(strain,group,gc_perc,measurement,lineage) %>%
  mutate(gc_perc = as.numeric(gc_perc)) %>%
  filter(measurement == "GC3") %>%
  select(!measurement) %>%
  distinct() %>%
  #replace na with a value?
  
  
  filter(!is.na(gc_perc)) %>%
  group_by(strain,group) %>%
  #calculate weight mean:
  summarise(GC3_weighted = weighted.mean(gc_perc,length)) %>%
  #now join in the average:
  left_join(average_entomophila %>% rename(av_gc3 = gc_perc), by = "group") %>%
  mutate(distance_from_gc3 = GC3_weighted - av_gc3) %>%
  # select then pivot wider to make heatmap
  select(strain,group,distance_from_gc3) %>%
  tidyr::pivot_wider(names_from = group, values_from = distance_from_gc3, values_fill = NA) %>%
  tibble::column_to_rownames(var = "strain")

##plot:

#change the colnames?
GC3_entomophila_genes_data <- GC3_entomophila_genes_data[sort(colnames(GC3_entomophila_genes_data))]

#now change the colnames:

colnames(GC3_entomophila_genes_data) <- colnames(selected$entomophila_heatmap)

gc3_1 <- new_tree %>% ggtree::gheatmap(GC3_entomophila_genes_data, color = NULL,
                              colnames_position = "top",
                              colnames_angle = 45,
                              colnames_offset_y = 5,
                              hjust = 0,
                              width = 3,
                              offset = 0.1,
                              font.size = 0.8) +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red",na.value = "white") +
  theme(legend.position = "none",
        text = element_text(size = 0.5),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) +
  ylim(0,max(t2$data$y)+300) +
  xlim(0,(max(t2$data$x) + (max(t2$data$x)*4)))

ggsave(plot = gc3_1,
       filename ="entomophila_GC3.svg", dpi = 600, height = (15*2)/2.54, width = (24*2)/2.54, units = c("cm"))


### get the legend:
library(cowplot)

l1 <- get_legend(gc3_1 + theme(legend.position = "right", text = element_text(size = 20)) +
                   guides(color = "none"))

ggsave(plot = ggdraw(l1),
       file = "entomophila_legend.svg")

## repeat this - but for liquefaciens example (2,3,5,6,7,8) and the proteamaculans example (11,10,14):

#re-pull out liq complex genes
prot_grim_liq <- plot_intersection_on_tree(ggtree_object = new_tree,
                                           intersection_members = c("3","2","5","6","7","8"),
                                           upset_df = upset_dfs$Core,
                                           eggnog_annotation = filtered_annotation,
                                           pangenome_df_transposed = pangenome_df_transposed)


liq_complex_names_comparison <- data.frame(original_names = colnames(prot_grim_liq),
                                           new_names = colnames(selected$prot_grim_liq))



average_liq_complex <- prot_grim_liq %>%
  tibble::rownames_to_column(var = "strain") %>%
  tidyr::pivot_longer(cols = !strain, names_to = "group", values_to = "presence") %>%
  filter(presence > 0 ) %>%
  #now join in the GC data:
  left_join(gc_details, by = c("strain","group")) %>%
  #add gals classifications:
  left_join(twilight_analysis %>% rename(group = gene_name) %>% select(group, specific_class), by = "group") %>%
  #merge in fastbaps cluster (lineage)
  left_join(fastbaps_l3 %>% rename(lineage = cluster), by = "strain") %>%
  #collapse the gc values into a single column:
  tidyr::pivot_longer(cols = c(GC1,GC2,GC3), names_to = "measurement", values_to = "gc_perc") %>%
  select(strain,group,gc_perc,measurement,lineage) %>%
  mutate(gc_perc = as.numeric(gc_perc)) %>%
  filter(measurement == "GC3") %>%
  filter(!is.na(gc_perc)) %>%
  filter(lineage %in% c(2,3,4,6,7,8)) %>%
  group_by(group) %>%
  summarise(gc_perc = mean(gc_perc))





liq_complex_amelioration_heatmap <- prot_grim_liq %>%
  tibble::rownames_to_column(var = "strain") %>%
  tidyr::pivot_longer(cols = !strain, names_to = "group", values_to = "presence") %>%
  filter(presence > 0 ) %>%
  #now join in the GC data:
  left_join(gc_details, by = c("strain","group")) %>%
  #collapse the gc values into a single column:
  tidyr::pivot_longer(cols = c(GC1,GC2,GC3), names_to = "measurement", values_to = "gc_perc") %>%
  #select(strain,group,gc_perc,measurement,lineage) %>%
  mutate(gc_perc = as.numeric(gc_perc)) %>%
  filter(measurement == "GC3") %>%
  distinct() %>%
  filter(!is.na(gc_perc)) %>%
  group_by(strain,group) %>%
  mutate(n = n()) %>%
  ungroup %>% group_by(strain, group) %>%
  #calculate weight mean:
  summarise(GC3_weighted = weighted.mean(gc_perc,length)) %>%
  # now compare to average value
  left_join(average_liq_complex %>% rename(av_gc3 = gc_perc), by = "group") %>%
  mutate(distance_from_gc3 = GC3_weighted - av_gc3) %>%
  # select then pivot wider to make heatmap
  select(strain,group,distance_from_gc3) %>%
  tidyr::pivot_wider(names_from = group, values_from = distance_from_gc3, values_fill = NA) %>%
  tibble::column_to_rownames(var = "strain")


#change the colnames?
liq_complex_amelioration_heatmap <- liq_complex_amelioration_heatmap[sort(colnames(liq_complex_amelioration_heatmap))]

colnames(liq_complex_amelioration_heatmap) <- colnames(selected$prot_grim_liq)

#plot:
gc3_2 <- new_tree %>% ggtree::gheatmap(liq_complex_amelioration_heatmap, color = NULL,
                              colnames_position = "top",
                              colnames_angle = 45,
                              colnames_offset_y = 5,
                              hjust = 0,
                              width = 3,
                              offset = 0.1,
                              font.size = 0.8) +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red",na.value = "white") +
  theme(legend.position = "none",
        text = element_text(size = 0.5),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) +
  ylim(0,max(t2$data$y)+300) +
  xlim(0,(max(t2$data$x) + (max(t2$data$x)*4))) 

ggsave(plot = gc3_2,
       filename = "liq_complex_GC3.svg", 
       dpi = 600, height = (15*2)/2.54, width = (24*2)/2.54, units = c("cm"))


#get the legend?
l2 <- get_legend(gc3_2 + theme(legend.position = "right", text = element_text(size = 20)) +
                   guides(color = "none"))

ggsave(plot = ggdraw(l2),
       file = "liq_complex_legend.svg")





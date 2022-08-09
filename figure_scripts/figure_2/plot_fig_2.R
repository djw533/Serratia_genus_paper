
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


#### 2 - load in gene_presence_absence data and plot pan genome curves ####
  
  
panaroo_data <- read.table("../../figshare_data/pangenome/gene_presence_absence.Rtab",
                   header = T,
                   sep="\t",
                   quote = "",
                   comment.char = "",
                   stringsAsFactors = F,
                   check.names = F)


colnames(panaroo_data) <- gsub("X26968", "26968", colnames(panaroo_data))


##### add in the eggnog data:

annotation <- read.csv("../../figshare_data/pangenome/eggnog_annotations.csv",
                       header = T,
                       quote = "",
                       stringsAsFactors = F,
                       comment.char = "",
                       check.names = F)

#### just take the query name and the cog letter:

letters <- subset.data.frame(annotation, select = c("query_name","COG_letter","description"))###also take the description
letters_table <- data.frame(table(unlist(letters$COG_letter)))
colnames(letters_table) <- c("COG_letter","Number")

### merge panaroo data and the annotation:
panaroo_data <- merge(panaroo_data,letters,by.x ="Gene",by.y = "query_name",all.x = T)



# ### 5 - plot pan-genome accumulation curves ###
panaroo_strains <- as.data.frame(colnames(panaroo_data)[2:(ncol(panaroo_data)-2)])
colnames(panaroo_strains) <- c("strain")


strains_and_phylogroups_copy <- fastani_data %>%
  select(Name,X95) %>%
  rename(strain = Name, phylogroup = X95) %>%
  mutate(strain = gsub("#","_",strain)) %>%
  left_join(panaroo_strains, by = "strain")
  

df_phylogroups = data.frame(cluster = character(0),
                            richness = numeric(0),
                            genomes = numeric(0),
                            sd = numeric(0))

for (group in unique(strains_and_phylogroups_copy$phylogroup)){
  print(group)
  ## get strains for this phylogroup:
  group_strains <- subset.data.frame(strains_and_phylogroups_copy, phylogroup == group)$strain

  #subset out the diff gene_presence_absence_data:
  temp_df <- panaroo_data[group_strains]
  ## set the rownames as genes:
  rownames(temp_df) <- panaroo_data$Gene


  if (ncol(temp_df) > 2 ) { # this will only work if more than 2 strains in the poppunk cluster
    no_values <- data.frame(t(temp_df[apply(temp_df[,-1], 1, function(x) !all(x==0)),]))
  } else {
    next
  }

  ##add to the accum dataset:
  sp <- vegan::specaccum(no_values, "random", permutations=100)

  df_phylogroups = rbind(df_phylogroups, data.frame(cluster = rep(group, length(sp$sites)),
                                                    richness = sp$richness,
                                                    genomes = sp$sites,
                                                    sd = sp$sd))
}
  
# write.table(df_phylogroups, "phylogroup_accumulation_n>2_curves.csv", col.names = T, row.names = F, quote = F, sep = ',')

df_phylogroups = cbind(df_phylogroups, min= df_phylogroups$richness-df_phylogroups$sd, max = df_phylogroups$richness+df_phylogroups$sd )
df_phylogroups = df_phylogroups[which(df_phylogroups$genomes %% 1 == 0),] ## to visualise more clearly
  
# # add species phylogroup:
# 
# phylogroups_and_poppunk_cluster <- unique(subset.data.frame(merge(strains_and_phylogroups, fastani_95, by = "strain"),
#                                                             select = c("phylogroup","Cluster")))
# 
# df_w_phylogroup <- merge(df, phylogroups_and_poppunk_cluster, by.x = "cluster",by.y = "Cluster")

##### plot:

ggplot(df_phylogroups, aes(x = genomes, y = richness, color = as.factor(cluster)))+ geom_line(alpha = 0.6, size = 1) +#+  geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 16) + xlab("Number of genomes") +
  ylab("Unique CDSs") +
  geom_ribbon(aes(ymin = min, ymax = max, fill = as.factor(cluster)), alpha  = 0.1, color = NA) +
  #facet_wrap(vars(phylogroup)) +
  scale_color_manual(name = "Species",values = c(tree_cols)) +
  scale_fill_manual(name = "Species",values = c(tree_cols)) +#+ scale_shape_manual(values = graphics$Shape, guide = F)  +
  theme(legend.position = "right")
ggsave("accum_plots.svg")




#save this dataframe for source data

write.csv(file = "figure_2c.csv",
          x = df_phylogroups %>%
            rename(phylogroup = cluster,
                   num.genomes = genomes,
                   unique.CDSs = richness),
          quote = F,
          row.names = F)



#### 3 - plot data like phandango: ####

#read in the twilight analysis:
classification <- read.table("../../figshare_data/pangenome/classification.tab",
                             sep = "\t",
                             header = T,
                             stringsAsFactors = F) %>%
  select(gene_name,specific_class) %>%
  rename(Gene = gene_name)

##phandango style plot:

#first transpose data 
pangenome_w_classification <- panaroo_data %>%
  select(!c(COG_letter,description)) %>%
  tidyr::pivot_longer(!Gene,names_to="strain",values_to = "presence") %>%
  left_join(classification, by = "Gene") 

### now order according to twilight groups in a specific order
twilight_order <- c("Collection core", "Multi-lineage core","Lineage specific core", 
                    #"Collection intermediate",
                    "Multi-lineage intermediate","Lineage specific intermediate",
                    #"Collection rare",
                    "Multi-lineage rare",
                    "Lineage specific rare","Core and intermediate","Core and rare","Core, intermediate and rare",
                    "Intermediate and rare")#,
                    #"Absent in large lineages")
  

#plot each heatmap: - there will be warnings that geom_text aesthetics are removed, this is due to custom set y-axes
heatmaps <- list()

for (i in 1:length(twilight_order)) {
  
  class <- twilight_order[i]
  print(class)
  
  temp_pangenome <- pangenome_w_classification %>%
    #filter(presence > 0) %>%
    filter(specific_class == class) %>%
    mutate(presence = as.factor(presence)) %>%
    select(!specific_class) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = presence) %>%
    tibble::column_to_rownames(var = "strain")
  
  if (nrow(temp_pangenome) < 664) {
    print(glue::glue("{class} has none"))
  } else if (nrow(temp_pangenome) == 664) {
    heatmaps[[i]] <- temp_pangenome
  }
  
  plot_width <- (ncol(temp_pangenome) / 50000) * 3
  
  #plots[[i]] <- 
  
  # png(glue::glue("plot_{i}.png"), res = 200, height = 7.7, width = 5 + (5*plot_width), units = "cm",
  #     type = "cairo-png")
  new_tree %>% ggtree::gheatmap(temp_pangenome, color = NULL,
                                colnames_offset_y = -100,
                                offset = 0.1,
                                width = plot_width) +
    ylim(0,710) +
    #scale_fill_gradient(low = "white", high = "black") + 
    scale_fill_manual(values = c("0" = "white","1" = "black")) +
    theme(legend.position = "none") +
    theme(plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "#bbbbbb", colour = NA))
  # dev.off()
  
  ggsave(glue::glue("plot_{class}.png"), dpi = 200, height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)
  
}



#### 4 - plot upsetR plot: ####

#create upset dataframes
upset_dfs <- micro.gen.extra::create_upset_dfs(gene_presence_absence_file = "../../figshare_data/pangenome/gene_presence_absence.Rtab",
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


# now add colours for the different twilight classifications

core_upset_with_twilight <- upset_dfs$Core %>% left_join(classification, by = "Gene")
shell_upset_with_twilight <- upset_dfs$Shell %>% left_join(classification, by = "Gene")
cloud_upset_with_twilight <- upset_dfs$Cloud %>% left_join(classification, by = "Gene")

core_genes_required <- core_upset_with_twilight %>% 
  #tibble::column_to_rownames(var = "Gene") %>%
  #select(Gene,`1`:`23`,specific_class)  %>% 
  rowwise() %>% 
  mutate(total_number = sum(c_across(`1`:`23`))) %>%
  filter(total_number > 0)



svg("UpSetR_plot_with_classification_colour.svg", width = 10, height = 8)
UpSetR::upset(core_upset_with_twilight, sets = c(lineage_order), order.by = "freq", keep.order = T,
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


#now draw a legend:

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

l1 <- cowplot::ggdraw(cowplot::get_legend(g1))

ggsave("legend.svg")

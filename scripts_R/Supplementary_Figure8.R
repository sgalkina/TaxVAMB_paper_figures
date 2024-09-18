# install.packages(setdiff(packages, rownames(installed.packages())), repos = "https://cloud.r-project.org/")

# install.packages("BiocManager", repos = "https://cloud.r-project.org/")
# BiocManager::install("ggtree")
# BiocManager::install("ggtreeExtra")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("treeio")

# install.packages("gridExtra", repos = "https://cloud.r-project.org/")
# install.packages("cowplot", repos = "https://cloud.r-project.org/")
# install.packages("patchwork", repos = "https://cloud.r-project.org/")
# install.packages("pheatmap", repos = "https://cloud.r-project.org/")

library(dplyr)
library(readr)
library(tidyr)
library(ggthemes)
library(ggtree)
library(tidytree) 
library(treeio)
library(ape)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(RColorBrewer)
library(grid)
library(ComplexHeatmap)
library(ape)
library(gridExtra)
library(patchwork)
library(pheatmap)
library(circlize)
library(cowplot)

pdf(file = "figures/supfigure8.pdf", height = 15, width = 26)

tree1 <- read.newick("data/MATRIX_MQ.tree")

summary_data <- read.table("data/MATRIX_MQ_summary.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pivot_table_abundance <- read.table("data/MATRIX_combined_table_sorted_MQ.tsv", sep = "\t", stringsAsFactors = FALSE)

samples_ann <- read.table('data/cultivars_ann_MQ.tsv', sep = "\t", header=1)


summary_data <- summary_data %>%
    select(user_genome, classification)
summary_data <- summary_data %>%
    separate(classification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), remove=FALSE, sep = ";")

data <- summary_data %>%
  mutate(
    pattern1 = str_extract(user_genome, "^TaxVAMB_.*"),
  )

data <- summary_data %>% mutate (
    id_source = case_when(
      !is.na(data$pattern1) ~ "TaxVAMB",
      TRUE ~ "Unknown"
    )
  )

df <- subset(data, id_source == 'TaxVAMB')
df_agg <- df %>% group_by(classification, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(TaxVAMB=n())
df2 <- df %>% distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, classification, .keep_all = TRUE)
vec <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "classification")
df2$gs <- apply(df2[, vec], 1, paste, collapse = " ")
df_agg$gs <- apply(df_agg[, vec], 1, paste, collapse = " ")
df_agg <- df_agg[order(match(df_agg$gs, df2$gs)), ]
df2$TaxVAMB <- df_agg$TaxVAMB

df_total = df2

df_total$TaxVAMB <- ifelse(is.na(df_total$TaxVAMB), 0, df_total$TaxVAMB)

df3 <- subset(df_total, select = c(user_genome, classification, TaxVAMB, Phylum))
summary_data_prune <- select(summary_data, c(user_genome, classification, Phylum, Class, Order, Family, Genus, Species))

valid_tips <- df3$user_genome

summary_data_prune <- summary_data_prune[summary_data_prune$user_genome %in% valid_tips, ]

summary_data_prune <- summary_data_prune %>% mutate(across('Species', str_replace, 's__', ''))
summary_data_prune$is_novel <- ifelse(summary_data_prune$Species == "", "*", NA)

summary_data_prune <- summary_data_prune %>% mutate(across('Genus', str_replace, 'g__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Genus, "Genus - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Family', str_replace, 'f__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Family, "Family - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Order', str_replace, 'o__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Order, "Order - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Class', str_replace, 'c__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Class, "Class - All Novel Species"), summary_data_prune$Species)



tree_phylo <- as.phylo(tree1)
tree_tidy_pruned <- drop.tip(tree_phylo, tree_phylo$tip.label[!(tree_phylo$tip.label %in% valid_tips)]) %>% as.treedata
x <- as_tibble(tree_tidy_pruned)

summary_data_prune$Specieslabel = label_pad(summary_data_prune$Species, justify = "right")

tree_tidy <- full_join(x, summary_data_prune, by = c("label" = "user_genome")) %>% as.treedata


p_tree <- ggtree(tree_tidy, layout="rectangular") + 
  xlim(0, 8) + 
  geom_tippoint(aes(shape = is_novel), size=2.5, color="#B30000") + 
  geom_tiplab(aes(label=Specieslabel, color=Family), align=TRUE, family='mono', size=4, linetype='dotted') + theme(legend.position = 'left') + guides(color = guide_legend(ncol = 1))

tip_labels <- get_taxa_name(p_tree)

data_matrix_ordered <- pivot_table_abundance[tip_labels, ]
data_matrix_ordered[is.na(data_matrix_ordered)] <- 0

col_fun = colorRamp2(c(min(data_matrix_ordered), max(data_matrix_ordered)), c("white", "black"))

N_samples = ncol(data_matrix_ordered)

col_date=c("2022-06-07"="#0000FFFF", 
"2022-06-10"= "#3F3FBFFF", 
"2022-06-14"= "#7F7F7FFF", 
"2022-06-17"= "#BFBF3FFF",
"2022-06-21"= "#FFFF00FF",
"2022-06-28"= "#FFBF00FF",
"2022-07-04"= "#FF7F00FF",
"2022-07-07"= "#FF3F00FF",
"2022-07-14"="#FF0000FF")

col_cultivar = c(
  "Heerup"="#a6d854",
  "Kvium"="#ffd92f",
  "Rembrandt"="#8da0cb"
)

col_fungi = c(
  "Sprayed"="#cab2d6",
  "Unsprayed"="#6a3d9a"
)

ht <- Heatmap(data_matrix_ordered,
              cluster_rows = FALSE,  # Rows are not clustered since we already have the tree
              cluster_columns = TRUE,  # Columns are clustered
              show_row_names = FALSE,
              show_column_names = FALSE,
              top_annotation = HeatmapAnnotation(text = anno_text(gsub("^X", "", colnames(data_matrix_ordered)),
                                                                  rot = 90, 
                                                                  just = "right"),
                                                
                                                Cultivar = samples_ann$Cultivar, 
                                                Fungicide = samples_ann$Fungicide,
                                                Unmapped_reads=samples_ann$unmapped_reads,
                                                Sampling_date = samples_ann$Sampling_date,
                                                
                                                col = list(Sampling_date=col_date, Cultivar=col_cultivar, Fungicide=col_fungi)
                                                ),
              col=col_fun,
              rect_gp = gpar(col = "grey", lwd = 1)
             )



heatmap_grob <- grid.grabExpr(draw(ht, heatmap_legend_side = "right"))

space <- ggdraw()
plot_grid(plot_grid(space, p_tree, rel_heights = c(.128, 1), ncol = 1), heatmap_grob, rel_widths = c(0.315, 0.7), align = 'h', axis = 'bt')

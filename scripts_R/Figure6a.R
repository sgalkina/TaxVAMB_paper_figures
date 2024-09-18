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


pdf(file = "figures/figure6a.pdf", height = 8, width = 10)

tree1 <- read.newick("data/MATRIX_MQ.tree")

summary_data <- read.table("data/MATRIX_MQ_summary.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pivot_table_abundance <- read.table("data/MATRIX_pivot_MQ_prevalence.tsv", sep = "\t", stringsAsFactors = FALSE)


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
df_agg <- df_agg[order(match(df_agg$gs, df2$gs)), ] # sort Tax according to the tree
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
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Genus, " Genus - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Family', str_replace, 'f__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Family, " Family - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Order', str_replace, 'o__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Order, " Order - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Class', str_replace, 'c__', ''))
summary_data_prune$Species <- ifelse(summary_data_prune$Species == "", paste(summary_data_prune$Class, " Class - All Novel Species"), summary_data_prune$Species)

summary_data_prune <- summary_data_prune %>% mutate(across('Phylum', str_replace, 'p__', ''))


condition <- pivot_table_abundance > 0.2
satisfying_rows <- rownames(pivot_table_abundance)[condition]


summary_data_prune$show_label <- ifelse(summary_data_prune$user_genome %in% satisfying_rows, summary_data_prune$Species, '')


tree_phylo <- as.phylo(tree1)
tree_tidy_pruned <- drop.tip(tree_phylo, tree_phylo$tip.label[!(tree_phylo$tip.label %in% valid_tips)]) %>% as.treedata
x <- as_tibble(tree_tidy_pruned)
tree_tidy <- full_join(x, summary_data_prune, by = c("label" = "user_genome")) %>% as.treedata



p <- ggtree(tree_tidy, layout="circular") + 
  xlim(0, 4) + 
  geom_tippoint(aes(shape = is_novel), size=2.5, color="#B30000") + 
  geom_tiplab(aes(label=show_label, color=Phylum), align=TRUE, size=2.8, linetype='solid', offset=0.1)

p2 <- gheatmap(
  p, 
  pivot_table_abundance, 
  low="white", 
  high="black", 
  width=0.1, 
  font.size=2, 
  colnames_position = "top", 
  colnames_angle = 90, 
  colnames_offset_y = 1.3, 
  offset=0.01,
  )

plot(p2)


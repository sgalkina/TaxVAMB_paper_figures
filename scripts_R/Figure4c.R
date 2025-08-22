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
# install.packages("optparse", repos = "https://cloud.r-project.org/")
library(optparse)



option_list <- list(
  make_option(c("-t", "--tree"), type="character", default=NULL,
              help="Newick tree file", metavar="character"),
  make_option(c("-s", "--summary"), type="character", default=NULL,
              help="Summary data file (TSV)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output PDF file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tree) | is.null(opt$summary) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("You must provide the tree file, summary data file, and output PDF file.\n", call.=FALSE)
}

tree1 <- read.newick(opt$tree)
tree1 <- drop.tip(tree1, grep("^(VAMB|SemiBin|TaxVAMB)", tree1$tip.label, invert = TRUE, value = TRUE))

summary_data <- read.table(opt$summary, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pdf(file = opt$output, width = 8, height = 8)

summary_data <- summary_data %>%
    select(user_genome, classification)
summary_data <- summary_data %>%
    separate(classification, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

print(summary_data)
data <- summary_data %>%
  mutate(
    pattern1 = str_extract(user_genome, "^(VAMB)"),
    pattern2 = str_extract(user_genome, "^(SemiBin)"),
  )

data <- summary_data %>% mutate (
    id_source = case_when(
      !is.na(data$pattern1) ~ "VAMB",
      !is.na(data$pattern2) ~ "SemiBin",
      TRUE ~ "TaxVAMB"
    )
  )

df <- subset(data, id_source == 'TaxVAMB')
df_agg <- df %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(TaxVAMB=n())
df2 <- df %>% distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all = TRUE)
vec <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
df2$gs <- apply(df2[, vec], 1, paste, collapse = " ")
df_agg$gs <- apply(df_agg[, vec], 1, paste, collapse = " ")
df_agg <- df_agg[order(match(df_agg$gs, df2$gs)), ] # sort Tax according to the tree
df2$TaxVAMB <- df_agg$TaxVAMB

df <- subset(data, id_source == 'VAMB')
df_agg <- df %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(VAMB=n())
df2_vamb <- df %>% distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all = TRUE)
vec <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
df2_vamb$gs <- apply(df2_vamb[, vec], 1, paste, collapse = " ")
df_agg$gs <- apply(df_agg[, vec], 1, paste, collapse = " ")
df_agg <- df_agg[order(match(df_agg$gs, df2_vamb$gs)), ] # sort Tax according to the tree
df2_vamb$VAMB <- df_agg$VAMB

df_total <- merge(df2, df2_vamb, by=vec, all.x=TRUE, all.y=TRUE)

df <- subset(data, id_source == 'SemiBin')
df_agg <- df %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(SemiBin2=n())
df2_vamb <- df %>% distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all = TRUE)
vec <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
df2_vamb$gs <- apply(df2_vamb[, vec], 1, paste, collapse = " ")
df_agg$gs <- apply(df_agg[, vec], 1, paste, collapse = " ")
df_agg <- df_agg[order(match(df_agg$gs, df2_vamb$gs)), ] # sort Tax according to the tree
df2_vamb$SemiBin2 <- df_agg$SemiBin2

df_total <- merge(df_total, df2_vamb, by=vec, all.x=TRUE, all.y=TRUE)

df_total$VAMB <- ifelse(is.na(df_total$VAMB), 0, df_total$VAMB)
df_total$TaxVAMB <- ifelse(is.na(df_total$TaxVAMB), 0, df_total$TaxVAMB)
df_total$SemiBin2 <- ifelse(is.na(df_total$TaxVAMB), 0, df_total$TaxVAMB)

df3 <- subset(df_total, select = c(user_genome, VAMB, TaxVAMB, SemiBin2, Phylum))
summary_data_prune <- select(summary_data, c(user_genome, Species))

valid_tips <- df3$user_genome

summary_data_prune <- summary_data_prune[summary_data_prune$user_genome %in% valid_tips, ]

tree_phylo <- as.phylo(tree1)
tree_tidy_pruned <- drop.tip(tree_phylo, tree_phylo$tip.label[!(tree_phylo$tip.label %in% valid_tips)]) %>% as.treedata
x <- as_tibble(tree_tidy_pruned)
tree_tidy <- full_join(x, summary_data_prune, by = c("label" = "user_genome")) %>% as.treedata

pwidth <- 0.015

ggtree(tree_tidy, layout="circular") + geom_fruit(
         data=df3,
         geom=geom_tile,
         mapping = aes(
                    x=1,
                    y=user_genome, 
                    fill=Phylum,
                    alpha=VAMB,
                    width=0.1,
                   ),
         color = "grey50", 
         pwidth=pwidth,
     ) + geom_fruit(
         data=df3,
         geom=geom_tile,
         mapping = aes(
                    x=1,
                    y=user_genome, 
                    fill=Phylum,
                    alpha=TaxVAMB,
                    width=0.1,
                   ),
         color = "grey50", 
         pwidth=pwidth,
     ) + geom_fruit(
         data=df3,
         geom=geom_tile,
         mapping = aes(
                    x=1,
                    y=user_genome, 
                    fill=Phylum,
                    alpha=SemiBin2,
                    width=0.1,
                   ),
         color = "grey50", 
         pwidth=pwidth,
     ) + geom_tippoint(aes(shape = ifelse(Species == "s__", "*", NA)), size=2, color="#B30000")

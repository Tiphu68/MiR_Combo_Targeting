###This file produces all the figures and tables for the paper.

#######top######
#####Import packages#####
#Data munging
library(reshape2)
library(readr)
library(purrr)
library(gtools)
library(tibble)
library(tidyr)
library(stringr)
library(forcats)

#Plotting
library(ggplot2)
library(ggpubr)
library(xtable)
library(gridExtra)
library(ggrepel)
library(ggdendro)
library(pheatmap)
library(cowplot)
library(lemon)
library(ggExtra)
library(ggvenn)

#custom
library(disruptr)
#load this last because of merge conflicts
library(dplyr)
source('./miR_Analysis_Functions.R')

#load mRNA-protein correlations from the cell paper. 
corr_df <- read_csv("./data_files/Protein_RNA_correlation.csv")

#Load the miRNA seq data
load('./data_files/averaged_miR_values.Rda')
mir_df <- clean_mir_seq(averaged_miR_values)
#Load the cleaned CCLE file
load('./data_files/CCLE_df.Rda')

#Import output from main Network Potential Pipeline
load("./data_files/final_df.Rda")

#Load output from main np pipeline including st. jude samples
load("./data_files/final_df_all.Rda")

#Load processed miR_targeting data file
load("./data_files/targeting_mat.Rda")

#Load evaluated miR cocktails
load("./data_files/MiRCocktail_output.Rda")
cocktail_df <- output[[2]]
cocktail_df$numtargets <- 10
cocktail_df$fraction_target <- 0.2
targets_df <- output[[1]]
rm(output)

load("./data_files/MiRCocktail5_output.Rda")
cocktail5_df <- output[[2]]
cocktail5_df$numtargets <- 5
cocktail5_df$fraction_target <- 0.2
targets5_df <- output[[1]]
rm(output)

load("./data_files/MiRCocktail15_output.Rda")
cocktail15_df <- output[[2]]
cocktail15_df$numtargets <- 15
cocktail15_df$fraction_target <- 0.2
targets15_df <- output[[1]]
rm(output)

# load("./data_files/MiRCocktailtargethalf_output.Rda")
# cocktailhalf_df <- output[[2]]
# cocktailhalf_df$numtargets <- 10
# cocktailhalf_df$fraction_target <- 0.5

load("./data_files/MiRCocktailtargettenth_output.Rda")
cocktailtenth_df <- output[[2]]
cocktailtenth_df$numtargets <- 10
cocktailtenth_df$fraction_target <- 0.1
cocktail_df <- bind_rows(cocktail_df, cocktail5_df,cocktail15_df, 
                         cocktailtenth_df)

#Load miR candidates matrix
load("./data_files/MiRcandidates_matrix.Rda")

#load Null dist info
load("./data_files/EWS_NPdnpNull.Rda")



###################Expression and Network Potential Overview Fig################

#This code creates a dataframe with total np calculated for each cell line -experiment pairing
#Commented out code aggregates even farther to just have average np for each ES cell line and 
#Then even farther to just have the global average np/standard deviation
final_df_fig3 <- final_df_all %>%
  filter(!is.na(np), !is.infinite(np)) %>% 
  group_by(cell_line, experiment_num) %>% 
  summarise(np_sum = sum(np), 
            expression_sum = sum(expression), 
            n_genes = n(), 
            np_normalized = np_sum / n_genes, 
            expression_normalized = expression_sum / n_genes) %>% ungroup() %>%
  mutate(cell_line = fct_reorder(cell_line, np_normalized, .desc = TRUE)) %>% 
  rename(replicate = experiment_num)

# final_df_fig3$experiment_num <- c("one","two","three","one","two","three",
#                                   "one","two","three","one","two","three",
#                                   "one","two","three","one","two","three")

#Also need to prep data to visualize description of G_i
#Goin to try this filtering out the observations with zero expression
#zero network potential - that should be one eto one. 
hist_df_allgenes <- final_df_all %>% 
  filter(expression > 0, np < 0) %>% 
  group_by(gene_name) %>% summarise(np = mean(np), 
                                    expression = mean(expression))

##Make the plots
g1 <- ggplot(data = hist_df_allgenes, aes(expression)) + theme_classic() +
  geom_histogram(bins = 70, color = "red3", fill = "darkred") +
  xlab("mRNA expression (normalized count)") + ylab("") + 
  labs(tag = "A") + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(face = "bold"))

g2 <- ggplot(data = hist_df_allgenes, aes(np)) + theme_classic() +
  geom_histogram(bins = 70, color = "blue", fill = "steelblue") + 
  xlab("Network potential") + ylab("") + labs(tag = "B") + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(face = "bold"))  

#code for boxplot of total expression for each cell line
g3 <- ggplot(data = final_df_fig3, aes(cell_line, expression_normalized))  + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        plot.tag = element_text(face = "bold"), 
        plot.margin = unit(c(0,2.5,0,0), "cm")) +
  ylab("mRNA expression \n(average per gene)") + 
  xlab("") + 
  labs(tag = "C")
# scale_y_continuous(limits = c(106800, 108600), 
#                    breaks = seq(106800, 108600, by = 400)) +

# Code for boxplot of total network potential for each cell line
g4 <- ggplot(data = final_df_fig3, aes(cell_line, np_normalized))  + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        plot.tag = element_text(face = "bold"), 
        plot.margin = unit(c(0,2.5,0,0), "cm")) +
  #scale_y_continuous(limits = c(-343000, -337000), 
  #                  breaks = seq(-343000 -337000, by = 1000)) +
  ylab("Network potential \n(average per gene)") + xlab("") +
  labs(tag = "D")


#Make a histogram of the distribution of G_i and also mRNA expression
png(filename = "./figures/histogram_allgenes_boxplot.png",
    width = 900, height = 700)

(g_final <- grid.arrange(g1, g2, g3, g4, nrow = 2))
dev.off()

ggsave(plot = g_final, filename = "./figures/histogram_allgenes_boxplot.eps", 
       device = "eps", width = 12, height = 7, units = "in")


######################Scatterplot with Marginal histograms###########
g<- ggplot(hist_df_allgenes, aes(x = np, y = expression)) + 
  geom_point(size = 1, alpha = 0.6) + 
  theme_classic() + 
  theme(text = element_text(size = 24)) +
  ylab("mRNA expression (normalized count)") + 
  xlab("network potential \n (average per gene)") 

png(filename = "./figures/jointplot.png",
    width = 900, height = 700)

ggMarginal(p = g, data = hist_df_allgenes, type = "histogram", 
           fill = "slateblue")

dev.off()

######################mRNA overview Figure###########################

#Calculate average dnp and average of several other data objects for each 
#cell_line-Gene_name pair #also can specify what we want the dnp cutoff 
#to be and whether we want to limit to just genes implicated causally in cancer
final_df_mean <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 50, justcancer = FALSE, 
                     keep_cell_line = FALSE) 

final_df_mean_corr <- final_df %>% 
  mean_cellline_calc(topx = FALSE, numgenes = 10, justcancer = FALSE, 
                     keep_cell_line = FALSE) 

#Okay lets prep data to be graphed
corr_df1 <- filter(corr_df, `Gene Symbol` %in% final_df_mean$gene_name)
colnames(corr_df1) <- c("gene_name", "pearson_corr", "spearman_corr")
final_df_mean <- final_df_mean %>% left_join(corr_df1)


#create a vector that tells GGPlot to color code the labels genes associated with cancer
final_df_mean$text_format <- "black"
final_df_mean$text_format[final_df_mean$housekeeping == "yes"] <- "blue"
final_df_mean$text_format[final_df_mean$cancer_associated == "yes"] <- "red"
color_vec_cells <- final_df_mean$text_format
final_df_mean <- final_df_mean %>% 
  filter(!is.na(se_dnp)) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_dnp, .desc = TRUE),
         gene_type = fct_recode(text_format, uncoded = "black", 
                                housekeep = "blue",
                                cancer = "red"))

#Get a corr_df for the plot that is overlaid on the main figure (top 1000 genes by network potential). 
#CCLE_df2 <- filter(CCLE_df, gene_name %in% final_df_mean_corr$gene_name)

corr_df_hist <- filter(corr_df, `Gene Symbol` %in% final_df_mean_corr$gene_name)
colnames(corr_df_hist) <- c("gene_name", "pearson_corr", "spearman_corr")
final_df_mean_hist <- left_join(final_df_mean_corr, corr_df_hist)
#Prep cancer only data to be graphed
final_df_mean_cancer <- final_df %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 50, justcancer = TRUE, 
                     keep_cell_line = FALSE)

corr_df2 <- filter(corr_df, `Gene Symbol` %in% final_df_mean_cancer$gene_name)
colnames(corr_df2) <- c("gene_name", "pearson_corr", "spearman_corr")
final_df_mean_cancer <- final_df_mean_cancer %>% left_join(corr_df2)

final_df_mean_cancer <- final_df_mean_cancer %>% 
  filter(!is.na(se_dnp)) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_dnp, .desc = TRUE))

#create a named vector for the color scale
cols <- as.character(unique(final_df_mean$gene_type))
names(cols) <- c("black", "blue", "red")

#Now we prep the Sj data for their plot
final_df_sj <- final_df_all %>% filter(condition == "PATIENT")%>% 
  group_by(sample_name) %>% 
  group_modify(~data.frame(dnp = .x$dnp, gene_name = .x$gene_name, 
                           dnp_rank = rank(.x$dnp)), .keep = TRUE) %>% 
  filter(dnp_rank <= 50) %>% 
  group_by(gene_name) %>% 
  summarise(num_top_50 = n(), 
            mean_dnp = mean(dnp)) %>% 
  arrange(desc(num_top_50)) %>% 
  slice_max(n = 30, order_by = num_top_50) %>% 
  mutate(gene_name = fct_reorder(gene_name, num_top_50, .desc = TRUE))

#add back the cosmic vars
final_df_sj <- add_cosmic_vars(final_df_sj)
final_df_sj$text_format <- "black"
final_df_sj$text_format[final_df_sj$housekeeping == "yes"] <- "blue"
final_df_sj$text_format[final_df_sj$cancer_associated == "yes"] <- "red"
color_vec_sj <- final_df_sj$text_format


#re-name null columns and join with final df
colnames(null_df) <- c("gene_name", "sample_name", "mean_dnp_null", "sd_dnp_null", "n_null")
null_df <- null_df %>% 
  group_by(gene_name) %>% 
  summarise(mean_dnp_null = mean(mean_dnp_null),
            sd_dnp_null = mean(sd_dnp_null), 
            n_null= sum(n_null))

final_df_mean <- final_df_mean %>% left_join(null_df) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_dnp, .desc = TRUE))
#Make the plots - this is a freaking nightmare because of the competing legends so bear with me 

#3.891 is 99.99% confidence interval - accounts for multiple hypothesis testing of 50 hypothesis
#first we make the entire first part of the plot with one legend
g1a <- ggplot(data = final_df_mean, 
              aes(gene_name, mean_dnp)) + 
  theme_classic() + 
  geom_point(size = 0.9, color = color_vec_cells) + 
  geom_errorbar(aes(x = gene_name, 
                    ymin = mean_dnp - 2*se_dnp,
                    ymax = mean_dnp + 2*se_dnp,
                    color = gene_type), show.legend = TRUE) + 
  geom_linerange(aes(ymin = mean_dnp_null - 3.891*sd_dnp_null,
                  ymax = mean_dnp_null + 3.891*sd_dnp_null)) + 
  scale_color_manual("Gene type", values = c("black", "blue", "red"), 
                     labels = c("", "housekeeping", "cancer-associated"),
                     aesthetics = c("color", "fill")) + 
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0, size = 12),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  ylab(expression(paste(Delta, "Network Potential"))) + 
  theme(axis.text.x = element_text(color = color_vec_cells)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(tag = "A")

leg1 <- g_legend(g1a)

#Next we make a completely different plot for the second part to get the second legend
g1b <- ggplot(data = final_df_mean, 
              aes(gene_name, mean_dnp))+ 
  geom_point() +
  geom_rect(aes(xmin = as.numeric(gene_name) - 0.5, 
                xmax = as.numeric(gene_name) + 0.5,
                ymin = -100, ymax = 0, fill = pearson_corr),
            show.legend = TRUE) +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_fill_viridis_c() + 
  labs(fill = "Pearson Coefficient")

leg2 <- g_legend(g1b)

#now we make the complete plot, with no legend
g1 <- ggplot(data = final_df_mean, 
             aes(gene_name, mean_dnp)) + 
  theme_classic() + 
  geom_point(size = 0.9, color = color_vec_cells) + 
  geom_errorbar(aes(x = gene_name, 
                    ymin = mean_dnp - 2*se_dnp,
                    ymax = mean_dnp + 2*se_dnp,
                    color = gene_type), show.legend = TRUE) +
  geom_linerange(aes(ymin = mean_dnp_null - 3.891*sd_dnp_null,
                     ymax = mean_dnp_null + 3.891*sd_dnp_null)) + 
  scale_color_manual("Gene type", values = c("black", "blue", "red"), 
                     labels = c("", "housekeeping", "cancer-associated"),
                     aesthetics = c("color")) + 
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0, size = 12),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none",
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.2,8,0.2,0.2), "cm")) +
  ylab(expression(paste(Delta, "Network Potential"))) + 
  theme(axis.text.x = element_text(color = color_vec_cells)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(tag = "A") + 
  geom_rect(aes(xmin = as.numeric(gene_name) - 0.5, 
                xmax = as.numeric(gene_name) + 0.5,
                ymin = -100, ymax = 0, fill = pearson_corr),
            show.legend = TRUE) +
  scale_fill_viridis_c() + 
  labs(fill = "Pearson Coefficient")


g2 <- ggplot(final_df_mean_hist, aes(x = pearson_corr)) +
  geom_histogram(fill = 'steelblue2', color = 'navy',bins = 35) +
  labs(x = "Pearson Correlation Coefficient", tag = "B") +
  theme_classic() +
  geom_label_repel(aes(
    label = ifelse(mean_dnp > 500, gene_name, ""), y = 200
  ), max.overlaps = 100000) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(text = element_text(size = 18),
        plot.tag = element_text(size = 20, face = "bold"))

g3 <- ggplot(final_df_sj, aes(x = gene_name, y = num_top_50)) + 
  geom_col(color = "lightgrey", aes(fill = text_format)) + 
  scale_fill_manual("Gene type", values = c("black", "blue", "red"), 
                    labels = c("", "housekeeping", "cancer-associated")) +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0, size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        plot.tag = element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) + 
  geom_hline(yintercept = 5) +
  geom_hline(yintercept = 10) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(y = "Number of times a gene \n was among the top \n50 projected targets", 
       tag = "C")  

#convert legends into ggplot2 objects and clear out some of the extra whitespace
leg1 <- plot_grid(leg1) + 
  theme(plot.margin = unit(c(-3,-3,-3,-3), "cm"))
leg2 <- plot_grid(leg2) + 
  theme(plot.margin = unit(c(-3,-3,-3,-3), "cm"))

gA <- g1 + annotation_custom(
  ggplotGrob(g2), xmin = 20, xmax = 50, ymin = 450, ymax =1850) + 
  annotation_custom(ggplotGrob(leg1), xmin = 54, xmax = 59, ymin = 600) + 
  annotation_custom(ggplotGrob(leg2), xmin = 54, xmax = 59, ymax = 950)

###Need to do the same for the SJ samples
#Draw the plots.
png(filename = "./figures/mRNA_Targetoverviewplot.png",
    width = 900, height = 700)

## I think a panel of this figure should be some kind of number of times 
## in the top 10 etc for x protein among the patient samples
#Here we'll draw the legends as annotations on the main plot
g_final <- plot_grid(gA, g3, nrow = 2, rel_heights = c(1,0.4), align = "v")
g_final

dev.off()
ggsave(plot = g_final, filename = "./figures/mRNA_Targtoverviewplot.eps", 
       device = "eps", width = 12, height = 9, units = "in")

##Make an mRNA target overview limited to genes associated in cancer
png(filename = "./figures/mRNA_Targetoverviewplot_cancer.png",
    width = 900, height = 600)
g <- ggplot(data = final_df_mean_cancer, 
            aes(gene_name, mean_dnp)) 
g <- g + geom_point(size = 0.9) + 
  geom_errorbar(aes(x = gene_name, 
                    ymin = mean_dnp - 2*se_dnp,
                    ymax = mean_dnp + 2*se_dnp)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  geom_rect(aes(xmin = as.numeric(gene_name) - 0.5, 
                xmax = as.numeric(gene_name) + 0.5,
                ymin = -50, ymax = 0, fill = pearson_corr),
            show.legend = TRUE) + 
  ylab(expression(paste(Delta, "Network Potential"))) + 
  xlab("Gene Name")  + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_viridis_c()
g
dev.off()
ggsave(plot = g, filename = "./figures/mRNA_Targtoverviewplot_cancer.eps", 
       device = "eps", width = 10.5, height = 6, units = "in")

########stringdb results##############
#Import output from main Network Potential Pipeline
load("./data_files/EWS_NP_string2dnp.Rda")

#load null dist info
load("./data_files/EWS_NP_string2dnpNull.Rda")
df_dnp <- add_cosmic_vars(df = df_dnp)
df_dnp <- add_entrez_id(df_dnp) %>% 
  mutate(entrez_id = as.character(entrez_id))
###################Expression and Network Potential Overview Fig - string################

#This code creates a dataframe with total np calculated for each cell line -experiment pairing
#Commented out code aggregates even farther to just have average np for each ES cell line and 
#Then even farther to just have the global average np/standard deviation
final_df_fig3 <- df_dnp %>%
  filter(!is.na(np), !is.infinite(np)) %>% 
  group_by(cell_line, experiment_num) %>% 
  summarise(np_sum = sum(np), 
            expression_sum = sum(expression), 
            n_genes = n(), 
            np_normalized = np_sum / n_genes, 
            expression_normalized = expression_sum / n_genes) %>% ungroup() %>%
  mutate(cell_line = fct_reorder(cell_line, np_normalized, .desc = TRUE)) %>% 
  dplyr::rename(replicate = experiment_num)

# final_df_fig3$experiment_num <- c("one","two","three","one","two","three",
#                                   "one","two","three","one","two","three",
#                                   "one","two","three","one","two","three")

#Also need to prep data to visualize description of G_i
#Goin to try this filtering out the observations with zero expression
#zero network potential - that should be one eto one. 
hist_df_allgenes <- df_dnp %>% 
  filter(expression > 0, np < 0) %>% 
  group_by(gene_name) %>% summarise(np = mean(np), 
                                    expression = mean(expression))

##Make the plots
g1 <- ggplot(data = hist_df_allgenes, aes(expression)) + theme_classic() +
  geom_histogram(bins = 70, color = "red3", fill = "darkred") +
  xlab("mRNA expression (normalized count)") + ylab("") + 
  labs(tag = "A") + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(face = "bold"))

g2 <- ggplot(data = hist_df_allgenes, aes(np)) + theme_classic() +
  geom_histogram(bins = 70, color = "blue", fill = "steelblue") + 
  xlab("Network potential") + ylab("") + labs(tag = "B") + 
  theme(text = element_text(size = 18),
        plot.tag = element_text(face = "bold"))  

#code for boxplot of total expression for each cell line
g3 <- ggplot(data = final_df_fig3, aes(cell_line, expression_normalized))  + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        plot.tag = element_text(face = "bold"), 
        plot.margin = unit(c(0,2.5,0,0), "cm")) +
  ylab("mRNA expression \n(average per gene)") + 
  xlab("") + 
  labs(tag = "C")
# scale_y_continuous(limits = c(106800, 108600), 
#                    breaks = seq(106800, 108600, by = 400)) +

# Code for boxplot of total network potential for each cell line
g4 <- ggplot(data = final_df_fig3, aes(cell_line, np_normalized))  + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        plot.tag = element_text(face = "bold"), 
        plot.margin = unit(c(0,2.5,0,0), "cm")) +
  #scale_y_continuous(limits = c(-343000, -337000), 
  #                  breaks = seq(-343000 -337000, by = 1000)) +
  ylab("Network potential \n(average per gene)") + xlab("") +
  labs(tag = "D")


#Make a histogram of the distribution of G_i and also mRNA expression
png(filename = "./figures/histogram_allgenes_boxplot_string2.png",
    width = 900, height = 700)

(g_final <- grid.arrange(g1, g2, g3, g4, nrow = 2))
dev.off()

ggsave(plot = g_final, filename = "./figures/histogram_allgenes_boxplot_string2.eps", 
       device = "eps", width = 12, height = 7, units = "in")


######################Scatterplot with Marginal histograms -string###########
g<- ggplot(hist_df_allgenes, aes(x = np, y = expression)) + 
  geom_point(size = 1, alpha = 0.6) + 
  theme_classic() + 
  theme(text = element_text(size = 24)) +
  ylab("mRNA expression \n (normalized count)") + 
  xlab("network potential \n (average per gene)") 

g_marg <- ggMarginal(p = g, data = hist_df_allgenes, type = "histogram", 
           fill = "slateblue")


######################mRNA overview Figure- string###########################

#Calculate average dnp and average of several other data objects for each 
#cell_line-Gene_name pair #also can specify what we want the dnp cutoff 
#to be and whether we want to limit to just genes implicated causally in cancer
final_df_mean_biogrid <- final_df_mean #save the previous final_df_mean
final_df_mean <- df_dnp %>% 
  mean_cellline_calc(topx = TRUE, numgenes = 50, justcancer = FALSE, 
                     keep_cell_line = FALSE) 



#create a vector that tells GGPlot to color code the labels genes associated with cancer
final_df_mean$text_format <- "black"
final_df_mean$text_format[final_df_mean$housekeeping == "yes"] <- "blue"
final_df_mean$text_format[final_df_mean$cancer_associated == "yes"] <- "red"
color_vec_cells <- final_df_mean$text_format
final_df_mean <- final_df_mean %>% 
  filter(!is.na(se_dnp)) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_dnp, .desc = TRUE),
         gene_type = fct_recode(text_format, uncoded = "black", 
                                housekeep = "blue",
                                cancer = "red"))

#create a named vector for the color scale
cols <- as.character(unique(final_df_mean$gene_type))
names(cols) <- c("black", "blue", "red")

#re-name null columns and join with final df
colnames(null_df) <- c("gene_name", "sample_name", "mean_dnp_null", "sd_dnp_null", "n_null")
null_df <- null_df %>% 
  group_by(gene_name) %>% 
  summarise(mean_dnp_null = mean(mean_dnp_null),
            sd_dnp_null = mean(sd_dnp_null), 
            n_null= sum(n_null))

final_df_mean <- final_df_mean %>% left_join(null_df) %>% 
  mutate(gene_name = fct_reorder(gene_name, mean_dnp, .desc = TRUE))
#Make the plots - this is a freaking nightmare because of the competing legends so bear with me 

#3.891 is 99.99% confidence interval - accounts for multiple hypothesis testing of 50 hypothesis
#first we make the entire first part of the plot with one legend
g1a <- ggplot(data = final_df_mean, 
              aes(gene_name, mean_dnp)) + 
  theme_classic() + 
  geom_point(size = 0.9, color = color_vec_cells) + 
  geom_errorbar(aes(x = gene_name, 
                    ymin = mean_dnp - 2*se_dnp,
                    ymax = mean_dnp + 2*se_dnp,
                    color = gene_type), show.legend = TRUE) + 
  geom_linerange(aes(ymin = mean_dnp_null - 3.891*sd_dnp_null,
                     ymax = mean_dnp_null + 3.891*sd_dnp_null)) + 
  scale_color_manual("Gene type", values = c("black", "blue", "red"), 
                     labels = c("", "housekeeping", "cancer-associated"),
                     aesthetics = c("color", "fill")) + 
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0, size = 11),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm")) +
  ylab(expression(paste(Delta, "Network Potential"))) + 
  theme(axis.text.x = element_text(color = color_vec_cells)) + 
  scale_y_continuous(expand = c(0,0))

#venn_diagram
x <- list(
  biogrid = as.character(final_df_mean_biogrid$gene_name),
  stringdb = as.character(final_df_mean$gene_name)
)
g_venn <- ggvenn(x, text_size = 7, set_name_size = 7, show_percentage = FALSE)

gA <- g1a + annotation_custom(
  ggplotGrob(g_venn), xmin = 12, xmax = 50, ymin = 350, ymax = 700) 

###Need to do the same bullshit for the SJ samples
#Draw the plots.
png(filename = "./figures/string_overview_plot.png",
    width = 900, height = 900)

## I think a panel of this figure should be some kind of number of times 
## in the top 10 etc for x protein among the patient samples
#Here we'll draw the legends as annotations on the main plot
(g_final <- plot_grid(g_marg, gA, labels = "AUTO", nrow =2, label_size = 20))

dev.off()
ggsave(plot = g_final, filename = "./figures/string_overview_plot.eps", 
       device = "eps", width = 12, height = 9, units = "in")


##########################miR Cocktail Fig######################################

#Just going to present the example top cocktail from one cell line + the worst cocktail from the same cell line
#Prepping the dataframe for plotting will use some of the code from the make cocktail function
cell_line_plot <- "A673"
main_targetnum <- 10
main_fractiontarget <- 0.2

targets <- unique(unlist(targets15_df))
#cell_line_targets <- unlist(targets_df[,cell_line_plot])
#pull out the non-target vector
non_targets_df <- filter(final_df, housekeeping == "yes")

#For plotting purposes need to pair down to 30 housekeeping genes
#Bit of a challenge because there are actually too many 
non_targets <- unique(non_targets_df$gene_name)
targeting <- filter(overall_targeting_mat.df, gene_name %in% c(non_targets, targets)) %>% 
  dplyr::rename(mean_dnp = mean_deltaGibbs, 
         weighted_dnp = weighted_deltaGibbs)
targeting$targets_value <- 0

#Set up the targeting matrix for the loss function
#This section allows the user to fiddle with both the priority of hitting targets vs. non-targets and the assumed 
#repression achieved from a given miR. Defaults are in the function call.
targeting$targets_value[
  targeting$housekeeping == "yes" & targeting$targets_logical ==1] <- 1
targeting$targets_value[
  targeting$housekeeping == "no" & targeting$targets_logical ==1] <- -1
#Create a new value for weighted delta np
targeting <- targeting %>% 
  mutate(weighted_target = targets_value*mean_dnp)

#fix the housekeeping variable here..
targeting$housekeeping[targeting$housekeeping == "yes"] <- 1
targeting$housekeeping <- as.numeric(targeting$housekeeping)
targeting$housekeeping[is.na(targeting$housekeeping)] <- 0

#Actually I just want to rank every cocktail by cell line/condition
#split-apply-combine
cocktail_df2 <- cocktail_df %>% 
  group_by(cell_line, numtargets, fraction_target) %>% 
  group_modify(~tibble(rank = rank(.x$loss, ties.method = "first"), 
                       miR1 = .x$V1, 
                       miR2 = .x$V2, 
                       miR3 = .x$V3, 
                       loss = .x$loss), 
               .keep = TRUE)


#Isolate top cocktail for base case condition
top_cocktail <- cocktail_df2 %>% 
  filter(rank == 1, cell_line == cell_line_plot,
         numtargets == main_targetnum,
         fraction_target == main_fractiontarget) %>% 
  ungroup() %>%
  dplyr::select(contains("miR")) %>% unlist(.)

#need to convert housekeeping to logical to rank gene name by that condition
top_cocktail_df <- targeting %>% filter(miR %in% top_cocktail, 
                                        cell_line == cell_line_plot) %>%
  mutate(gene_name = fct_reorder(gene_name, housekeeping, .desc = FALSE))

#removing extra white space in the fig by getting rid of all non-targeted genes.
nottargeted <- top_cocktail_df %>% group_by(gene_name) %>%
  summarise(sum_target = sum(targets_logical)) %>% 
  filter(sum_target == 0)
top_cocktail_df <- filter(top_cocktail_df, !(gene_name %in% nottargeted$gene_name))

top_cocktail_dfhist <- top_cocktail_df %>% group_by(gene_name) %>% 
  filter(housekeeping == TRUE) %>% 
  summarise(num_hits = sum(targets_logical))
top_cocktail_dfhist$status <- "best cocktail"

#Prep the bottom cocktail dfs for plotting
#Need max_rank to be in a different statement because different numbers of cocktails were 
#evaluated for the different conditions. 
bottom_cocktail <- cocktail_df2 %>% 
  filter(cell_line == cell_line_plot,
         numtargets == main_targetnum, 
         fraction_target == main_fractiontarget) %>% 
  filter(rank == max(rank)) %>%
  ungroup %>%
  dplyr::select(contains("miR")) %>% unlist(.)

bottom_cocktail_df <- targeting %>% 
  filter(miR %in% bottom_cocktail, cell_line == cell_line_plot) %>% 
  mutate(gene_name = fct_reorder(gene_name, housekeeping, .desc = FALSE))

#same story as above. 
nottargeted <- bottom_cocktail_df %>% group_by(gene_name) %>%
  summarise(sum_target = sum(targets_logical)) %>% 
  filter(sum_target == 0)
bottom_cocktail_df <- filter(bottom_cocktail_df, !(gene_name %in% nottargeted$gene_name))

bottom_cocktail_dfhist <- bottom_cocktail_df %>% group_by(gene_name) %>% 
  filter(housekeeping == 1) %>% 
  summarise(num_hits = sum(targets_logical))
bottom_cocktail_dfhist$status <- "worst cocktail"
cocktail_dfhist <- bind_rows(top_cocktail_dfhist, bottom_cocktail_dfhist)

#Lets plot the number of times a given miR shows up in the bottom 10 or top 10 cocktails 
#(averaged across cell lines and conditions tested)
mir_rank_df <- cocktail_df2 %>% 
  filter(numtargets == 10) %>% 
  group_by(miR1,miR2,miR3) %>% 
  dplyr::summarise(mean_rank = mean(rank)) %>% 
  tidyr::pivot_longer(cols = c(miR1,miR2,miR3), values_to = "miR", 
               names_to = "cocktail_position") %>% 
  dplyr::arrange(-dplyr::desc(mean_rank))

mir_rank_dftop <- mir_rank_df %>% slice_head(n = 30)
mir_rank_dftop$cocktail_cat <- "top 10"

mir_rank_dfbottom <- mir_rank_df %>% slice_tail(n = 30)
mir_rank_dfbottom$cocktail_cat <- "bottom 10"

mir_rank_df <- bind_rows(mir_rank_dftop, mir_rank_dfbottom)

mir_rank_df <- mir_rank_df %>% 
  group_by(miR, cocktail_cat) %>% 
  summarise(n = n()) %>% ungroup() %>%
  mutate(miR = fct_reorder2(.f = miR, .x = cocktail_cat, .y = n, .desc = TRUE))

#Start of plotting code for main figure
g1 <- ggplot(data = top_cocktail_df, 
             aes(x = gene_name, y = miR, fill = targets_value)) + 
  geom_tile() + 
  scale_fill_gradient2(name = "", 
                       low = "red", mid = "white", high = "blue", 
                       guide = "legend", 
                       breaks = c(-1,0,1),
                       labels = c("intended target", "not targeted", 
                                  "housekeeping gene")) + 
  theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 10), 
        axis.title.x = element_blank(),
        text = element_text(size = 16),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm"),
        legend.position = "bottom",
        legend.background = element_rect(fill = "lightgrey", colour = "black")) + 
  labs(y = "miRNA", 
       title = "Best cocktail") 
g2 <- ggplot(data = top_cocktail_dfhist, aes(num_hits)) + 
  geom_histogram(bins =4, fill = "steelblue") + 
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(1,3)) + 
  theme_classic() + 
  theme(text = element_text(size = 16)) + 
  xlab("Number of hits")
g3 <- ggplot(data = bottom_cocktail_df, 
             aes(x = gene_name, y = miR, fill = targets_value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") + 
  theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 7.5),
        axis.title.x = element_blank(),
        legend.position = "none", 
        text = element_text(size = 16),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm")) + 
  labs(y = "miRNA", 
       title = "Worst cocktail") 
g4 <- ggplot(data = bottom_cocktail_dfhist, aes(num_hits)) + 
  geom_histogram(bins =4, fill = "steelblue") + 
  scale_x_continuous(breaks = c(0,1,2,3)) + 
  theme_classic() + 
  theme(text = element_text(size = 16)) + 
  xlab("Number of hits") 
g5 <- ggplot(data = mir_rank_df, aes(miR, n, fill =cocktail_cat)) + 
  geom_col() + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle = -60, vjust = 0.5, hjust = 0, size = 10),
        text = element_text(size = 16),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm"),
        legend.title = element_blank()) + 
  labs(y = "Number of appearances", x = "miRNA")
##Start plot - main text
png(filename = "./figures/miRCocktailPlot.png", height = 800, width = 1000) 
#g3 is s histogram showing the distribution of hits.
g_final <- plot_grid(g1,g2,g3, g4, g5, nrow = 3, ncol = 2, rel_widths = c(1, 0.5),
          align = "hv", axis = "btl", labels = "AUTO", label_size = 20, 
          label_x = c(0,0.1,0,0.1,0))
g_final
dev.off()

ggsave(plot = g_final, filename = "./figures/miR_CocktailPlot.eps", 
       device = "eps", width = 12, height = 9, units = "in")

#######################miR Overview Figure#########################################

#Set up miR_table to be graphed
#need to aggregate down to gene-miR level (average away cell lines.)
#The data is only processed to include cell lines because it was easier to do the miR
#cocktail dplyr::selection with that data structure. 
#This takes some time.

summary_tbl <- overall_targeting_mat.df %>% 
  group_by(miR, gene_name, cell_line) %>% 
  summarise(dnp_mean = mean(mean_deltaGibbs), 
            targets_logical = unique(targets_logical), 
            weighted_dnp = mean(weighted_deltaGibbs), 
            sd_dnp = sd(mean_deltaGibbs),
            housekeeping = unique(housekeeping)) %>% 
  mutate(miR = paste0("hsa-", miR))  

#We are going to replace NAs in standard deviation with 0. The NAs arise because our
#process of dplyr::selecting the top 5% of genes for each cell lines results in a slightly 
#different set of genes for each cell line. In some cases, a given gene will not be in the top 
#5 % for each cell line. Our treatment of SD implicitly imputes the gibb_diff to be the mean
#of the dnp for the cell_lines that were represented for a given gene. For genes that 
#only showed up in the top 5% for one cell line, the SD is 0 because we are implicitly 
#imputing the dnp for the other 5 cell lines to be the value for the one cell line that we have.
summary_tbl$sd_dnp[is.na(summary_tbl$sd_dnp)] <- 0
summary_tbl <- mutate(summary_tbl, 
                      weighted_sd_dnp = sd_dnp * targets_logical)
#Now we can calculate our top miRs
summary_tbl2 <- summary_tbl %>% group_by(miR, cell_line) %>%
  summarise(num_target = sum(targets_logical), 
            np_target = sum(weighted_dnp),
            sd_target = sum(weighted_sd_dnp),
            num_housekeeping = sum(housekeeping == "yes" & targets_logical ==TRUE)) %>% 
  filter(num_target > 0) %>% arrange(desc(np_target))

miR_table2 <- summary_tbl2 %>% 
  filter(np_target > 10000) %>% ungroup() %>%
  mutate(miR = fct_reorder(miR, np_target, .desc = TRUE))

summary_tbl2 <- summary_tbl2 %>%
  group_by(miR) %>% 
  summarise(num_target = unique(num_target),
            np_target = mean(np_target)) %>% 
  distinct(miR, .keep_all = TRUE) 



#Isolate vector of top 20% of miRs
top_miR_vec <- summary_tbl2 %>% ungroup() %>%
  filter(miR %in% rownames(averaged_miR_values)) %>%
  slice_max(order_by = np_target, prop = 0.2) %>% 
  dplyr::select(miR) %>% unlist()

toptop_miR_vec <- summary_tbl2 %>% ungroup() %>%
  slice_max(order_by = np_target, n = 0) %>% 
  dplyr::select(miR) %>% unlist()

bottom_miR_vec <- summary_tbl2 %>% ungroup() %>%
  filter(miR %in% rownames(averaged_miR_values)) %>%
  slice_min(order_by = np_target, prop = 0.2) %>% 
  dplyr::select(miR) %>% unlist()

mir_rank_df <- filter(mir_rank_df, n>1) %>% 
  mutate(miR = paste0("hsa-", miR))  

#put together annotation data
averaged_miR_values <- scale(averaged_miR_values)

#filter so we only have miRNA for which we have projected disruption in network potential (and lose that random extra cell line)
averaged_miR_values <- averaged_miR_values[rownames(averaged_miR_values) %in% summary_tbl2$miR,]
averaged_miR_values <- averaged_miR_values[,colnames(averaged_miR_values) %in% unique(final_df$cell_line)]

# 
# annotation_df <- data.frame(
#   dnp_rank = case_when(
#     rownames(averaged_miR_values) %in% top_miR_vec ~ "top 20",
#     rownames(averaged_miR_values) %in% bottom_miR_vec ~ "bottom 20", 
#     TRUE ~ "middle 60"
#     ), 
#   row.names = rownames(averaged_miR_values)) %>% 
#   arrange(desc(dnp_rank)) 

###try annotation again w/ continuous
annotation_df <- data.frame(dnp = log(summary_tbl2$np_target), 
                            row.names = summary_tbl2$miR) %>% 
  arrange(desc(dnp))

#arrange the expression matrix so our discrete groups "cluster" together
new_index = match(rownames(averaged_miR_values), rownames(annotation_df))
averaged_miR_values <- averaged_miR_values[order(new_index), ]



g1 <- ggplot(data = miR_table2, aes(miR, (np_target))) +
  geom_boxplot(color = "steelblue", fill = "steelblue") + 
  geom_rect(aes(xmin = as.numeric(miR) - 0.5, 
                xmax = as.numeric(miR) + 0.5,
                ymin = 9000, ymax = 9750, fill = num_housekeeping),
            show.legend = TRUE) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        plot.margin = unit(c(0.2,1.5,0.2,0.2), "cm"),
        plot.tag = element_text(face = "bold")) +
  scale_y_continuous(limits = c(9000, 19000),
                     expand = c(0,0), 
                     breaks = seq(10000,19000, by = 1000)) +
  scale_fill_viridis_c() +
  ylab("Projected disruption in network potential") + 
  xlab("miRNA") + 
  labs(tag = "A",
       fill = "Number of\nhousekeeping\ngenes\ntargeted")

g2 <-ggscatter(summary_tbl2, x = 'num_target', y = 'np_target', add = "reg.line",
               size = 1.75, color = "steelblue") +
  stat_cor(label.y = 17000, 
           aes(label = ..rr.label..), size = 5.5) +
  geom_label_repel(aes(
    label = ifelse(miR %in% toptop_miR_vec, miR,"")
  ), 
  box.padding = 1, size = 4.5, 
  nudge_x = 15, 
  max.overlaps = 10000) + 
  scale_x_continuous(breaks = seq(0,150,by = 25)) +
  xlab("Number of putative \n mRNA targets") + 
  ylab("") + 
  theme(text = element_text(size = 18),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        plot.tag = element_text(face = "bold"),
        legend.position = "none") + 
  labs(tag = "B")

g3 <- pheatmap(averaged_miR_values, 
               show_rownames = FALSE, 
               cluster_rows = FALSE, 
               annotation_row = annotation_df,
               annotation_legend = FALSE,
               annotation_names_row = FALSE,
               annotation_colors = list(dnp = c("purple", "yellow")),
               fontsize = 18)

g3 <- grid.arrange(g3[[4]], 
                   left = grid::textGrob(paste0("\U0394", "G"),
                                         gp = grid::gpar(fontsize = 20)))
g_blank <- ggplot() + theme_void()
g3 <- plot_grid(g3, g_blank, ncol = 2, rel_widths = c(0.9,0.1), 
                labels = c("C", ""), label_size = 26)
#add arrows


gA <- g1 + annotation_custom(
  ggplotGrob(g2), xmin = 2, xmax = 15, ymin = 12250, ymax = 19750)
png(filename = "./figures/miROverviewPlot.png", width = 750, height = 1000)
g_final <- plot_grid(gA , g3,nrow = 2)
g_final
dev.off()
ggsave(plot = g_final, filename = "./figures/miROverviewPlot.eps", 
       device = "eps", width = 10, height = 9, units = "in")

########################Tables#####################

#The xtable code generates Latex tables
xtable(targets_matrix)
xtable(miR_cocktail_matrix)

#Need to make a table of the 144 distinct transcripts repressed
#need to limit to just the "target set" of the top 5% of mRNA by network potential
top_mir <- overall_targeting_mat.df %>% 
  filter(miR == "miR-3613-3p", targets_logical == 1)
top_mir_targets <- unique(as.character(top_mir$gene_name))
top_mir_targets_df <- enframe(top_mir_targets, value = "miR-3613-3p targets") %>% 
  dplyr::select(-name)
xtable(top_mir_targets_df)

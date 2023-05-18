############################# Plot Percentage Bar Plot ############################
################################# 2023.05.15 ######################################
################################### Yanwen ########################################
###################################################################################

library(Seurat)
options(Seurat.object.assay.version = "v5")
library(tidyverse)
######################### Set Path #########################

data_path <- "/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen"
plot_path <- file.path(data_path, "plot/cell_type_vs_region_per_histological_region_plot_20230515")
inter_data_path <- file.path(data_path, "intermediate_data")
# read data

UCI604.Xenium <- readRDS(file = file.path(data_path, "UCI604.Xenium_with_hist_region.rds"))

plot_df <- UCI604.Xenium@meta.data %>% select(region, big_hist_region, new_names) %>%
  dplyr::rename("cell_types" = "new_names") %>%
  filter(cell_types != "mixed epithelial")
plot_df$cell_types <- factor(plot_df$cell_types, 
                             levels = c("Adipocytes", "B cells", "Basal",
                                        "Fibroblasts", "LumHR", "LumSec",
                                        "Lymphatic", "Mast cells", "Myeloid",
                                        "Pericytes", "T cells", "Vascular"))

pdf(file.path(plot_path, "duct_perc_bar.pdf"),
    width = 12, height = 5)
plot_df %>% filter(big_hist_region == "duct") %>%
  group_by(region, cell_types) %>%
  count() %>%
  group_by(region) %>%
  mutate(percentage = n/sum(n)) %>%
  ggplot(aes(x = region, y = percentage, fill = cell_types)) +
  geom_bar(position = 'fill', stat = 'identity')+ 
  #scale_fill_manual(values=col_vector[1:n])+
  #scale_fill_manual(values = selected_color)+
  theme_bw() +
  ggtitle("Duct") + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
dev.off()
plot_df %>% filter(big_hist_region == "duct") %>%
  group_by(region, cell_types) %>%
  count()%>%
  as.data.frame()%>%
  write.csv(file.path(inter_data_path, "duct_cell_count_by_region_celltype.csv"))




pdf(file.path(plot_path, "adipose_perc_bar.pdf"),
    width = 12, height = 5)
plot_df %>% filter(big_hist_region == "adipose") %>%
  group_by(region, cell_types) %>%
  count() %>%
  group_by(region) %>%
  mutate(percentage = n/sum(n)) %>%
  ggplot(aes(x = region, y = percentage, fill = cell_types)) +
  geom_bar(position = 'fill', stat = 'identity')+ 
  #scale_fill_manual(values=col_vector[1:n])+
  #scale_fill_manual(values = selected_color)+
  theme_bw() +
  ggtitle("Adipose") + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
dev.off()
plot_df %>% filter(big_hist_region == "adipose") %>%
  group_by(region, cell_types) %>%
  count()%>%
  as.data.frame()%>%
  write.csv(file.path(inter_data_path, "adipose_cell_count_by_region_celltype.csv"))





pdf(file.path(plot_path, "connective_tissue_perc_bar.pdf"),
    width = 12, height = 5)
plot_df %>% filter(big_hist_region == "connective_tissue") %>%
  group_by(region, cell_types) %>%
  count() %>%
  group_by(region) %>%
  mutate(percentage = n/sum(n)) %>%
  ggplot(aes(x = region, y = percentage, fill = cell_types)) +
  geom_bar(position = 'fill', stat = 'identity')+ 
  #scale_fill_manual(values=col_vector[1:n])+
  #scale_fill_manual(values = selected_color)+
  theme_bw() +
  ggtitle("Connective Tissue") + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
dev.off()
plot_df %>% filter(big_hist_region == "connective_tissue") %>%
  group_by(region, cell_types) %>%
  count()%>%
  as.data.frame()%>%
  write.csv(file.path(inter_data_path, "connective_tissue_cell_count_by_region_celltype.csv"))





pdf(file.path(plot_path, "lobule_perc_bar.pdf"),
    width = 12, height = 5)
plot_df %>% filter(big_hist_region == "lobule") %>%
  group_by(region, cell_types) %>%
  count() %>%
  group_by(region) %>%
  mutate(percentage = n/sum(n)) %>%
  ggplot(aes(x = region, y = percentage, fill = cell_types)) +
  geom_bar(position = 'fill', stat = 'identity')+ 
  #scale_fill_manual(values=col_vector[1:n])+
  #scale_fill_manual(values = selected_color)+
  theme_bw() +
  ggtitle("Lobule") + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
dev.off()

plot_df %>% filter(big_hist_region == "lobule") %>%
  group_by(region, cell_types) %>%
  count()%>%
  as.data.frame()%>%
  write.csv(file.path(inter_data_path, "lobule_cell_count_by_region_celltype.csv"))

# ## output meta data for neighborhood analysis 
# UCI604.Xenium@meta.data %>% #select(region, big_hist_region, new_names) %>%
#   dplyr::rename("cell_types" = "new_names") %>%
#   filter(cell_types != "mixed epithelial")%>%head()

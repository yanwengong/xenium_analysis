###################################################################################
############################# Add Regions to Object Meta ##########################
################################# 2023.05.10 ######################################
################################### Yanwen ########################################
###################################################################################


######################### install packages #########################
install.packages("remotes")

# remotes::install_github("satijalab/seurat")
install.packages("Matrix")
install.packages("vctrs")
install.packages("sctransform")
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
## updated seurat_object package 

# remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
# remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)

# remotes::install_github("bnprks/BPCells", quiet = TRUE)


library(Seurat)
options(Seurat.object.assay.version = "v5")
library(tidyverse)


######################### Set Path #########################

data_path <- "/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen"

region_path <- file.path(data_path, "Cell_coordinates_with_region")

######################### Read in Data #########################

UCI604.Xenium <- readRDS(file.path(data_path, "UCI604.Xenium.rds"))

cell_file <- list.files(region_path, "*.csv")
cell_file_name <- strsplit(cell_file, "_cell.csv") %>% unlist()
cell_file_name <- gsub("_", "", cell_file_name)

## read in cell vs hist_region

rm(all_sample_region)

for (i in seq(length(cell_file))){
  file_name <- cell_file[i]
  sample_name <- cell_file_name[i]
  df <- read.csv(file.path(region_path, file_name)) %>%
    dplyr::rename(hist_region = region,  big_hist_region = big_region) %>%
    select(-X) %>%
    mutate(region = sample_name)%>%
    mutate(cell_ident = paste(region, cell_id, sep = "_"))
  
  if (exists("all_sample_region")){
    all_sample_region <- rbind(all_sample_region, df)
  }
  else{
    all_sample_region <- df
  }
}


################### Merge hist region info to meta ######################
meta_cell_id <- UCI604.Xenium@meta.data %>% rownames_to_column(var = "cell_ident") %>% select(cell_ident)
comb <- dplyr::left_join(meta_cell_id, all_sample_region, by = "cell_ident")

################### Add hist region info to meta ######################

## add hist_region
hist_region <- comb$hist_region
names(hist_region) <- comb$cell_ident
UCI604.Xenium <- AddMetaData(
  object = UCI604.Xenium,
  metadata = hist_region,
  col.name = 'hist_region'
)

## add big_hist_region
big_hist_region <- comb$big_hist_region
names(big_hist_region) <- comb$cell_ident
UCI604.Xenium <- AddMetaData(
  object = UCI604.Xenium,
  metadata = big_hist_region,
  col.name = 'big_hist_region'
)



#############################  save data #################

saveRDS(UCI604.Xenium, file = file.path(data_path, "UCI604.Xenium_with_hist_region.rds"), compress = FALSE)

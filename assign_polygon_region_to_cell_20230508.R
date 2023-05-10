###################################################################################
######################### Xenium assign cell to polygon region ####################
################################# 2023.05.08 ######################################
################################### Yanwen ########################################
###################################################################################

library(sp)
library(tidyverse)
input_path <- "/share/crsp/lab/dalawson/share/HBCA_Phase2/DeepDive_MP_JIR/Jupyter/Yanwen"

cell_coor_path <- file.path(input_path, "Cell_coordinates")
polygon_coor_path <- file.path(input_path, "Polygon_coordinates")
output_path <- file.path(input_path, "Cell_coordinates_with_region")
plot_path <- file.path(input_path, "Cell_region_plot")

####################### functions #######################
check_file_is_match <- function(cell_file_name, poly_file_name){
  if (length(cell_file_name) != length(poly_file_name)) {
    stop("File Number Is Different!")
  }
  
  for (i in length(cell_file_name)){
    if (cell_file_name[i] != poly_file_name[i]){
      stop("File Name Does Not Match!")
    } 
  }
  print("Cell Files Match With Polygon Files.")
  
}

assign_region<- function(poly, cell, r, cell_update){
  pol_x = poly %>% filter(Selection == r) %>% pull(X)
  pol_y = poly %>% filter(Selection == r) %>% pull(Y)
  
  point_x = cell%>%pull(x_centroid)
  point_y = cell%>%pull(y_centroid)
  
  # check if point in polygon: 0 no, 1 yes
  res = point.in.polygon(point.x = point_x, point.y = point_y, pol.x = pol_x, pol.y = pol_y)
  
  in_index = which(res == 1)
  
  cell_update <- cell_update %>% mutate(region = ifelse(cell_id %in% in_index, r, region))
  
}



####################### get list of the sample name #######################

cell_file <- list.files(cell_coor_path, "*.csv")
poly_file <- list.files(polygon_coor_path, "*.csv")


cell_file_name <- strsplit(cell_file, "_cells.csv") %>% unlist()
poly_file_name <- strsplit(poly_file, "_coordinates.csv") %>% unlist()

check_file_is_match(cell_file_name, poly_file_name)


####################### generate cell to region df and plot #######################
for (sample_name in cell_file_name){
  print(sample_name)
  
  # read file
  cell <- read.csv(file.path(cell_coor_path, paste(sample_name, "cells.csv", sep = "_")))
  poly <- read.csv(file.path(polygon_coor_path, 
                             paste(sample_name, "coordinates.csv", sep = "_")), 
                   skip = 2)
  
  # get unique region name list
  region_names <- unique(poly$Selection)
  region_names <- region_names[region_names != "Tissue"]
  
  cell_update <- cell %>% mutate(region = "temp")
  
  # iterate over region
  for (r in region_names){
    print(r)
    cell_update <- assign_region(poly, cell, r, cell_update)
  }
  cell_update$region <- gsub("temp", "Tissue", cell_update$region)
  ## give high level regions: 
  
  cell_update <- cell_update %>% 
    mutate(big_region = case_when(
      str_starts(region, "A") ~ "adipose",
      str_starts(region, "D") ~ "duct",
      str_starts(region, "L") ~ "lobule",
      TRUE ~ "connective_tissue"
    ))
  
  write.csv(cell_update, file.path(output_path, paste(sample_name, "cell.csv", sep = "_")))
  
  pdf(file.path(plot_path, paste(sample_name, "cell_colored_by_region.pdf", sep = "_")),
      width = 8, height = 5)
  print(ggplot(cell_update, aes(x=x_centroid, y=y_centroid, color=big_region)) +
          geom_point(size=0.2)+theme_bw())
  dev.off()
  
}




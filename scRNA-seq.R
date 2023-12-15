library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(tidyverse)
library(gridExtra)

setwd("C:/Users/dimit/Desktop/12342")
getwd()

directories <- list.dirs(path = 'data/', recursive = F, full.names = F)

# Create a for loop and in each iteration
# it should go in each folder and fetch
# the files to create counts matrix
# and then using the count matrix
# we create a Seurat object and assign it
# a variable that has both the patient 
# information as well as the type of the tissue.
### Remember that it takes some time to loop all the folders ###

for(x in directories){
  name <- gsub('_filtered_feature_bc_matrix','', x)

  counts <- ReadMtx(mtx = paste0('data/',x,'/matrix.mtx.gz'),
                 features = paste0('data/',x,'/features.tsv.gz'),
                 cells = paste0('data/',x,'/barcodes.tsv.gz'))

  # create seurat objects
  assign(name, CreateSeuratObject(counts = counts))
}


# Merge these objects together into one object
# Why we merge them?
# Because there is a big chance to 
# make mistakes by proceeding to the 
# quality control immediately after
# creating the Seurat object we prefer
# to merge all the objects into one!!
# very important step

# Under in the console type ls()
# and check the names of all
# the objects that have created.
# Based on them merge the files

merge_seurat <- merge(HB17_background, y=c(HB17_PDX, HB17_tumor,
                                           HB30_PDX, HB30_tumor,
                                           HB53_background, HB53_tumor),
                      # add.cell.ids --> we use it because when we merge we want to know
                      # what
                      add.cell.ids = ls()[3:9],
                      project = "HB")

merge_seurat
# the result is 33538 features (genes)
# across 77936 samples(cells)

# Preprocesing the data for QC & Trimming

view(merge_seurat@meta.data)

# Based on the information now 
# I want to create a column with what patient and
# tissue is the cell barcode origin from.
# Create a sample column


merge_seurat$sample <- rownames(merge_seurat@meta.data)
view(merge_seurat@meta.data)


# Split the column (Sample)
# that we just created.
# ??separate()

merge_seurat@meta.data <- separate(merge_seurat@meta.data,
                                   col = sample,
                                   into= c("Patient_ID","Tissue","Barcodes"),
                                   sep = '_')
view(merge_seurat@meta.data)

# Check that everything has merged perfectly
unique(merge_seurat@meta.data$Patient_ID)
unique(merge_seurat@meta.data$Tissue)


## QC & Trimming

# How to filter low quality cells?
# The seurat function PercentageFeatureSet()
# enables you to calculate the mitochondrial transcript

merge_seurat$mitoPer <- PercentageFeatureSet(merge_seurat, pattern = "^MT-")
view(merge_seurat$mitoPer)

## EXPLORE THE DATA##

# Filtering 

# Based on the paper:
# Cells with fewer than 500 
# expressed genes or 800 UMIs, 
# or greater than 10% 
# mitochondrial counts were removed

merge_seurat_filtered <- subset(merge_seurat, subset = nCount_RNA > 800 &
       nFeature_RNA > 500 &
         mitoPer < 10)
merge_seurat_filtered

# How can I find how many genes I had before filtering?

merge_seurat


# Now we want to figure out if we have any batch effects
# Normalization


merge_seurat_filtered <- NormalizeData(object = merge_seurat_filtered)
merge_seurat_filtered <- FindVariableFeatures(object = merge_seurat_filtered)
merge_seurat_filtered <- ScaleData(object = merge_seurat_filtered)
merge_seurat_filtered <- RunPCA(object = merge_seurat_filtered)
# find the dimension of the dataset
ElbowPlot(merge_seurat_filtered)
merge_seurat_filtered <-FindNeighbors(object = merge_seurat_filtered)

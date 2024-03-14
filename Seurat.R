#Limpia informaci?n actual
rm(list=ls(all=TRUE))

#Seurat 
install.packages("cpp11")
install.packages("systemfonts")
install.packages("textshaping")
install.packages("ragg")
install.packages("tidyverse")
install.packages("igraph")
install.packages("Matrix")
library(Matrix)
library(Seurat)
library(data.table)
library(tidyverse)
library(magrittr)
library(gridExtra)
library(DoubletFinder)
#Control de calidad
#UMI por celda, igual o mayor a 500
#Minimo 200 genes y maximo 35000 detectados por celula 
#% mitoncdrial igual o mayor a 20%
#Eliminar dobletes
#Minimo 500 celulas
# Directorio que contiene los archivos de datos
data_dir <- '/media/paca/Acer/Datos_tesis/CONTROLES/Control_362/outs/filtered_feature_bc_matrix' #PC
#data_dir <- '//media/datos/paca/features_matrix/filtered_feature_bc_matrix3838' #Servidor

# Lista los archivos en el directorio de datos
list.files(data_dir)

# Cargar datos de expresión génica y de captura de anticuerpos
data <- Read10X(data.dir = data_dir)

# Recargar o recrear el objeto Seurat
seurat_object <- CreateSeuratObject(counts = data)

# Asignar un nombre al proyecto utilizando el slot project.name
seurat_object@project.name <- "C362"

# Obtener el nombre del proyecto
project_name <- seurat_object@project.name

class(seurat_object)

# Si hay datos de captura de anticuerpos disponibles, agrégalos al objeto Seurat
if ("Antibody Capture" %in% names(seurat_object@assays)) {
  seurat_object[["Protein"]] <- CreateAssayObject(counts = seurat_object[["Antibody Capture"]])
}

# Obtener los datos de expresión génica del Assay RNA
rna_counts <- GetAssayData(seurat_object, assay = "RNA")

# Calcular el número de genes detectados por célula
nFeature_RNA <- colSums(rna_counts > 0)

# Agregar los metadatos al objeto Seurat
seurat_object <- AddMetaData(seurat_object, metadata = data.frame(nFeature_RNA = nFeature_RNA))

# Calcular la proporción de expresión génica relacionada con mitocondrias
seurat_object [["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Filtrar células basadas en el porcentaje de expresión génica relacionada con mitocondrias
C362 <- subset(seurat_object, 
                  subset = nFeature_RNA > 200 & nFeature_RNA <= 35000 &
                    percent.mt < 20 & 
                    seurat_object@meta.data$nCount_RNA > 500)


# Escalar y normalizar los datos
seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object)


# Extraer los metadatos de las células filtradas
filtered_metadata <- as.data.frame(C362@meta.data)

# Guardar los metadatos filtrados como un archivo CSV
write.csv(filtered_metadata, file = "C362_filtered_metadata.csv", row.names = TRUE)


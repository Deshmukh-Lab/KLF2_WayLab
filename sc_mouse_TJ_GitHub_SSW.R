#Full analysis for murine single-cell data seen in Fig. 1 heatmap#
#Shao et al., unpublished#
#compiled by: Tony Jiang#
#Last updated: 10/26/2021#
#For GitHub#

#Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(msigdbr)
library(presto)
library(devtools)
library(fgsea)
library(tibble)
library(ggplot2)
library(SingleR)
library(celldex)
library(gprofiler2)
library(SingleR)
library(monocle)
library(corrplot)
library(scales)
library(ggrepel)

setwd('/pathtowd')

# Load the dataset- can be found in subfolder "Selected Data"
data = read.csv("/pathtodata/exp.filtermerge.csv", header = T, row.names = 1)
ensb = rownames(data)

# convert from ensembl to gene names
geneInfo.full = gconvert(query = ensb, organism = "mmusculus", target = "ENSG", filter_na = T)
geneNames = data.frame(geneInfo.full[,"name"]) 
geneNames.upper <- mutate_each(geneNames, funs(toupper))

# remove NaN's and duplicates
nan.indices = which(geneNames.upper[,1] %in% "NAN")
dupes.indices = which(duplicated(geneNames.upper) %in% "TRUE")
merge_indices = unique(c(nan.indices, dupes.indices))

filtered.data = data[-merge_indices,]
rownames(filtered.data) = geneNames.upper[-merge_indices,]
#fwrite(t(filtered.data), row.names = T, file = "/pathtosave/scRNA_filtered.csv")


# separate into T helper lineages
thresh = 0.5
th1 = which(filtered.data["TBX21",] > thresh)
th2 = which(filtered.data["IRF4",] > thresh)
tfh = which(filtered.data["CXCR5",] > thresh)
treg = which(filtered.data["FOXP3",] > thresh)
th17 = which(filtered.data["RORC",] > thresh)

th.labels = rep("Tmys",ncol(filtered.data))
th.labels[th1] = "Th1"
th.labels[tfh] = "Tfh"
th.labels[treg] = "Treg"
th.labels[th2] = "Th2"
th.labels[th17] = "Th17"

# Initialize the Seurat object with the raw data.
seurat = CreateSeuratObject(counts = filtered.data, project = "klf", min.cells = 0, min.features = 0)
seurat$th = th.labels
seurat$th = factor(seurat$th, c("Th17","Th2","Treg","Tfh","Tmys"))

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >2% mitochondrial counts
seurat <- subset(seurat, subset = nFeature_RNA > 3000 & nFeature_RNA < 6500 & percent.mt < 2)

# Feature selection. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# scale data
all.genes = rownames(seurat)
seurat = ScaleData(seurat, features = all.genes)

# perform PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

#compare across lineages
geneNum = 50
Idents(seurat) = seurat$th
seurat.markers <- FindAllMarkers(seurat, only.pos = T, min.pct = 0.2, logfc.threshold = 0.5)
seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
topN = seurat.markers %>% group_by(cluster) %>% top_n(n = geneNum, wt = avg_log2FC)

textSize = 12
DoHeatmap(seurat, features = topN$gene, size = textSize/2.5, 
          lines.width = 1, group.bar = F, label = F) + 
  theme(axis.text.y = element_blank()) +
  theme(text = element_text(size = textSize*1.5)) + #+ NoLegend()
  theme(legend.position = c(0.5,-0.07)) +
  theme(legend.direction = "horizontal") +
  theme(legend.title = element_text(size = textSize)) +
  theme(legend.text = element_text(size = textSize)) +
  theme(aspect.ratio = 0.68) +
  theme(plot.margin=unit(c(0,0,0,20),"mm")) + # top, right, bottom, left
  guides(color = "none") +
  guides(fill = guide_colorbar(title = expression(paste("Expression ", italic("Z"), "-score")), 
                               frame.colour = "black", frame.linewidth = 1, 
                               ticks.color = "black", ticks.linewidth = 1,
                               title.vjust=0.9, barwidth = unit(3, "cm"), barheight = unit(0.5, "cm"))) 


#ggsave('heatmap_lineages.png', width = 22, height = 16, bg="transparent", dpi = 600, units = c("cm")) 




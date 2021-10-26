#Full analysis for CD and UC human datasets#
#Shao et al., unpublished#
#compiled by: Jake Stevens#
#Last updated: 10/26/2021#
#For GitHub#

#Load main libraries#

library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(magrittr)
library(ggplot2)

#Set directory

setwd("/path_to_working_directory")

#Code used to make CD object from Teichmann dataset
library(SeuratDisk)

#Convert H5ad Crohn's Disease Dataset to a Seurat dataset#
#This is the dataset from Teichmann that can be downloaded directly at gutcellatlas.com
Convert("/Users/steqz7/Desktop/pediatric_RAWCOUNTS_cellxgene.h5ad", dest = "h5seurat", overwrite = TRUE)

#Load converted h5Seurat data into Seurat object
CD <- LoadH5Seurat("/Users/steqz7/Desktop/pediatric_RAWCOUNTS_cellxgene.h5seurat")

CD
head(CD[[]])

#Check QC Data#
VlnPlot(CD, features = c("n_genes", "n_counts", "percent_mito"), ncol = 3)

plot1 <- FeatureScatter(CD, feature1 = "n_counts", feature2 = "percent_mito")
plot2 <- FeatureScatter(CD, feature1 = "n_counts", feature2 = "n_genes")
plot1 + plot2

#Normalize Data#
CD <- NormalizeData(CD, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify highly variable features#
CD <- FindVariableFeatures(CD, selection.method = "vst", nfeatures = 2000)

#Identify 10 most highly variable genes#
top10 <- head(VariableFeatures(CD), 10)

#Plot variable features
plot1 <- VariableFeaturePlot(CD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data#
all.genes <- rownames(CD)
CD <- ScaleData(CD, features = all.genes)

#Perform linear dimensional reduction#
CD <- RunPCA(CD, features = VariableFeatures(object = CD))

#Assess dimensionality of dataset
ElbowPlot(CD)
#Will use 15 dimensions for downstream

#Cluster the cells
CD <- FindNeighbors(CD, dims = 1:15)
CD <- FindClusters(CD, resolution = 0.5)

CD <- RunUMAP(CD, dims = 1:15)

#Visualize
Idents(CD) <- "annotation_V2"
DimPlot(CD)

##LOAD OBJECTS##
#CD- Teichmann dataset from Elmentaite et al., Dev Cell, 2020

#Load if needed
CD <- readRDS("CD.RDS") #read in RDS we made from this study
CD #22502 samples 
head(CD[[]])
Idents(CD) <- "Sample.name"
levels(CD)
table(CD$Sample.name, CD$Diagnosis)
#8 CD (T017, T019, T160, T176, T189, T197, T202, and T203) 
#7 Healthy (T024, T036, T44, T057, T110POS, T161, T182)
#There should be 8 healthy controls ans 7 CD samples according to Dr. Teichmann
#See below

#UC- Regev dataset from Smillie et al., Cell, 2019
UC <- readRDS("UC_full.RDS") #Note: this is FULL dataset from the study above that includes training and discovery sets
UC #210614 cells
head(UC[[]])
Idents(UC) <- "Subject"
table(UC$Subject, UC$Health)
#18 UC (N106, N110, N111, N12, N14, N19, N23, N24, N26, N44, N49, N50, N52, N539, N58, N661, N7, N9) 
#12 Healthy (N10, N11, N13, N15, N16, N17, N18, N20, N21, N46, N51, N8)

#Take steps for both datasets to add necessary metadata and refine#
#CD, this is what was mentioned above
#T160 is actually healthy, so make new metadata incorporating that
Idents(CD) <- "Sample.name"
levels(CD)

CD <- RenameIdents(CD, "T017" = "Crohn Disease", "T019" = "Crohn Disease", "T024" = "Healthy control", "T036" = "Healthy control", "T44" = "Healthy control", "T057" = "Healthy control", "T110POS" = "Healthy control", "T160" = "Healthy control", "T161" = "Healthy control", "T176" = "Crohn Disease", "T182" = "Healthy control", "T189" = "Crohn Disease", "T197" = "Crohn Disease", "T202" = "Crohn Disease", "T203" = "Crohn Disease") 
levels(CD)
CD[["final.health"]] <- Idents(CD)
#Check
table(CD$Sample.name, CD$final.health)
#7 CD (T017, T019, T176, T189, T197, T202, and T203) 
#8 Control (T024, T036, T44, T057, T110POS, T161, T182, T160)
#Perfect

#UC#
Idents(UC) <- "Health"
levels(UC)
#We only want to look at healthy versus inflammed, so we will remove non-inflamed
UC.final <- subset(UC, idents = c("Healthy", "Inflamed"))
UC.final #143202 samples
levels(UC.final)

Idents(UC.final) <- "Cluster"
levels(UC.final)

#This is just to better identify and group our clusters of interest
#We ended up using label transfer, so this did not end up being in final analysis
UC.final <- RenameIdents(UC.final, "CD4+ Activated Fos-hi" = "CD4 T cells", "CD4+ Activated Fos-lo" = "CD4 T cells", "CD4+ Memory" = "CD4 T cells", "CD4+ PD1+" = "CD4 T cells", "CD8+ IL17+" = "CD8 T cells", "CD8+ IELs" = "CD8 T cells", "CD8+ LP" = "CD8 T cells") 
levels(UC.final)
UC.final[["Final Clusters"]] <- Idents(UC.final)

head(UC.final[[]])
Idents(UC.final) <- "Health"
DimPlot(UC.final) + RestoreLegend()
table(UC.final$Subject, UC.final$Health)

#Cool#

##BEGIN ANALYSIS##

#Start here with CD-specific data#

# Plot tSNE
Idents(CD) <- "annotation_V2"
DimPlot(CD)+RestoreLegend()

#Save the RDS with the correct patient metadata if desired
saveRDS(CD, file = "CD_updated.rds")

#Collapse different cell states into common types, again to get
#our particular clusters of interest
Idents(CD) <- "annotation_V2"
levels(CD)
cell.state <- CD$annotation_V2
CD <- AddMetaData(CD, metadata = cell.state, col.name = "cell.state")
Idents(CD) <-"cell.state"
levels(CD)

CD <- RenameIdents(CD, "crypt" = "non-immune", "TA" = "non-immune", "early enterocyte" = "non-immune","enterocyte" = "non-immune", "enteroendocrine" = "non-immune", "BEST4 enterocyte" = "non-immune", "Goblet cell" = "non-immune", "IL2RG+ enterocyte (M cell)" = "non-immune",
                   "Paneth cell" = "non-immune",  "Tuft" = "non-immune", "S1 fibroblasts" = "non-immune",  "S2 fibroblasts" = "non-immune","S4 fibroblasts" = "non-immune", "myofibroblast" = "non-immune",  
                   "Glial cell" = "non-immune", "pericyte" = "non-immune", "Arterial endothelial cell" = "non-immune",
                   "Venous endothelial cell"= "non-immune",  "Lymphatic endothelial cell" = "non-immune", "B cell" = "B cell", "FCER2 B cell" = "B cell",  "Memory B cell" = "Memory B cell",
                   "Activated B cell"= "B cell",  "Cycling B cell" = "B cell", "Cycling plasma cell" = "Plasma cell", "IgA plasma cell" = "Plasma cell",  "IgG plasma cell" = "Plasma cell",
                   "CD8 T cell"= "CD8 T cell",  "CD4 T cell" = "CD4 T cell", "Tfh" = "Tfh", "Activated T" = "Activated T",  "Treg" = "Treg",
                   "gd T/NK cell"= "gd T/NK cell",  "mast cells" = "mast cells", "Monocyte" = "Monocyte", "Cycling myeloid cells" = "Cycling myeloid cells",  "Macrophage" = "Macrophage",
                   "cDC1"= "cDC1",  "cDC2" = "cDC2", "activated DC" = "activated DC", "pDC" = "pDC")

CD@meta.data$"cell.id" <- as.factor(Idents(CD))
table(CD$cell.id, CD$final.health)

#This is how we will designate cells on KLF2 expression
#We do it this way for easier counting and for maintaining this expression
#as a separate metadata column

# Add metadata column corresponding to different expression levels of KLF2
poscells.0.1 <- WhichCells(CD, expression = KLF2 >= 0.1)
poscells.1 <- WhichCells(CD, expression = KLF2 >= 1)
CD$KLF2.exp.0.1 <- ifelse(colnames(CD) %in% poscells.0.1, "KLF2.Pos.0.1", "KLF2.Neg.0.1")
CD$KLF2.exp.1 <- ifelse(colnames(CD) %in% poscells.1, "KLF2.Pos.1", "KLF2.Neg.1")
head(CD[[]])
#Now we should have two columns, one for cells with KLF2 expression > 0.1 (which)
#we didn't end up using and one for cells with KLF2 > 1, with each cell being
#either postive or negative based on expression of KLF2

#If you want you can also combine with the cell type identification, but 
#we didn't actually end up using this
# Create new metadata column to hold cell.id and KLF2.status
cell.id.KLF2.stats.0.1 <- paste(CD@meta.data$cell.id,CD@meta.data$KLF2.exp.0.1)
cell.id.KLF2.stats.1 <- paste(CD@meta.data$cell.id,CD@meta.data$KLF2.exp.1)
CD <- AddMetaData(CD, metadata = cell.id.KLF2.stats.0.1, col.name = "id.exp.0.1")
CD <- AddMetaData(CD, metadata = cell.id.KLF2.stats.1, col.name = "id.exp.1")
head(CD[[]])

#View KLF2 expression
VlnPlot(CD, features = "KLF2", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0.5, split.by = "Sample.name")
FeatureScatter(CD, feature1 = "KLF2", feature2 = "IL10")

#These were initial analyses looking at the T cells in the dataset... 
#We didn't end up using this, but it is retained here for completion
#Create T cells group#
Idents(CD) <- "cell.id"
levels(CD)
CD.Tcell <- subset(CD, idents = c("CD8 T cell", "CD4 T cell", "Tfh", "Activated T", "Treg"))
levels(CD.Tcell)
CD.Tcell #2268 cells

#Looking at various features in the KLF2+ T cells
Idents(CD.Tcell) <- "KLF2.exp.1"
lM<- subset(CD.Tcell, idents = "KLF2.Pos.1")
head(lM[[]])
Idents(lM) <-"final.health"
a <- FeatureScatter(lM, feature1 = "KLF2", feature2 = "IL10") #They really have no IL10 expression#
b <- FeatureScatter(lM, feature1 = "KLF2", feature2 = "TGFB1")
c <- FeatureScatter(lM, feature1 = "KLF2", feature2 = "PDCD1")
a+b+c
a
VlnPlot(subset(CD.Tcell, idents = "KLF2.Pos.1"), features = "IL10", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0.5)


# Subset CD4 T cells.
Idents(CD) <-"cell.id"
CD.CD4Tcells <- subset(CD, idents = "CD4 T cell")
CD.CD4Tcells #1258 cells#

head(CD.CD4Tcells[[]])
Idents(CD.CD4Tcells) <-"KLF2.exp.0.1"
levels(CD.CD4Tcells)
FeaturePlot(CD.CD4Tcells, features = "KLF2", split.by = "final.health", min.cutoff = "q1") + RestoreLegend()
FeaturePlot(CD.CD4Tcells, features = "IL10", split.by = "final.health", min.cutoff = "q1") + RestoreLegend()
View(rownames(CD.CD4Tcells))

# Order according to health status
Idents(CD.CD4Tcells) <- "Sample.name"
levels(CD.CD4Tcells)
my.levels = c("T017","T019","T024","T036","T057","T160","T161","T176","T182","T189","T197","T202","T203")
#Note: there are only 13 samples above because 2 just don't have any CD4 T cells#
CD.CD4Tcells$Sample.name <- factor(CD.CD4Tcells$Sample.name, levels = my.levels)

##A first attempt at plotting frequencies and looking at these KLF2+ T cells
#Again, did not end up using this, but good to see

# Now Plot KLF2 expression for CD4+ T cells

plot1 <- VlnPlot(CD.CD4Tcells, features = "KLF2", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0)
plot1 + geom_jitter(data = plot1$data, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.1), alpha =0.2) + geom_boxplot(width=0.5) + ggtitle("All CD4 T cells") + NoLegend() + xlab("") +
  ylab("KLF2 expression")

# Now Plot KLF2 expression for KLF2 expressing cells first at thereshold of KLF2>0.1 and KLF2>1
Idents(CD.CD4Tcells) <-"KLF2.exp.0.1"
levels(CD.CD4Tcells)
VlnPlot(subset(CD.CD4Tcells, idents = "KLF2.Pos.0.1"), features = "KLF2", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0.5, split.by = "Sample.name")
plot1 <- VlnPlot(subset(CD.CD4Tcells, idents = "KLF2.Pos.0.1"), features = "KLF2", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0)
P1.1 <- plot1 + geom_jitter(data = plot1$data, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.1), alpha =0.2) + geom_boxplot(width=0.1) + ggtitle("KLF2 pos cells 0.1") + NoLegend() + xlab("") +
  ylab("KLF2 expression")
P1.1

Idents(CD.CD4Tcells) <-"KLF2.exp.1"
levels(CD.CD4Tcells)
VlnPlot(subset(CD.CD4Tcells, idents = "KLF2.Pos.1"), features = "KLF2", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0.5, split.by = "Sample.name")
plot1 <- VlnPlot(subset(CD.CD4Tcells, idents = "KLF2.Pos.1"), features = "KLF2", group.by ="final.health", cols = c("slateblue", "salmon"), pt.size = 0)
P1.2 <- plot1 + geom_jitter(data = plot1$data, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.1), alpha =0.2) + geom_boxplot(width=0.1) + ggtitle("KLF2 pos cells 1") + NoLegend() + xlab("") +
  ylab("KLF2 expression")

P1.1 + P1.2

# Get cell numbers for KLF2 > 0.1 and KLF2>1
head(CD.CD4Tcells[[]])
table <- table(CD.CD4Tcells$Sample.name, CD.CD4Tcells$final.health, CD.CD4Tcells$cell.id, CD.CD4Tcells$KLF2.exp.0.1, CD.CD4Tcells$KLF2.exp.1)
t1 <- as.data.frame(table(CD.CD4Tcells$Sample.name, CD.CD4Tcells$final.health, CD.CD4Tcells$cell.id, CD.CD4Tcells$KLF2.exp.0.1, CD.CD4Tcells$KLF2.exp.1))
head(t1)
write.csv(table, file = "data1.csv")
t2 <- t1[t1$Var3 %in% c("CD4 T cell"),] # To only plot CD4 Tcells.
head(t2)
t2$Var5 <- as.factor(t2$Var5)
t2$Var4 <- as.factor(t2$Var4)
t2$Var3 <- as.factor(t2$Var3)
t2$Var2 <- as.factor(t2$Var2)
t2$Var1 <- as.factor(t2$Var1)

# Rearrange by health status and only plot KLF2 positive cells
t2 <- t2 %>% arrange(t2$Var2)
t3 <- t2[t2$Var4 %in% c("KLF2.Pos.0.1"),]
t4 <- t2[t2$Var5 %in% c("KLF2.Pos.1"),]
head(t4)

#Different frequency plots now, again did not use these in final version
#Frequency of all CD4 T cells
A <- ggplot(t2, aes(x= Var3, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all T cells") + ggtitle("CD4 Tcells")

#Frequency of all KLF2 cells
B <- ggplot(t3, aes(x= Var4, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-0.1") + ggtitle("KLF2>0.1 CD4 Tcells")

#Frequency of all KLF2 cells
C <- ggplot(t4, aes(x= Var5, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-1") + ggtitle("KLF2> 1 CD4 Tcells")

A+B+C


#Frequency of all CD4 T cells, split by subjects
D <- ggplot(t2, aes(x= Var1, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all T cells") + ggtitle("CD4 Tcells")

#Frequency of all KLF2 cells
E <- ggplot(t3, aes(x= Var1, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-0.1") + ggtitle("KLF2>0.1 CD4 Tcells")

#Frequency of all KLF2 cells
G <- ggplot(t4, aes(x= Var1, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-1") + ggtitle("KLF2> 1 CD4 Tcells")

D+E+G

# Since there appears to be not a lot of difference in numbers, I will look at quantitative difference using cut off of KLF2>1
Idents(CD.CD4Tcells) <-"KLF2.exp.1"
KLF2.Tcells.CD <- subset(CD.CD4Tcells, idents ="KLF2.Pos.1")
KLF2.Tcells.CD # 408 samples within 1 assay

#First look at DEG, for quick pass look 
#Note: this DEG analysis was not used in final analysis- for that see section on making GO plots
Idents(KLF2.Tcells) <-"final.health"
levels(KLF2.Tcells)
MAST.DE <- FindMarkers(KLF2.Tcells, ident.1 = "Crohn Disease", ident.2 = "Healthy control")
View(MAST.DE)
write.csv(MAST.DE, file = "DE_genes_CD.csv", col.names = T)


##BEGIN UC ANALYSIS##
#This is the same as above, and again, not the final analysis

#Start here with UC-specific data#
head(UC.final[[]])
table(UC.final$Final.Clusters, UC.final$Health)
Table <- prop.table(table(UC.final$Final.Clusters, UC.final$Health, UC.final$Subject))
View(Table)

# Plot tSNE
Idents(UC.final) <- "Final.Clusters"
DimPlot(UC.final)+RestoreLegend()

table(UC.final$Final.Clusters, UC.final$Health)

# Add metadata column corresponding to different expression levels of KLF2
poscells.0.1 <- WhichCells(UC.final, expression = KLF2 >= 0.1)
poscells.1 <- WhichCells(UC.final, expression = KLF2 >= 1)
UC.final$KLF2.exp.0.1 <- ifelse(colnames(UC.final) %in% poscells.0.1, "KLF2.Pos.0.1", "KLF2.Neg.0.1")
UC.final$KLF2.exp.1 <- ifelse(colnames(UC.final) %in% poscells.1, "KLF2.Pos.1", "KLF2.Neg.1")
head(UC.final[[]])

#Create T cells group
Idents(UC.final) <- "Final.Clusters"
levels(UC.final)
UC.Tcell <- subset(UC.final, idents = c("CD4 T cells", "CD8 T cells", "Tregs")) #Leaving Cycling T out because its only ~400 cells total
levels(UC.Tcell)
UC.Tcell #46672 cells

Idents(UC.Tcell) <- "KLF2.exp.1"
lM<- subset(UC.Tcell, idents = "KLF2.Pos.1")
head(lM[[]])
Idents(lM) <-"Health"
a <- FeatureScatter(lM, feature1 = "KLF2", feature2 = "IL10") #They really have no IL10 expression#
b <- FeatureScatter(lM, feature1 = "KLF2", feature2 = "TGFB1")
c <- FeatureScatter(lM, feature1 = "KLF2", feature2 = "PDCD1")
a+b+c
a
VlnPlot(subset(UC.Tcell, idents = "KLF2.Pos.1"), features = "KLF2", group.by ="Health", cols = c("slateblue", "salmon"), pt.size = 0.5)

# Subset CD4 T cells
Idents(UC.final) <-"Final.Clusters"
UC.CD4Tcells <- subset(UC.final, idents = "CD4 T cells")
UC.CD4Tcells #30294 cells

#Look at features in KLF2+ CD4 T cells
head(UC.CD4Tcells[[]])
Idents(UC.CD4Tcells) <-"KLF2.exp.0.1"
levels(UC.CD4Tcells)
FeaturePlot(UC.CD4Tcells, features = "KLF2", split.by = "Health", min.cutoff = "q1") + RestoreLegend()
FeaturePlot(UC.CD4Tcells, features = "IL10", split.by = "Health", min.cutoff = "q1") + RestoreLegend()
View(rownames(UC.CD4Tcells))

# Order according to health status
Idents(UC.CD4Tcells) <- "Subject"
levels(UC.CD4Tcells)
my.levels = c("N7","N9","N10","N8","N11","N12","N13","N14","N15","N16","N17","N18","N19","N20","N21","N23","N24","N26","N51","N52","N58","N111","N661","N44","N46","N49","N50","N106","N539","N110")

UC.CD4Tcells$Subject <- factor(UC.CD4Tcells$Subject, levels = my.levels)

# Now Plot KLF2 expression for CD4 T cells
plot1 <- VlnPlot(UC.CD4Tcells, features = "KLF2", group.by ="Health", cols = c("slateblue", "salmon"), pt.size = 0)
plot1 + geom_jitter(data = plot1$data, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.1), alpha =0.2) + geom_boxplot(width=0.5) + ggtitle("All CD4 T cells") + NoLegend() + xlab("") +
  ylab("KLF2 expression")

# Now Plot KLF2 expression for KLF2 expressing cells first at threshold of KLF2>0.1 and KLF2>1
Idents(UC.CD4Tcells) <-"KLF2.exp.0.1"
levels(UC.CD4Tcells)
VlnPlot(subset(UC.CD4Tcells, idents = "KLF2.Pos.0.1"), features = "KLF2", group.by ="Health", cols = c("slateblue", "salmon"), pt.size = 0.5, split.by = "Subject")
plot1 <- VlnPlot(subset(UC.CD4Tcells, idents = "KLF2.Pos.0.1"), features = "KLF2", group.by ="Health", cols = c("slateblue", "salmon"), pt.size = 0)
P1.1 <- plot1 + geom_jitter(data = plot1$data, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.1), alpha =0.2) + geom_boxplot(width=0.1) + ggtitle("KLF2 pos cells 0.1") + NoLegend() + xlab("") +
  ylab("KLF2 expression")
P1.1

Idents(UC.CD4Tcells) <-"KLF2.exp.1"
levels(UC.CD4Tcells)
VlnPlot(subset(UC.CD4Tcells, idents = "KLF2.Pos.1"), features = "KLF2", group.by ="Health", cols = c("slateblue", "salmon"), pt.size = 0.5, split.by = "Subject")
plot1 <- VlnPlot(subset(UC.CD4Tcells, idents = "KLF2.Pos.1"), features = "KLF2", group.by ="Health", cols = c("slateblue", "salmon"), pt.size = 0)
P1.2 <- plot1 + geom_jitter(data = plot1$data, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.1), alpha =0.2) + geom_boxplot(width=0.1) + ggtitle("KLF2 pos cells 1") + NoLegend() + xlab("") +
  ylab("KLF2 expression")
P1.2

P1.1 + P1.2

# Get cell numbers for KLF2 > 0.1 and KLF2>1 ####
head(UC.CD4Tcells[[]])
table <- table(UC.CD4Tcells$Subject, UC.CD4Tcells$Health, UC.CD4Tcells$Final.Clusters, UC.CD4Tcells$KLF2.exp.0.1, UC.CD4Tcells$KLF2.exp.1)
t1 <- as.data.frame(table(UC.CD4Tcells$Subject, UC.CD4Tcells$Health, UC.CD4Tcells$Final.Clusters, UC.CD4Tcells$KLF2.exp.0.1, UC.CD4Tcells$KLF2.exp.1))
head(t1)
write.csv(table, file = "UCdata1.csv")
t2 <- t1[t1$Var3 %in% c("CD4 T cells"),] # To only plot CD4 Tcells.
head(t2)
t2$Var5 <- as.factor(t2$Var5)
t2$Var4 <- as.factor(t2$Var4)
t2$Var3 <- as.factor(t2$Var3)
t2$Var2 <- as.factor(t2$Var2)
t2$Var1 <- as.factor(t2$Var1)

# Rearrange by health status and only plot KLF2 positive cells
t2 <- t2 %>% arrange(t2$Var2)
t3 <- t2[t2$Var4 %in% c("KLF2.Pos.0.1"),]
t4 <- t2[t2$Var5 %in% c("KLF2.Pos.1"),]
head(t4)

#Frequency of all CD4 T cells
A <- ggplot(t2, aes(x= Var3, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all T cells") + ggtitle("CD4 Tcells")

#Frequency of all KLF2 cells
B <- ggplot(t3, aes(x= Var4, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-0.1") + ggtitle("KLF2>0.1 CD4 Tcells")

#Frequency of all KLF2 cells
C <- ggplot(t4, aes(x= Var5, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-1") + ggtitle("KLF2> 1 CD4 Tcells")

A+B+C


#Frequency of all CD4 T cells, split by subject s
D <- ggplot(t2, aes(x= Var1, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all T cells") + ggtitle("CD4 Tcells")

#Frequency of all KLF2 cells
E <- ggplot(t3, aes(x= Var1, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-0.1") + ggtitle("KLF2>0.1 CD4 Tcells")

#Frequency of all KLF2 cells
G <- ggplot(t4, aes(x= Var1, y = Freq, color = Var2)) + geom_boxplot() + scale_color_manual(values = c("slateblue", "salmon")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Number of all KLF2 pos T cells(exp-1") + ggtitle("KLF2> 1 CD4 Tcells")

D+E+G

#Since there appears to be no difference in numbers  of KLF2+ T cells, I will look at quantitative difference using cut off of KLF2>1
Idents(UC.CD4Tcells) <-"KLF2.exp.1"
KLF2.Tcells.UC <- subset(UC.CD4Tcells, idents ="KLF2.Pos.1")
KLF2.Tcells.UC #5577 samples within 1 assay

#Again, first look at DEG in UC. Not used in final analysis.
Idents(KLF2.Tcells) <-"Health"
levels(KLF2.Tcells)
MAST.DE <- FindMarkers(KLF2.Tcells, ident.1 = "Inflamed", ident.2 = "Healthy")
View(MAST.DE)
write.csv(MAST.DE, file = "DE_genes_UC.csv", col.names = T)


##KEY##
#This was where we first found the best way to assess frequencies of
#KLF2+CD4+ T cells between conditions
#Ended up doing by patient too to be able to copy to Prism for statistical analysis

table(CD.CD4Tcells$Sample.name, CD.CD4Tcells$KLF2.exp.1)
prop.table(table(CD.CD4Tcells$final.health, CD.CD4Tcells$KLF2.exp.1), margin = 1) #Use margin 1 to make the individual patients comparable

##FOR PRISM##
#This was the way to look at expression of KLF2 in cells in Prism
#i.e. make violin plots in Prism so we could make them look nicer and utilize statistics
#Did not end up using in final analysis

#Set the ident you want to pull- will need to pull Health status, KLF2 expression > 1, and Sample number
#Will need to run 3 times and make a combined Excel from printed .csv#
#Do for both CD and UC on the CD4 object

#For UC set
Idents(UC.CD4Tcells) <- "Health"
Idents(UC.CD4Tcells) <- "KLF2.exp.1"
Idents(UC.CD4Tcells) <- "Subject"

#This is the base expression to write
FetchData(UC.CD4Tcells, vars = c("KLF2", "ident"))

#This is it fully
write.csv(FetchData(UC.CD4Tcells, vars = c("KLF2", "ident")), file = "UC.Health.KLF2.csv", sep = ",")
write.csv(FetchData(UC.CD4Tcells, vars = c("KLF2", "ident")), file = "UC.KLF2.KLF2.csv", sep = ",")
write.csv(FetchData(UC.CD4Tcells, vars = c("KLF2", "ident")), file = "UC.Subject.KLF2.csv", sep = ",")


#For CD set
Idents(CD.CD4Tcells) <- "final.health"
Idents(CD.CD4Tcells) <- "KLF2.exp.1"
Idents(CD.CD4Tcells) <- "Sample.name"

#This is it
write.csv(FetchData(CD.CD4Tcells, vars = c("KLF2", "ident")), file = "CD.Health.KLF2.csv", sep = ",")
write.csv(FetchData(CD.CD4Tcells, vars = c("KLF2", "ident")), file = "CD.KLF2.KLF2.csv", sep = ",")
write.csv(FetchData(CD.CD4Tcells, vars = c("KLF2", "ident")), file = "CD.Subject.KLF2.csv", sep = ",")

#OK, so now that we have a handle on our basic processes for the Teichmann CD dataset
#and Regev UC dataset, we wanted to assess the other large, human CD datasets
#This includes one from the Colonna group (Jaeger et al., Nat Comm, 2021) and
#one that only included CD samples (inflammed and non-inflammed; Martin et al., 
#Cell, 2019) but from which we just used the inflammed samples

#First
#Bring in Colonna dataset

#Load all libraries that may be necessary to create Seurat object if not already loaded
library(Seurat)
library(dplyr)
library(tidyverse)
library(scater)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

#This is using the matrices downloaded from the GEO database for Jaeger et al., Nat Comm, 2021
#In this study, there were 2 CD patients and 2 Healthy controls, each with an IEL sample
#And an LP samples (4 total from each group)
#Non-inflammed samples were not included 

# Read count matrices and add metadata 
CD1 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/CD/filtered_gene_bc_matrices/GRCh38")
CD1 <-CreateSeuratObject(counts = CD1, min.cells = 1, min.features = 200) 
CD1@meta.data[, "id"] <- "CD_col_1"

HC1 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/HC/filtered_gene_bc_matrices/GRCh38")
HC1 <-CreateSeuratObject(counts = HC1, min.cells = 1, min.features = 200) 
HC1@meta.data[, "id"] <- "HC_col_1"

CD2 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/CD/filtered_gene_bc_matrices 2/GRCh38")
CD2 <-CreateSeuratObject(counts = CD2, min.cells = 1, min.features = 200) 
CD2@meta.data[, "id"] <- "CD_col_2"

HC2 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/HC/filtered_gene_bc_matrices 2/GRCh38")
HC2 <-CreateSeuratObject(counts = HC2, min.cells = 1, min.features = 200) 
HC2@meta.data[, "id"] <- "HC_col_2"

CD3 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/CD/filtered_gene_bc_matrices 3/GRCh38")
CD3 <-CreateSeuratObject(counts = CD3, min.cells = 1, min.features = 200) 
CD3@meta.data[, "id"] <- "CD_col_3"

HC3 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/HC/filtered_gene_bc_matrices 3/GRCh38")
HC3 <-CreateSeuratObject(counts = HC3, min.cells = 1, min.features = 200) 
HC3@meta.data[, "id"] <- "HC_col_3"

CD4 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/CD/filtered_gene_bc_matrices 4/GRCh38")
CD4 <-CreateSeuratObject(counts = CD4, min.cells = 1, min.features = 200) 
CD4@meta.data[, "id"] <- "CD_col_4"

HC4 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE157477_RAW/HC/filtered_gene_bc_matrices 4/GRCh38")
HC4 <-CreateSeuratObject(counts = HC4, min.cells = 1, min.features = 200) 
HC4@meta.data[, "id"] <- "HC_col_4"


# Merge all matrices and add metadata #
CD_col = merge(x = CD1, y = c(CD2, CD3, CD4))
CD_col@meta.data[, "health"] <- "CD"
HC_col = merge(x = HC1, y = c(HC2, HC3, HC4))
HC_col@meta.data[, "health"] <- "HC"
CD_colonna = merge(x = CD_col, y = HC_col)
CD_col
head(CD_col[[]])
Idents(CD_col) <- "id"
levels(CD_col)

CD_colonna
head(CD_colonna[[]])

#Now we have our CD Colonna Seurat object
#proceed to normal, basic analysis of object

CD_colonna[["percent.mt"]] <- PercentageFeatureSet(CD_colonna, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(CD_colonna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(CD_colonna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CD_colonna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

CD_colonna <- subset(CD_colonna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
CD_colonna #35019 samples 

CD_colonna <- NormalizeData(CD_colonna, normalization.method = "LogNormalize", scale.factor = 10000)

CD_colonna <- FindVariableFeatures(CD_colonna, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CD_colonna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(CD_colonna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#scale data
all.genes <- rownames(CD_colonna)
CD_colonna <- ScaleData(CD_colonna, features = all.genes)

#perform linear dimensional reduction
CD_colonna <- RunPCA(CD_colonna, features = VariableFeatures(object = CD_colonna))

ElbowPlot(CD_colonna) #Probably use 20 dimensions as the cutoff

CD_colonna <- FindNeighbors(CD_colonna, dims = 1:20)
CD_colonna <- FindClusters(CD_colonna, resolution = 0.5)

#Run UMAP
CD_colonna <- RunUMAP(CD_colonna, dims = 1:20)
DimPlot(CD_colonna, reduction = "umap")

#We now will transfer cell identities/labels from Teichmann dataset
#Note: we chose to use cell definitions from Elmentaite et al. because
#it is from the Gut cell atlas, which is part of the Human cell atlas
#Thus, having these labels apply to all objects for downstream analysis
#lends consistency

#Load libraries#
library(patchwork)

#Add better metadata to Colonna dataset right quick
#For identifying patients

Idents(CD_colonna) <- "id"
levels(CD_colonna)
#rename metadata to patient
CD_colonna <- RenameIdents(CD_colonna, "CD_col_1" = "CD1813", "CD_col_2" = "CD1818", "CD_col_3" = "CD1813", "CD_col_4" = "CD1818", "HC_col_1" = "HC7420", "HC_col_2" = "HC1425", "HC_col_3" = "HC7420", "HC_col_4" = "HC1425")
levels(CD_colonna)
CD_colonna[["Patient"]] <- Idents(CD_colonna)
head(CD_colonna[[]])

#Increase memory size- this label transfer is computationally expensive

options(future.globals.maxSize= 2500000000000)

#Cell type classification using Teichmann object as reference#
#Find transfer anchors that will use Teichmann (CD) as a reference for cell type and the Colonna set as a query#
CD.anchors <- FindTransferAnchors(reference = CD, query = CD_colonna, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = CD.anchors, refdata = CD$annotation_V2, dims = 1:30)

#Add the predicted cell types to the query, i.e. the IBD integrated object, metadata
CD_colonna <- AddMetaData(CD_colonna, metadata = predictions)

head(CD_colonna[[]])
DimPlot(CD_colonna, group.by = "predicted.id")

#Clean off some of the irrelevant metadata columns that accompany a label transfer

CD_colonna$prediction.score.CD4.T.cell <- NULL
CD_colonna$prediction.score.TA <- NULL
CD_colonna$prediction.score.IgA.plasma.cell <- NULL
CD_colonna$prediction.score.Goblet.cell <- NULL
CD_colonna$prediction.score.IgG.plasma.cell <- NULL
CD_colonna$prediction.score.Cycling.B.cell <- NULL
CD_colonna$prediction.score.S4.fibroblasts <- NULL
CD_colonna$prediction.score.crypt <- NULL
CD_colonna$prediction.score.early.enterocyte <- NULL
CD_colonna$prediction.score.Memory.B.cell <- NULL
CD_colonna$prediction.score.FCER2.B.cell <- NULL
CD_colonna$prediction.score.Tfh <- NULL
CD_colonna$prediction.score.cDC2 <- NULL
CD_colonna$prediction.score.Arterial.endothelial.cell <- NULL
CD_colonna$prediction.score.IL2RG..enterocyte..M.cell. <- NULL
CD_colonna$prediction.score.enterocyte <- NULL
CD_colonna$prediction.score.enteroendocrine <- NULL
CD_colonna$prediction.score.Venous.endothelial.cell <- NULL
CD_colonna$prediction.score.Activated.T <- NULL
CD_colonna$prediction.score.Monocyte <- NULL
CD_colonna$prediction.score.S1.fibroblasts <- NULL
CD_colonna$prediction.score.cDC1 <- NULL
CD_colonna$prediction.score.myofibroblast <- NULL
CD_colonna$prediction.score.pericyte <- NULL
CD_colonna$prediction.score.Treg <- NULL
CD_colonna$prediction.score.Lymphatic.endothelial.cell <- NULL
CD_colonna$prediction.score.B.cell <- NULL
CD_colonna$prediction.score.Activated.B.cell <- NULL
CD_colonna$prediction.score.pDC <- NULL
CD_colonna$prediction.score.Cycling.plasma.cell <- NULL
CD_colonna$prediction.score.gd.T.NK.cell <- NULL
CD_colonna$prediction.score.mast.cells <- NULL
CD_colonna$prediction.score.Paneth.cell <- NULL
CD_colonna$prediction.score.BEST4.enterocyte <- NULL
CD_colonna$prediction.score.activated.DC <- NULL
CD_colonna$prediction.score.S2.fibroblasts <- NULL
CD_colonna$prediction.score.CD8.T.cell <- NULL
CD_colonna$prediction.score.Macrophage <- NULL
CD_colonna$prediction.score.Glial.cell <- NULL
CD_colonna$prediction.score.Cycling.myeloid.cells <- NULL
CD_colonna$prediction.score.Tuft <- NULL
CD_colonna$prediction.score.max <- NULL
CD_colonna$prediction.match <- NULL
CD_colonna$health_meta <- NULL
CD_colonna$sample_meta <- NULL
CD_colonna$dataset <- NULL

#Note: we will assess the apparent validity of this label transfer later when we subset T cell groups
#But we trust the cell identities from Teichmann again, as it is published and part of the 
#Human cell atlas

head(CD_colonna[[]])
Idents(CD_colonna) <- "predicted.id"
levels(CD_colonna)

#save RDS now that it has predicted ids 
saveRDS(CD_colonna, file = "CD_colonna.RDS")

#Now we can do KLF2 analyses on the Colonna dataset
## Add metadata column corresponding to different expression levels of KLF2
poscells.0.1 <- WhichCells(CD_colonna, expression = KLF2 >= 0.1)
poscells.1 <- WhichCells(CD_colonna, expression = KLF2 >= 1)
CD_colonna$KLF2.exp.0.1 <- ifelse(colnames(CD_colonna) %in% poscells.0.1, "KLF2.Pos.0.1", "KLF2.Neg.0.1")
CD_colonna$KLF2.exp.1 <- ifelse(colnames(CD_colonna) %in% poscells.1, "KLF2.Pos.1", "KLF2.Neg.1")
head(CD_colonna[[]])

#Initial look in JUST CD4 T cells
# Subset CD4 T cells.
Idents(CD_colonna) <-"predicted.id"
CD_colonna.CD4Tcells <- subset(CD_colonna, idents = "CD4 T cell")
CD_colonna.CD4Tcells #9159 cells#

#for Prism#
prop.table(table(CD_colonna.CD4Tcells$health, CD_colonna.CD4Tcells$KLF2.exp.1), margin = 1)

#Now 3rd dataset
#Add inflammed cells from Martin et al., Cell, 2019, PMID 31474370#
#Again, this study only has CD samples with paired non-inflammed and inflammed samples
#For consistency, we just pulled out the Inflammed samples
#Now we have all the human CD datasets and enough of a number to analyze fully

#Again, load libraries if not already done
library(Seurat)
library(dplyr)
library(tidyverse)
library(scater)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

#Just pulled matrices from GEO for Martin et al., Cell, 2019
# Read count matrices and add metadata 
inf1 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_69")
inf1 <-CreateSeuratObject(counts = inf1, min.cells = 1, min.features = 200) 
inf1@meta.data[, "id"] <- "inflam_1"

inf2 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_122")
inf2 <-CreateSeuratObject(counts = inf2, min.cells = 1, min.features = 200) 
inf2@meta.data[, "id"] <- "inflam_2"

inf3 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_128")
inf3 <-CreateSeuratObject(counts = inf3, min.cells = 1, min.features = 200) 
inf3@meta.data[, "id"] <- "inflam_3"

inf4 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_138")
inf4 <-CreateSeuratObject(counts = inf4, min.cells = 1, min.features = 200) 
inf4@meta.data[, "id"] <- "inflam_4"

inf5 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_158")
inf5 <-CreateSeuratObject(counts = inf5, min.cells = 1, min.features = 200) 
inf5@meta.data[, "id"] <- "inflam_5"

inf6 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_181")
inf6 <-CreateSeuratObject(counts = inf6, min.cells = 1, min.features = 200) 
inf6@meta.data[, "id"] <- "inflam_6"

inf7 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_187")
inf7 <-CreateSeuratObject(counts = inf7, min.cells = 1, min.features = 200) 
inf7@meta.data[, "id"] <- "inflam_7"

inf8 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_190")
inf8 <-CreateSeuratObject(counts = inf8, min.cells = 1, min.features = 200) 
inf8@meta.data[, "id"] <- "inflam_8"

inf9 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_193")
inf9 <-CreateSeuratObject(counts = inf9, min.cells = 1, min.features = 200) 
inf9@meta.data[, "id"] <- "inflam_9"

inf10 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_196")
inf10 <-CreateSeuratObject(counts = inf10, min.cells = 1, min.features = 200) 
inf10@meta.data[, "id"] <- "inflam_10"

inf11 <- Read10X(data.dir = "/Users/steqz7/Desktop/SSW_analysis/GSE134809/Involved_209")
inf11 <-CreateSeuratObject(counts = inf11, min.cells = 1, min.features = 200) 
inf11@meta.data[, "id"] <- "inflam_11"

# Merge all matrices and add metadata #
inflamm = merge(x = inf1, y = c(inf2, inf3, inf4, inf5, inf6, inf7, inf8, inf9, inf10, inf11))
inflamm@meta.data[, "health"] <- "inflammed"

inflamm #58351 cells
head(inflamm[[]])
Idents(inflamm) <- "id"
levels(inflamm)

#Now we have our the Seurat object from Martin et al. study
#proceed to running clustering and basic analysis

inflamm[["percent.mt"]] <- PercentageFeatureSet(CD_colonna, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(inflamm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(inflamm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inflamm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

inflamm <- subset(inflamm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
inflamm #57175 samples 

inflamm <- NormalizeData(inflamm, normalization.method = "LogNormalize", scale.factor = 10000)

inflamm <- FindVariableFeatures(inflamm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(inflamm), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(inflamm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#scale data
all.genes <- rownames(inflamm)
inflamm <- ScaleData(inflamm, features = all.genes)

#perform linear dimensional reduction
inflamm <- RunPCA(inflamm, features = VariableFeatures(object = inflamm))

ElbowPlot(inflamm) #Probably use 18 dimensions as the cutoff

inflamm <- FindNeighbors(inflamm, dims = 1:18)
inflamm <- FindClusters(inflamm, resolution = 0.5)

#Run UMAP
inflamm <- RunUMAP(inflamm, dims = 1:18)
DimPlot(inflamm, reduction = "umap")

#Again, transfer labels from Teichmann dataset for consistent analysis 
#Will do an use the CD dataset from Teichmann (CD) as the reference object
#and the CD data set from Martin et al. as the query 

#First, add better patient metadata to dataset according to numbering from Martin et al. study

Idents(inflamm) <- "id"
levels(inflamm)
#rename metadata to patient
inflamm <- RenameIdents(inflamm, "inflam_1" = "inflamm69", "inflam_2" = "inflamm122", "inflam_3" = "inflamm128", "inflam_4" = "inflamm138", "inflam_5" = "inflamm158", "inflam_6" = "inflamm181", "inflam_7" = "inflamm187", "inflam_8" = "inflamm190", "inflam_9" = "inflamm193", "inflam_10" = "inflamm196", "inflam_11" = "inflamm209")
levels(inflamm)
inflamm[["Patient"]] <- Idents(inflamm)
head(inflamm[[]])


#Find transfer anchors that will use CD as a reference for cell type and the Martin et al. set as a query#
CD.anchors <- FindTransferAnchors(reference = CD, query = inflamm, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = CD.anchors, refdata = CD$annotation_V2, dims = 1:30)

#Add the predicted cell types to the query, i.e. the IBD integrated object, metadata
inflamm <- AddMetaData(inflamm, metadata = predictions)


#Clean off some of the irrelevant metadata columns
inflamm$prediction.score.CD4.T.cell <- NULL
inflamm$prediction.score.TA <- NULL
inflamm$prediction.score.IgA.plasma.cell <- NULL
inflamm$prediction.score.Goblet.cell <- NULL
inflamm$prediction.score.IgG.plasma.cell <- NULL
inflamm$prediction.score.Cycling.B.cell <- NULL
inflamm$prediction.score.S4.fibroblasts <- NULL
inflamm$prediction.score.crypt <- NULL
inflamm$prediction.score.early.enterocyte <- NULL
inflamm$prediction.score.Memory.B.cell <- NULL
inflamm$prediction.score.FCER2.B.cell <- NULL
inflamm$prediction.score.Tfh <- NULL
inflamm$prediction.score.cDC2 <- NULL
inflamm$prediction.score.Arterial.endothelial.cell <- NULL
inflamm$prediction.score.IL2RG..enterocyte..M.cell. <- NULL
inflamm$prediction.score.enterocyte <- NULL
inflamm$prediction.score.enteroendocrine <- NULL
inflamm$prediction.score.Venous.endothelial.cell <- NULL
inflamm$prediction.score.Activated.T <- NULL
inflamm$prediction.score.Monocyte <- NULL
inflamm$prediction.score.S1.fibroblasts <- NULL
inflamm$prediction.score.cDC1 <- NULL
inflamm$prediction.score.myofibroblast <- NULL
inflamm$prediction.score.pericyte <- NULL
inflamm$prediction.score.Treg <- NULL
inflamm$prediction.score.Lymphatic.endothelial.cell <- NULL
inflamm$prediction.score.B.cell <- NULL
inflamm$prediction.score.Activated.B.cell <- NULL
inflamm$prediction.score.pDC <- NULL
inflamm$prediction.score.Cycling.plasma.cell <- NULL
inflamm$prediction.score.gd.T.NK.cell <- NULL
inflamm$prediction.score.mast.cells <- NULL
inflamm$prediction.score.Paneth.cell <- NULL
inflamm$prediction.score.BEST4.enterocyte <- NULL
inflamm$prediction.score.activated.DC <- NULL
inflamm$prediction.score.S2.fibroblasts <- NULL
inflamm$prediction.score.CD8.T.cell <- NULL
inflamm$prediction.score.Macrophage <- NULL
inflamm$prediction.score.Glial.cell <- NULL
inflamm$prediction.score.Cycling.myeloid.cells <- NULL
inflamm$prediction.score.Tuft <- NULL
inflamm$prediction.score.max <- NULL
inflamm$prediction.match <- NULL
inflamm$health_meta <- NULL
inflamm$sample_meta <- NULL
inflamm$dataset <- NULL

head(inflamm[[]])
Idents(inflamm) <- "predicted.id"
levels(inflamm)

#save RDS
saveRDS(inflamm, file = "CD_martin_et_al.RDS")

# Add metadata column corresponding to allow for downstream KLF2 analyses
poscells.0.1 <- WhichCells(inflamm, expression = KLF2 >= 0.1)
poscells.1 <- WhichCells(inflamm, expression = KLF2 >= 1)
inflamm$KLF2.exp.0.1 <- ifelse(colnames(inflamm) %in% poscells.0.1, "KLF2.Pos.0.1", "KLF2.Neg.0.1")
inflamm$KLF2.exp.1 <- ifelse(colnames(inflamm) %in% poscells.1, "KLF2.Pos.1", "KLF2.Neg.1")
head(inflamm[[]])

# Subset CD4 T cells ONLY at first
Idents(inflamm) <-"predicted.id"
inflamm.CD4Tcells <- subset(inflamm, idents = "CD4 T cell")
inflamm.CD4Tcells #19229 cells#

#for Prism#
prop.table(table(inflamm$Patient, inflamm$KLF2.exp.1), margin = 1)

#So now we have all the right pieces for KLF2 analysis in the human CD datasets 
#(Elmentaite et al., Jaeger et al., Martin et al.)

#We went to also be able to look at FOXP3 and show whether differences between CD 
#and control is related or not to differences in FOXP3 too
#Same premise, just with FOXP3 instead of KLF2

#Note: this is first instance, but CD_colonna and colonna were both used at different points
#Due to re-loading the object over time, but they both obviously refer to the same dataset 
#from the Colonna group 

#Look at 3 CD4 object#
head(CD[[]])
head(colonna[[]])
head(inflamm[[]])
Idents(CD) <- "final.health"
Idents(colonna) <- "health"
Idents(inflamm) <- "health"

#Assess FOXP3 expression in sets
VlnPlot(CD, features = "FOXP3", split.by = "final.health") 
VlnPlot(colonna, features = "FOXP3", split.by = "health") 
VlnPlot(inflamm, features = "FOXP3", split.by = "health")

## Add metadata column corresponding to different expression levels of FOXP3
poscells.foxp3.0.1 <- WhichCells(CD, expression = FOXP3 >= 0.1)
poscells.foxp3.1 <- WhichCells(CD, expression = FOXP3 >= 1)
CD$FOXP3.exp.0.1 <- ifelse(colnames(CD) %in% poscells.foxp3.0.1, "FOXP3.Pos.0.1", "FOXP3.Neg.0.1")
CD$FOXP3.exp.1 <- ifelse(colnames(CD) %in% poscells.foxp3.1, "FOXP3.Pos.1", "FOXP3.Neg.1")
head(CD[[]])

poscells.foxp3.0.1 <- WhichCells(colonna, expression = FOXP3 >= 0.1)
poscells.foxp3.1 <- WhichCells(colonna, expression = FOXP3 >= 1)
colonna$FOXP3.exp.0.1 <- ifelse(colnames(colonna) %in% poscells.foxp3.0.1, "FOXP3.Pos.0.1", "FOXP3.Neg.0.1")
colonna$FOXP3.exp.1 <- ifelse(colnames(colonna) %in% poscells.foxp3.1, "FOXP3.Pos.1", "FOXP3.Neg.1")
head(colonna[[]])

poscells.foxp3.0.1 <- WhichCells(inflamm, expression = FOXP3 >= 0.1)
poscells.foxp3.1 <- WhichCells(inflamm, expression = FOXP3 >= 1)
inflamm$FOXP3.exp.0.1 <- ifelse(colnames(inflamm) %in% poscells.foxp3.0.1, "FOXP3.Pos.0.1", "FOXP3.Neg.0.1")
inflamm$FOXP3.exp.1 <- ifelse(colnames(inflamm) %in% poscells.foxp3.1, "FOXP3.Pos.1", "FOXP3.Neg.1")
head(inflamm[[]])

#Assess freq for prism
prop.table(table(CD$Sample.name, CD$FOXP3.exp.1), margin = 1)
prop.table(table(colonna$Patient, colonna$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm$Patient, inflamm$FOXP3.exp.1), margin = 1)

##Let's do pairwise comparisons of all CD4 T cell groups to complete up the analysis##
#This is where we start generating the data that is graphed in the final analyses
#Groups
#CD4, Treg, and Tfh
#Treg and Tfh
#Treg alone
#CD4 and Tfh

# Subset all CD4s for each group.
CD
head(CD[[]])
Idents(CD) <-"cell.id"
levels(CD)
CD.all.CD4 <- subset(CD, idents = c("CD4 T cell", "Treg", "Tfh"))
CD.all.CD4 #1627 cells#

#Note here: these feature plots assess the overall sense of 
#Both the label transfer and the subsets
FeaturePlot(CD.all.CD4, features = "CD4")
FeaturePlot(CD.all.CD4, features = "LCK")
FeaturePlot(CD.all.CD4, features = "CD3E")
FeaturePlot(CD.all.CD4, features = "PTPRC")
CD.all.CD4 <-  subset(CD.all.CD4, CD4 > 1) #Could do this, but is unnecessary here and in all other sets below, but keep for completeness; NOT DONE IN FINAL ANALYSES
CD.all.CD4 #1627 cells#
head(CD.all.CD4[[]])
#I won't comment anymore on it, but all subsets used in these analyses generally fit 
#what would be expected; not ever cell will always be expressing x gene at time of collection
#However, I trust the label transfer after looking at all these


colonna
head(colonna[[]])
Idents(colonna) <-"predicted.id"
levels(colonna)
colonna.all.CD4 <- subset(colonna, idents = c("CD4 T cell", "Treg", "Tfh"))
colonna.all.CD4 #15284 cells#

FeaturePlot(colonna.all.CD4, features = "CD4")
FeaturePlot(colonna.all.CD4, features = "LCK")
FeaturePlot(colonna.all.CD4, features = "CD3E")
FeaturePlot(colonna.all.CD4, features = "PTPRC")
colonna.all.CD4 #15284 cells#
head(colonna.all.CD4[[]])

inflamm
head(inflamm[[]])
Idents(inflamm) <-"predicted.id"
levels(inflamm)
inflamm.all.CD4 <- subset(inflamm, idents = c("CD4 T cell", "Treg", "Tfh"))
inflamm.all.CD4 #20440 cells#

FeaturePlot(inflamm.all.CD4, features = "CD4")
FeaturePlot(inflamm.all.CD4, features = "LCK")
FeaturePlot(inflamm.all.CD4, features = "CD3E")
FeaturePlot(inflamm.all.CD4, features = "PTPRC")
inflamm.all.CD4 #20440 cells#
head(inflamm.all.CD4[[]])

# Subset Treg and Tfh for each group.
CD
head(CD[[]])
Idents(CD) <-"cell.id"
levels(CD)
CD.treg.tfh <- subset(CD, idents = c("Treg", "Tfh"))
CD.treg.tfh #369 cells#

FeaturePlot(CD.treg.tfh, features = "CD4")
FeaturePlot(CD.treg.tfh, features = "LCK")
FeaturePlot(CD.treg.tfh, features = "CD3E")
FeaturePlot(CD.treg.tfh, features = "PTPRC")
CD.treg.tfh #369 cells#
head(CD.treg.tfh[[]])


colonna
head(colonna[[]])
Idents(colonna) <-"predicted.id"
levels(colonna)
colonna.treg.tfh <- subset(colonna, idents = c("Treg", "Tfh"))
colonna.treg.tfh #6125 cells#

FeaturePlot(colonna.treg.tfh, features = "CD4")
FeaturePlot(colonna.treg.tfh, features = "LCK")
FeaturePlot(colonna.treg.tfh, features = "CD3E")
FeaturePlot(colonna.treg.tfh, features = "PTPRC")
colonna.treg.tfh #6125 cells#
head(colonna.treg.tfh[[]])

inflamm
head(inflamm[[]])
Idents(inflamm) <-"predicted.id"
levels(inflamm)
inflamm.treg.tfh <- subset(inflamm, idents = c("Treg", "Tfh"))
inflamm.treg.tfh #1211 cells#

FeaturePlot(inflamm.treg.tfh, features = "CD4")
FeaturePlot(inflamm.treg.tfh, features = "LCK")
FeaturePlot(inflamm.treg.tfh, features = "CD3E")
FeaturePlot(inflamm.treg.tfh, features = "PTPRC")
inflamm.treg.tfh #1211 cells#
head(inflamm.treg.tfh[[]])

# Subset Treg for each group.
CD
head(CD[[]])
Idents(CD) <-"cell.id"
levels(CD)
CD.treg <- subset(CD, idents = c("Treg"))
CD.treg #137 cells#

FeaturePlot(CD.treg, features = "CD4")
FeaturePlot(CD.treg, features = "LCK")
FeaturePlot(CD.treg, features = "CD3E")
FeaturePlot(CD.treg, features = "PTPRC")
FeaturePlot(CD.treg, features = "FOXP3") #is much more here; good
CD.treg #137#
head(CD.treg[[]])


colonna
head(colonna[[]])
Idents(colonna) <-"predicted.id"
levels(colonna)
colonna.treg <- subset(colonna, idents = c("Treg"))
colonna.treg #1178 cells#

FeaturePlot(colonna.treg, features = "CD4")
FeaturePlot(colonna.treg, features = "LCK")
FeaturePlot(colonna.treg, features = "CD3E")
FeaturePlot(colonna.treg, features = "PTPRC")
FeaturePlot(colonna.treg, features = "FOXP3")
colonna.treg #1178 cells#
head(colonna.treg[[]])

inflamm
head(inflamm[[]])
Idents(inflamm) <-"predicted.id"
levels(inflamm)
inflamm.treg <- subset(inflamm, idents = c("Treg"))
inflamm.treg #747 cells#

FeaturePlot(inflamm.treg, features = "CD4")
FeaturePlot(inflamm.treg, features = "LCK")
FeaturePlot(inflamm.treg, features = "CD3E")
FeaturePlot(inflamm.treg, features = "PTPRC")
FeaturePlot(inflamm.treg, features = "FOXP3")
inflamm.treg #747 cells#
head(inflamm.treg[[]])

# Subset CD4 and Tfh for each group.
CD
head(CD[[]])
Idents(CD) <-"cell.id"
levels(CD)
CD.CD4.tfh <- subset(CD, idents = c("CD4 T cell", "Tfh"))
CD.CD4.tfh #1490 cells#

FeaturePlot(CD.CD4.tfh, features = "CD4")
FeaturePlot(CD.CD4.tfh, features = "LCK")
FeaturePlot(CD.CD4.tfh, features = "CD3E")
FeaturePlot(CD.CD4.tfh, features = "PTPRC")
CD.CD4.tfh #1490 cells#
head(CD.CD4.tfh[[]])


colonna
head(colonna[[]])
Idents(colonna) <-"predicted.id"
levels(colonna)
colonna.CD4.tfh <- subset(colonna, idents = c("CD4 T cell", "Tfh"))
colonna.CD4.tfh #14106 cells#

FeaturePlot(colonna.CD4.tfh, features = "CD4")
FeaturePlot(colonna.CD4.tfh, features = "LCK")
FeaturePlot(colonna.CD4.tfh, features = "CD3E")
FeaturePlot(colonna.CD4.tfh, features = "PTPRC")
colonna.CD4.tfh #14106 cells#
head(colonna.CD4.tfh[[]])

inflamm
head(inflamm[[]])
Idents(inflamm) <-"predicted.id"
levels(inflamm)
inflamm.CD4.tfh <- subset(inflamm, idents = c("CD4 T cell", "Tfh"))
inflamm.CD4.tfh #19693 cells#

FeaturePlot(inflamm.CD4.tfh, features = "CD4")
FeaturePlot(inflamm.CD4.tfh, features = "LCK")
FeaturePlot(inflamm.CD4.tfh, features = "CD3E")
FeaturePlot(inflamm.CD4.tfh, features = "PTPRC")
inflamm.CD4.tfh #19693 cells#
head(inflamm.CD4.tfh[[]])

#Subset activated T too for each group, just for completeness and to initially check
#We don't end up using this as it is not clear what these T cells actually include
#Regardless, there are not many anyway

CD
head(CD[[]])
Idents(CD) <-"cell.id"
levels(CD)
CD.activated <- subset(CD, idents = c("CD4 T cell", "Tfh", "Treg", "Activated T"))
CD.activated #1928 cells#

FeaturePlot(CD.activated, features = "CD4")
FeaturePlot(CD.activated, features = "LCK")
FeaturePlot(CD.activated, features = "CD3E")
FeaturePlot(CD.activated, features = "PTPRC")
CD.activated #1928 cells#
head(CD.activated[[]])


colonna
head(colonna[[]])
Idents(colonna) <-"predicted.id"
levels(colonna)
colonna.activated <- subset(colonna, idents = c("CD4 T cell", "Tfh", "Treg", "Activated T"))
colonna.activated #20978 cells#

FeaturePlot(colonna.activated, features = "CD4")
FeaturePlot(colonna.activated, features = "LCK")
FeaturePlot(colonna.activated, features = "CD3E")
FeaturePlot(colonna.activated, features = "PTPRC")
colonna.activated #20978 cells#
head(colonna.activated[[]])

inflamm
head(inflamm[[]])
Idents(inflamm) <-"predicted.id"
levels(inflamm)
inflamm.activated <- subset(inflamm, idents = c("CD4 T cell", "Tfh", "Treg", "Activated T"))
inflamm.activated #23141 cells#

FeaturePlot(inflamm.activated, features = "CD4")
FeaturePlot(inflamm.activated, features = "LCK")
FeaturePlot(inflamm.activated, features = "CD3E")
FeaturePlot(inflamm.activated, features = "PTPRC")
inflamm.activated #23141 cells#
head(inflamm.activated[[]])

#Start KLF2 and FoxP3 freq analyses in these populations
#CD4, Treg, and Tfh
#Treg and Tfh
#Treg alone
#CD4 and Tfh#

head(CD.all.CD4[[]])
head(colonna.all.CD4[[]])
head(inflamm.all.CD4[[]])

#Note: if you have already added the metadata columns for KLF2 and FOXP3 before subsetting, 
#they will carry to subsets. If not, you'll have to re-add them to each subset

#lets do freq tables for 4 cell groupings across 3 datasets
#lets also do number tables so we can now what numbers we're working with
#The raw numbers can be found in Supplemental data

#We end up using KLF2 and FOXP3 > 1 for all analyses for consistency and stringency 

#KLF2
#Assess all CD4 (CD4, Treg, and Tfh)
table(CD.all.CD4$Sample.name, CD.all.CD4$KLF2.exp.1)
table(colonna.all.CD4$Patient, colonna.all.CD4$KLF2.exp.1)
table(inflamm.all.CD4$Patient, inflamm.all.CD4$KLF2.exp.1)

prop.table(table(CD.all.CD4$Sample.name, CD.all.CD4$KLF2.exp.1), margin = 1)
prop.table(table(colonna.all.CD4$Patient, colonna.all.CD4$KLF2.exp.1), margin = 1)
prop.table(table(inflamm.all.CD4$Patient, inflamm.all.CD4$KLF2.exp.1), margin = 1)

#Assess other CD4 (Treg and Tfh)
table(CD.treg.tfh$Sample.name, CD.treg.tfh$KLF2.exp.1)
table(colonna.treg.tfh$Patient, colonna.treg.tfh$KLF2.exp.1)
table(inflamm.treg.tfh$Patient, inflamm.treg.tfh$KLF2.exp.1)

prop.table(table(CD.treg.tfh$Sample.name, CD.treg.tfh$KLF2.exp.1), margin = 1)
prop.table(table(colonna.treg.tfh$Patient, colonna.treg.tfh$KLF2.exp.1), margin = 1)
prop.table(table(inflamm.treg.tfh$Patient, inflamm.treg.tfh$KLF2.exp.1), margin = 1)

#Assess Treg alone (Treg)
table(CD.treg$Sample.name, CD.treg$KLF2.exp.1)
table(colonna.treg$Patient, colonna.treg$KLF2.exp.1)
table(inflamm.treg$Patient, inflamm.treg$KLF2.exp.1)

prop.table(table(CD.treg$Sample.name, CD.treg$KLF2.exp.1), margin = 1)
prop.table(table(colonna.treg$Patient, colonna.treg$KLF2.exp.1), margin = 1)
prop.table(table(inflamm.treg$Patient, inflamm.treg$KLF2.exp.1), margin = 1)

#Assess CD4 and Tfh together (CD4 and Tfh)
table(CD.CD4.tfh$Sample.name, CD.CD4.tfh$KLF2.exp.1)
table(colonna.CD4.tfh$Patient, colonna.CD4.tfh$KLF2.exp.1)
table(inflamm.CD4.tfh$Patient, inflamm.CD4.tfh$KLF2.exp.1)

prop.table(table(CD.CD4.tfh$Sample.name, CD.CD4.tfh$KLF2.exp.1), margin = 1)
prop.table(table(colonna.CD4.tfh$Patient, colonna.CD4.tfh$KLF2.exp.1), margin = 1)
prop.table(table(inflamm.CD4.tfh$Patient, inflamm.CD4.tfh$KLF2.exp.1), margin = 1)

#Assess ALL CD4 T cells, including activated (CD4, activated T, Tfh, and Treg)
table(CD.activated$Sample.name, CD.activated$KLF2.exp.1)
table(colonna.activated$Patient, colonna.activated$KLF2.exp.1)
table(inflamm.activated$Patient, inflamm.activated$KLF2.exp.1)

prop.table(table(CD.activated$Sample.name, CD.activated$KLF2.exp.1), margin = 1)
prop.table(table(colonna.activated$Patient, colonna.activated$KLF2.exp.1), margin = 1)
prop.table(table(inflamm.activated$Patient, inflamm.activated$KLF2.exp.1), margin = 1)

#FOXP3
#Assess all CD4 (CD4, Treg, and Tfh)
table(CD.all.CD4$Sample.name, CD.all.CD4$FOXP3.exp.1)
table(colonna.all.CD4$Patient, colonna.all.CD4$FOXP3.exp.1)
table(inflamm.all.CD4$Patient, inflamm.all.CD4$FOXP3.exp.1)

prop.table(table(CD.all.CD4$Sample.name, CD.all.CD4$FOXP3.exp.1), margin = 1)
prop.table(table(colonna.all.CD4$Patient, colonna.all.CD4$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm.all.CD4$Patient, inflamm.all.CD4$FOXP3.exp.1), margin = 1)

#Assess other CD4 (Treg and Tfh)
table(CD.treg.tfh$Sample.name, CD.treg.tfh$FOXP3.exp.1)
table(colonna.treg.tfh$Patient, colonna.treg.tfh$FOXP3.exp.1)
table(inflamm.treg.tfh$Patient, inflamm.treg.tfh$FOXP3.exp.1)

prop.table(table(CD.treg.tfh$Sample.name, CD.treg.tfh$FOXP3.exp.1), margin = 1)
prop.table(table(colonna.treg.tfh$Patient, colonna.treg.tfh$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm.treg.tfh$Patient, inflamm.treg.tfh$FOXP3.exp.1), margin = 1)

#Assess Treg alone (Treg)
table(CD.treg$Sample.name, CD.treg$FOXP3.exp.1)
table(colonna.treg$Patient, colonna.treg$FOXP3.exp.1)
table(inflamm.treg$Patient, inflamm.treg$FOXP3.exp.1)

prop.table(table(CD.treg$Sample.name, CD.treg$FOXP3.exp.1), margin = 1)
prop.table(table(colonna.treg$Patient, colonna.treg$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm.treg$Patient, inflamm.treg$FOXP3.exp.1), margin = 1)

#Assess CD4 and Tfh together (CD4 and Tfh)
table(CD.CD4.tfh$Sample.name, CD.CD4.tfh$FOXP3.exp.1)
table(colonna.CD4.tfh$Patient, colonna.CD4.tfh$FOXP3.exp.1)
table(inflamm.CD4.tfh$Patient, inflamm.CD4.tfh$FOXP3.exp.1)

prop.table(table(CD.CD4.tfh$Sample.name, CD.CD4.tfh$FOXP3.exp.1), margin = 1)
prop.table(table(colonna.CD4.tfh$Patient, colonna.CD4.tfh$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm.CD4.tfh$Patient, inflamm.CD4.tfh$FOXP3.exp.1), margin = 1)

#Assess ALL CD4 T cells, including activated (CD4, Activated T cell, Tfh, and Treg)
table(CD.activated$Sample.name, CD.activated$FOXP3.exp.1)
table(colonna.activated$Patient, colonna.activated$FOXP3.exp.1)
table(inflamm.activated$Patient, inflamm.activated$FOXP3.exp.1)

prop.table(table(CD.activated$Sample.name, CD.activated$FOXP3.exp.1), margin = 1)
prop.table(table(colonna.activated$Patient, colonna.activated$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm.activated$Patient, inflamm.activated$FOXP3.exp.1), margin = 1)


##Do these same analyses now for UC dataset##
#This will lend consistency to the final analyses 

Idents(UC) <- "Cluster"
levels(UC)

#Transfer labels from Teichmann for consistency 
##Cell type classification using an integrated reference##
#Find transfer anchors that will use Teichmann's CD as a reference for cell type and the Smillie et al. set as a query#
CD.anchors <- FindTransferAnchors(reference = CD, query = UC.final, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = CD.anchors, refdata = CD$annotation_V2, dims = 1:30)

#Add the predicted cell types to the query, i.e. the IBD integrated object, metadata
UC.final <- AddMetaData(UC.final, metadata = predictions)


#Clean off some of the irrelevant metadata columns
UC.final$prediction.score.CD4.T.cell <- NULL
UC.final$prediction.score.TA <- NULL
UC.final$prediction.score.IgA.plasma.cell <- NULL
UC.final$prediction.score.Goblet.cell <- NULL
UC.final$prediction.score.IgG.plasma.cell <- NULL
UC.final$prediction.score.Cycling.B.cell <- NULL
UC.final$prediction.score.S4.fibroblasts <- NULL
UC.final$prediction.score.crypt <- NULL
UC.final$prediction.score.early.enterocyte <- NULL
UC.final$prediction.score.Memory.B.cell <- NULL
UC.final$prediction.score.FCER2.B.cell <- NULL
UC.final$prediction.score.Tfh <- NULL
UC.final$prediction.score.cDC2 <- NULL
UC.final$prediction.score.Arterial.endothelial.cell <- NULL
UC.final$prediction.score.IL2RG..enterocyte..M.cell. <- NULL
UC.final$prediction.score.enterocyte <- NULL
UC.final$prediction.score.enteroendocrine <- NULL
UC.final$prediction.score.Venous.endothelial.cell <- NULL
UC.final$prediction.score.Activated.T <- NULL
UC.final$prediction.score.Monocyte <- NULL
UC.final$prediction.score.S1.fibroblasts <- NULL
UC.final$prediction.score.cDC1 <- NULL
UC.final$prediction.score.myofibroblast <- NULL
UC.final$prediction.score.pericyte <- NULL
UC.final$prediction.score.Treg <- NULL
UC.final$prediction.score.Lymphatic.endothelial.cell <- NULL
UC.final$prediction.score.B.cell <- NULL
UC.final$prediction.score.Activated.B.cell <- NULL
UC.final$prediction.score.pDC <- NULL
UC.final$prediction.score.Cycling.plasma.cell <- NULL
UC.final$prediction.score.gd.T.NK.cell <- NULL
UC.final$prediction.score.mast.cells <- NULL
UC.final$prediction.score.Paneth.cell <- NULL
UC.final$prediction.score.BEST4.enterocyte <- NULL
UC.final$prediction.score.activated.DC <- NULL
UC.final$prediction.score.S2.fibroblasts <- NULL
UC.final$prediction.score.CD8.T.cell <- NULL
UC.final$prediction.score.Macrophage <- NULL
UC.final$prediction.score.Glial.cell <- NULL
UC.final$prediction.score.Cycling.myeloid.cells <- NULL
UC.final$prediction.score.Tuft <- NULL
UC.final$prediction.score.max <- NULL
UC.final$prediction.match <- NULL
UC.final$health_meta <- NULL
UC.final$sample_meta <- NULL
UC.final$dataset <- NULL

head(UC.final[[]])
Idents(UC.final) <- "predicted.id"
levels(UC.final)

#save RDS
saveRDS(UC.final, file = "UC_final_predicted_ids.RDS")

#Add KLF2 and FOXP3 metadata for UC
#KLF2
## Add metadata column corresponding to different expression levels of KLF2
poscells.0.1 <- WhichCells(UC.final, expression = KLF2 >= 0.1)
poscells.1 <- WhichCells(UC.final, expression = KLF2 >= 1)
UC.final$KLF2.exp.0.1 <- ifelse(colnames(UC.final) %in% poscells.0.1, "KLF2.Pos.0.1", "KLF2.Neg.0.1")
UC.final$KLF2.exp.1 <- ifelse(colnames(UC.final) %in% poscells.1, "KLF2.Pos.1", "KLF2.Neg.1")
head(UC.final[[]])

#FOXP3
## Add metadata column corresponding to different expression levels of FOXP3
poscells.foxp3.0.1 <- WhichCells(UC.final, expression = FOXP3 >= 0.1)
poscells.foxp3.1 <- WhichCells(UC.final, expression = FOXP3 >= 1)
UC.final$FOXP3.exp.0.1 <- ifelse(colnames(UC.final) %in% poscells.foxp3.0.1, "FOXP3.Pos.0.1", "FOXP3.Neg.0.1")
UC.final$FOXP3.exp.1 <- ifelse(colnames(UC.final) %in% poscells.foxp3.1, "FOXP3.Pos.1", "FOXP3.Neg.1")
head(UC.final[[]])

#Create different subsets- Should be 5 groupings
#CD4, Tfh, and Treg
#Treg and Tfh
#Treg only
#CD4 and Tfh
#CD4, Tfh, Treg, and Activated T
#Oh also JUST the CD4s for completeness

# Subset CD4, Tfh, and Treg for each group
UC.final #143202 cells 
head(UC.final[[]])
Idents(UC.final) <-"predicted.id"
levels(UC.final)
UC.final.all.CD4 <- subset(UC.final, idents = c("CD4 T cell", "Treg", "Tfh"))
UC.final.all.CD4 #30915 cells#

FeaturePlot(UC.final.all.CD4, features = "CD4")
FeaturePlot(UC.final.all.CD4, features = "LCK")
FeaturePlot(UC.final.all.CD4, features = "CD3E")
FeaturePlot(UC.final.all.CD4, features = "PTPRC")
UC.final.all.CD4 #30915 cells#
head(UC.final.all.CD4[[]])


# Subset Treg and Tfh for each group.
UC.final
head(UC.final[[]])
Idents(UC.final) <-"predicted.id"
levels(UC.final)
UC.final.treg.tfh <- subset(UC.final, idents = c("Treg", "Tfh"))
UC.final.treg.tfh #5624 cells#

FeaturePlot(UC.final.treg.tfh, features = "CD4")
FeaturePlot(UC.final.treg.tfh, features = "LCK")
FeaturePlot(UC.final.treg.tfh, features = "CD3E")
FeaturePlot(UC.final.treg.tfh, features = "PTPRC")
UC.final.treg.tfh #5624 cells#
head(UC.final.treg.tfh[[]])

# Subset Treg for each group.
UC.final
head(UC.final[[]])
Idents(UC.final) <-"predicted.id"
levels(UC.final)
UC.final.treg <- subset(UC.final, idents = c("Treg"))
UC.final.treg #1558 cells#

FeaturePlot(UC.final.treg, features = "CD4")
FeaturePlot(UC.final.treg, features = "LCK")
FeaturePlot(UC.final.treg, features = "CD3E")
FeaturePlot(UC.final.treg, features = "PTPRC")
FeaturePlot(UC.final.treg, features = "FOXP3") #is much more here; good
UC.final.treg #1558#
head(UC.final.treg[[]])

# Subset CD4 and Tfh for each group.
UC.final
head(UC.final[[]])
Idents(UC.final) <-"predicted.id"
levels(UC.final)
UC.final.CD4.tfh <- subset(UC.final, idents = c("CD4 T cell", "Tfh"))
UC.final.CD4.tfh #29357 cells#

FeaturePlot(UC.final.CD4.tfh, features = "CD4")
FeaturePlot(UC.final.CD4.tfh, features = "LCK")
FeaturePlot(UC.final.CD4.tfh, features = "CD3E")
FeaturePlot(UC.final.CD4.tfh, features = "PTPRC")
UC.final.CD4.tfh #29357 cells#
head(UC.final.CD4.tfh[[]])


# Subset activated too for each group.
UC.final
head(UC.final[[]])
Idents(UC.final) <-"predicted.id"
levels(UC.final)
UC.final.activated <- subset(UC.final, idents = c("CD4 T cell", "Tfh", "Treg", "Activated T"))
UC.final.activated #39133 cells#

FeaturePlot(UC.final.activated, features = "CD4")
FeaturePlot(UC.final.activated, features = "LCK")
FeaturePlot(UC.final.activated, features = "CD3E")
FeaturePlot(UC.final.activated, features = "PTPRC")
UC.final.activated #39133 cells#
head(UC.final.activated[[]])


# Subset JUST CD4 for UC
UC.final
head(UC.final[[]])
Idents(UC.final) <-"predicted.id"
levels(UC.final)
UC.final.CD4 <- subset(UC.final, idents = c("CD4 T cell"))
UC.final.CD4 #25291 cells#

FeaturePlot(UC.final.CD4, features = "CD4")
FeaturePlot(UC.final.CD4, features = "LCK")
FeaturePlot(UC.final.CD4, features = "CD3E")
FeaturePlot(UC.final.CD4, features = "PTPRC")
UC.final.CD4 #25291 cells#
head(UC.final.CD4[[]])

##Begin prop analyses for prism for UC
head(UC.final[[]])
#KLF2
#Assess all CD4 (CD4, Treg, and Tfh)
table(UC.final.all.CD4$Subject, UC.final.all.CD4$KLF2.exp.1)
prop.table(table(UC.final.all.CD4$Subject, UC.final.all.CD4$KLF2.exp.1), margin = 1)

#Assess other CD4 (Treg and Tfh)
table(UC.final.treg.tfh$Subject, UC.final.treg.tfh$KLF2.exp.1)
prop.table(table(UC.final.treg.tfh$Subject, UC.final.treg.tfh$KLF2.exp.1), margin = 1)

#Assess Treg alone (Treg)
table(UC.final.treg$Subject, UC.final.treg$KLF2.exp.1)
prop.table(table(UC.final.treg$Subject, UC.final.treg$KLF2.exp.1), margin = 1)

#Assess CD4 and Tfh together (CD4 and Tfh)
table(UC.final.CD4.tfh$Subject, UC.final.CD4.tfh$KLF2.exp.1)
prop.table(table(UC.final.CD4.tfh$Subject, UC.final.CD4.tfh$KLF2.exp.1), margin = 1)

#Assess ALL CD4 T cells (CD4, activated T, Tfh, and Treg)
table(UC.final.activated$Subject, UC.final.activated$KLF2.exp.1)
prop.table(table(UC.final.activated$Subject, UC.final.activated$KLF2.exp.1), margin = 1)

#Assess CD4 T cells alone
table(UC.final.CD4$Subject, UC.final.CD4$KLF2.exp.1)
prop.table(table(UC.final.CD4$Subject, UC.final.CD4$KLF2.exp.1), margin = 1)

#FOXP3
#Assess all CD4 (CD4, Treg, and Tfh)
table(UC.final.all.CD4$Subject, UC.final.all.CD4$FOXP3.exp.1)
prop.table(table(UC.final.all.CD4$Subject, UC.final.all.CD4$FOXP3.exp.1), margin = 1)

#Assess other CD4 (Treg and Tfh)
table(UC.final.treg.tfh$Subject, UC.final.treg.tfh$FOXP3.exp.1)
prop.table(table(UC.final.treg.tfh$Subject, UC.final.treg.tfh$FOXP3.exp.1), margin = 1)

#Assess Treg alone (Treg)
table(UC.final.treg$Subject, UC.final.treg$FOXP3.exp.1)
prop.table(table(UC.final.treg$Subject, UC.final.treg$FOXP3.exp.1), margin = 1)

#Assess CD4 and Tfh together (CD4 and Tfh)
table(UC.final.CD4.tfh$Subject, UC.final.CD4.tfh$FOXP3.exp.1)
prop.table(table(UC.final.CD4.tfh$Subject, UC.final.CD4.tfh$FOXP3.exp.1), margin = 1)

#Assess ALL CD4 T cells (CD4, activated T, Tfh, and Treg)
table(UC.final.activated$Subject, UC.final.activated$FOXP3.exp.1)
prop.table(table(UC.final.activated$Subject, UC.final.activated$FOXP3.exp.1), margin = 1)

#Assess CD4 T cells alone
table(UC.final.CD4$Subject, UC.final.CD4$FOXP3.exp.1)
prop.table(table(UC.final.CD4$Subject, UC.final.CD4$FOXP3.exp.1), margin = 1)

#We then wanted to add an analysis just looking at Tfh because the mouse data
#originally indicated that the important KLF2+CD4+ cell was NOT
#within existing Tfh populations

#Do analysis on CD sets for Tfh only for KLF2 and FOXP3

#Subset out tfh
CD
head(CD[[]])
Idents(CD) <-"annotation_V2"
levels(CD)
CD.tfh <- subset(CD, idents = c("Tfh"))
CD.tfh #232 cells#

FeaturePlot(CD.tfh, features = "CD4")
FeaturePlot(CD.tfh, features = "LCK")
FeaturePlot(CD.tfh, features = "CD3E")
FeaturePlot(CD.tfh, features = "PTPRC")
CD.tfh #232 cells#


colonna
head(colonna[[]])
Idents(colonna) <-"predicted.id"
levels(colonna)
colonna.tfh <- subset(colonna, idents = c("Tfh"))
colonna.tfh #4947 cells#

FeaturePlot(colonna.tfh, features = "CD4")
FeaturePlot(colonna.tfh, features = "LCK")
FeaturePlot(colonna.tfh, features = "CD3E")
FeaturePlot(colonna.tfh, features = "PTPRC")
colonna.tfh #4947 cells#


inflamm
head(inflamm[[]])
Idents(inflamm) <-"predicted.id"
levels(inflamm)
inflamm.tfh <- subset(inflamm, idents = c("Tfh"))
inflamm.tfh #464 cells#

FeaturePlot(inflamm.tfh, features = "CD4")
FeaturePlot(inflamm.tfh, features = "LCK")
FeaturePlot(inflamm.tfh, features = "CD3E")
FeaturePlot(inflamm.tfh, features = "PTPRC")
inflamm.tfh #464 cells#

#prop tables for prism
#KLF2
prop.table(table(CD.tfh$Sample.name, CD.tfh$KLF2.exp.1), margin = 1)
prop.table(table(colonna.tfh$Patient, colonna.tfh$KLF2.exp.1), margin = 1)
prop.table(table(inflamm.tfh$Patient, inflamm.tfh$KLF2.exp.1), margin = 1)

#FOXP3
prop.table(table(CD.tfh$Sample.name, CD.tfh$FOXP3.exp.1), margin = 1)
prop.table(table(colonna.tfh$Patient, colonna.tfh$FOXP3.exp.1), margin = 1)
prop.table(table(inflamm.tfh$Patient, inflamm.tfh$FOXP3.exp.1), margin = 1)

##Run code to get raw numbers for 4 datasets (3 CD and 1 UC) of
#CD4, Tfh, and Treg, both total number and KLF2+##
#This was done to get the raw numbers of cells for the Supplemental data
#Was done after the fact, so this may not need to be run if already done, 
#but we retain this section for completeness

#Now run tables for numbers by patient
#CD4
#KLF2
table(CD.CD4$Sample.name, CD.CD4$KLF2.exp.1)
table(colonna.CD4$Patient, colonna.CD4$KLF2.exp.1)
table(inflamm.CD4$Patient, inflamm.CD4$KLF2.exp.1)
table(UC.final.CD4$Subject, UC.final.CD4$KLF2.exp.1)

#FOXP3
table(CD.CD4$Sample.name, CD.CD4$FOXP3.exp.1)
table(colonna.CD4$Patient, colonna.CD4$FOXP3.exp.1)
table(inflamm.CD4$Patient, inflamm.CD4$FOXP3.exp.1)
table(UC.final.CD4$Subject, UC.final.CD4$FOXP3.exp.1)

#Tfh
#KLF2
table(CD.tfh$Sample.name, CD.tfh$KLF2.exp.1)
table(colonna.tfh$Patient, colonna.tfh$KLF2.exp.1)
table(inflamm.tfh$Patient, inflamm.tfh$KLF2.exp.1)
table(UC.final.tfh$Subject, UC.final.tfh$KLF2.exp.1)

#FOXP3
table(CD.tfh$Sample.name, CD.tfh$FOXP3.exp.1)
table(colonna.tfh$Patient, colonna.tfh$FOXP3.exp.1)
table(inflamm.tfh$Patient, inflamm.tfh$FOXP3.exp.1)
table(UC.final.tfh$Subject, UC.final.tfh$FOXP3.exp.1)

#Treg
#KLF2
table(CD.treg$Sample.name, CD.treg$KLF2.exp.1)
table(colonna.treg$Patient, colonna.treg$KLF2.exp.1)
table(inflamm.treg$Patient, inflamm.treg$KLF2.exp.1)
table(UC.final.treg$Subject, UC.final.treg$KLF2.exp.1)

#FOXP3
table(CD.treg$Sample.name, CD.treg$FOXP3.exp.1)
table(colonna.treg$Patient, colonna.treg$FOXP3.exp.1)
table(inflamm.treg$Patient, inflamm.treg$FOXP3.exp.1)
table(UC.final.treg$Subject, UC.final.treg$FOXP3.exp.1)

#Finally, to assess the integration of these 3 CD datasets, we will
#Integrate them using the standard Seurat workflow and then
#visualize with UMAP embeddings seen in figures

#Note: this is NOT what was used for the final UMAP figures
#This was an initial process where we integrated the datasets into one
#For the final figures, we instead used UMAP projection as a more intuitive
#way of looking at the UMAP

###Create integrated object with all 3 CD datasets to make UMAPs###
head(CD[[]])
head(colonna[[]])
head(inflamm[[]])

#Create column in original CD set called predicted.id, so they will all mesh in 
#combined object

Idents(CD) <- "annotation_V2"
levels(CD)
cell.state <- CD$annotation_V2
CD <- AddMetaData(CD, metadata = cell.state, col.name = "predicted.id")
Idents(CD) <-"predicted.id"
levels(CD)
head(CD[[]])

#Add other metadata before integrating!
CD@meta.data[, "health_meta"] <- paste0("Teichmann- ", CD$Diagnosis )
CD@meta.data[, "sample_meta"] <- paste0("Teichmann- ", CD$Sample.name)
colonna@meta.data[, "health_meta"] <- paste0("Colonna- ", colonna$health)
colonna@meta.data[, "sample_meta"] <- paste0("Colonna- ", colonna$Patient)
inflamm@meta.data[, "health_meta"] <- paste0("Martin- ", inflamm$health )
inflamm@meta.data[, "sample_meta"] <- paste0("Martin- ", inflamm$Patient)

CD@meta.data[, "dataset"] <- "Teichmann"
colonna@meta.data[, "dataset"] <- "Colonna"
inflamm@meta.data[, "dataset"] <- "Martin"


#check
Idents(CD) <- "health_meta"
levels(CD)
Idents(CD) <- "sample_meta"
levels(CD)
Idents(colonna) <- "health_meta"
levels(colonna)
Idents(colonna) <- "sample_meta"
levels(colonna)
Idents(inflamm) <- "health_meta"
levels(inflamm)
Idents(inflamm) <- "sample_meta"
levels(inflamm)

Idents(CD) <- "dataset"
levels(CD)
Idents(colonna) <- "dataset"
levels(colonna)
Idents(inflamm) <- "dataset"
levels(inflamm)

#merge three datasets
combined <- merge(x = CD, y = c(colonna, inflamm))
combined 
head(combined[[]])

#Create list of objects#

ibd.list <- SplitObject(combined, split.by = "dataset")
ibd.list

#Integrate dataset#

#again, make sure this is set as this is a computationally expensive process
options(future.globals.maxSize= 2500000000000)

#normalize and identify variables for each dataset independently
ibd.list <- lapply(X = ibd.list, FUN = function(x) {
 x <- NormalizeData(x)
 x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
 })

#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ibd.list)

#Perform integration#
#Will use modified integration technique for large datasets to get FindIntegrationAnchors to work
#Again, this is as described in Seurat

#Run PCA, which is needed for running the alternative reciprocal PCA workflow
ibd.list <- lapply(X = ibd.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x<- RunPCA(x, features = features, verbose = FALSE)
 })

ibd.anchors <- FindIntegrationAnchors(object.list = ibd.list, reduction = "rpca", dims = 1:50, anchor.features = features)
ibd.integrated <- IntegrateData(anchorset = ibd.anchors, dims = 1:50)

#Run workflow for visualization and clustering#

#set assay to integrated#
ibd.integrated
DefaultAssay(ibd.integrated) <- "integrated"

#workflow#
ibd.integrated <- ScaleData(ibd.integrated, verbose = FALSE)
ibd.integrated <- RunPCA(ibd.integrated, npcs = 30, verbose = FALSE)
ibd.integrated <- RunUMAP(ibd.integrated, reduction = "pca", dims = 1:30)
ibd.integrated <- FindNeighbors(ibd.integrated, reduction = "pca", dims = 1:30)
ibd.integrated <- FindClusters(ibd.integrated, resolution = 0.5)

#SAVE RDS FOR INTEGRATED OBJECT
saveRDS(ibd.integrated, file = "integrated_CD_all3_sets.RDS")

#Run UMAP visualization
#Can do many ways, but these are what we visualized
Idents(ibd.integrated) <- "dataset"
DimPlot(ibd.integrated)

#Subset out T cell groups
#To see how they integrate
ibd.tcell <- subset(ibd.integrated, idents = c("CD4 T cell", "Tfh", "Treg"))
ibd.tcell
Idents(ibd.tcell) <- "predicted.id"
DimPlot(ibd.tcell, split.by = "dataset")
FeaturePlot(ibd.tcell, features = "FOXP3", split.by = "predicted.id", min.cutoff = "q10") + RestoreLegend()


#Making GO plots for KLF2 pos and KLF2 neg CD4 T cells 
#These object names will be different due to the times they were made 

imm.seur.Tcells
head(imm.seur.Tcells[[]])
Idents(imm.seur.Tcells) <-"cell.id"
levels(imm.seur.Tcells)
Idents(imm.seur.Tcells) <-"Health"

DE.all.T.cells <- FindMarkers(imm.seur.Tcells, ident.1 = "Healthy", ident.2 = "Inflamed", min.pct = 0.25, logfc.threshold = 0.5)
DE.all.T.cells <- DE.all.T.cells %>% rownames_to_column(var = "gene")
View(DE.all.T.cells)
write.table(DE.all.T.cells, file = "/volumetowriteto", sep = ",", quote = FALSE, row.names = T)

# DEG for CD4 KLF2 pos cells
Idents(imm.seur.Tcells) <-"id.exp.1"
levels(imm.seur.Tcells)
KLF2.pos.CD4.Tcells <- subset(imm.seur.Tcells, idents = "CD4 KLF2.Pos.1")
Idents(KLF2.pos.CD4.Tcells) <-"Health"
levels(KLF2.pos.CD4.Tcells)
DE.KLF2.pos.CD4.Tcells <- FindMarkers(KLF2.pos.CD4.Tcells, ident.1 = "Healthy", ident.2 = "Inflamed", min.pct = 0.25, logfc.threshold = 0.5)
DE.KLF2.pos.CD4.Tcells <- DE.KLF2.pos.CD4.Tcells %>% rownames_to_column(var = "gene")
View(DE.KLF2.pos.CD4.Tcells)
write.table(DE.KLF2.pos.CD4.Tcells, file = "/volumetowriteto", sep = ",", quote = FALSE, row.names = T)

# DEG for CD4 KLF2 neg cells
Idents(imm.seur.Tcells) <-"id.exp.1"
levels(imm.seur.Tcells)
KLF2.neg.CD4.Tcells <- subset(imm.seur.Tcells, idents = "CD4 KLF2.Neg.1")
Idents(KLF2.neg.CD4.Tcells) <-"Health"
levels(KLF2.neg.CD4.Tcells)
DE.KLF2.neg.CD4.Tcells <- FindMarkers(KLF2.neg.CD4.Tcells, ident.1 = "Healthy", ident.2 = "Inflamed", min.pct = 0.25, logfc.threshold = 0.5)
DE.KLF2.neg.CD4.Tcells <- DE.KLF2.neg.CD4.Tcells %>% rownames_to_column(var = "gene")
View(DE.KLF2.neg.CD4.Tcells)
write.table(DE.KLF2.neg.CD4.Tcells, file = "/volumetowriteto", sep = ",", quote = FALSE, row.names = T)

# Make GO Dotplot 
GO <- read_delim("/volumetoreadfrom/GO.import.csv")
GO 
library(forcats)
ggplot(data = GO, aes(x = Type, y = fct_reorder(GO, FE), color = FDR, size = FE)) + geom_point() + scale_color_gradient(limits=c(0, 0.10), low="red") + theme_classic() + facet_grid(.~ID)


GO <- read_delim("volumetoreadfrom/GO.import.UC.CD.csv")
GO 
library(forcats)
ggplot(data = GO, aes(x = Type, y = fct_reorder(GO, FE), color = FDR, size = FE)) + geom_point() + scale_color_gradient(limits=c(0, 0.10), low="red") + theme_classic()
ggplot(data = GO, aes(x = Type, y = fct_reorder(GO, Type), color = FDR, size = FE)) + geom_point() + scale_color_gradient(limits=c(0, 0.10), low="red") + theme_classic()


##Code for plotting UMAPs as seen in Extended Data 8 and 9##
#This code is used after creating the Reference Plotted UMAPs of the data
#The input for this section of code can be found in GitHub subfolder called "Selected data"
#This code created/compiled by Tony Jiang

#Load necessary libraries if needed 
library(dplyr)
library(Seurat)
library(patchwork)
library(msigdbr)
library(presto)
library(devtools)
library(tibble)
library(ggplot2)
library(SingleR)
library(celldex)
library(gprofiler2)
library(SingleR)
library(cowplot)
library(scales)
library(ggrepel)

setwd('/pathtowd')

textSize = 8
lineSize = 0.5
dotSize = 0.25

# Load the dataset
data = read.csv("pathtodata/UMAP_metadata_Extended_data_8_9.csv", header = T)

# Setup groupings
data$isKLF2 = "KLF2-"
data$isKLF2[data$KLF2.exp.1 %in% "KLF2.Pos"] = "KLF2+"

data$isFOXP3 = "FOXP3-"
data$isFOXP3[data$FOXP3.poscells.1 %in% "FOXP3.Pos"] = "FOXP3+"

# split into UC and CD dataframes
cd.df = subset(data, dataset != "UC")
hc.idx = c(grep("HC", cd.df$health_meta), grep("control", cd.df$health_meta))
cd.df$ptID = "CD"
cd.df$ptID[hc.idx] = "HC"

uc.df = subset(data, dataset %in% "UC")

# plot KLF2 expression among Crohns set for subsets CD4 T/Tfh/Treg
cd.hc.df = subset(cd.df, ptID %in% "HC")
a = ggplot(cd.hc.df) +
  geom_point(data = subset(cd.hc.df, isKLF2 %in% "KLF2-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  geom_point(data = subset(cd.hc.df, isKLF2 %in% "KLF2+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("KLF2-", "KLF2+"),
                     values = c("gray","orange")) +
  labs(color = "Healthy controls") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-13, 1), breaks = seq(-12,0,4))

cd.cd.df = subset(cd.df, ptID %in% "CD")
b = ggplot(cd.cd.df) +
  geom_point(data = subset(cd.cd.df, isKLF2 %in% "KLF2-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  geom_point(data = subset(cd.cd.df, isKLF2 %in% "KLF2+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("KLF2-", "KLF2+"),
                     values = c("gray","orange")) +
  labs(color = "Crohn's disease") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-13, 1), breaks = seq(-12,0,4))

wrap_plots(a, b, nrow = 2)
ggsave("CD_KLF2.png", width = 18, height = 12, unit = "cm", dpi = 600, bg="transparent")

#plot FOXP3 expression among Crohns set for subsets CD4 T/Tfh/Treg
c = ggplot(cd.hc.df) +
  geom_point(data = subset(cd.hc.df, isFOXP3 %in% "FOXP3-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  geom_point(data = subset(cd.hc.df, isFOXP3 %in% "FOXP3+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("FOXP3-", "FOXP3+"),
                     values = c("gray","darkgreen")) +
  labs(color = "Healthy controls") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-13, 1), breaks = seq(-12,0,4))

d = ggplot(cd.cd.df) +
  geom_point(data = subset(cd.cd.df, isFOXP3 %in% "FOXP3-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  geom_point(data = subset(cd.cd.df, isFOXP3 %in% "FOXP3+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("FOXP3-", "FOXP3+"),
                     values = c("gray","darkgreen")) +
  labs(color = "Crohn's disease") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-13, 1), breaks = seq(-12,0,4))

wrap_plots(c, d, nrow = 2)
ggsave("CD_foxp3.png", width = 18, height = 12, unit = "cm", dpi = 600, bg="transparent")

# plot KLF2 expression in UC set for subsets CD4 T/Tfh/Treg
uc.hc.df = subset(uc.df, health_meta %in% "UC-Healthy")
a2 = ggplot(uc.hc.df) +
  geom_point(data = subset(uc.hc.df, isKLF2 %in% "KLF2-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  geom_point(data = subset(uc.hc.df, isKLF2 %in% "KLF2+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("KLF2-", "KLF2+"),
                     values = c("gray","orange")) +
  labs(color = "Healthy controls") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-11, 3), breaks = seq(-10,0,5))

uc.inf.df = subset(uc.df, health_meta %in% "UC-Inflamed")
b2 = ggplot(uc.inf.df) +
  geom_point(data = subset(uc.inf.df, isKLF2 %in% "KLF2-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  geom_point(data = subset(uc.inf.df, isKLF2 %in% "KLF2+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isKLF2), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("KLF2-", "KLF2+"),
                     values = c("gray","orange")) +
  labs(color = "Ulcerative colitis") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-11, 3), breaks = seq(-10,0,5))

wrap_plots(a2, b2, nrow = 2)
ggsave("UC_KLF2.png", width = 18, height = 14, unit = "cm", dpi = 600, bg="transparent")

#plot FOXP3 expression among in UC set for subsets CD4 T/Tfh/Treg
c2 = ggplot(uc.hc.df) +
  geom_point(data = subset(uc.hc.df, isFOXP3 %in% "FOXP3-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  geom_point(data = subset(uc.hc.df, isFOXP3 %in% "FOXP3+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("FOXP3-", "FOXP3+"),
                     values = c("gray","darkgreen")) +
  labs(color = "Healthy controls") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-11, 3), breaks = seq(-10,0,5))

d2 = ggplot(uc.inf.df) +
  geom_point(data = subset(uc.inf.df, isFOXP3 %in% "FOXP3-"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  geom_point(data = subset(uc.inf.df, isFOXP3 %in% "FOXP3+"),
             aes(x = refUMAP_1, y = refUMAP_2, color = isFOXP3), size = dotSize) +
  facet_wrap(~predicted.celltype) +
  scale_color_manual(breaks = c("FOXP3-", "FOXP3+"),
                     values = c("gray","darkgreen")) +
  labs(color = "Ulcerative colitis") +
  theme_cowplot() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = lineSize) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = lineSize) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-11, 3), breaks = seq(-10,0,5))

wrap_plots(c2, d2, nrow = 2)
ggsave("UC_foxp3.png", width = 18, height = 14, unit = "cm", dpi = 600, bg="transparent")


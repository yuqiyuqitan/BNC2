library(Seurat)
library(clustree)
library(dplyr)
library(plyr)
library(ggtree)
library(MetaNeighbor)
library(Matrix)
library(dendextend)
library(SingleCellExperiment)

flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}


setwd('./')

HSatlasNeuron_mt10 = readRDS('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HSatlasNeuro_mt10_integrated.rds')

genes = c('LEPR','BNC2','AGRP','NPY','POMC')

for(i in genes){

HSatlasNeuron_mt10@meta.data[,i] = 'N'

HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10@assays$RNA[i,]>1),i] = 'Y'

}


table(HSatlasNeuron_mt10$LEPR,HSatlasNeuron_mt10$BNC2)

table(HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10$LEPR=='Y' & HSatlasNeuron_mt10$BNC2=='Y'),'SampleID'])

table(HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10$LEPR=='Y' & HSatlasNeuron_mt10$BNC2=='Y'),'Tissue'])


K560_asn = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/Adult_CellTypeAssigns.csv')

treeTrim = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/tree80_L15_treeTrim.rds'))


HSatlasNeuron_mt10@meta.data = cbind(HSatlasNeuron_mt10@meta.data,treeTrim[colnames(HSatlasNeuron_mt10),])


HSatlasNeuron_mt10@meta.data$Nuclei_K560 = K560_asn[match(HSatlasNeuron_mt10$K560,K560_asn$K560) ,'Nuclei_K560']


table(HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10$LEPR=='Y' & HSatlasNeuron_mt10$BNC2=='Y'),'Nuclei_K560'])

#Lepr/Pomc, Lepr/Agrp, Lepr/Npy 

table(HSatlasNeuron_mt10$LEPR,HSatlasNeuron_mt10$POMC)

table(HSatlasNeuron_mt10$LEPR,HSatlasNeuron_mt10$AGRP)

table(HSatlasNeuron_mt10$LEPR,HSatlasNeuron_mt10$NPY)


table(HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10$LEPR=='Y' & HSatlasNeuron_mt10$POMC=='Y'),'BNC2'])

table(HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10$LEPR=='Y' & HSatlasNeuron_mt10$AGRP=='Y'),'BNC2'])

table(HSatlasNeuron_mt10@meta.data[which(HSatlasNeuron_mt10$LEPR=='Y' & HSatlasNeuron_mt10$NPY=='Y'),'BNC2'])


## try clustering and nuclei from H108

humanMeta = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/FIG2Meta_OCT23.csv',row.names=1)


HSatlasNeuron_mt10@meta.data = cbind(HSatlasNeuron_mt10@meta.data,humanMeta[colnames(HSatlasNeuron_mt10),c('H108_Class', 'H108_Nuclei')])

HSarc = subset(HSatlasNeuron_mt10,cells=colnames(HSatlasNeuron_mt10)[which(HSatlasNeuron_mt10$H108_Nuclei=='ARC')])

HSarc@meta.data = cbind(HSarc@meta.data,humanMeta[colnames(HSarc),'H108'])

colnames(HSarc@meta.data)[54]="H108"

DefaultAssay(HSarc) = 'RNA'

HSarc@active.ident=as.factor(HSarc$H108)

dir.create('./Tan')

HSarc.test = HSarc

HSarc.test <- ScaleData(HSarc.test, features=rownames(HSarc.test))


colnames(HSarc.test@meta.data)[31:35] = c("LEPR_Bi","BNC2_Bi","AGRP_Bi","NPY_Bi","POMC_Bi")

colnames(HSarc@meta.data)[31:35] = c("LEPR_Bi","BNC2_Bi","AGRP_Bi","NPY_Bi","POMC_Bi")

HSarc$LEPR_BNC2 = 'N'
HSarc@meta.data[which(HSarc$LEPR_Bi=='Y' & HSarc$BNC2_Bi=='Y'),'LEPR_BNC2'] = 'Y'

## Plotting LEPR and BNC2 exp


DefaultAssay(HSarc) = 'RNA'


pdf('./Tan/test_Arc_BCN2_LEPR_Exp.pdf',width=12,height=8)
FeaturePlot(HSarc,features = 'BNC2')
FeaturePlot(HSarc,features = 'LEPR')
DimPlot(HSarc,group.by = 'BNC2_Bi',order='Y')
DimPlot(HSarc,group.by = 'LEPR_Bi',order='Y')
DimPlot(HSarc,group.by = 'LEPR_BNC2',order='Y')
dev.off()



HSarc$LEPR_POMC = 'N'
HSarc@meta.data[which(HSarc$LEPR_Bi=='Y' & HSarc$POMC_Bi=='Y'),'LEPR_POMC'] = 'Y'

HSarc$LEPR_NPY = 'N'
HSarc@meta.data[which(HSarc$LEPR_Bi=='Y' & HSarc$NPY_Bi=='Y'),'LEPR_NPY'] = 'Y'


pdf('./Tan/test_Arc_BCN2_POMC_NPY_LEPR_Exp.pdf',width=12,height=8)
FeaturePlot(HSarc,features = 'BNC2')
FeaturePlot(HSarc,features = 'LEPR')
FeaturePlot(HSarc,features = 'POMC')
FeaturePlot(HSarc,features = 'NPY')
DimPlot(HSarc,group.by = 'BNC2_Bi',order='Y')
DimPlot(HSarc,group.by = 'LEPR_Bi',order='Y')
DimPlot(HSarc,group.by = 'POMC_Bi',order='Y')
DimPlot(HSarc,group.by = 'NPY_Bi',order='Y')
DimPlot(HSarc,group.by = 'LEPR_BNC2',order='Y')
DimPlot(HSarc,group.by = 'LEPR_POMC',order='Y')
DimPlot(HSarc,group.by = 'LEPR_NPY',order='Y')
dev.off()


## https://divingintogeneticsandgenomics.com/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

## Vln plots for visualization
pdf('./Tan/test_Arc_VlnPlt_Heatmap.pdf',width=12,height=8)
VlnPlot(HSarc.test,features = rev(c('NOS1','BNC2','POMC','NPY','LEPR','GAD2','GAD1','SLC32A1','SLC17A6')),group.by='H108')
StackedVlnPlot(obj = HSarc.test, features = rev(c('NOS1','BNC2','POMC','NPY','LEPR','GAD2','GAD1','SLC32A1','SLC17A6')))


DoHeatmap(HSarc.test,features = rev(c('NOS1','BNC2','POMC','NPY','LEPR','GAD2','GAD1','SLC32A1','SLC17A6')),group.by='H108',label=FALSE,raster=FALSE,assay='RNA')
dev.off()


pdf('./Tan/AllTri2_Top_NonNeuDEG_heatmap.pdf',width=12,height=8)
DoHeatmap(AllTri2.integrated.test, features = topPerDif,group.by="Maj_Clust_From_Ref",label=FALSE,raster=FALSE,assay='RNA') + NoLegend()
dev.off()


for(i in c('K16','K28','K40','K55','K116','K169','K199','K229')){

HSarc.test@active.ident=as.factor(HSarc.test@meta.data[,i])

pdf(paste0('./Tan/test_Res_',i,'Arc_VlnPlt_Heatmap.pdf'),width=12,height=8)
print(VlnPlot(HSarc.test,features = rev(c('NOS1','BNC2','POMC','NPY','LEPR','GAD2','GAD1','SLC32A1','SLC17A6')),group.by=i))

print(DoHeatmap(HSarc.test,features = rev(c('NOS1','BNC2','POMC','NPY','LEPR','GAD2','GAD1','SLC32A1','SLC17A6')),group.by=i,label=FALSE,raster=FALSE,assay='RNA'))
dev.off()
}

test=HSarc@meta.data[,c('LEPR','BNC2','AGRP','NPY','POMC')]

test[which(test$LEPR=='Y'),'LEPR'] = 'LEPR'
test[which(test$BNC2=='Y'),'BNC2'] = 'BNC2'
test[which(test$AGRP=='Y'),'AGRP'] = 'AGRP'
test[which(test$NPY=='Y'),'NPY'] = 'NPY'
test[which(test$POMC=='Y'),'POMC'] = 'POMC'

test$comb = apply( test[ , 1:5 ] , 1 , paste , collapse = "_" )


## out of ARC, - 16819 

## lepr - 2231
## BNC2 - 137
## both -28 

## fishers test - p=0.0218


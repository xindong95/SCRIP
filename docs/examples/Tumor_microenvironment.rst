BCC Tumor Microenvironment 
======================================

We applied SCRIP to a Basal Cell Carcinoma (BCC) tumor microenvironment (TME) dataset `(Satpathy et al., Nat Biotechnol, 2019) <https://doi.org/10.1038/s41587-019-0206-z>`_ to investigate how TRs and their target genes were changed in different cell states under disease status. Data was downloaded from GEO with accession `GSE129785 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785>`_.

.. code:: shell

    SCRIP enrich -i example/TME/data/BCC_peak_count.h5 -s hs -p example/TME/BCC_SCRIP -t 32
    SCRIP enrich -i example/TME/data/Tcell_peak_count.h5 -s hs -p example/TME/Tcell_SCRIP -t 32

.. code:: r

    library(Seurat)
    library(plyr)
    library(patchwork)
    library(dplyr)
    library(pheatmap)
    library(stringr)

    se <- readRDS("example/TME/data/scATAC_TME_All_SummarizedExperiment.final.rds")
    all_tumor_table <- read.table("example/TME/BCC_SCRIP/enrichment/SCRIP_enrichment.txt",comment.char = "")

    cluster <- mapvalues(rownames(all_tumor_table),se@colData@listData[["Group_Barcode"]],se@colData@listData[["Clusters"]])
    rownames(all_tumor_table) <- paste(cluster,1:36657,"Total_Post#TTACG",sep = "_")
    all_tumor_table_norm <- apply(all_tumor_table, 2, function (x) (x-min(x))/(max(x)-min(x)))
    all_tumor_seurat <- CreateSeuratObject(counts = t(all_tumor_table_norm),project = "t_cell")
                            
    all_tumor_seurat <- NormalizeData(all_tumor_seurat)
    all.genes <- rownames(all_tumor_seurat)
    all_tumor_seurat <- ScaleData(all_tumor_seurat, features = all.genes)

    cluster_all <- c("Cluster17","Cluster18","Cluster19","Cluster20")
    cluster_normal <- paste("Cluster",c(1:16),sep = "")
    out <- vector("list", length(cluster_all))
    for (i in seq_along(cluster_all)) {
        cluster_markers <- FindMarkers(all_tumor_seurat, 
                                        ident.1 = cluster_all[i], 
                                        ident.2 = cluster_normal,
                                        logfc.threshold = 1)
        out[[i]] <- rownames(cluster_markers)
    }
    out <- out[out %in% colnames(all_tumor_table)]

    CLU <- paste("Cluster",1:20,sep = "")
    CLU_type <- c("Naive CD4 T","TH17","Tfh","Treg","Naive CD8 T","Th1","Memory CD8 T","CD8 TEx","Effector CD8 T","NK1","NK2","B","Plasma B","Myeloid","Endothelial","Fibroblasts","Tumor 1","Tumor 2","Tumor 3","Tumor 4")

    tumor_type <- cbind(cluster = se@colData@listData[["Clusters"]],barcode = se@colData@listData[["Group_Barcode"]])
    tumor_type <- as.data.frame(tumor_type)
    
    tumor_cell_f <- all_tumor_table[,out]
    tumor_cell_f$barcode <- rownames(tumor_cell_f)
    tem <- inner_join(tumor_cell_f,tumor_type,by = "barcode")
    TF_mean <- as.data.frame(group_by(tem, cluster) %>% summarise_each(funs = mean))
    rownames(TF_mean) <- mapvalues(TF_mean$cluster, CLU, CLU_type)
    TF_mean <- TF_mean[,-c(1,length(TF_mean))]

    TF_mean <- apply(TF_mean, 2, function (x) (x-min(x))/(max(x)-min(x)))
    all_tumor <- pheatmap(t(TF_mean), border_color= NA, fontsize_row = 13, fontsize_col = 17, treeheight_col = 0, treeheight_row = 0)

.. image:: ../_static/img/Tumors/Tumor_heatmap.png
    :alt: tumor heatmap
    :width: 50%
    :align: center


.. code:: r

    tftable <- read.table("example/TME/Tcell_SCRIP/enrichment/SCRIP_enrichment.txt")

.. image:: ../_static/img/Tumors/Tcell_na_ex.png
    :alt: t cell heatmap
    :width: 30%
    :align: center


.. code:: shell

    SCRIP impute -i example/TME/data/Tcell_peak_count.h5 -s hs -p example/TME/Tcell_SCRIP/ -f h5 --factor JUNB
    SCRIP target -i example/TME/Tcell_SCRIP/imputation/imputed_JUNB.h5ad -s hs -o JUNB_target.h5ad

.. code:: r

    cluster_markers_all <- FindMarkers(test_seurat, 
                                     ident.1 = "Cluster13", 
                                     ident.2 = "Cluster17",
                                     logfc.threshold = 0.1)
    volcano <- cluster_markers_all[,c(1,2)]
    colnames(volcano) <- c("Pvalue","Foldchange")
    threshold<-as.factor((volcano$Foldchange>0.25|volcano$Foldchange<(-0.25))&volcano$Pvalue<0.01)
    r = ggplot(volcano,aes(Foldchange,-log10(Pvalue),colour=threshold)) +
        geom_point() +
        labs(title="Volcanoplot") +
        theme(plot.title = element_text(hjust = 0.25)) +
        xlim(-1,1) +
        theme_classic() +
        geom_vline(xintercept=c(-0.25,0.25),linetype="dotted",size=1) +
        geom_hline(yintercept=-log10(0.01),col="blue")
    r

.. image:: ../_static/img/Tumors/JUNB_target_vol.png
    :alt: JUNB target vol
    :width: 50%
    :align: center

.. code:: r

    library(plyr)
    library(ggplot2)
    library(Seurat)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)



.. image:: ../_static/img/Tumors/junb_go.png
    :alt: JUNB target heatmap
    :width: 50%
    :align: center
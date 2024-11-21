## run pyscenic
{
    library(Seurat)
    library(stringr)
    library(readr)
    library(dplyr)
    library(loomR)
    library(pheatmap)

    seurat2loom <- function(inRDS, outloom){
        dat.loom <- Seurat::as.loom(inRDS, filename = outloom, assay = "RNA", verbose = FALSE, overwrite = TRUE)
        dat.loom$close_all()
    }

    Tcells <- readRDS("./result/Tcells/Tcells_integrated2.RDS")
    DefaultAssay(Tcells) <- "RNA"
    Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
    seurat2loom(Tcells, "./result/Tcells/pyscenic/Tcells.loom")
    # dat <- t(as.matrix(Tcells_integrated@assays$RNA@data))
    # write.csv(dat, "./result/Tcells/pyscenic/Tcells.csv")    
}



```
{
    ####@ parameters
    DATABASE="/public/workspace/zhumy/ref/SCENIC"
    OUTDIR="/public/workspace/zhumy/AITL/result/Tcells/pyscenic/"
    # grnboost2
    expression_mtx_fname="${OUTDIR}/Tcells.loom"
    tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
    grnboost2_output="${OUTDIR}/Tcells.adjacencies.tsv"
    # ctx
    database_fname1=${DATABASE}/hg19-500bp-upstream-7species.mc9nr.feather
    database_fname2=${DATABASE}/hg19-tss-centered-10kb-7species.mc9nr.feather
    annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
    ctx_output="${OUTDIR}/Tcells.regulons.tsv"
    # aucell
    aucell_output="${OUTDIR}/Tcells.auc_mtx.csv"

    python ${DATABASE}/arboreto_with_multiprocessing.py \
        $expression_mtx_fname $tfs_fname --method grnboost2 \
        --output $grnboost2_output \
        --num_workers 8 --seed 12345
    pyscenic ctx \
        $grnboost2_output \
        $database_fname1 $database_fname2 \
        --annotations_fname $annotations_fname \
        --expression_mtx_fname $expression_mtx_fname \
        --output $ctx_output \
        --num_workers 6 --mode "custom_multiprocessing" 
    pyscenic aucell \
        $expression_mtx_fname $ctx_output \
        --output $aucell_output \
        --num_workers 8 --seed 12345    
}


import pandas as pd
import os
from pyscenic.cli.utils import load_signatures
from pyscenic.cli.utils import save_enriched_motifs
regulon_file="./result/Tcells/pyscenic/Tcells.regulons.tsv"
regulons = load_signatures(regulon_file)
#save_enriched_motifs(regulons, './result/Tcells/pyscenic/Tcells.regulons.gmt')

with open("./result/Tcells/pyscenic/Tcells.regulons.txt", "w") as f:
    for i in regulons:
        if len(i.genes) >= 20:
            info_ = i.name + "\t" + "NA" + "\t" + ",".join(i.genes) + "\n"
            f.write(info_)
with open("./result/Tcells/pyscenic/Tcells.regulons_weight.txt", "w") as f:
    for i in regulons:
        if len(i.genes) >= 20:
            info_ = i.name + "\t" + "NA" + "\t" + ",".join([ str(x) for x in i.weights ]) + "\n"
            f.write(info_)


```

{
    Tcells <- readRDS("./result/Tcells/Tcells_integrated2.RDS")
    Tcells1 <- subset(Tcells, cells = colnames(Tcells)[Tcells$subtype != "Unclassified"])
    aucmtx <- read_csv("/public/workspace/zhumy/AITL/result/Tcells/pyscenic/Tcells.auc_mtx.csv")
    aucmtx <- as.data.frame(aucmtx)
    rownames(aucmtx) <- aucmtx$Cell
    aucmtx <- aucmtx[,-1]
    aucmtx1 <- aucmtx[colnames(Tcells1), ]
    colnames(aucmtx1) <- str_replace_all(colnames(aucmtx1), "\\(\\+\\)", "")
    genes <- intersect(colnames(aucmtx1), rownames(Tcells1))
}

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
  #eval = FALSE
)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages({library(scCloneVar)
  library(Seurat)})

data(scCloneVar_test_demo)

seurat_obj <- scCloneVar_test_demo

## -----------------------------------------------------------------------------
seurat_obj


## ----fig.width=4, fig.height=2.5, out.width="70%", dpi=300--------------------
comparison_list <- list(
  list(seurat_obj = seurat_obj,samples = grep("^Y_vitro", unique(seurat_obj$sampleName), value = TRUE),label_counts = data.frame(
      Sample = grep("^Y_vitro", unique(seurat_obj$sampleName), value = TRUE),
      Replicate = unique(seurat_obj$Rep[grepl("^Y_vitro", seurat_obj$sampleName)]),
      SampleLabel = grep("^Y_vitro", unique(seurat_obj$sampleName), value = TRUE)),
    threshold = 0.05,
    title = "Relative Clone Size Distribution: Y_vitro",
    name = "Y_vitro"),
  list(
    seurat_obj = seurat_obj,
    samples = grep("^O_vitro", unique(seurat_obj$sampleName), value = TRUE),
    label_counts = data.frame(
      Sample = grep("^O_vitro", unique(seurat_obj$sampleName), value = TRUE),
      Replicate = unique(seurat_obj$Rep[grepl("^O_vitro", seurat_obj$sampleName)]),
      SampleLabel = grep("^O_vitro", unique(seurat_obj$sampleName), value = TRUE)),
    threshold = 0.05,
    title = "Relative Clone Size Distribution: O_vitro",
    name = "O_vitro"))

dist_out <- run_clone_distribution_engine(comparison_list)

dist_out$plots$Y_vitro
dist_out$plots$O_vitro

## -----------------------------------------------------------------------------
dist_out$descriptive_stats$Y_vitro

## -----------------------------------------------------------------------------
# Define low and high output marker 

low_output_markers <- c(
  "Mpl","Ifitm3","Ifitm1","Tgm2","H2-K1","Socs2","Mycn","Nupr1","Hacd4",
  "Mllt3","Gda","B2m","Procr","Txnip","Clu","Sult1a1","S100a6","Rpl36a",
  "Tsc22d1","Gng11","Ccnd2","mt-Co3","Tmem176b","Lmo2","Ly6a","Mmrn1",
  "Trim47","St3gal1","Mecom","Pik3r1","Adgrg1","Esam","Ryk","Rps21","Hoxb8",
  "Ccnd1","Uba7","Rps28","Serpina3g","Rpl37a","Fkbp1a","Pdzk1ip1","Selp",
  "Eif4a2","Tmem176a","Bex4","Grb10","Jam3","Ptk2b","Gimap8","Csgalnact1",
  "H2-Q7","App","Mylk","Casp12","Rhof","Laptm4a","Tcf15","Gabarapl1","Ctsl",
  "Bex1","Rpl21","D630039A03Rik","Myof","Glul","Ppic","Tpt1","Rpl5","Cdkn1c",
  "mt-Nd5","Abcg3","Aplp2","Art4","Tbxas1","Cish","Vwf","Scarf1","mt-Atp6",
  "H2-Eb1","H2-D1","Cd74","Iigp1","Aldh1a1","Gstm1","mt-Nd4","mt-Cytb","mt-Co1"
)
high_output_markers <- c(
  "Plac8","H2afy","Cdk6","Cd34","Nkg7","Ptma","Mpo","Cd48","Stmn1","Fam117a",
  "Slc22a3","Adgrg3","Ppia","Car1","Ctsg","Flt3","Lgals1","Muc4","Gpx1",
  "Hmgb2","Ndufa4","Serpinb1a","Ccl9","Oaz1","H2afz","Crip1","Mcm7","Cpa3",
  "Vim","Ybx1","Sell","Sh3bgrl3","H3f3a","Dut","Atpif1","Ran","Hnrnpa2b1",
  "Hdgf","Mcm4","Elane","Rabgap1l","Cmtm7","Rpsa","Mcm6","Plek","Set","Atp5g3",
  "Myc","Taldo1","Tuba1b","Cks2","Slc25a5","Fam65a","Cebpa","Tmsb10","Cd52",
  "Klf1","Anp32b","Hn1","Parvg","Ffar2","Bex6","Emilin2","Itgal","Cox5a",
  "Hnrnpab","Ighm","BC035044","Lmnb1","Golm1","Bin1","Igfbp4","Tyrobp","E2f8",
  "Banf1","Cst7","Sh2d5","Aqp1","Ptbp3","Snrpf","Rfc2","Ramp1","Hmgn5","Nme1"
)
# Can be set to NULL if no comparison needed, for simplicity we show a few genes here
ref_deg_genes_list <- c("Olfm3","Scin","Tox","Pbx3","Mecom","Fry","Slc24a3","Pcdh7" )



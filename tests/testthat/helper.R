# Helper objects reused across test files

make_fake_deg_list <- function() {
  c("TP53", "MYC", "BRCA1", "EGFR", "VEGFA", "CDK2", "MDM2",
    "BCL2", "PTEN", "RB1", "AKT1", "STAT3", "CCND1", "IRF1")
}

make_fake_ora_df <- function() {
  data.frame(
    ID          = c("hsa04110", "hsa04115"),
    Description = c("Cell cycle", "p53 signaling pathway"),
    pvalue      = c(0.001, 0.005),
    qvalue      = c(0.01, 0.04),
    Count       = c(8L, 5L),
    Gene        = c("TP53,CDK2,CCND1", "TP53,MDM2,PTEN"),
    stringsAsFactors = FALSE
  )
}

make_fake_gsea_df <- function() {
  data.frame(
    ID              = c("hsa04110", "hsa04151"),
    Description     = c("Cell cycle", "PI3K-Akt signaling"),
    NES             = c(1.85, -1.62),
    enrichmentScore = c(0.72, -0.58),
    pvalue          = c(0.002, 0.008),
    Gene            = c("CDK2,CCND1,TP53", "AKT1,PTEN,MYC"),
    stringsAsFactors = FALSE
  )
}

make_fake_gsva_summary <- function() {
  data.frame(
    GeneSet   = c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_APOPTOSIS"),
    MeanScore = c(0.45, -0.38),
    SDScore   = c(0.12, 0.20),
    Genes     = c("CDK2,CCND1", "TP53,BCL2"),
    stringsAsFactors = FALSE
  )
}

make_fake_roast_df <- function() {
  data.frame(
    GeneSet   = c("CELL_CYCLE_KEGG", "APOPTOSIS_KEGG"),
    PValue    = c(0.001, 0.042),
    Direction = c("Up", "Down"),
    Genes     = c("CDK2,CCND1,MYC", "TP53,BCL2,PTEN"),
    stringsAsFactors = FALSE
  )
}

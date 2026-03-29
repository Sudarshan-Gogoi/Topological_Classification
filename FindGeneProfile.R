FindGeneProfile <- function(
    cancer_type = "Breast",
    recurrence_threshold = 0.5,
    excel_file = "Significant_Genes.xlsx",
    h5_filename = NULL,
    geo_tar_filename = NULL,
    barcodes_filename = NULL,
    features_filename = NULL,
    matrix_filename = NULL,
    project_name = "Cancer",
    base_path = NULL
){
  
  # ------------------------------
  # Install and Load Packages
  # ------------------------------
  
  install_if_needed <- function(pkg){
    if(!require(pkg, character.only = TRUE)){
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  
  install_if_needed("Seurat")
  install_if_needed("Matrix")
  install_if_needed("data.table")
  install_if_needed("readxl")
  
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
  
  if(!requireNamespace("AnnotationDbi", quietly = TRUE))
    BiocManager::install("AnnotationDbi")
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  if(is.null(base_path)){
    stop("Please provide base_path")
  }
  
  # ------------------------------
  # 1 Load Dataset
  # ------------------------------
  
  if(!is.null(h5_filename)){
    
    filepath <- file.path(base_path, h5_filename)
    
    data <- Read10X_h5(filepath)
    
    if(is.list(data)){
      
      if("Gene Expression" %in% names(data)){
        cts <- data$`Gene Expression`
      } else {
        cts <- data[[1]]
      }
      
    } else {
      cts <- data
    }
    
  }
  
  else if(!is.null(geo_tar_filename)){
    
    tar_path <- file.path(base_path, geo_tar_filename)
    
    untar(tar_path, exdir = base_path)
    
    files <- list.files(base_path, full.names = TRUE)
    
    barcodes <- files[grep("barcodes", files)]
    features <- files[grep("features", files)]
    matrix <- files[grep("matrix", files)]
    
    cts <- ReadMtx(
      mtx = matrix[1],
      features = features[1],
      cells = barcodes[1]
    )
  }
  
  else if(!is.null(barcodes_filename)){
    
    cts <- ReadMtx(
      mtx = file.path(base_path, matrix_filename),
      features = file.path(base_path, features_filename),
      cells = file.path(base_path, barcodes_filename)
    )
  }
  
  else{
    stop("Provide dataset input")
  }
  
  # ------------------------------
  # 2 Convert to Seurat Object
  # ------------------------------
  
  Cancer.seurat.obj <- CreateSeuratObject(
    counts = cts,
    project = project_name,
    min.cells = 10,
    min.features = 300
  )
  
  dataset_genes <- rownames(Cancer.seurat.obj)
  
  # Remove version numbers
  dataset_genes_clean <- sub("\\..*", "", dataset_genes)
  
  # ------------------------------
  # 3 Load Significant Genes Excel
  # ------------------------------
  
  excel_path <- file.path(base_path, excel_file)
  
  sig_df <- readxl::read_excel(excel_path, .name_repair = "unique")
  
  cancer_cols <- grep(cancer_type, colnames(sig_df), value = TRUE)
  
  if(length(cancer_cols) == 0){
    stop("No columns found for cancer type")
  }
  
  cancer_data <- sig_df[, cancer_cols]
  
  # ------------------------------
  # 4 Calculate Recurring Genes
  # ------------------------------
  
  total_samples <- ncol(cancer_data)
  
  threshold <- ceiling(recurrence_threshold * total_samples)
  
  gene_counts <- list()
  
  for(i in 1:ncol(cancer_data)){
    
    genes <- unique(na.omit(cancer_data[[i]]))
    
    for(g in genes){
      
      if(!is.null(gene_counts[[g]])){
        gene_counts[[g]] <- gene_counts[[g]] + 1
      } else {
        gene_counts[[g]] <- 1
      }
      
    }
  }
  
  gene_counts_vec <- unlist(gene_counts)
  
  selected_genes <- names(gene_counts_vec[gene_counts_vec >= threshold])
  
  selected_genes_clean <- sub("\\..*", "", selected_genes)
  
  # ------------------------------
  # 5 Detect Gene Type
  # ------------------------------
  
  dataset_is_ensembl <- mean(grepl("^ENSG", dataset_genes_clean)) > 0.5
  selected_is_ensembl <- mean(grepl("^ENSG", selected_genes_clean)) > 0.5
  
  # ------------------------------
  # 6 Convert Gene IDs if Needed
  # ------------------------------
  
  if(dataset_is_ensembl == selected_is_ensembl){
    
    converted_genes <- selected_genes_clean
    
  }
  
  else if(selected_is_ensembl & !dataset_is_ensembl){
    
    converted_genes <- mapIds(
      org.Hs.eg.db,
      keys = selected_genes_clean,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    converted_genes <- na.omit(converted_genes)
    
  }
  
  else if(!selected_is_ensembl & dataset_is_ensembl){
    
    converted_genes <- mapIds(
      org.Hs.eg.db,
      keys = selected_genes_clean,
      column = "ENSEMBL",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    converted_genes <- na.omit(converted_genes)
    
  }
  
  selected_genes <- unique(converted_genes)
  
  # ------------------------------
  # 7 Save Recurring Genes
  # ------------------------------
  
  recurring_file <- file.path(
    base_path,
    paste0(project_name,"_",cancer_type,"_recurring_genes.csv")
  )
  
  data.table::fwrite(
    data.frame(Gene = selected_genes),
    recurring_file
  )
  
  # ------------------------------
  # 8 Extract Expression Matrix
  # ------------------------------
  
  expr_matrix <- as.matrix(
    GetAssayData(Cancer.seurat.obj, layer = "counts")
  )
  
  # Clean matrix gene IDs
  expr_genes_clean <- sub("\\..*", "", rownames(expr_matrix))
  
  rownames(expr_matrix) <- expr_genes_clean
  
  # ------------------------------
  # 9 Find Common Genes
  # ------------------------------
  
  common_genes <- intersect(selected_genes, rownames(expr_matrix))
  
  if(length(common_genes) == 0){
    stop("No overlapping genes between dataset and recurring genes")
  }
  
  filtered_matrix <- expr_matrix[common_genes, , drop = FALSE]
  
  # ------------------------------
  # 10 Calculate Gene Expression Profile
  # ------------------------------
  
  gene_expression <- Matrix::rowMeans(filtered_matrix)
  
  gene_profile <- data.frame(
    Gene = names(gene_expression),
    Expression = gene_expression
  )
  
  # ------------------------------
  # 11 Save Gene Expression Profile
  # ------------------------------
  
  profile_file <- file.path(
    base_path,
    paste0(project_name,"_",cancer_type,"_gene_expression_profile.csv")
  )
  
  data.table::fwrite(gene_profile, profile_file)
  
  message("Recurring genes saved at: ", recurring_file)
  message("Gene expression profile saved at: ", profile_file)
  
  return(gene_profile)
  
}

# Example Files

# Run .h5 Files (Example)
FindGeneProfile(
  
  cancer_type = "Breast",
  
  recurrence_threshold = 0.5,
  
  excel_file = "Significant_Genes.xlsx",
  
  h5_filename = "V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5",
  
  base_path = "D:/Topologicalcancer-signatures",
  
  project_name = "BreastCancer"
  
)

# Run .tar files (Example)
FindGeneProfile(
  
  cancer_type = "Breast",
  
  recurrence_threshold = 0.5,
  
  excel_file = "Significant_Genes.xlsx",
  
  geo_tar_filename = "GSE279219_RAW.tar",
  
  base_path = "D:/Topologicalcancer-signatures",
  
  project_name = "BreastCancer"
  
)

# Run feature barcode matrix files (Example)
FindGeneProfile(
  
  cancer_type = "Breast",
  
  recurrence_threshold = 0.5,
  
  excel_file = "Significant_Genes.xlsx",
  
  barcodes_filename = "GSM9102020_Sample6_barcodes.tsv.gz",
  
  features_filename = "GSM9102020_Sample6_features.tsv.gz",
  
  matrix_filename = "GSM9102020_Sample6_matrix.mtx.gz",
  
  base_path = "D:/Topologicalcancer-signatures",
  
  project_name = "Healthy-Breast"
  
)




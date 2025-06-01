FindCorrMat <- function(h5_filename = NULL, 
                        geo_tar_filename = NULL, 
                        project_name = "Cancer", 
                        num_pcs = 13, 
                        n_top_genes = 300, 
                        alpha = 0.05,
                        base_path = NULL) {
  
  # Install and load necessary libraries
  install_if_needed <- function(package) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
  
  install_if_needed("Seurat")
  install_if_needed("tidyverse")
  install_if_needed("data.table")
  install_if_needed("Hmisc")
  install_if_needed("Matrix")
  
  # Check if base path is provided
  if (is.null(base_path)) {
    stop("Please provide a base path for input/output files.")
  }
  
  # Define full paths for input and output files
  h5_filepath <- if (!is.null(h5_filename)) file.path(base_path, h5_filename) else NULL
  geo_tar_filepath <- if (!is.null(geo_tar_filename)) file.path(base_path, geo_tar_filename) else NULL
  output_path <- base_path
  
  # Step 1: Load the Data
  if (!is.null(h5_filename)) {
    message("Loading 10x Genomics data from ", h5_filepath)
    Cancer.sparse.m <- Read10X_h5(filename = h5_filepath)
    cts <- Cancer.sparse.m
  } else if (!is.null(geo_tar_filename)) {
    message("Extracting GEO data from ", geo_tar_filepath)
    untar(geo_tar_filepath, exdir = base_path)
    
    # List the extracted files
    extracted_files <- list.files(base_path, full.names = TRUE)
    message("Extracted files: ", paste(extracted_files, collapse = ", "))
    
    # Try to map the extracted files to expected names
    barcodes_path <- NULL
    features_path <- NULL
    matrix_path <- NULL
    
    for (file in extracted_files) {
      if (grepl("barcodes", file) && grepl(".tsv.gz", file)) {
        barcodes_path <- file
      } else if (grepl("features", file) && grepl(".tsv.gz", file)) {
        features_path <- file
      } else if (grepl("matrix", file) && grepl(".mtx.gz", file)) {
        matrix_path <- file
      }
    }
    
    # Check if all required files are found
    if (is.null(barcodes_path) || is.null(features_path) || is.null(matrix_path)) {
      stop("Could not find all required files (barcodes, features, matrix) in the extracted tar file.")
    }
    
    # Read the extracted GEO files
    barcodes <- read.table(barcodes_path, sep = "\t", header = FALSE)
    features <- read.table(features_path, sep = "\t", header = FALSE)
    matrix_data <- readMM(matrix_path)
    matrix_sparse <- as(matrix_data, "CsparseMatrix")
    count_matrix <- as.matrix(matrix_sparse)
    rownames(count_matrix) <- features$V1
    colnames(count_matrix) <- barcodes$V1
    cts <- count_matrix
  } else {
    stop("Please provide either an H5 filename or a GEO tar filename.")
  }
  
  # Step 2: Initialize Seurat Object
  Cancer.seurat.obj <- CreateSeuratObject(counts = cts, project = project_name, min.cells = 10, min.features = 300)
  
  # Step 3: Quality Control Plots
  Cancer.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(Cancer.seurat.obj, pattern = "^MT-")
  
  # Create violin plot
  graphics.off()  # Close all devices before starting
  violin_plot_path <- file.path(output_path, "violin_plot.png")
  png(violin_plot_path, width = 1600, height = 1600, units = "px", res = 300)
  print(VlnPlot(Cancer.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  Sys.sleep(2)
  dev.off()
  
  # Create scatter plot
  graphics.off()  # Close all devices before starting
  scatter_plot_path <- file.path(output_path, "scatter_plot.png")
  png(scatter_plot_path, width = 1600, height = 1600, units = "px", res = 300)
  print(FeatureScatter(Cancer.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm'))
  Sys.sleep(2)
  dev.off()
  
  # Step 4: Normalize Data
  Cancer.seurat.obj <- NormalizeData(Cancer.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Step 5: Identify Variable Features
  Cancer.seurat.obj <- FindVariableFeatures(Cancer.seurat.obj, selection.method = "vst", nfeatures = 4000)
  top2 <- head(VariableFeatures(Cancer.seurat.obj), 2)
  
  # Create variable features plot
  graphics.off()  # Close all devices before starting
  variable_features_plot_path <- file.path(output_path, "variable_features_plot.png")
  png(variable_features_plot_path, width = 4000, height = 2000, units = "px", res = 600)
  plot1 <- VariableFeaturePlot(Cancer.seurat.obj)
  print(LabelPoints(plot = plot1, points = top2))
  Sys.sleep(2)
  dev.off()
  
  # Step 6: Scale Data
  all.genes <- rownames(Cancer.seurat.obj)
  Cancer.seurat.obj <- ScaleData(Cancer.seurat.obj, features = all.genes)
  
  # Step 7: PCA and Elbow Plot
  Cancer.seurat.obj <- RunPCA(Cancer.seurat.obj, features = VariableFeatures(object = Cancer.seurat.obj))
  
  # Create elbow plot
  graphics.off()  # Close all devices before starting
  elbow_plot_path <- file.path(output_path, "elbow_plot.png")
  png(elbow_plot_path, width = 1600, height = 1600, units = "px", res = 300)
  print(ElbowPlot(Cancer.seurat.obj))
  Sys.sleep(2)
  dev.off()
  
  # Step 8: Feature Selection
  gene_loadings <- Loadings(Cancer.seurat.obj)
  selected_genes <- list()
  
  for (pc in 1:num_pcs) {
    ranked_genes <- rownames(gene_loadings)[order(abs(gene_loadings[, pc]), decreasing = TRUE)]
    selected_genes[[pc]] <- ranked_genes[1:n_top_genes]
  }
  
  unique_genes <- unique(unlist(selected_genes))
  
  # Step 9: Construct Gene Correlation Network
  seurat_df <- as.data.frame(GetAssayData(Cancer.seurat.obj, layer = "data"))
  filtered_df <- seurat_df[rownames(seurat_df) %in% unique_genes, ]
  
  m <- filtered_df[, -c(1, ncol(filtered_df))]
  m1 <- t(m)
  expression_matrix <- as.matrix(m1)
  
  res <- rcorr(m1)
  n1 <- as.data.frame(res$P)
  n2 <- as.data.frame(res$r)
  
  diag(n1) <- 0
  filtered_correlation <- n2
  filtered_correlation[n1 >= alpha] <- 0
  
  rownames(filtered_correlation) <- rownames(filtered_df)
  colnames(filtered_correlation) <- rownames(filtered_df)
  
  # Save the correlation matrix as a CSV file
  correlation_matrix_path <- file.path(output_path, "CancerCorrMat.csv")
  fwrite(x = filtered_correlation, row.names = TRUE, file = correlation_matrix_path)
  
  # Print messages about saved files
  message("Violin plot saved at: ", violin_plot_path)
  message("Scatter plot saved at: ", scatter_plot_path)
  message("Variable features plot saved at: ", variable_features_plot_path)
  message("Elbow plot saved at: ", elbow_plot_path)
  message("Cancer correlation matrix saved at: ", correlation_matrix_path)
}

# Example usage:
# FindCorrMat(h5_filename = "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5", base_path = "D:/My PhD information folder second paper/Paper2DatasetsAnalysis/FunctionTest")
FindCorrMat(geo_tar_filename = "GSE235168_RAW.tar", base_path = "D:/My PhD information folder second paper/Paper2DatasetsAnalysis/FunctionTest")

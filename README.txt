## FindCorrMat and Topological_Classification Function Documentation:

# FindCorrMat Function:
# Run this Function first

Overview:
The FindCorrMat function processes single-cell RNA sequencing data (from 10x Genomics or GEO datasets) to create a gene correlation network (CancerCorrMat.csv). It performs quality control, normalization, feature selection, and correlation analysis, and generates diagnostic plots and a correlation matrix.

Requirements:
- R (version 3.6 or higher recommended)
- Required R packages: Seurat, tidyverse, data.table, Hmisc, Matrix
- These packages will install automatically if missing

Input Parameters:
1. h5_filename: Name of 10x Genomics H5 file (optional)
2. geo_tar_filename: Name of GEO tar file (optional)
3. project_name: Project name (default: "Cancer")
4. num_pcs: Number of principal components (User Defined)
5. n_top_genes: Top genes per PC (default: 300)
6. alpha: Significance threshold (default: 0.05)
7. base_path: Directory path for files (required)

Note: Provide either h5_filename OR geo_tar_filename, not both.

Output Files (saved in base_path):
1. violin_plot.png - Quality control plots
2. scatter_plot.png - Feature scatter plot
3. variable_features_plot.png - Variable features
4. elbow_plot.png - PCA elbow plot
5. CancerCorrMat.csv - Final correlation matrix

Usage Examples:

For 10x Genomics data:
FindCorrMat(h5_filename = "filtered_feature_bc_matrix.h5", 
           base_path = "D:/MyData/Analysis")

For GEO data:
FindCorrMat(geo_tar_filename = "GSE123456_RAW.tar", 
           base_path = "D:/MyData/Analysis")

Processing Steps:
1. Load data (H5 or GEO tar)
2. Create Seurat object and QC
3. Normalize and scale data
4. Identify variable features
5. PCA analysis
6. Select genes based on PCA
7. Build correlation matrix
8. Generate outputs (plots and CSV)

Notes:
- Required packages auto-install if missing
- Must specify base_path containing input file
- GEO tar needs standard 10x files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
- Processing time increases with dataset size
- Close all graphic devices before running if errors occur
----------------------------------------------------------------------------------------------------------
# Topological_Classification Function:
# Run this Function after running FindCorrMat Function

Overview:
This function performs topological analysis on gene correlation networks (CancerCorrMat.csv) to classify cancer types. It uses persistent homology and simplicial complexes to extract topological features from correlation matrices.

Requirements:
- Python 3.6 or higher
- Required packages: pandas, numpy, matplotlib, gudhi, persim, scipy, networkx
- These packages will install automatically if missing

Input Parameters:
1. input_csv_path: Path to input correlation matrix CSV file (CancerCorrMat.csv)
2. output_directory: Directory to save all output files (required)
3. desired_num_holes: Desired number of 0-dimensional holes (default: 600)
4. max_lambda: Maximum lambda values for persistence landscape (default: 10)

Output Files:
1. 0-th Betti Numbers Plot.png - Shows 0-dimensional holes
2. FilteredDistanceMatrix.csv - Filtered correlation matrix
3. SignificantGenes.csv - List of significant genes
4. Persistence Diagram.png - Persistence diagram visualization
5. Persistence Barcode.png - Persistence barcode visualization
6. simplicial_complex_plot.png - Network visualization
7. simplicial_complex_data.csv - Network statistics
8. feature_vector1.csv - Extracted topological features
9. landscape_plot.png - Persistence landscape
10. LandscapeAreas.csv - Calculated landscape areas
11. TrendlinePlot.png - Area trendline
12. TrendlineDetails.csv - Trendline statistics
13. ClassificationResults.xlsx - Final classification results

Processing Steps:
1. Reads input correlation matrix
2. Computes persistent homology
3. Filters the distance matrix
4. Visualizes topological features
5. Computes persistence landscapes
6. Extracts topological features
7. Performs cancer classification
8. Generates all output files

Classification Models:
1. Model 1: Uses 9 Persistent Homology features for cancer type prediction
   - Compares against breast, lung, colorectal, ovarian, prostate references
   - Returns top 2 predictions with similarity scores

2. Model 2: Uses Euler characteristic for prostate cancer detection
   - Returns Positive/Negative prediction based on threshold

Example Usage:
results = Topological_Classification(
    input_csv_path='D:/Path/To/CancerCorrMat.csv',
    output_directory='D:/Path/To/Output',
    desired_num_holes=600,
    max_lambda=10
)

Important Notes:
1. Input CSV must be a square correlation matrix with gene names
2. Output directory must exist before running
3. Processing time depends on matrix size
4. Visualizations may take several seconds to render

Output Interpretation:
1. Top 1 Prediction: Most likely cancer type
2. Top 1 Similarity (%): Confidence score (0-100)
3. Top 2 Prediction: Second most likely cancer type
4. Most_Similar_Cancer: Best match based on radius_life_range
5. Prostate Cancer: Positive/Negative prediction

Troubleshooting:
1. If no simplices exist at threshold:
   - Try reducing desired_num_holes
   - Check input matrix quality
2. If visualizations fail:
   - Close all matplotlib windows before running
   - Update matplotlib package
3. If installation fails:
   - Run pip install manually for each package
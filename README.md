# markdown
# # Inferring Ancient Relatedness: A Statistical Framework
#
# This repository provides a statistical framework to infer ancient relatedness between individuals or populations by estimating the genetic distances between them. Our approach bridges two widely used methods in population genetics and archaeogenetics: **Principal Component Analysis (PCA)** and **F-statistics**.
#
# ## Overview
#
# Both PCA and F-statistics are essential tools in genetic studies:
#
# - **PCA**: A dimensionality reduction technique used for visualizing population structure and forming hypotheses about historical admixtures.
# - **F-statistics**: A set of statistics used for statistical testing of admixture.
#
# Although these methods are often used together, they make slightly different assumptions, use different models, and rely on different software. Our framework **combines PCA and F-statistics** into a joint analysis, ensuring consistent assumptions across both methods.
#
# ### Key Advantages
# - Unveils the effect of modeling assumptions.
# - Provides recommendations for interpreting PCA-based analyses.
# - Offers a better understanding of how different PCA algorithms emphasize population structure versus sampling noise.
#
# ## Implementations
#
# This framework includes different PCA implementations:
# - **Probabilistic PCA (PPCA)**
# - **Latent Subspace Estimation (LSE)**
# - **Classical PCA**
#
# We find that F-statistics are more naturally interpreted in a PPCA or LSE-based framework. Our method also allows for accurate F-statistics estimation from PPCA, even in the presence of large amounts of missing data.
#
# ## Getting Started
#
# ### Clone the Repository
# To get started, clone the repository to your current directory:
#
# ```bash
# git clone https://github.com/DivyaratanPopli/pca_paper
# ```
#
# ### Running the Pipeline
#
# To run the pipeline on the toy example data, follow these steps:
#
# 1. Navigate to the real data directory:
#     ```bash
#     cd Snakemake_pipelines/real_data
#     ```
#
# 2. Run the pipeline using Snakemake with the following command:
#     ```bash
#     snakemake all --cores 50 --config npcs=2
#     ```
#
#     - `--cores 50`: Use 50 cores for computation (you can adjust this based on your system).
#     - `npcs=2`: Set the number of principal components for probabilistic PCA (you can adjust this as needed).
#
# ### Output
#
# The output files will be stored in the following directories:
#
# - `data_files/f2mat_mean_ppca_scale{npcs}.csv`: Mean matrix for estimated F2s from PPCA for all individuals.
# - `data_files/f2mat_std_ppca_scale{npcs}.csv`: Standard deviation matrix for estimated F2s from PPCA.
# - `data_files/admixtools2/f2mat_mean.csv`: Mean matrix for estimated F2s using ADMIXTOOLS 2 for all populations.
# - `data_files/admixtools2/f2mat_std.csv`: Standard deviation matrix for estimated F2s using ADMIXTOOLS 2.
#
# ## Calculating F3 and F4 Statistics
#
# You can calculate **F3** and **F4** statistics as linear combinations of F2s using the following equations (Patterson, 2012):
#
# - **F3(X1; X2, X3)**:
#     ```
#     F3(X1; X2, X3) = ( F2(X1, X2) + F2(X1, X3) - F2(X2, X3) ) / 2
#     ```
#
# - **F4(X1, X2; X3, X4)**:
#     ```
#     F4(X1, X2; X3, X4) = ( F2(X1, X3) + F2(X2, X4) - F2(X1, X4) - F2(X2, X3) ) / 2
#     ```
#
# ### Standard Deviations
# You can calculate standard deviations for **F3** and **F4** using the standard deviations of F2s, applying standard statistical methods.
#
# ## Example Datasets
#
# We have provided the following datasets and pipelines for testing:
#
# - Example input files with Neandertal data are available in:
#     - `Snakemake_pipelines/real_data/toy_example`
#
# - Pipelines for simulating genetic data and testing our method are available in:
#     - `Snakemake_pipelines/simulations`
#
# ## Contributing
#
# Feel free to open issues, suggest enhancements, or submit pull requests if you'd like to contribute to this project.
#
# ## License
#
# This project is licensed under the [MIT License](LICENSE).
#
# ---
#
# ### References
# - Patterson, N., et al. (2012). "Ancient admixture in human history."

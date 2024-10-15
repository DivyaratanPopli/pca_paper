This is a statistical framework to infer ancient relatedness between individuals or populations by estimating the genetic distances between them. Our framework provides the link between two widely used methods, principal component analysis (PCA) and F-statistics, which are routinely used in population genetic and archaeogenetic studies. PCA is a dimentionality reduction technique that is used for visualizing structure in populations, and for forming hypotheses about past admixtures. In many studies, PCA is followed by statistical testing for admixture using methods based on a set of statistics called F-statistics. It has been shown that PCA and F-statistics are closely related analyses and reveal the same biological signal. Many ancient genetic studies use both of these tools, but make slightly different assumptions, using slightly different models and different software for them.

Thus, we combine PCA and F-statistics into a joint analysis. This framework allows us to ensure that the assumptions for PCA and F-statistics are consistent throughout the analysis. The key advantage is that the effect of modelling assumptions becomes apparent, and this also allows us to make novel recommendations about how PCA-based analyses should be performed and interpreted. The connections between F-statistics and PCA allow us to provide a better understanding on how different PCA algorithms emphasize population structure versus sampling noise.

Here, we provide different implementations of PCA - probabilistic PCA (PPCA), Latent Subspace Estimation (LSE) and classical PCA - and find that F-statistics are more naturally interpreted in a PPCA or LSE-based framework. Using this framework, F-statistics can be accurately estimated from PPCA in the presence of large amounts of missing data.

We have provided a pipeline to implement this method in the repository Snakemake_pipelines/real_data. We also have example input files with Neandertal data in the repository Snakemake_pipelines/real_data/toy_example. The pipelines that we used to simulate genetic data and to test our method are available in the repository Snakemake_pipelines/simulations.

To run the pipeline on the data files in toy example, follow these steps:

git clone https://github.com/DivyaratanPopli/pca_paper #clone the repository to your current directory
cd Snakemake_pipelines/real_data
snakemake all --cores 50 --config npcs=2 #For example, use 50 cores, and the number of principle componenets for probabilistic PCA is 2

Output:
data_files/f2mat_mean_ppca_scale{npcs}.csv , data_files/f2mat_std_ppca_scale{npcs}.csv contain matrix with mean and standard deviations respectively for estimated F2s from PPCA for all individuals.
data_file/admixtools2/f2mat_mean.csv , data_files/admixtools2/f2mat_std.csv contain matrix with mean and standard deviations respectively for estimated F2s using ADMIXTOOLS 2 for all populations.

F3s and F4s can then be calculated as linear combinations of F2s.

To calculate mean F3 and F4 values using mean F2s, use the following equations (from Patterson, 2012):
F3(X1;X2,X3) = [ F2(X1,X2) + F2(X1,X3) - F2(X2,X3) ]/2
F4(X1,X2;X3,X4) = [ F2(X1,X3) + F2(X2,X4) - F2(X1,X4) - F2(X2,X3) ]/2

One can also calculate standard deviations for F3 and F4 using standard deviations of F2s using the standard statistics.

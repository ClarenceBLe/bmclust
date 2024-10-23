# set CRAN mirror--network of ftp and web servers that store identical and updated versions of code and R documentation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install packages
install.packages("cli")
library(cli)
install.packages("NbClust")
library(NbClust)
install.packages("cluster")
library(cluster)
install.packages("kohonen")
library(kohonen)
install.packages("stats")
library(stats)
install.packages("fpc")
library(fpc)
install.packages("dbscan")
library(dbscan)
install.packages("clusterCrit")
library(clusterCrit)

# define functions for command-line interface (CLI)
create_directories <- function(base_dir) {
  if(!dir.exists(base_dir)){
    dir.create(base_dir, recursive=TRUE)
    cli_alert_success("base directory created: {base_dir}.")
  }
  subdirs <- c("cluster_ensemble", "cluster_quality", "cluster_stats", "cluster_plots")
  for(subdir in subdirs){
    if(!dir.exists(file.path(base_dir,subdir))){
      dir.create(file.path(base_dir,subdir), recursive=TRUE)
    }
  }
}

read_scale_files <- function(filename) {
  cli_alert_info("reading in {filename}.")
  raw_df <- read.csv(filename, header=T, sep=,)
  dropcat_df <- unlist(lapply(raw_df, is.numeric), use.names=FALSE)
  numerical_df <- raw_df[,dropcat_df]
  df <- as.data.frame(scale(numerical_df)) 
  cli_alert_info("categorical variables removed and data scaled.")
}

optimal_cluster_count <- function(x,d) {
  nbclust_wardD <- NbClust(x, distance="{d}", method="ward.D", min.nc=2, max.nc=6, index="silhouette")
  cli_alert_success("optimal clusters based on 'ward.D' method and '{d}' distance--saved as 'nbclust_wardD$Best.nc'.")
  cli_alert_success("{nbclust_wardD$Best.nc}")
}

distance_dissimilarity_matrices <- function(x,d) {
  distance_matrix <- dist(x, method="{d}")
  cli_alert_success("distance matrix computed and saved as 'distance_matrix'")
  
  dissim_matrix <- daisy(x, metric="{d}")
  cli_alert_success("dissimilarity matrix computed and saved as 'dissim_matrix'")
}

ensemble_clustering <- function(x,k,d,cluster_outdir) {
  #k-means
  cluster_kmeans <- kmeans(x, centers=k) 
  kmeans_append_df <- data.frame(cluster_kmeans$cluster, x)
  write.csv(kmeans_append_df,file="{cluster_outdir}/kmeans_appended_df.csv")
  cli_alert_success("kmeans cluster-appended dataframe saved to '{cluster_outdir}/kmeans_appended_df.csv'.")
  cli_alert_success("{k} clusters identied with sizes: {cluster_kmeans$size}")
  cli_alert_success("{k} clusters centered at: {cluster_kmeans$centers}")
  
  #partitioning around medoids (PAM)
  cluster_pam <- pam(x, k, metric="{d}")
  pam_append_df <- data.frame(cluster_pam$clustering, x)
  write.csv(pam_append_df,file="{cluster_outdir}/pam_appended_{d}_df.csv")
  cli_alert_success("PAM cluster-appended dataframe saved to '{cluster_outdir}/pam_appended_{d}_df.csv'.")
  cli_alert_success("{k} clusters identied with sizes: {cluster_pam$clusinfo}")
  cli_alert_success("{k} clusters centered at: {cluster_pam$medoids}")
  
  #clustering large applications (CLARA)
  cluster_clara <- clara(x, k, metric="{d}", stand=TRUE)
  clara_append_df <- data.frame(cluster_clara$clustering, x)
  write.csv(clara_append_df,file="{cluster_outdir}/clara_appended_{d}_df.csv")
  cli_alert_success("CLARA cluster-appended dataframe saved to '{cluster_outdir}/clara_appended_{d}_df.csv'.")
  cli_alert_success("{k} clusters identied with sizes: {cluster_clara$clusinfo}")
  cli_alert_success("{k} clusters centered at: {cluster_clara$medoids}")
  
  #hierarchical clustering (HCLUST)
  cluster_hclust <- hclust(distance_matrix, method="ward.D")
  hclust_cut <- cutree(cluster_hclust, k=k) #cut tree into 'k' clusters
  hclust_rect <- rect.hclust(cluster_hclust, k=k, border="blue") #draw boxes around clusters
  hclust_append_df <- data.frame(hclust_cut, x) 
  write.csv(hclust_append_df,file="{cluster_outdir}/hclust_appended_df.csv")
  cli_alert_success("HCLUST cluster-appended dataframe saved to '{cluster_outdir}/hclust_appended_df.csv'.")
  
  #agglomerative nesting hierarchical clustering (AGNES)
  cluster_agnes <- agnes(dissim_matrix, diss=TRUE, method="ward")
  agnes_cut <- cutree(as.hclust(cluster_agnes), k=k) 
  agnes_rect <- rect.hclust(cluster_agnes, k=k, border="blue") 
  agnes_append_df <- data.frame(agnes_cut, x)
  write.csv(agnes_append_df,file="{cluster_outdir}/agnes_appended_df.csv")
  cli_alert_success("AGNES cluster-appended dataframe saved to '{cluster_outdir}/agnes_appended_df.csv'.")
  
  #divisive analysis clustering (DIANA) 
  cluster_diana <- diana(x, metric="{d}", stand=TRUE) 
  diana_cut <- cutree(as.hclust(cluster_diana), k=k) 
  diana_rect <- rect.hclust(cluster_diana, k=k, border="blue") 
  diana_append_df <- data.frame(diana_cut, x)
  write.csv(diana_append_df,file="{cluster_outdir}/diana_appended_df.csv")
  cli_alert_success("DIANA cluster-appended dataframe saved to '{cluster_outdir}/diana_appended_df.csv'.")
  
  #fuzzy analysis clustering (FANNY)
  cluster_fanny <- fanny(x=dissim_matrix, metric="{d}", k=k, diss=TRUE, stand=TRUE)
  fanny_append_df <- data.frame(cluster_fanny$clustering, x)
  write.csv(fanny_append_df,file="{cluster_outdir}/fanny_appended_{d}_df.csv")
  cli_alert_success("FANNY cluster-appended dataframe saved to '{cluster_outdir}/fanny_appended_{d}_df.csv'.")
  cli_alert_success("additional cluster information accessible at 'cluster_fanny$silinfo'.")
}

cluster_quality_assess <- function(x, cluster_kmeans, cluster_pam, cluster_clara, hclust_cut, agnes_cut, diana_cut, cluster_fanny) {
  for(col_name in names(x)){
    dunn_kmeans <- dunn.test(x[[col_name]], as.factor(cluster_kmeans$clusters), method ="Bonferroni", kw=TRUE, label=TRUE, list=TRUE)
    write.table(output, file="/Users/clarencele/Desktop/dunn_kmeans.csv", row.names=FALSE, sep="\t")
  }
}

#command-line arguments parsing using base R
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  cli_alert_danger("no arguments provided.")
  quit(status = 1)
}

#function dispatcher based on the first argument
command <- args[1]
additional_arguments <- as.list(args[-1])

#call function dispatcher with command and arguments
do.call(dispatch, c(command, additional_arguments))

#use case in terminal
Rscript bmclust.R 

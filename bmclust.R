options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("cli")
library(cli)

# define functions
import_packages <- function() {
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
  install.packages("calinhara")
  library(calinhara)
}

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
  return(df)
}

optimal_cluster_count <- function(df,d) {
  nbclust_wardD <- NbClust(df, distance="{d}", method="ward.D", min.nc=2, max.nc=6, index="silhouette")
  cli_alert_success("optimal clusters based on 'ward.D' method and '{d}' distance--saved as 'nbclust_wardD$Best.nc'.")
  cli_alert_success("{nbclust_wardD$Best.nc}")
}

distance_dissimilarity_matrices <- function(df,d) {
  distance_matrix <- dist(df, method="{d}")
  cli_alert_success("distance matrix computed and saved as 'distance_matrix'")
  
  dissim_matrix <- daisy(df, metric="{d}")
  cli_alert_success("dissimilarity matrix computed and saved as 'dissim_matrix'")
  
  return(list(distance=distance_matrix, dissimilarity=dissim_matrix))
}

ensemble_clustering <- function(df,k,d,base_dir,matrices) {
  # access distance and dissimilarity matrices 
  distance_matrix <- matrices$distance
  dissim_matrix <- matrices$dissimilarity
  
  # k-means
  cluster_kmeans <- kmeans(df, centers=k) 
  kmeans_append_df <- data.frame(cluster_kmeans$cluster, df)
  file_path <- file.path(base_dir, "cluster_ensemble", "kmeans_appended_df.csv")
  write.csv(kmeans_append_df,file=file_path)
  cli_alert_success("kmeans cluster-appended dataframe saved to '{file_path}'.")
  cli_alert_success("{k} clusters identied with sizes: {cluster_kmeans$size}")
  cli_alert_success("{k} clusters centered at: {cluster_kmeans$centers}")
  
  # partitioning around medoids (PAM)
  cluster_pam <- pam(df, k, metric="{d}")
  pam_append_df <- data.frame(cluster_pam$clustering, df)
  file_path <- file.path(base_dir, "cluster_ensemble", "pam_appended_{d}_df.csv")
  write.csv(pam_append_df,file=file_path)
  cli_alert_success("PAM cluster-appended dataframe saved to '{file_path}'.")
  cli_alert_success("{k} clusters identied with sizes: {cluster_pam$clusinfo}")
  cli_alert_success("{k} clusters centered at: {cluster_pam$medoids}")
  
  # clustering large applications (CLARA)
  cluster_clara <- clara(df, k, metric="{d}", stand=TRUE)
  clara_append_df <- data.frame(cluster_clara$clustering, df)
  file_path <- file.path(base_dir, "cluster_ensemble", "clara_appended_{d}_df.csv")
  write.csv(clara_append_df,file=file_path)
  cli_alert_success("CLARA cluster-appended dataframe saved to '{file_path}'.")
  cli_alert_success("{k} clusters identied with sizes: {cluster_clara$clusinfo}")
  cli_alert_success("{k} clusters centered at: {cluster_clara$medoids}")
  
  # hierarchical clustering (HCLUST)
  cluster_hclust <- hclust(distance_matrix, method="ward.D")
  hclust_cut <- cutree(cluster_hclust, k=k) 
  hclust_rect <- rect.hclust(cluster_hclust, k=k, border="blue") 
  hclust_append_df <- data.frame(hclust_cut, df) 
  file_path <- file.path(base_dir, "cluster_ensemble", "hclust_appended_df.csv")
  write.csv(hclust_append_df,file=file_path)
  cli_alert_success("HCLUST cluster-appended dataframe saved to '{file_path}'.")
  
  # agglomerative nesting hierarchical clustering (AGNES)
  cluster_agnes <- agnes(dissim_matrix, diss=TRUE, method="ward")
  agnes_cut <- cutree(as.hclust(cluster_agnes), k=k) 
  agnes_rect <- rect.hclust(cluster_agnes, k=k, border="blue") 
  agnes_append_df <- data.frame(agnes_cut, df)
  file_path <- file.path(base_dir, "cluster_ensemble", "agnes_appended_df.csv")
  write.csv(agnes_append_df,file=file_path)
  cli_alert_success("AGNES cluster-appended dataframe saved to '{file_path}'.")
  
  # divisive analysis clustering (DIANA) 
  cluster_diana <- diana(df, metric="{d}", stand=TRUE) 
  diana_cut <- cutree(as.hclust(cluster_diana), k=k) 
  diana_rect <- rect.hclust(cluster_diana, k=k, border="blue") 
  diana_append_df <- data.frame(diana_cut, df)
  file_path <- file.path(base_dir, "cluster_ensemble", "diana_appended_df.csv")
  write.csv(diana_append_df,file=file_path)
  cli_alert_success("DIANA cluster-appended dataframe saved to '{file_path}'.")
  
  # fuzzy analysis clustering (FANNY)
  cluster_fanny <- fanny(x=dissim_matrix, metric="{d}", k=k, diss=TRUE, stand=TRUE)
  fanny_append_df <- data.frame(cluster_fanny$clustering, df)
  file_path <- file.path(base_dir, "cluster_ensemble", "fanny_appended_{d}_df.csv")
  write.csv(fanny_append_df,file=file_path)
  cli_alert_success("FANNY cluster-appended dataframe saved to '{file_path}'.")
  cli_alert_success("additional cluster information accessible at 'cluster_fanny$silinfo'.")

  #return cluster variables
  return(list(c_kmeans=cluster_kmeans,c_pam=cluster_pam,c_clara=cluster_clara,c_hclust=cluster_hclust,c_agnes=cluster_agnes,c_diana=cluster_diana,c_fanny=cluster_fanny))
}

plot_clusters <- function(mlist,k,ec) {
  for(m in mlist){
    if(m=="pam"){
      file_path <- file.path(base_dir, "cluster_plots", "pam_multi_cluster_silhouette.png")
      png(file_path, width = 1200, height = 600)
      layout(matrix(1:2, nrow = 1)) 
      plot(ec$c_pam, which.plot = 1, main = "pam cluster")
      plot(ec$c_pam, which.plot = 2, main = "pam silhouette")
      dev.off()
      cli_alert_success("PAM cluster plot saved to {file_path}.")
      }
    elif(m=="clara"){
      file_path <- file.path(base_dir, "cluster_plots", "clara_multi_cluster_silhouette.png")
      png(file_path, width = 1200, height = 600)
      layout(matrix(1:2, nrow = 1))  
      plot(ec$c_clara, which.plot = 1, main = "clara cluster")
      plot(ec$c_clara, which.plot = 2, main = "clara silhouette")
      dev.off()
      cli_alert_success("CLARA cluster plot saved to {file_path}.")
      }
    elif(m=="hclust"){
      file_path <- file.path(base_dir, "cluster_plots", "hclust_dendrogram.png")
      png(file_path, width = 1200, height = 600)
      plot(ec$c_hclust, main="hclust dendrogram",xlab="observations",ylab="height")
      rect.hclust(ec$c_hclust,k=k,border="blue")
      dev.off()
      cli_alert_success("HCLUST cluster plot saved to {file_path}.")
      }
    elif(m=="agnes"){
      file_path <- file.path(base_dir, "cluster_plots", "agnes_multi_banner_dendrogram.png")
      png(file_path, width = 1200, height = 600)
      layout(matrix(1:2, nrow = 1))  
      plot(ec$c_agnes, which.plot = 1, main = "agnes banner")
      plot(ec$c_agnes, which.plot = 2, main = "agnes dendrogram")
      rect.hclust(ec$c_agnes,k=k,border="blue")
      dev.off()
      cli_alert_success("AGNES cluster plot saved to {file_path}.")
      }
    elif(m=="diana"){
      file_path <- file.path(base_dir, "cluster_plots", "diana_multi_banner_dendrogram.png")
      png(file_path, width = 1200, height = 600)
      layout(matrix(1:2, nrow = 1))
      plot(ec$c_diana, which.plot = 1, main = "diana banner")
      plot(ec$c_diana, which.plot = 2, main = "diana dendrogram")
      rect.hclust(ec$c_diana,k=k,border="blue")
      dev.off()
      cli_alert_success("DIANA cluster plot saved to {file_path}.")
      }
    elif(m=="fanny"){
      file_path <- file.path(base_dir, "cluster_plots", "fanny_silhouette.png")
      png(file_path, width = 1200, height = 600)
      plot(ec$c_fanny, main = "fanny silhouette")
      dev.off()
      cli_alert_success("FANNY cluster plot saved to {file_path}.")
      }
  }
}

# command-line (CLI) argument parsing
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  cli_alert_danger("incorrect input.")
  cli_alert_danger("correct format: Rscript bmclust.R <base_dir, str> <filename, str> <distance metric, str> <number of clusters, int>")
  cli_alert_danger("example: Rscript bmclust.R . data.csv manhattan 2")
  quit(status = 1)
}

# define parameters based on user input
base_dir <- args[1] # complete base directory where output files created/stored
filename <- args[2] # raw data-frame to be assessed
d <- args[3] # distance metric--i.e."euclidean", "manhattan"
k <- as.integer(args[4]) # number of clusters 
mlist <- c("pam","clara","hclust","agnes","diana","fanny") # subset list for specific methods to be visualized

# function dispatcher
import_packages()
create_directories(base_dir)
df <- read_scale_files(filename)
optimal_cluster_count(df,d)
matrices <- distance_dissimilarity_matrices(df,d) 
ec <- ensemble_clustering(df,k,d,base_dir,matrices)
plot_clusters(mlist,k,ec)

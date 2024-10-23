# bmclust
* Perform ensemble clustering on numerical dataset--kmeans, pam, clara, hclust, agnes, diana, fanny
## How to run it
```
cd bmclust 
```
* Run bmclust.R
```
Rscript bmclust.R <base_dir, str> <filename, str> <distance metric, str> <number of clusters, int>
i.e. Rscript bmclust.R . /path/to/data.csv manhattan 2
```
## Summary of the pipeline
* Necessary packages are imported
* Output directories are created within user-input base directory
* Raw data imported in *.csv format
* Categorical columns are removed and numerical columns are scaled
* Optimal number of clusters determined using 'NbClust' package--(select from 30 indices)
* Distance and dissimilarity matrices computed from scaled data
* Ensemble clustering performed based on scaled data, distance matrix, and dissimilarity matrix
* Cluster assignments appended onto scaled dataframe for downstream statistical analysis

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
* Imports necessary packages
* Creates output sub-directories based on user-input base directory
* Reads in raw data in *.csv format
* Removes categorical columns and scales numerical data columns
* Determine optimal number of clusters based on 'NbClust' package--(select from 30 indices)
* Calculate distance and dissimilarity matrices based on scaled data
* Performs ensemble clustering based on scaled data, distance matrix, and dissimilarity matrix
* Cluster assignments appended onto scaled dataframe for downstream statistical analysis

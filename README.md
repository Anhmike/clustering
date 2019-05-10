# Clustering

## Introduction

This repository contains numerous clustering methods, including hierarchical clustering, CURE (Clustering Using REpresentatives), k-means and DBSCAN (Density-Based Spatial Clustering of Applications with Noise). Clustering is a useful technique in data mining and statistical data analysis used to group similar data together and identify patterns in distributions.

The clustering algorithms above have been separated into two scripts - k-means can be found in `kmeans.q`, while the other algorithms can be found in `clust.q`. Additionally, example notebooks have been provided to show how the algorithms perform on a variety of datasets.

A k-dimensional tree (k-d tree) is used by the single and centroid hierarchical algorithms, as well as for CURE which can use both q and C implementations of the k-d tree.

## Requirements

- embedPy

The python packages required to allow successful exectution of all functions within the machine learning toolkit can be installed via:

pip:
```bash
pip install -r requirements.txt
```

or via conda:
```bash
conda install --file requirements.txt
```

*Running of the notebook examples contained within the FRESH section of this library will require the installation of JupyterQ however this is not a dependancy for the running of functions at an individual level.*

## Status

The clustering library is still in development, further improvements will be made to the library in the coming months.

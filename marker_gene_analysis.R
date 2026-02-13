# install.packages('Seurat')
# install.packages("hdf5r")
# install.packages("arrow") #
# install.packages("sf")
# install.packages("Matrix")
# install.packages('BiocManager')
# BiocManager::install('glmGamPoi')
# install.packages("foreach")
# install.packages("doParallel")

library(Seurat)
library(ggplot2)
library(arrow)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)
library(glmGamPoi)
library(foreach)
library(doParallel)

## Set a node_number
node_number <- 2

## Set path to RDA of interest:
path <- "F:/Spacial_transcriptomics/336/336DataObject_newSeurat.rda"

## Set a series of cutoff values to generate data for:
cutoff_values <- c(4,5,6,7,8,9,10,11,12,15,20,30,40,50)

## Load the rda:
load(path)

# Get a list of all object names in the global environment
all_objects <- ls()

# Calculate the size of each object
object_sizes <- sapply(all_objects, function(obj_name) {
    # Use get() to access the object by its name
    object.size(get(obj_name))
})

# Sort the objects by size in descending order
sorted_sizes <- sort(object_sizes, decreasing = TRUE)

dataset_name <- names(sorted_sizes)[1]

# First remove the smaller dataset to save as much space as possible.
rm(list = c(as.character(names(sorted_sizes)[2])))

## Assign the dataset to a placeholder, then delete it, but keep the old name.
placeholder <- get(dataset_name)

rm(list = c(as.character(names(sorted_sizes)[1])))

## Do a PCA of the SCTransformed data.
placeholder <- RunPCA(placeholder)

## Generate an elbow plot to examine where the "elbow" is to determine the
## number of PCs to plot on. This seems to suggest maybe 10-11?
DatasetElbowPlot <- ElbowPlot(placeholder, ndims = 50)

if(
    !file.exists(
        paste0(
            dirname(path),
            "/",
            names(sorted_sizes)[1],
            "_ElbowPlot.png"
        )
    )
) {
    png(
        paste0(
            dirname(path),
            "/",
            names(sorted_sizes)[1],
            "_ElbowPlot.png"
        ),
        width = 1170,
        height = 783
    )
    
    plot(DatasetElbowPlot)
    
    Sys.sleep(30)
    
    dev.off()
}

plot(DatasetElbowPlot)

for(i in cutoff_values) {
    
    # Construct nearest-neighbor graph using top PCs
    placeholder <- FindNeighbors(placeholder, dims = seq_len(i))
    
    # Identify cell clusters
    placeholder <- FindClusters(
        placeholder, 
        resolution = 0.2
    )
    
    ## Add the clusters and mitochondrial gene percentage to the dataset
    ## metadata.
    placeholder@meta.data["percent_mitochondrial"] <- PercentageFeatureSet(
        placeholder, 
        pattern = "^mt-"
    )
    
    placeholder@meta.data[
        as.character(i)
    ] <- placeholder@meta.data$seurat_clusters
    
    ## Create a directory to store results for the given PC in the directory
    ## containing the original .rda.
    cluster_directory <- paste0(
        dirname(path),
        "/PC_",
        i,
        "_clusters/"
    )
    
    if(!dir.exists(cluster_directory)) {
        dir.create(cluster_directory)
    }
    
    ## Identify the counts of cells per cluster.
    if(
        !file.exists(
            paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_cells_by_cluster_data.tsv"
            )
        )
    ) {
        
        ## Identify the sequencing depth per cluster.
        cluster_data <- dplyr::select(
            placeholder@meta.data, 
            seurat_clusters, 
            nCount_Spatial.Polygons,
            percent_mitochondrial
        )
        
        summary_by_cluster <- dplyr::summarize(
            dplyr::group_by(
                cluster_data, 
                seurat_clusters
            ),
            median_nCount_RNA = median(nCount_Spatial.Polygons),
            mean_nCount_RNA = mean(nCount_Spatial.Polygons),
            n_cells = n(),
            mean_mitochondrial_percentage = mean(
                percent_mitochondrial,
                na.rm = TRUE
            )
        )
        
        write.table(
            summary_by_cluster,
            file =  paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_cells_by_cluster_data.tsv"
            ),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t"
        )
    }
    
    # Project cells into 2D using UMAP for visualization
    placeholder <- RunUMAP(
        placeholder, 
        reduction = "pca", 
        dims = seq_len(i)
    )
    
    # Visualize the UMAP plot
    
    ## Create a png for the UMAP:
    if(
        !file.exists(
            paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_UMAP_Plot.png"
            )
        )
    ) {
        UmapPlot <- DimPlot(
            placeholder, 
            reduction = "umap", 
            label = TRUE
        )
        
        png(
            paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_UMAP_Plot.png"
            ),
            width = 1170,
            height = 783
        )
        
        ## Plot the UMAP:
        plot(UmapPlot)
        
        ## Wait:
        Sys.sleep(5)
        
        ## clear plot:
        dev.off()
    }
    
    ## Visualize the clusters:
    if(
        !file.exists(
            paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_Overall_Cluster_Plot.png"
            )
        )
    ) {
        
        ClusterPlot <- SpatialDimPlot(
            object = placeholder,
            group.by = "seurat_clusters",
            images = "slice.polygons",
            plot_segmentations = FALSE,
            pt.size.factor = 0.5
        ) + guides(
            fill = guide_legend(override.aes = list(size=10))
        ) + theme(
            legend.key.size = unit(1, "cm")
        )
        
        png(
            paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_Overall_Cluster_Plot.png"
            ),
            width = 1170,
            height = 783
        )
        
        plot(ClusterPlot)
        
        Sys.sleep(5)
        
        dev.off()
    }
    
    ## Get top cluster markers:
    markers <- FindAllMarkers(
        placeholder, 
        min.pct = 0.25, 
        logfc.threshold = 0.25
    )
    
    if(
        !file.exists(
            paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_top_cluster_markers.tsv"
            )
        )
    ) {
        
        write.table(
            markers,
            file =  paste0(
                cluster_directory,
                dataset_name,
                "_",
                i,
                "_PCAs_top_cluster_markers.tsv"
            ),
            quote = FALSE,
            sep = "\t"
        )
    }
    
    ## Generate plots for each cluster
    
    ## Before doing cluster-level analyses, first let's generate a table with
    ## data on the centroids of each cluster in the UMAP data, which we can
    ## use to identify the closest clusters to each other.
    
    ## Get the cluster IDs and UMAP coordinates of each cell.
    cell_clusters <- Idents(placeholder)
    
    umap_coords <- Embeddings(placeholder[["umap"]])
    
    ## Create a table noting the UMAP coordinates and cluster of each cell.
    cell_data <- data.frame(
        UMAP1 = umap_coords[ ,1], 
        UMAP2 = umap_coords[ ,2], 
        Cluster = cell_clusters
    )
    
    ## Create a final table with the mean UMAP coordinates by each cluster. 
    centroids <- aggregate(
        cbind(UMAP1, UMAP2) ~ Cluster, 
        data = cell_data, 
        FUN = mean
    )
    
    ## Then create a euclidean distance matrix noting the distance each cluster
    ## is from the others.
    distance_matrix <- as.data.frame(
        as.matrix(
            stats::dist(
                centroids[, -1], 
                method = "euclidean"
            )
        )
    )
    
    ## Initialize a holding value for the cluster_number, which will be used
    ## shortly when setting up parallelization using foreach().
    GSMIndices <- NULL
    
    ## Create a list of the nodes that will be used, and then register them for
    ## use with the foreach package.
    makeClusterList <- parallel::makeCluster(node_number)
    
    doParallel::registerDoParallel(makeClusterList)
    
    ## Use foreach() to iterate across the clusters and generate a plot for each
    foreach::foreach(
        cluster_number = seq_along(
            levels(placeholder@meta.data$seurat_clusters)
        )
    ) %dopar% {
        
        ## Show where cells from the given cluster are in the sample.
        if (
            !file.exists(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Location_Plot.png"
                )
            )
        ) {
            
            nthClusterPlot <- Seurat::SpatialDimPlot(
                placeholder, 
                images = "slice.polygons", 
                pt.size.factor = 0.5, 
                cells.highlight = Seurat::CellsByIdentities(
                    object = placeholder, 
                    idents = as.numeric(
                        levels(placeholder@meta.data$seurat_clusters)[
                            cluster_number
                        ]
                    )
                ), 
                facet.highlight = TRUE, 
                ncol = 1
            )
            
            png(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Location_Plot.png"
                ),
                width = 1170,
                height = 783
            )
            
            plot(nthClusterPlot)
            
            Sys.sleep(5)
            
            dev.off()
        }
        
        ## Get the marker genes associated with the given cluster.
        overall_cluster_marker_genes_df <- markers[
            markers$cluster == levels(placeholder@meta.data$seurat_clusters)[
                cluster_number
            ],
        ]
        
        ## Next get the genes that identify the given cluster compared to all
        ## others.
        individual_cluster_marker_genes_df <- Seurat::FindMarkers(
            placeholder, 
            ident.1 = levels(placeholder@meta.data$seurat_clusters)[
                cluster_number
            ]
        )
        
        ## Finally, get genes that identify the cluster compared to its nearest
        ## neighbors.
        
        ## First, set the number of closest neighbors to be relative to the
        ## total clusters identified (25%, rounded up).
        n_nearest_neighbors <- ceiling(
            length(levels(placeholder@meta.data$seurat_clusters))/4
        )
        
        ## Get the distances for the other clusters to the cluster of interest
        ## from the distance matrix.
        distances_to_cluster_of_int <- distance_matrix[
            , cluster_number
        ]
        
        names(distances_to_cluster_of_int) <- levels(
            placeholder@meta.data$seurat_clusters
        )
        
        ## Order the values from least to greatest, then grab the 2nd through 
        ## number of closest neighbors +1 (since the closest will always be 0
        ## for the cluster itself).
        distances_to_cluster_of_int <- sort(
            distances_to_cluster_of_int, 
            decreasing = FALSE
        )
        
        closest_clusters_to_cluster_of_int <- names(
            distances_to_cluster_of_int
        )[2:(n_nearest_neighbors+1)]
        
        ## Now get the marker genes for the cluster of interest relative to
        ## its nearest clusters.
        nearest_clusters_marker_genes_df <- Seurat::FindMarkers(
            placeholder, 
            ident.1 = levels(placeholder@meta.data$seurat_clusters)[
                cluster_number
            ], 
            ident.2 = closest_clusters_to_cluster_of_int
        )
        
        ## Now take the genes found as markers for this cluster in the overall
        ## analysis, and for each set of 9, save a violin plot showing the
        ## expression of these genes.
        overall_cluster_marker_genes_split <- split(
            overall_cluster_marker_genes_df$gene, 
            ceiling(seq_along(overall_cluster_marker_genes_df$gene) / 9)
        )
        
        if(length(overall_cluster_marker_genes_df$gene)>0) {
        
            for(j in seq_len(length(overall_cluster_marker_genes_split))) {
                
                ## Get the marker genes for the given split:
                split_genes <- overall_cluster_marker_genes_split[[j]]
                
                ## Determine what the range of genes being plotted is.
                lower_range_bound <- (((j-1)*9)+1)
                upper_range_bound <- (length(split_genes) + ((j-1)*9))
                
                ## Check to see if a plot has already been created, if it hasn't,
                ## create it.
                if(
                    !file.exists(
                        paste0(
                            cluster_directory,
                            dataset_name,
                            "_",
                            i,
                            "_PCAs_Cluster_",
                            levels(placeholder@meta.data$seurat_clusters)[
                                cluster_number
                            ],
                            "_Overall_Marker_Genes_",
                            lower_range_bound,
                            "_",
                            upper_range_bound,
                            "_Counts.png"
                        )
                    )
                ) {
                    
                    ## Create the violin plot.
                    overall_violin_plot <- Seurat::VlnPlot(
                        placeholder,
                        features = split_genes,
                        layer = "counts"
                    )
                    
                    ## Open the PNG.
                    png(
                        paste0(
                            cluster_directory,
                            dataset_name,
                            "_",
                            i,
                            "_PCAs_Cluster_",
                            levels(placeholder@meta.data$seurat_clusters)[
                                cluster_number
                            ],
                            "_Overall_Marker_Genes_",
                            lower_range_bound,
                            "_",
                            upper_range_bound,
                            "_Counts.png"
                        ),
                        width = 1170,
                        height = 783
                    )
                    
                    ## Plot the violin plot.
                    plot(overall_violin_plot)
                    
                    Sys.sleep(5)
                    
                    ## Close the plot:
                    dev.off()
                }
            }
        } else {
            
            ## Open an empty PNG.
            png(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Overall_Marker_Genes_",
                    0,
                    "_",
                    0,
                    "_Counts.png"
                ),
                width = 1170,
                height = 783,
                bg = "transparent"
            )
            
            ## Close it.
            dev.off()
        }
        
        ## Next, let's identify the nine most overexpressed genes which mark the
        ## cluster of interest compared to all others, which also weren't
        ## already included in the overall marker genes.
        individual_cluster_marker_genes <- rownames(
            individual_cluster_marker_genes_df[
                !(
                    rownames(
                        individual_cluster_marker_genes_df
                    ) %in% overall_cluster_marker_genes_df$gene
                ) &
                individual_cluster_marker_genes_df$avg_log2FC > 0,
            ]
        )[
            seq_len(
                min(
                    9,
                    nrow(
                        individual_cluster_marker_genes_df[
                            !(
                                rownames(
                                    individual_cluster_marker_genes_df
                                ) %in% overall_cluster_marker_genes_df$gene
                            ) &
                                individual_cluster_marker_genes_df$avg_log2FC > 0,
                        ]
                    )
                )
            )
        ]
        
        ## First write a small file that says which clusters are closest to the
        ## cluster of interest.
        if(
            !file.exists(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Closest_Clusters.txt"
                )
            )
        ) {
            writeLines(
                closest_clusters_to_cluster_of_int, 
                con = paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Closest_Clusters.txt"
                )
            )
        }
        
        ## Check to see if a plot for the marker genes for the individual
        ## cluster has already been created, if it hasn't, create it.
        if(
            !file.exists(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Top_9_Individual_Cluster_Marker_Genes_Counts.png"
                )
            )
        ) {
            
            if(length(individual_cluster_marker_genes)>0) {
                ## Create the violin plot.
                individual_cluster_violin_plot <- Seurat::VlnPlot(
                    placeholder,
                    features = individual_cluster_marker_genes,
                    layer = "counts"
                )
                
                ## Open the PNG.
                png(
                    paste0(
                        cluster_directory,
                        dataset_name,
                        "_",
                        i,
                        "_PCAs_Cluster_",
                        levels(placeholder@meta.data$seurat_clusters)[
                            cluster_number
                        ],
                        "_Top_9_Upregulated_Individual_Cluster_Marker_Genes_Counts",
                        ".png"
                    ),
                    width = 1170,
                    height = 783
                )
                
                ## Plot the violin plot.
                plot(individual_cluster_violin_plot)
                
                Sys.sleep(5)
                
                ## Close the plot:
                dev.off()
            } else{
                
                ## Open the blank PNG.
                png(
                    paste0(
                        cluster_directory,
                        dataset_name,
                        "_",
                        i,
                        "_PCAs_Cluster_",
                        levels(placeholder@meta.data$seurat_clusters)[
                            cluster_number
                        ],
                        "_Top_9_Upregulated_Individual_Cluster_Marker_Genes_Counts",
                        ".png"
                    ),
                    width = 1170,
                    height = 783,
                    bg = "transparent"
                )
                
                ## Close the plot:
                dev.off()
            }
        }
        
        ## Now get the top 9 genes which differentiate the given cluster from
        ## its nearest clusters
        nearest_clusters_marker_genes <- rownames(
            nearest_clusters_marker_genes_df[
                !(
                    rownames(
                        nearest_clusters_marker_genes_df
                    ) %in% c(
                        overall_cluster_marker_genes_df$gene,
                        individual_cluster_marker_genes
                    )
                ), 
            ]
        )[
            seq_len(
                min(
                    9,
                    nrow(
                        nearest_clusters_marker_genes_df[
                            !(
                                rownames(
                                    nearest_clusters_marker_genes_df
                                ) %in% c(
                                    overall_cluster_marker_genes_df$gene,
                                    individual_cluster_marker_genes
                                )
                            ), 
                        ]
                    )
                )
            )
        ]
        
        ## Check to see if a plot for the marker genes for the cluster relative
        ## to its nearest neighbors has already been created, if it hasn't,
        ## create it.
        if(
            !file.exists(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Top_9_Nearest_Clusters_Marker_Genes_Counts.png"
                )
            )
        ) {
            
            if(length(nearest_clusters_marker_genes)>0) {
            
                ## Create the violin plot.
                nearest_clusters_violin_plot <- Seurat::VlnPlot(
                    placeholder,
                    features = nearest_clusters_marker_genes,
                    layer = "counts"
                )
                
                ## Open the PNG.
                png(
                    paste0(
                        cluster_directory,
                        dataset_name,
                        "_",
                        i,
                        "_PCAs_Cluster_",
                        levels(placeholder@meta.data$seurat_clusters)[
                            cluster_number
                        ],
                        "_Top_9_Nearest_Clusters_Marker_Genes_Counts.png"
                    ),
                    width = 1170,
                    height = 783
                )
                
                ## Plot the violin plot.
                plot(nearest_clusters_violin_plot)
                
                Sys.sleep(5)
                
                ## Close the plot:
                dev.off()
            } else{
                ## Open the PNG.
                png(
                    paste0(
                        cluster_directory,
                        dataset_name,
                        "_",
                        i,
                        "_PCAs_Cluster_",
                        levels(placeholder@meta.data$seurat_clusters)[
                            cluster_number
                        ],
                        "_Top_9_Nearest_Clusters_Marker_Genes_Counts.png"
                    ),
                    width = 1170,
                    height = 783,
                    bg = "transparent"
                )
                
                ## Close the plot:
                dev.off()
            }
        }
        
        ## Now for all the genes found either from either the overall marker
        ## analysis, the individual cluster marker analysis, or the marker analysis
        ## of the cluster's nearest neighbors, calculate the mean expression across
        ## the clusters, as well as the proportion of cells expressing each gene
        ## in each cluster.
        
        ## Organize the data first.
        genes_by_cell_data <- SeuratObject::FetchData(
            object = placeholder, 
            vars = c(
                overall_cluster_marker_genes_df$gene,
                individual_cluster_marker_genes,
                nearest_clusters_marker_genes, 
                "ident"
            ),
            layer = 'counts'
        )
        
        ## Write a little function to calculate proportion greater than 0.
        proportion_more_than_zero_func <- function(x) {
            return(mean(x>0))
        }
        
        ## Calculate the average expression of each gene per cluster.
        cluster_averages <- data.frame(
            t(
                dplyr::summarize(
                    dplyr::group_by(genes_by_cell_data, ident), 
                    dplyr::across(dplyr::everything(), mean)
                )[,-1]
            )
        )
        
        colnames(cluster_averages) <- paste0(
            "cluster_",
            levels(placeholder@meta.data$seurat_clusters)
        )
        
        ## Calculate the proportion of cells expressing each gene per cluster.
        cluster_expressed_proportion <- data.frame(
            t(
                dplyr::summarize(
                    dplyr::group_by(genes_by_cell_data, ident), 
                    dplyr::across(dplyr::everything(), proportion_more_than_zero_func)
                )[,-1]
            )
        )
        
        colnames(cluster_expressed_proportion) <- paste0(
            "cluster_",
            levels(placeholder@meta.data$seurat_clusters)
        )
        
        ## For each table, note what type of gene each is.
        cluster_averages$overall_cluster_marker_gene <- c(
            rownames(cluster_averages) %in% overall_cluster_marker_genes_df$gene
        )
        
        cluster_averages$individual_cluster_marker_gene <- c(
            rownames(cluster_averages) %in% individual_cluster_marker_genes
        )
        
        cluster_averages$marker_gene_from_nearest_clusters <- c(
            rownames(cluster_averages) %in% nearest_clusters_marker_genes
        )
        
        cluster_expressed_proportion$overall_cluster_marker_gene <- c(
            rownames(cluster_expressed_proportion) %in% overall_cluster_marker_genes_df$gene
        )
        
        cluster_expressed_proportion$individual_cluster_marker_gene <- c(
            rownames(cluster_expressed_proportion) %in% individual_cluster_marker_genes
        )
        
        cluster_expressed_proportion$marker_gene_from_nearest_clusters <- c(
            rownames(cluster_expressed_proportion) %in% nearest_clusters_marker_genes
        )
        
        ## Save the plots if they are not already saved.
        if(
            !file.exists(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Marker_Genes_Average_Expression_Counts_By_Cluster.tsv"
                )
            )
        ) {
            write.table(
                cluster_averages,
                file = paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Marker_Genes_Average_Expression_Counts_By_Cluster.tsv"
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE
            )
        }
        
        if(
            !file.exists(
                paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Marker_Genes_Expression_Proportion_By_Cluster.tsv"
                )
            )
        ) {
            write.table(
                cluster_expressed_proportion,
                file = paste0(
                    cluster_directory,
                    dataset_name,
                    "_",
                    i,
                    "_PCAs_Cluster_",
                    levels(placeholder@meta.data$seurat_clusters)[
                        cluster_number
                    ],
                    "_Marker_Genes_Expression_Proportion_By_Cluster.tsv"
                ),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE
            )
        }
    }
    
    ## Clear the global environment and save the remaining variables:
    rm(
        list = setdiff(
            ls(),
            c(
                "placeholder",
                "path",
                "sorted_sizes",
                "cutoff_values",
                "i",
                "dataset_name"
            )
        ),
        envir = .GlobalEnv
    )
    
    ## rename the placeholder object based on the original dataset name.
    assign(
        dataset_name,
        placeholder
    )
    
    rm(placeholder)
    
    save.image(
        paste0(
            dirname(path),
            "/",
            names(sorted_sizes)[1],
            "_",
            i,
            "_PCAs_data.rda"
        )
    )
}

## Save cluster identities across the PC cutoffs selected.
if(
    !file.exists(
        paste0(
            dirname(path),
            "/",
            names(sorted_sizes)[1],
            "_",
            paste(
                cutoff_values,
                collapse = "_"
            ),
            "_PCAs_identities.tsv"
        )
    )
) {
    
    write.table(
        placeholder@meta.data[,as.character(cutoff_values)],
        file =  paste0(
            dirname(path),
            "/",
            names(sorted_sizes)[1],
            "_",
            paste(
                cutoff_values,
                collapse = "_"
            ),
            "_PCAs_identities.tsv"
        ),
        quote = FALSE,
        row.names = FALSE,
        sep = "\t"
    )
}
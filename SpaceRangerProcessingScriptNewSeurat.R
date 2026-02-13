library(arrow, lib.loc = "/project2/ilaird_1342/Spacial_Transcriptomics/r_env_spacial/lib/R/library")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)
library(glmGamPoi)

## List out a vector of directories to the SpaceRanger outputs. These should
## be the directories that contain the 'spatial', 'segmented_outputs', and
## 'binned_outputs' subdirectories, along with various other files output by
## Spaceranger.
dataDirs <- c(
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/v3d_example",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/206/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/215/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/244/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/336/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/Visium_HD_Spatial_Gene_Expression_Library_Human_Pancreas_FFPE",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/Visium_HD_Spatial_Gene_Expression_Library_Mouse_Kidney_FFPE",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/A2C1",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/209/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/238/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/327_dan_HE_swap_352/outs",
    "/project2/ilaird_1342/Spacial_Transcriptomics/samples/352_dan_HE_swap_327/outs"
)

## Set a vector of sample names to name the datasets after. These should be
## paired on a 1:1 basis with the directories listed in the dataDirs vector,
## i.e. the first of the sampleNames should match with the first of the
## dataDirs entries.
sampleNames <- c(
    "v3d",
    "206",
    "215",
    "244",
    "336",
    "FFPE_Human_Pancreas",
    "FFPE_Mouse_Kidney",
    "Jonathan_A2C1",
    "209",
    "238",
    "327",
    "352"
)

## Finally, set a vector of genes to use for the one example plot this will
## generate. Again, these will be matched 1:1 with the elements in the
## dataDirs and sampleNames.
geneNames <- c(
    "rho",
    "Sftpc",
    "Sftpc",
    "Sftpc",
    "Sftpc",
    "CTRB1",
    "Hnf4a",
    "SFTPC",
    "Sftpc",
    "Sftpc",
    "Sftpc",
    "Sftpc"
)

## Write a function which will load the relevant space ranger data from the
## directories, subset and normalize them, and perform some initial checks
## (for now) to ensure my analysis has been performed correctly.
spaceRangerProcessorFunction <- function(
    internalDataDir,
    internalSampleName,
    internalGeneName
) {

    ## First, load the data in the specified internalDataDir,
    internalObjPlaceholder <- Load10X_Spatial(
        data.dir = internalDataDir,
        slice= "slice",
        bin.size= c("polygons")
    )

    ## Subset the object to include just genes that have at least one read
    ## across all the cells.
    internalObjPlaceholderSub <- subset(
        internalObjPlaceholder ,
        nCount_Spatial.Polygons>0
    )

    ## Now do the SCTranform on the subsetted dataset.
    internalObjPlaceholderSub <- SCTransform(
        internalObjPlaceholderSub,
        assay          = "Spatial.Polygons",
        new.assay.name = "Polygon"
    )

    DefaultAssay(internalObjPlaceholderSub) <- "Polygon"

    ## Let's perform and save a pair of plots just to confirm they work.
    ## This will save the plots in each internalDataDir
    genePlotBig <- SpatialFeaturePlot(
        object = internalObjPlaceholderSub,
        features = internalGeneName,
        images = "slice.polygons",
        plot_segmentations = TRUE,
        pt.size.factor = 0.5
    )

    genePlotZoom <- SpatialFeaturePlot(
        object = internalObjPlaceholderSub,
        features = internalGeneName,
        images = "slice.polygons",
        plot_segmentations = TRUE,
        pt.size.factor = 0.5,
        crop = TRUE,
        stroke = 0.2
    ) + coord_sf(
        xlim = c(200, 300),
        ylim = c(200, 300),
        expand = FALSE
    )

    # Open a PNG device for the larger plot
    png(
        paste0(
            internalDataDir,
            "/",
            "bigPlot_newSeurat.png"
        ),
        width = 1000,
        height = 1000
    )

    # Plot the recorded plot to the device again
    plot(genePlotBig)

    # Close the device to save the file
    dev.off()

    # Open a PNG device for the zoom plot
    png(
        paste0(
            internalDataDir,
            "/",
            "zoomPlot_newSeurat.png"
        ),
        width = 1000,
        height = 1000
    )

    # Plot the recorded plot to the device again
    plot(genePlotZoom)

    # Close the device to save the file
    dev.off()

    ## Now assign the overall data objects to new variables based on the
    ## internalSampleName.
    internalSampleNameFinal <- paste0(
        internalSampleName,
        "_obj"
    )

    internalSampleNameSubsetFinal <- paste0(
        internalSampleName,
        "_obj_sub"
    )

    assign(
        internalSampleNameFinal,
        internalObjPlaceholder,
        envir = .GlobalEnv
    )

    assign(
        internalSampleNameSubsetFinal,
        internalObjPlaceholderSub,
        envir = .GlobalEnv
    )

    ## Save the object.
    save.image(
        paste0(
            internalDataDir,
            "/",
            internalSampleName,
            "DataObject_newSeurat.rda"
        )
    )

    ## Clear memory
    gc()

    ## Remove objects.
    rm(list = setdiff(ls(envir = .GlobalEnv), c("dataDirs", "sampleNames", "geneNames", "spaceRangerProcessorFunction")), envir = .GlobalEnv)
}

## Use mapply with the function to run on all the supplied samples.
mapply(
    FUN = spaceRangerProcessorFunction,
    internalDataDir = dataDirs,
    internalSampleName = sampleNames,
    internalGeneName = geneNames
)

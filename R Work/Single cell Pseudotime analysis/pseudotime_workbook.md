```R
library(Seurat)
library(monocle3)
library(dplyr)
library(slingshot)
library(SingleCellExperiment)
library(phateR)
```

Welcome to the workbook on pseudotime/Trajectory Inference

Throughout this task we are going to be applying two Trajectory Inference (TI) to two different scRNA-seq datasets: one simulated and one real. At each step you will be asked to process the datasets and perform a variety of dimensionality reduction datasets before applying two different TI methods. Following that you will analyse the trajectories; interpreting how accurately they reflect the underlying process

Let's start of by carrying out TI on a scRNA-seq dataset simulated with Dyngen. Simulated datasets have several advantages, for one we can tell it what type of trajectory to make (branching, linear, cyclical) and it calculates the actual pseudotime values of the cells and returns their gene expression.

Using Dyngen we've simulated a dataset which has a linear trajectory, therefore the results of the Trajectory Inference tools should be linear trajectories and the pseudotime values they generate should correlate with the actual pseudotime values.

Let's load in the simulated dataset:


```R
simulated_data <- readRDS("start_shared_wt.rds")
```

Perform a standard Seurat analysis pipeline on the dataset (excluding quality control/filtering steps) up until you generate the UMAP.


```R

```

Just from the UMAP you can see structure in the data, it is one long continious line of cells that goes from one end to the other or, to put it shortly, it's a linear process; which is what we expect.

As this is a simulated dataset, we can also show the actual pseudotime values of the cells, which we do in the code below


```R
FeaturePlot(simulated_data, features = "sim_time")
```

You can see that we can see a continious time gradient that seems to be captured in the UMAP space, i.e. the cells with similar pseudotime values are near each other and are far away from those with different pseudotime values.

Now we have our reduced dimension space let's use some TI methods.

The first method we will use is Monocle, which is the first TI method ever published. It has changed quite a bit over the three versions that have been out so far, with one of the biggest new additions being the ability to split the dataset up into different trajectories if it thinks there are more than one biological process occuring. <b> (foreshadowing incoming)  </b> The problem is that Monocle could end up splitting trajectories up into seperates ones incorrectly.

Compared to Slingshot, which we'll arrive at after this, Monocle is restricted to only calculating pseudotime from UMAP space (you can make it run on other spaces, but it doesn't like it).

You can find the Monocle tutorial at this link: https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/. Work your way through it, adapting the code to fit the simulated dataset, until you reach the 'order_cells' function and plot your cells coloured by pseudotime, and then progress no further.

You'll notice that Monocle has it's own functions for reducing dimensions and clustering the data. Perform them as written in the tutorial, but for the 'preprocess_cds' function, run it first with 50 dims, plot the elbow plot then rerun 'preprocess_cds' with your specific number of dimensions.

Furthermore, when it gets to running the 'order_cells' function, some users may find it returns an error to do with an 'interactive node'. If you get this error, you will have to specify the root state (the cells with 0 sim_time, i.e. the start of the process) through other means, which the tutorial tells you about.


```R

```

Depending on your dimension setting, a lot of you probably have very different looking trajectories. Perhaps you can begin to see some of the issues with using Monocle on your datasets.

Let's see how well the pseudotime values from Monocle correlate with the actual pseudotime values of the simulated cells.

If you have some cells without pseudotime values then you first need to remove them from the actual pseudotime vector and monocle pseudotime vector before calculating correlation.


```R
#Get Monocle pseudotime values of cells 
monocle_pseudo <- pseudotime(simulated_cds) 


```

Even if Monocle did not assign pseudotime values to some of the cells, you should hopefully be able to see that Monocle pseudotime and actual pseudotime are highly correlated, meaning it has captured the underlying biological process well.

Let's test out another pseudotime method, namely Slingshot which lacks a lot of in-built plotting functions, interactivity and simplicity that Monocle has, but what it makes up for is how customisable it is.

With Slingshot you can use any reduced dimension space to build your trajectory from (It will not break up trajectories into smaller ones like Monocle does however) which is useful as many people have issues with UMAP (See the Biorxv paper by Lior Pachter and colleagues (2021) titled 'The Specious Art of Single-Cell Genomics' for more information on that).

Slingshot runs on SingleCellExperiment (sce) objects, so the first step is converting our Seurat dataset to sce which we've given you the code for.

The tutorial for Slingshot can be found here: https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html, but as Slingshot is highly customisable, we don't need to perform normalization, clustering and dimensionality reduction again, we can just use the ones we generated from Seurat. Shortly put, just follow step 3 of the tutorial and do nothing else. We will get to differential expression later.

Follow step 3 of the tutorial and run Slingshot across the UMAP reduced dimension space, using the Seurat clusters as the cluster information.


```R
simulated_sce <- as.SingleCellExperiment(simulated_data, assay = "RNA") # give them this

```

Hopefully you should see a nice smooth trajectory line over the UMAP plot connecting all the clusters together.

It may look nice, but let's see how it compares to the actual pseudotime values being finding the correlation between the Slingshot and actual pseudotime values.


```R

```

You should hopefully see that like Monocle, Slingshot seems to accurately capture the pseudotime values of the simulated dataset.

So far, we've only done TI on the UMAP space, but there are other dimensionality reduction methods that we can calculate pseudotime on using Slingshot. Rerun the Slingshot analysis, but this time give it the PCA embeddings rather than the UMAP embeddings.

See how the pseudotime values generated by Slingshot on the PCA space correlate with the actual pseudotime values and then see how well the UMAP and PCA pseudotime values correlate.


```R

```

You should be able to see that the pseudotime values generated by Slingshot from the PCA and UMAP space are very similar to each other - but not exactly the same. Therefore what dimensionality reduction method you use has an effect on what results you get.

Simulated datasets represent the best case scenario a dataset could be. We control exactly what the shape of the trajectory should be, we control how many cells and genes it simulates and most importantly of all, it is no where near as noisy as real data is.

Let's now apply these TI methods to real data and see how they perform. We've chosen a dataset of the <i> Plasmodium Berghei </i> lifecycle (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9045723/) as we already know what the overall lifecycle stage transitions should be, and thus we can make judgements on how well the TI methods have captured the process.

In order to save you time and computing power we are going to subset the dataset down to just one condition i.e. plasmodium cells found only in the blood (denotated in the Organ metadata slot as: B) as well as remove some of the stages: Outlier, Schizonts, Males and Females.

After you've done that you should be left with two blood replicate sequencing runs. Run the Seurat integration pipeline and stop once you've generated the UMAP for the data



```R
pb_full <- readRDS("Pb.combined.rds") 


```

Now we've integrated the datasets let's put them into Monocle and Slingshot and see what trajectories they return.

Let's start off with Monocle, which has it's own way of remove batch effects so make sure to give it the raw counts from the object and use the 'align_cds' command to remove the batch effects between the two replicates.


```R

```

We don't have real pseudotime values to compare the Monocle trajectory against, but what we do know is the ordering of the cell stages. Follow the Monocle trajectories from your starting point and compare it against the Stage progression and discuss whether Monocle has captured the process or not.

Next we'll move onto Slingshot, but before we do that we're going to introduce another dimensionality reduction method: PHATE (https://www.nature.com/articles/s41587-019-0336-3).

Unlike the other methods, PHATE is the first dimensionality reduction method purpose built for scRNA-seq data. It aims to identify and embed the biological transitions that occur across scRNA-seq data, making them ideal for TI methods to be performed on.

You can find a short tutorial to using PHATE at this link: https://github.com/KrishnaswamyLab/phateR. Create PHATE embeddings for your data (think carefully about what expression matrix you give to phate) and then put the embeddings into your integrated Seurat object and perform Slingshot using PCA, UMAP and PHATE space embeddings (make sure to save each of them into different objects/variables). Look at the trajectories drawn by Slingshot and see whether they accurately reflect the lifecycle stage transitions.

PHATE is quite an intense method for your computer, so it may take a while to run and cause your computer fan to ramp up (if you are running it locally).


```R

```

Now you've generated the three trajectories, see if you can identify any major differences between them and see how well they capture the transition through the lifecycle.

It can be hard to see where branching points are on a 2D plot. Fortunately, we can create a Slingshot object, which shows the progression of clusters across the trajectory.

Look at the order of clusters and see if it matches the biology of the system.


```R
#Give them all of this 
pb_sling_umap <- SlingshotDataSet(pb_sce_umap)
pb_sling_pca <- SlingshotDataSet(pb_sce_pca)
pb_sling_phate <- SlingshotDataSet(pb_sce_phate)


pb_sling_umap@lineages
pb_sling_pca@lineages
pb_sling_phate@lineages
```

Now let's look at the correlation between the different pseudotime values. Perhaps we can use corrgram to easily plot the correlation between the four vectors.

Identify where the biggest differences in correlation are between the different pseudotime calculations and try and think why they might be different.


```R
library(corrgram)

```

While pseudotime can help us identify branching points in datasets and new end states, it is mostly applied on systems where the overall progression is already known.

In both cases, it is followed by performing differential expression (DE); identifying genes that are differentially expressed across the pseudotime axis.

Monocle has it's own in-built DE method which we will not cover today, but it is worth a look if you are interested in pseudotime DE methods.

For today we are going to work with tradeSeq, a DE method which is made to be used with (but not exclusive to) Slingshot (https://www.nature.com/articles/s41467-020-14766-3).

We are going to follow the tradeSeq workflow (https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html) using the sce object where we have pseudotime generated by Slingshot from PHATE embeddings. Because of this we are going to start from the section with the title 'Differential Expression'.

TradeSeq works by fitting a model to each of the genes over pseudotime, essentially creating a line which reflects how the gene expression changes over pseudotime (i.e. goes down or goes up).

It takes a long time to fit a model for all the cells across all the genes, so what we are going to do is subsample our dataset down in terms of cells (only keep 33% of the cells) and remove genes whose total expression is less than 50 across all the cells.

Find the best knot value to use for the data and fit the GAMs to your object. After that, perform the Association test to see what genes are significantly different over the pseudotime axis and plot the top 5 genes (in terms of p-value) with the PlotSmoothers function for both the sexual and asexual trajectories.

<b> Note </b>

The tutorial for tradeSeq is quite old so some of the functions don't work exactly as it does in the tutorial. If something doesn't work, look up the documentation with the ? function.


```R
getwd()
```


'/home1/2254704l/pseudotime_workbook'


We would encourage you to put the genes into plasmoDB and see if their function lines up with them being important in the lifecycle of the parasite. You could even colour the cells in plotSmoother by the different stages and see which stage the genes are upregulated/downregulated in.

With that we've reached the end of this workbook. You may come away with the impression that pseudotime is a very accurate process, and even with use of different dimensionality reduction methods and TI methods they all arrive at a similar conclusion. Before you jump to that conclusion let's focus on two datasets we analysed:

1. The first is a simulated dataset which lacks any significant noise or bias where almost all of the genes contribute to the underlying process

2. The second is a real dataset, but it represents quite a homogenous process. Compare this with a sample from mammalian tissue, where the cells are all undergoing their own individual processes. By contrast the parasite is commited to other sexual or asexual lifecycle and most of the genes will thus be focused on this and this alone.

If you are interested, try applying Monocle and Slingshot to mammalian processes and you will see a much bigger fluctuation in results. 

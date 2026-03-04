# WGCNA - North vs South 

#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#49_Explore_our_WGCNA_results

##About this analysis 

#using a weighted gene co-expression anaylsis (WGCNA) to id co expressed gene modules --> a series of correlations to id sets of 
#genes that are expressed together in our dataset 

#output - wgcna gives groups of co-expressed genes 
#and eigengene x sample matrix (where the values for each eigengene represent the summarized expression for a group of co-expressed genes)



#set working dir
setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA")

## Install Pakages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("tximportData")
#BiocManager::install("DESeq2")
#BiocManager::install("vsn")
#BiocManager::install("genefilter")
#BiocManager::install("ComplexHeatmap")

## Load libraries 
library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library("vsn")
library("genefilter")
library("tidyverse")
library("magrittr") 
library("WGCNA") 
library("ggplot2")
library("ggforce")
library("ComplexHeatmap")

dir <- getwd()
list.files()

#### import and set up data ####

# set up metadata file (libraries_to_stages file)
samples <- read.table(file.path("libraries_to_stages_2.txt"), header=T)
samples$condition <- factor(rep(c("A","B","C","D","E","F","G","H","I","J","K","L","M","N"),each=3))
samples$time_point <- factor(samples$time_point, levels = c("T0","T14"))
samples$location <- factor(samples$location, levels = c("NS","ME","NH","FIELD","MA","SC","FL"))
samples$North_South <- factor(samples$North_South, levels = c("North","Field","South"))
rownames(samples) <- samples$sample_loc_rep
samples[,c("condition","time_point","location","sample_loc_rep","loc_time_point","North_South")]

# identify the files to import and add full header info for genes (from stella uk transcriptome)
files <- file.path(dir,"mapping",samples$sample_loc_rep,"quant.sf")
names(files) <- samples$sample_loc_rep
tx2gene <- read_csv(file.path("uk_transcriptome_headers.csv"))

## Import files
txi.salmon <- tximport(files, type="salmon", tx2gene=tx2gene)

# Create a DESeq2 dataset 
#ddsMF <- DESeqDataSetFromTximport(txi.salmon,colData = samples,design = ~ time_point + location) #multifactor analysis for DE 
dds <- DESeqDataSetFromTximport(txi.salmon,colData = samples,design = ~1) #not specifyng a model - not a DE analysis 



#### Perform DESeq2 normalization and transformation #### 

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)
#library(ggplot2)
# view the data
meanSdPlot(assay(dds_norm))
head(assay(dds_norm))
hist(assay(dds_norm))
plotPCA(dds_norm, intgroup="loc_time_point", ntop=2000)
plotPCA(dds_norm, intgroup="North_South",ntop=2000)
plotPCA(dds_norm, intgroup="location",ntop=2000)
plotPCA(dds_norm, intgroup="time_point",ntop=2000)




#### Format normalized data for WGCNA #### 

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data


#### Determine parameters for WGCNA ####

#To identify which genes are in the same modules, WGCNA first creates a weighted network to define which genes are near each other. 
#The measure of “adjacency” it uses is based on the correlation matrix, but requires the definition of a threshold value, which in 
#turn depends on a “power” parameter that defines the exponent used when transforming the correlation values. The choice of power parameter 
#will affect the number of modules identified, and the WGCNA modules provides the pickSoftThreshold() function to help identify good 
#choices for this parameter.

sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

# We want to plot some of this info to figure out our power soft-threshold - 
#first calculate a measure of the model fit (r^2) and make it a new variable
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

# Plot the model so we can decide on a soft-threshold for power
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

#Using the plot = I'm picking 10 --> above 0.80 power


#### Run WGCNA ####

#We will use the blockwiseModules() function to find gene co-expression modules in WGCNA, using 10 for the power argument like we determined above.

## ONLY NEED TO RUN ONCE - LOAD PREVIOUS

#bwnet <- blockwiseModules(normalized_counts,
#maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
#TOMType = "signed", # topological overlap matrix
#power = 10, # soft threshold for network construction
#numericLabels = TRUE, # Let's use numbers instead of colors for module labels
#randomSeed = 1234, # there's some randomness associated with this calculation
# so we should set a seed
#)

#The TOMtype argument specifies what kind of topological overlap matrix (TOM) should be used to make gene modules. 
#You can safely assume for most situations a signed network represents what you want – we want WGCNA to pay attention 
#to directionality. However if you suspect you may benefit from an unsigned network, where positive/negative is ignored 
#see this article to help you figure that out (Langfelder 2018).


#### Write main wgcna results object to file ####

# We will save our whole results object to an RDS file in case we want to return to our original WGCNA results.
#readr::write_rds(bwnet,file = file.path(dir, "stella_wgcna_results.RDS")) #only need to do once

#load original bwnet file <-- START FROM HERE*********
bwnet<-readRDS("stella_wgcna_results.RDS")



#### Explore WGCNA results ####

#In bwnet we have a data frame of eigengene module data for each sample in the MEs slot. 
#These represent the collapsed, combined, and normalized expression of the genes that make up each module.

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)


#### Which modules have biggest differences across treatment groups ####

#first see if samples are in order
all.equal(samples$sample_loc_rep, rownames(module_eigengenes))

# Create the design matrix from the `time_point` variable ## Maybe have this as two variables next time? 
#des_mat <- model.matrix(~ samples$time_point)
des_mat_2<- model.matrix(~ samples$North_South)

# Run linear model on each module. Limma wants our tests to be per row, so we need to transpose so the eigengenes are rows
# lmFit() needs a transposed version of the matrix
#fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit_2 <- limma::lmFit(t(module_eigengenes), design = des_mat_2)

# Apply empirical Bayes to smooth standard errors
#fit <- limma::eBayes(fit)
fit_2 <- limma::eBayes(fit_2)

# Apply multiple testing correction and obtain stats
#stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
#  tibble::rownames_to_column("module")

stats_df_2 <- limma::topTable(fit_2, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

# Look at the results - the most signficant results are at the top 
head(stats_df)
head(stats_df_2)

# Write out statsdf as a tsv 
#readr::write_tsv(stats_df,
#                 file = file.path(dir, "stats_wgcna_gene_to_module.tsv")
#)

readr::write_tsv(stats_df_2,
                 file = file.path(dir, "stats_2-northVSsouth_wgcna_gene_to_module.tsv")
)

## Now investigate the modules with the most differentially express across your treatment groups 




#### Lets make plots of significant modules ####

# First see what the module eigengene looks like between treatment groups

module_35_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample_loc_rep, time_point,location,loc_time_point,North_South),
                    by = c("sample_name" = "sample_loc_rep")
  )

module_10_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample_loc_rep, time_point,location,loc_time_point,North_South),
                    by = c("sample_name" = "sample_loc_rep")
  )

module_0_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample_loc_rep, time_point,location,loc_time_point,North_South),
                    by = c("sample_name" = "sample_loc_rep")
  )

module_4_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample_loc_rep, time_point,location,loc_time_point,North_South),
                    by = c("sample_name" = "sample_loc_rep")
  )

module_83_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample_loc_rep, time_point,location,loc_time_point,North_South),
                    by = c("sample_name" = "sample_loc_rep")
  )


module_2_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample_name") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(samples %>%
                      dplyr::select(sample_loc_rep, time_point,location,loc_time_point,North_South),
                    by = c("sample_name" = "sample_loc_rep")
  )




# Plot box plot - Module 35 
ggplot(
  module_35_df,
  aes(
    x = location,
    y = ME35,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


ggplot(
  module_35_df,
  aes(
    x = North_South,
    y = ME35,
    color = location
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


# Plot box plot - Module 10 
ggplot(
  module_10_df,
  aes(
    x = location,
    y = ME10,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


ggplot(
  module_10_df,
  aes(
    x = North_South,
    y = ME10,
    color = location
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


# Plot box plot - Module 0 
ggplot(
  module_0_df,
  aes(
    x = location,
    y = ME0,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


ggplot(
  module_0_df,
  aes(
    x = North_South,
    y = ME0,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


# Plot box plot - Module 4 
ggplot(
  module_4_df,
  aes(
    x = location,
    y = ME4,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


ggplot(
  module_4_df,
  aes(
    x = North_South,
    y = ME4,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


# Plot box plot - Module 83 
ggplot(
  module_83_df,
  aes(
    x = location,
    y = ME83,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


ggplot(
  module_83_df,
  aes(
    x = North_South,
    y = ME83,
    color = time_point
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()






#### What genes are a part of sifnificant Modules ####

# If you want to know which of your genes make up a modules, you can look at the $colors slot. This is a 
#named list which associates the genes with the module they are a part of. We can turn this into a data frame for handy use.

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

# Now we can find what genes are apart of current module *****(107 genes)******
gene_module_key_35 <- gene_module_key %>%
  dplyr::filter(module == "ME35")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME35-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(387 genes)******
gene_module_key_10 <- gene_module_key %>%
  dplyr::filter(module == "ME10")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME10-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(19609 genes)******
gene_module_key_0 <- gene_module_key %>%
  dplyr::filter(module == "ME0")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME0-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(710 genes)******
gene_module_key_4 <- gene_module_key %>%
  dplyr::filter(module == "ME4")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME4-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(47 genes)******
gene_module_key_83 <- gene_module_key %>%
  dplyr::filter(module == "ME83")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME83-NvS_wgcna_gene_to_module.tsv")
)


# Now we can find what genes are apart of current module *****(253 genes)******
gene_module_key_15 <- gene_module_key %>%
  dplyr::filter(module == "ME15")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME15-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(697 genes)******
gene_module_key_5 <- gene_module_key %>%
  dplyr::filter(module == "ME5")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME5-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(74 genes)******
gene_module_key_55 <- gene_module_key %>%
  dplyr::filter(module == "ME55")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME55-NvS_wgcna_gene_to_module.tsv")
)

# Now we can find what genes are apart of current module *****(117 genes)******
gene_module_key_32 <- gene_module_key %>%
  dplyr::filter(module == "ME32")

# Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(dir, "ME32-NvS_wgcna_gene_to_module.tsv")
)







#### Make custom heatmap ####

make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = samples,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its sample_loc_rep
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample_loc_rep")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(sample_loc_rep, time_point, location,North_South) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample_loc_rep") %>%
    # Arrange by patient and time point
    dplyr::arrange(time_point, location,North_South) %>%
    # Store sample
    tibble::column_to_rownames("sample_loc_rep")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    North_South = col_annot_df$North_South,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point #(New colors: "North" = '#6B0504', "Field" = "#D64045","South" = "#F19C79") #Original colors: "North" = "#f1a340", "Field" = "#2A9D8F","South" = "#998ec3"
    col = list(North_South = c("North" = '#6B0504', "Field" = "#D64045","South" = "#F19C79"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  
  # Create a color function based on standardized scale #(original colors: "#67a9cf", "#f7f7f7", "#ef8a62") #viridis: "#440154","#21918c","#fde725"
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#440154","#21918c","#fde725")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = T,
                                     show_column_names = T
  )
  
  # Return heatmap
  return(heatmap)
}

# Print out heatmaps 
mod_19_heatmap <- make_module_heatmap(module_name = "ME19") #tutorial example one 

#Top 10 in order of most significant P values
mod_35_heatmap <- make_module_heatmap(module_name = "ME35") # 1 **

mod_10_heatmap <- make_module_heatmap(module_name = "ME10") # 2
mod_0_heatmap <- make_module_heatmap(module_name = "ME0") # 3
mod_4_heatmap <- make_module_heatmap(module_name = "ME4") # 4
mod_83_heatmap <- make_module_heatmap(module_name = "ME83") # 5
mod_2_heatmap <- make_module_heatmap(module_name = "ME2") # 6
mod_77_heatmap <- make_module_heatmap(module_name = "ME77") # 7
mod_15_heatmap <- make_module_heatmap(module_name = "ME15") # 8 **
mod_82_heatmap <- make_module_heatmap(module_name = "ME82") # 9
mod_3_heatmap <- make_module_heatmap(module_name = "ME3") # 10
mod_7_heatmap <- make_module_heatmap(module_name = "ME7") # 11
mod_23_heatmap <- make_module_heatmap(module_name = "ME23") # 12
mod_5_heatmap <- make_module_heatmap(module_name = "ME5") # 13 **
mod_20_heatmap <- make_module_heatmap(module_name = "ME20") # 14
mod_55_heatmap <- make_module_heatmap(module_name = "ME55") # 15 **
mod_27_heatmap <- make_module_heatmap(module_name = "ME27") # 16
mod_32_heatmap <- make_module_heatmap(module_name = "ME32") # 17 **


dev.off()

#### Write out matrix ####

make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = samples,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its sample_loc_rep
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample_loc_rep")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(sample_loc_rep, time_point, location) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample_loc_rep") %>%
    # Arrange by patient and time point
    dplyr::arrange(time_point, location) %>%
    # Store sample
    tibble::column_to_rownames("sample_loc_rep")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    North_South = col_annot_df$North_South,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(North_South = c("North" = "#f1a340", "Field" = "#2A9D8F","South" = "#998ec3"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  #write out the norm matrix
  return(mod_mat)
  
  
}

#Top 5 in order of most significant P values - pull out your matrix from the heatmap code, turn it into a df, add row names
mod_35_matrix <- make_module_heatmap(module_name = "ME35") # 1
mod_35_df <- as.data.frame(mod_35_matrix)
mod_35_df$gene_name=c(rownames(mod_35_df))

mod_21_matrix <- make_module_heatmap(module_name = "ME21") # 2
mod_21_df <- as.data.frame(mod_21_matrix)
mod_21_df$gene_name=c(rownames(mod_21_df))

mod_60_matrix <- make_module_heatmap(module_name = "ME60") # 3
mod_60_df <- as.data.frame(mod_60_matrix)
mod_60_df$gene_name=c(rownames(mod_60_df))

mod_17_matrix <- make_module_heatmap(module_name = "ME17") # 4
mod_17_df <- as.data.frame(mod_17_matrix)
mod_17_df$gene_name=c(rownames(mod_17_df))

mod_96_matrix <- make_module_heatmap(module_name = "ME96") # 5
mod_96_df <- as.data.frame(mod_96_matrix)
mod_96_df$gene_name=c(rownames(mod_96_df))

# Write out the matricies of normalized counts used to make the heatmap
readr::write_tsv(mod_35_df,
                 file = file.path(dir, "ME35_NvS_norm_matrix.tsv"))

readr::write_tsv(mod_21_df,
                 file = file.path(dir, "ME21_norm_matrix.tsv"))

readr::write_tsv(mod_60_df,
                 file = file.path(dir, "ME60_norm_matrix.tsv"))

readr::write_tsv(mod_17_df,
                 file = file.path(dir, "ME17_norm_matrix.tsv"))

readr::write_tsv(mod_96_df,
                 file = file.path(dir, "ME96_norm_matrix.tsv"))


00\_MN\_Final\_Project
================
MN
4/19/2023

# Loading in all of the .broadpeak files for the 480 DBPs

``` r
# Loading in ChIP peak files of each protein and their replicates
broadpeak_path <- "/scratch/Shares/rinnclass/CLASS_2023/data/data/peaks"

# Importing peaks to get peak list (Grange)
peak_list <- import_peaks(consensus_file_path = broadpeak_path)

# Identifying how many peaks are in each file 
peak_numbers <- sapply(peak_list, length) %>% as.data.frame()

# Labeling column as peak_numbers
names(peak_numbers) <- c("peak_numbers")

# Separating DBP name and replicate number
peak_numbers <- peak_numbers %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")

# Saving peak_numbers data frame
write_csv(peak_numbers, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_numbers_df.csv")
```

# Creating consensus peaks

``` r
# Generating list of unique DBPs from peak_list (need for consensus_from_reduced function)
dbps <- unique(sapply(names(peak_list), function(x) {
  unlist(strsplit(x, "_"))[1]
  }))

# Using the consensus_from_reduced function, which finds peaks that overlap between replicates and combines them
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)

# Adding names to dbps
names(consensus_list) <- dbps 

# Making a data frame of the consensus peaks for each DBP
num_consensus_peaks_df <- sapply(consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(number_consensus_peaks = ".")

# For plotting
num_consensus_peaks_plotting <- sapply(consensus_list, length)

# Making a histogram of the number of consensus peaks to identify and remove outliers

Histogram_peak_numbers <- hist(num_consensus_peaks_plotting, breaks = 1000)
```

![](00_MN_Final_Project_files/figure-gfm/Consensus%20Peaks-1.png)<!-- -->

``` r
Histogram_peak_numbers <- hist(num_consensus_peaks_plotting, breaks = 1000, xlim = c(0,50000))
```

![](00_MN_Final_Project_files/figure-gfm/Consensus%20Peaks-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/Histogram_peak_numbers.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Merging the consensus peaks into the peak_numbers data frame
peak_numbers <- left_join(peak_numbers, num_consensus_peaks_df)
```

    ## Joining with `by = join_by(dbp)`

``` r
# Saving updates to data frame 
write_csv(peak_numbers, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_numbers_df.csv")
```

# Filtering consensus\_list to dbps with &gt; 1000 peaks

``` r
# Filtering to 1000 peaks 
filtered_consensus_list <- consensus_list[sapply(consensus_list, length) > 1000] 
# 486 dbps to 430, 56 dpbs had < 1000 consensus peaks

# Identifying which dbps were lost
lost_dbps <- names(consensus_list[sapply(consensus_list, length) < 1000]) %>% as.data.frame()

# Saving 
save(filtered_consensus_list, file = "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/filtered_consensus_list.RData")
write.table(lost_dbps, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/lost_dbps.csv")
```

# Exporting filtered\_consensus\_list

``` r
# Exporting filtered consensus peaks to results folder in .bed format

# Setting file path to export
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/mina8856"
filtered_consensus_path <- "CLASS_2023/CLASSES/06_final_going_global/results/filtered_consensus_peaks/"
exportpath <- file.path(basepath, filtered_consensus_path)

# Exporting each DBP consensus peak as .bed file
for(i in 1:length(filtered_consensus_list)) {rtracklayer::export(filtered_consensus_list[[i]], paste0(exportpath, names(filtered_consensus_list)[i], "_filtered_consensus.bed") )}
```

# Loading in the human genome from gencode

``` r
# Loading in the human genome
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

# gencode genes
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 

# Creating files for mRNA and lncRNA genes
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"]
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

# Combining the two
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]

# lncrna_mrna_promoters
lncrna_mrna_promoters <- promoters(mrna_lncrna_genes, upstream = 1000, downstream = 1000)

# Exporting files
rtracklayer::export(gencode_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/gencode_genes.gtf")
rtracklayer::export(mrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/mrna_genes.gtf")
rtracklayer::export(lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/lncrna_genes.gtf")
rtracklayer::export(mrna_lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/mrna_lncrna_genes.gtf")
rtracklayer::export(lncrna_mrna_promoters, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/lncrna_mrna_promoters.gtf")


# Importing files shortcuts for future use

# Importing annotation files 
mrna_lncrna_genes <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/mrna_lncrna_genes.gtf")
lncrna_mrna_promoters <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/gene_annotations/lncrna_mrna_promoters.gtf")

# For subsetting
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]
```

# Using human transcription factor (TF) annotations from a cell paper to identify which of my ChIP proteins are TF’s, their DNA binding domain, and their ensemble ID identification

``` r
# Adding column names to .bed consensus files 
# lapply (for loop) across consensus file list to add colnames
# Column names for .broadPeak are: chr, start, end, name, score, strand

# Importing filtered consensus list (430 dbps)
filtered_consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/filtered_consensus_peaks", 
                            pattern = "*.bed",
                            full.names = TRUE)
filtered_consensus_list <- lapply(filtered_consensus_file_list, rtracklayer::import)
names(filtered_consensus_list) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/filtered_consensus_peaks/|_filtered_consensus.bed","", filtered_consensus_file_list)

# Making simplified data frame
# Removing replicates and their peak numbers from "peak_numbers" df
# Adding the total_peak_length

num_peaks_df <- data.frame("dbp" = names(filtered_consensus_list),
                           "num_consensus_peaks" = sapply(filtered_consensus_list, length),
                           "total_peak_length" = sapply(filtered_consensus_list, function(x) sum(width(x))))


# Downloading TF annotations from cell paper
#redx1::read_excel to import
suppressWarnings({human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/mmc2.xlsx",
                                sheet = 2, skip = 1)})
```

    ## New names:
    ## • `` -> `...4`

``` r
# Renaming the 4th column of human_tfs to indicate if it is a TF (yes/no).
names(human_tfs)[4] <- "is_tf"

# Intersecting gene names that are in our ChIP data and has TF identity.
length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 407

``` r
# Filtering and grabbing the first 4 columns that match DBPs in num_peaks_df
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]


# Adding new column names
names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

# Adding human_tf's column to num_peaks_df
num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T)

# Exporting new data frame
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/num_peaks_df.csv")
```

# ChIP Proteins Data Analysis

# Determining how my proteins consensus peaks overlap mRNA and lncRNA genes and promoters

``` r
# Counting promoter overlaps using consensus peaks (columns are promoters, rows are DBPs)
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "counts")

# Summing each row for each DBP to identify the amount of promoter overlaps  
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

# Using gene_id objects to index into "lncrna" and "mrna" and separate them into two groups
# lncRNA promoter overlaps
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

# mrna promoter overlaps
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

# Exporting updates data frame
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/num_peaks_df.csv")
```

# Identifying amount of protein binding over gene bodies and comparing with promoters

``` r
# Creating genebody variable using count_peaks_per_feature function
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                filtered_consensus_list, 
                                                type = "counts")

# Extracting the overlaps the same way we did for promoters above
# All gene bodies
num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

# lncRNA gene bodies 
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])

# Exporting data frame
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/num_peaks_df.csv")
```

# Identifying how many of my proteins bind specific promoters

``` r
# Creating promoter peak occurence variable using count_peaks_per_feature function and type = occurrence
promoter_peak_occurrence_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, 
                                               type = "occurrence")

# Let's double check that all lncrna & mrna genes are accounted for:
stopifnot(all(colnames(promoter_peak_occurrence_matrix) == lncrna_mrna_promoters$gene_id))

# Saving
write.table(promoter_peak_occurrence_matrix, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/promoter_peak_occurrence_matrix.tsv")

# Now let's use the 'data.frame()' function. 
peak_occurrence_df <- data.frame("gene_id" = colnames(promoter_peak_occurrence_matrix),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurrence_matrix))

# Exporting data frame
write_csv(peak_occurrence_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_occurrence_dataframe.csv")
```

``` r
save(filtered_consensus_list, gencode_genes, lncrna_gene_ids, mrna_gene_ids, num_peaks_df, peak_occurrence_df, promoter_peak_occurrence_matrix, lncrna_mrna_promoters, mrna_lncrna_genes, file = "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData")

load("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData", verbose = T)
```

    ## Loading objects:
    ##   filtered_consensus_list
    ##   gencode_genes
    ##   lncrna_gene_ids
    ##   mrna_gene_ids
    ##   num_peaks_df
    ##   peak_occurrence_df
    ##   promoter_peak_occurrence_matrix
    ##   lncrna_mrna_promoters
    ##   mrna_lncrna_genes

# Making a promoter data frame that tells which dbps are bound

``` r
# dbps on promoters object
DBPs_on_promoter <- lncrna_mrna_promoters %>%
                    as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

# creating promoter dbps by pivot longer of promoter_peak_occurrence_matrix
promoter_dbps <- promoter_peak_occurrence_matrix %>%
  as.data.frame() %>%
  rownames_to_column("dbp") %>%
pivot_longer(2:ncol(.), names_to = "gene_id", values_to = "occurrence") %>%
  filter(occurrence == 1) %>%
  dplyr::select(-occurrence) %>%
  left_join(DBPs_on_promoter)
```

    ## Joining with `by = join_by(gene_id)`

``` r
# checking Firre promoter
firre_promoter <- promoter_dbps %>%
  filter(gene_name == "FIRRE")

# XIST promoter (should be off since male cells)
XIST_promoter <- promoter_dbps %>%
  filter(gene_name == "XIST")
# GAPDH
GAPDH_promoter <- promoter_dbps %>%
  filter(gene_name == "GAPDH")

# saving
promoter_dbps_df <- promoter_dbps %>% as.data.frame()
write.csv(promoter_dbps, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/promoter_dbps.csv")
```

# Plotting results for num\_peaks\_df (ggplot)

``` r
# Environment re-load
load("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData", verbose = T)
```

    ## Loading objects:
    ##   filtered_consensus_list
    ##   gencode_genes
    ##   lncrna_gene_ids
    ##   mrna_gene_ids
    ##   num_peaks_df
    ##   peak_occurrence_df
    ##   promoter_peak_occurrence_matrix
    ##   lncrna_mrna_promoters
    ##   mrna_lncrna_genes

``` r
# Number of peaks -vs- total peak length
consensus_peaks_vs_total_peak_length <-ggplot(num_peaks_df, aes(x = num_consensus_peaks, 
                                      y = total_peak_length)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
consensus_peaks_vs_total_peak_length
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/consensus_peaks_vs_total_peak_length.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Plotting number of consensus peaks on promoters

peak_num_vs_promoter_coverage.pdf <- ggplot(num_peaks_df,
                                     aes(x = num_consensus_peaks, y = peaks_overlapping_promoters)) +
                                    xlab("Peaks per DBP") +
                                    ylab("Number of peaks overlapping promoters") +
                                    ggtitle("Relationship Between Number of DBP Peaks and Promoter Overlaps")+
                                    geom_point() +
                                    geom_abline(slope = 1, linetype="dashed") +
                                    geom_smooth(method = "lm", se=FALSE, formula = 'y ~ x',
                                                color = "#a8404c") +
                                    stat_regline_equation(label.x = 35000, label.y = 18000) +
                                    ylim(0,60100) +
                                    xlim(0,60100)
peak_num_vs_promoter_coverage.pdf
```

    ## Warning: Removed 25 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 25 rows containing non-finite values
    ## (`stat_regline_equation()`).

    ## Warning: Removed 25 rows containing missing values (`geom_point()`).

![](00_MN_Final_Project_files/figure-gfm/ggplot-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/peak_num_vs_promoter_coverage.pdf")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 25 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 25 rows containing non-finite values
    ## (`stat_regline_equation()`).

    ## Warning: Removed 25 rows containing missing values (`geom_point()`).

``` r
# Plotting number of consensus peaks on gene bodies
peak_num_vs_gene_body_coverage <- ggplot(num_peaks_df,
                                 aes(x = num_consensus_peaks, y = peaks_overlapping_genebody)) +
                                xlab("Peaks per DBP") +
                                ylab("Number of peaks overlapping genes") +
                                ggtitle("Relationship Between Number of DBP Peaks and Gene Body Overlaps")+
                                geom_point() +
                                geom_abline(slope = 1, linetype="dashed") +
                                geom_smooth(method = "lm", se=F, formula = 'y ~ x',
                                            color = "#a8404c") +
                                stat_regline_equation(label.x = 35000, label.y = 18000) +
                                ylim(0,60100) +
                                xlim(0,60100)

peak_num_vs_gene_body_coverage
```

    ## Warning: Removed 25 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 25 rows containing non-finite values
    ## (`stat_regline_equation()`).

    ## Warning: Removed 25 rows containing missing values (`geom_point()`).

![](00_MN_Final_Project_files/figure-gfm/ggplot-3.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/peak_num_vs_gene_body_coverage.pdf")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 25 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 25 rows containing non-finite values
    ## (`stat_regline_equation()`).

    ## Warning: Removed 25 rows containing missing values (`geom_point()`).

``` r
# Relationship between promoter and genebody peaks 
promoter_peaks_vs_genebody_peaks <-ggplot(num_peaks_df, aes(x = peaks_overlapping_promoters, 
                                      y = peaks_overlapping_genebody)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
promoter_peaks_vs_genebody_peaks
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-4.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/promoter_peaks_vs_genebody_peaks.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Density plot of binding events

num_binding_events_per_promoter <- ggplot(peak_occurrence_df, aes(x = number_of_dbp)) +
                                  geom_density(alpha = 0.2, color = "#424242", fill = "#424242") +
                                    theme_paperwhite() +
                                    xlab(expression("Number of DBPs")) +
                                    ylab(expression("Density")) +
                                    ggtitle("Promoter binding events",
                                            subtitle = "mRNA and lncRNA genes")
          
num_binding_events_per_promoter
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-5.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/num_binding_events_per_promoter.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Relationship between mRNA and lncRNA promoters 
mRNA_promoter_peaks_vs_lncRNA_promoter_peaks <- ggplot(num_peaks_df, aes(x = peaks_overlapping_mrna_promoters, 
                                      y = peaks_overlapping_lncrna_genebody)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
mRNA_promoter_peaks_vs_lncRNA_promoter_peaks
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-6.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mRNA_promoter_peaks_vs_lncRNA_promoter_peaks.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Relationship between mRNA and lncRNA genebodies 
mRNA_genebody_peaks_vs_lncRNA_genebody_peaks <- ggplot(num_peaks_df, aes(x = peaks_overlapping_mrna_genebody, 
                                      y = peaks_overlapping_lncrna_promoters)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
mRNA_genebody_peaks_vs_lncRNA_genebody_peaks
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-7.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mRNA_genebody_peaks_vs_lncRNA_genebody_peaks.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Coloring plot based on if the protein is a TF or not.
consensus_peaks_vs_total_peak_length_Is_TF <- ggplot(num_peaks_df, aes(x = log2(num_consensus_peaks/1e3), 
                                      y = total_peak_length/1e6,
                                       color = tf == "Yes")) +
                                       geom_point()
consensus_peaks_vs_total_peak_length_Is_TF
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-8.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/consensus_peaks_vs_total_peak_length_Is_TF.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Making a histogram of the number of peaks for each of my proteins
# Histogram of peak#/DBP
Histogram_num_consensus_peaks <- ggplot(num_peaks_df, aes(x = num_consensus_peaks)) + 
                                geom_histogram(bins = 70)
Histogram_num_consensus_peaks
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-9.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/Histogram_num_consensus_peaks.pdf")
```

    ## Saving 7 x 5 in image

``` r
Histogram_TF_consensus_peaks_vs_count <- ggplot(num_peaks_df, aes(x = num_consensus_peaks, fill = tf)) +
                                       geom_histogram(bins = 30, position = "dodge")
Histogram_TF_consensus_peaks_vs_count
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-10.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/Histogram_TF_consensus_peaks_vs_count.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Now I want to facet this by the type of DNA binding domain my protein has.
DBP_facet_wrap <- ggplot(num_peaks_df, aes(x = num_consensus_peaks, 
                 y = total_peak_length )) +
                facet_wrap(dbd ~ .) +
                geom_point (shape = 'circle',
                color =  'red')
DBP_facet_wrap
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-11.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/DBP_facet_wrap.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Creating distance matrix & dendrogram

# creating distance matrix
peak_occurrence_dist <- dist(promoter_peak_occurrence_matrix, method = "binary")

# clustering distance values
bin_hier <- hclust(peak_occurrence_dist, method = "complete")

# Dendrogram of binding profiles by promoter (not binding profile - below)
promoter_overlap_dendrogram <- ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3,
                                 theme_dendro = TRUE) +
                                 coord_flip() +
                                 scale_y_continuous() +
                                 scale_x_continuous(position = "top") +
                                 scale_x_continuous(breaks = seq_along(bin_hier$labels[bin_hier$order]),
                                           labels = bin_hier$labels[bin_hier$order], position = "top",
                                           expand = c(0,0)) +
                                 theme(axis.text.x = element_text(angle = 90, hjust  = 1)) +
                                 theme(axis.text.y = element_text(angle = 0,hjust = 1)) +
                                 scale_y_reverse(expand = c(0.01, 0)) +
                                 theme(
                                   plot.background = element_blank(),
                                   panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                   panel.border = element_blank()
                                 )
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
promoter_overlap_dendrogram
```

![](00_MN_Final_Project_files/figure-gfm/ggplot-12.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/promoter_overlap_dendrogram.pdf")
```

    ## Saving 7 x 5 in image

\#RNA-seq \#Data is from the cell line HEPG2 that is fractionated by
cellular compartment(nuc, cyto, total etc)

``` r
# Reading in sample sheet from RNA-seq NF-Core pipeline
samplesheet <- read_rds("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/05_RNAseq/01_differential_expression/results/final_samplesheet.rds")
```

# Reading in TPM (transcripts per million) values from Salmon for our analyses

``` r
# Reading in salmon tpm
salmon_tpm <- read.csv("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/05_RNAseq/00_RNAseq_download_NF_core_pipeline/00_NF_CORE_RNAseq_Pipeline_run/results/salmon/salmon_merged_gene_tpm.csv")

# TPM table is in same order as samplesheet
tpm <- salmon_tpm %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  merge(samplesheet) %>%
  group_by(gene_id, condition) %>%
  summarize(tpm = mean(tpm, na.rm = T)) %>%
  pivot_wider(names_from = condition, values_from = tpm, names_prefix = "tpm_")
```

    ## `summarise()` has grouped output by 'gene_id'. You can override using the
    ## `.groups` argument.

# Determining how many DBPs are bound at each promoter

``` r
load("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData", verbose = T)
```

    ## Loading objects:
    ##   filtered_consensus_list
    ##   gencode_genes
    ##   lncrna_gene_ids
    ##   mrna_gene_ids
    ##   num_peaks_df
    ##   peak_occurrence_df
    ##   promoter_peak_occurrence_matrix
    ##   lncrna_mrna_promoters
    ##   mrna_lncrna_genes

``` r
# Merging peak_occurence_df (contains promoters where dbps are bound) and tpm
promoter_features_df <- merge(peak_occurrence_df, tpm)

# DBP density plot
DBP_binding_density_plot <- ggplot(promoter_features_df, aes(x = number_of_dbp)) +
                            geom_density() 
DBP_binding_density_plot 
```

![](00_MN_Final_Project_files/figure-gfm/Loading%20in%20promoter_features_df-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/DBP_binding_density_plot.pdf")
```

    ## Saving 7 x 5 in image

# Abundance of genes in each cellular fraction

``` r
# Making tpm_matrix
tpm_matrix <- tpm %>% 
  column_to_rownames("gene_id") %>%
  as.matrix()
tpm_scaled <- t(scale(t(tpm_matrix)))
tpm_scaled <- tpm_scaled[complete.cases(tpm_scaled),]

# Plotting
new.env()
```

    ## <environment: 0x24544c948>

``` r
pheatmap::pheatmap(tpm_scaled, show_rownames = FALSE)
```

![](00_MN_Final_Project_files/figure-gfm/TPM%20of%20genes%20in%20each%20fraction-1.png)<!-- -->

``` r
pdf("figures/heatmap_expression.pdf", height =49, width = 12)
```

# Plotting binding versus expression

``` r
# Plotting binding vs total RNA expression
binding_vs_expression_total_rna <- ggplot(promoter_features_df, 
                                  aes(y = log2(tpm_homo_sapiens_hepg2 + 0.001), x = number_of_dbp, color = gene_type)) + 
                                  geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_hepg2 > 0.001),
                                               shape = 17, alpha = 0.7) +
                                  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
                                  stat_cor() +
                                  geom_smooth(method = "lm") +
                                  scale_x_continuous(expand = c(0,0)) +
                                  scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
                                  ggtitle("Expression vs. promoter binding events") + 
                                  xlab(expression('Number of TFs')) +
                                  ylab(expression(log[2](TPM))) 
binding_vs_expression_total_rna
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](00_MN_Final_Project_files/figure-gfm/DBP%20promoter%20binding%20versus%20total%20RNA%20expression-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/binding_vs_expression_total_rna.pdf")
```

    ## Saving 7 x 5 in image
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
#linear trend with number of DBPS and expression levels, as expected
```

# Binding versus nuclear expression

``` r
# Now let's make a similar plot for nuclear RNA abundance versus #DBPs bound to their promoter
nuclear_expression_vs_promoter_binding <- ggplot(promoter_features_df, 
                                          aes(y = log2(tpm_homo_sapiens_nuclear_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
                                          geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_nuclear_fraction > 0.001),
                                                     shape = 17, alpha = 0.7) +
                                          geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
                                          stat_cor() +
                                          scale_x_continuous(expand = c(0,0)) +
                                          scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
                                          ggtitle("Nuclear Expression vs. promoter binding events") + 
                                          xlab(expression('Number of TFs')) +
                                          ylab(expression(log[2](TPM))) 
nuclear_expression_vs_promoter_binding
```

![](00_MN_Final_Project_files/figure-gfm/Binding%20versus%20nuclear%20expression-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/nuclear_expression_vs_promoter_binding.pdf")
```

    ## Saving 7 x 5 in image

# Binding versus cytoplasmic expression

``` r
# Same thing just seeing if there is a difference of cyto RNAs versus DBPs on promoter
cytoplasmic_expression_vs_promoter_binding <- ggplot(promoter_features_df, 
                                            aes(y = log2(tpm_homo_sapiens_cytosolic_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
                                            geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_cytosolic_fraction > 0.001),
                                                       shape = 17, alpha = 0.7) +
                                            geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
                                            stat_cor() +
                                            scale_x_continuous(expand = c(0,0)) +
                                            scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
                                            ggtitle("Cytoplasmic Expression vs. promoter binding events") + 
                                            xlab(expression('Number of TFs')) +
                                            ylab(expression(log[2](TPM))) 
cytoplasmic_expression_vs_promoter_binding
```

![](00_MN_Final_Project_files/figure-gfm/Binding%20versus%20cytoplasmic%20expression-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/cytoplasmic_expression_vs_promoter_binding.pdf")
```

    ## Saving 7 x 5 in image

# lncRNA versus mRNA expression in total RNA

``` r
# Determining lncRNA and mRNA expression levels in total RNA}
mrna_lncrna_tpm_total_rna <- ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_hepg2 + 0.01), color = gene_type))+
                              geom_density()
mrna_lncrna_tpm_total_rna
```

![](00_MN_Final_Project_files/figure-gfm/Determining%20lncRNA%20and%20mRNA%20expression%20levels%20in%20total%20RNA-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mrna_lncrna_tpm_total_rna.pdf")
```

    ## Saving 7 x 5 in image

``` r
#  Determining lncRNA and mRNA expression levels in nuclear RNA
mrna_lncrna_tpm_nuclear <- ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_nuclear_fraction + 0.01), color = gene_type))+
                            geom_density()
mrna_lncrna_tpm_nuclear
```

![](00_MN_Final_Project_files/figure-gfm/Determining%20lncRNA%20and%20mRNA%20expression%20levels%20in%20total%20RNA-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mrna_lncrna_tpm_nuclear.pdf")
```

    ## Saving 7 x 5 in image

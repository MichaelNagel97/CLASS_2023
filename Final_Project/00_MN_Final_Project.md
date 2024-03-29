00\_MN\_Final\_Project
================
MN
4/19/2023

# Loading in all of the .broadpeak files for the 486 DBPs

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

# Determining how the DBPs consensus peaks overlap mRNA and lncRNA genes and promoters

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

# Identifying what proteins bind at specific promoters

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

# Making a promoter data frame that tells which dbps are bound

``` r
# dbps on promoters object
DBPs_on_promoter <- lncrna_mrna_promoters %>%
                    as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

# Creating promoter dbps by pivot longer of promoter_peak_occurrence_matrix
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
# GAPDH
GAPDH_promoter <- promoter_dbps %>%
  filter(gene_name == "GAPDH")

# saving
promoter_dbps_df <- promoter_dbps %>% as.data.frame()
write.csv(promoter_dbps, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/promoter_dbps.csv")
```

# Using profile\_tss for all 430 DBPs

``` r
# establishing DF
metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())

# for loop to populate DF 
for(i in 1:length(filtered_consensus_list)) {
  # print(names(filtered_consensus_list)[[i]])
  tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters)
  tmp_df$dbp <- names(filtered_consensus_list)[[i]]
  metaplot_df <- bind_rows(metaplot_df, tmp_df)
}

# saving
write_rds(metaplot_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/metaplot_df_final.rds")
```

``` r
save(filtered_consensus_list, gencode_genes, lncrna_gene_ids, mrna_gene_ids, num_peaks_df, peak_occurrence_df, promoter_dbps_df, metaplot_df, promoter_peak_occurrence_matrix, lncrna_mrna_promoters, mrna_lncrna_genes, file = "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData")

load("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData", verbose = T)
```

    ## Loading objects:
    ##   filtered_consensus_list
    ##   gencode_genes
    ##   lncrna_gene_ids
    ##   mrna_gene_ids
    ##   num_peaks_df
    ##   peak_occurrence_df
    ##   promoter_dbps_df
    ##   metaplot_df
    ##   promoter_peak_occurrence_matrix
    ##   lncrna_mrna_promoters
    ##   mrna_lncrna_genes

``` r
# Main data frames for plotting: num_peaks_df, peak_occurrence_df, promoter_dbps_df, metaplot_df
```

# Plotting results using ggplot

# Plotting number of consensus peaks on promoters and gene bodies

``` r
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

![](00_MN_Final_Project_files/figure-gfm/consensus%20peaks-1.png)<!-- -->

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

![](00_MN_Final_Project_files/figure-gfm/consensus%20peaks-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/peak_num_vs_gene_body_coverage.pdf")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 25 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 25 rows containing non-finite values
    ## (`stat_regline_equation()`).

    ## Warning: Removed 25 rows containing missing values (`geom_point()`).

``` r
# Peak number vs genome coverage (total_peak_length) 
```

# Results

# 1) There is a saturation point of peaks overlapping the promoter (approx. 20,000)

# 2) There is a linear relationship between DBP peaks and peaks overlapping the genebody, indicating that DBPs bind the genebody more than elsewhere in the genome

``` r
# Number of peaks -vs- total peak length
num_consensus_peaks_vs_total_peak_length <- ggplot(num_peaks_df, aes(x = num_consensus_peaks, 
                                      y = total_peak_length)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
num_consensus_peaks_vs_total_peak_length
```

![](00_MN_Final_Project_files/figure-gfm/Peak%20number%20vs%20peak%20length%20(genome%20coverage)-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/num_consensus_peaks_vs_total_peak_length.pdf")
```

    ## Saving 7 x 5 in image

# Results

# 1) Generally, the more consensus peaks a given DBP has, the greater the genome coverage

# 2) Most DBPs have under 25,000 consensus peaks, and have relatively specific binding

# 3) The top 3 DBPs with the highest number of consensus peaks are all unknown proteins (ZNF280B, ZNF407, and ZBTB38)

# Distribution of promoter overlaps versus gene-bodies

``` r
# Relationship between promoter and genebody peaks 
promoter_peaks_vs_genebody_peaks <-ggplot(num_peaks_df, aes(x = peaks_overlapping_promoters, 
                                      y = peaks_overlapping_genebody)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
promoter_peaks_vs_genebody_peaks
```

![](00_MN_Final_Project_files/figure-gfm/Promoter%20vs%20gene%20body%20peaks-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/promoter_peaks_vs_genebody_peaks.pdf")
```

    ## Saving 7 x 5 in image

``` r
promoter_peaks_vs_genebody_peaks_isTF <-ggplot(num_peaks_df, aes(x = peaks_overlapping_promoters, 
                                      y = peaks_overlapping_genebody,
                                      color = tf == "Yes")) +
                                       geom_point()
promoter_peaks_vs_genebody_peaks_isTF
```

![](00_MN_Final_Project_files/figure-gfm/Promoter%20vs%20gene%20body%20peaks-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/promoter_peaks_vs_genebody_peaks_isTF.pdf")
```

    ## Saving 7 x 5 in image

# Results

# 1) Plotting promoter vs genebody overlaps resembles an exponential curve

# 2) This graph indicates that the DBPs that bind to a vast number of promoters also bind to a vast number of gene body sequences, indicating that they are permissive binders

# 3) An interesting outlier in the graph is H3K4me2, with 22865 peaks overlapping promoters, and 49225 peaks overlapping the gene body. Indicating that this histone mark is mainly associated with the 5’ end of a gene.

# 4) Many of the outliers are not annotated as TF’s (False or True), and thus are likley histone marks or phosphorylation marks on the CTD of POLII

# Density plot of promoter binding events

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

![](00_MN_Final_Project_files/figure-gfm/DBP%20binding%20events%20per%20promoter-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/num_binding_events_per_promoter.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Making superbinders data frame, superbinders = 200-400 dbps bound to promoter

superbinders_peak_occurrence_df <- filter(peak_occurrence_df,number_of_dbp %in% c(200:400))
length(superbinders_peak_occurrence_df$gene_type) #11712
```

    ## [1] 11712

# Results

# 1) Most promoters have &lt; 100 DBPs bound

# 2) A subset of promoters have, on average, 300 proteins bound. These promoters are termed superbinders

# 3) The majority of super binders are present on protein coding genes, and not lncRNA, which may indicate that these genes require more intricate regulation of their expression

# Looking at gene ontology of data set

``` r
# Looking at all of the DBPs
# Coloring plot based on if the protein is a TF or not.
consensus_peaks_vs_total_peak_length_Is_TF <- ggplot(num_peaks_df, aes(x = log2(num_consensus_peaks/1e3), 
                                      y = total_peak_length/1e6,
                                       color = tf == "Yes")) +
                                       geom_point()
consensus_peaks_vs_total_peak_length_Is_TF
```

![](00_MN_Final_Project_files/figure-gfm/Gene%20Ontology%20All%20DBPs-1.png)<!-- -->

``` r
length(which(num_peaks_df$tf=="Yes")) #344 (80%)
```

    ## [1] 344

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/consensus_peaks_vs_total_peak_length_Is_TF.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Looking at number of consensus peaks for TF's vs not
Histogram_TF_consensus_peaks_vs_count <- ggplot(num_peaks_df, aes(x = num_consensus_peaks, fill = tf)) +
                                       geom_histogram(bins = 30, position = "dodge")
Histogram_TF_consensus_peaks_vs_count
```

![](00_MN_Final_Project_files/figure-gfm/Gene%20Ontology%20All%20DBPs-2.png)<!-- -->

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

![](00_MN_Final_Project_files/figure-gfm/Gene%20Ontology%20All%20DBPs-3.png)<!-- -->

``` r
length(which(num_peaks_df$dbd=="C2H2 ZF")) #169 (39%)
```

    ## [1] 169

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/DBP_facet_wrap.pdf")
```

    ## Saving 7 x 5 in image

# Results

# 1) Most DBPs transcription factor or not, contain less than 50,000 consensus peaks

# 2) This data set contains many more TF’s than not (80%)

# 3) The main DNA binding domain is a C2H2 ZF (39%)

# Looking at what gene types are most prevalent on three categories of promoters: nonbinders (0 DBPs), normalbinders (1-100 DBPs), and superbinders (200-400 DBPs)

``` r
length(superbinders_peak_occurrence_df$gene_type) 
```

    ## [1] 11712

``` r
# Looking at only the promoters which are superbinders

num_binding_events_superbinder_promoters <- ggplot(superbinders_peak_occurrence_df, aes(x = number_of_dbp)) +
                                  geom_density(alpha = 0.2, color = "#424242", fill = "#424242") +
                                    theme_paperwhite() +
                                    xlab(expression("Number of DBPs")) +
                                    ylab(expression("Density")) +
                                    ggtitle("Superbinders",
                                            subtitle = "mRNA and lncRNA genes")
num_binding_events_superbinder_promoters
```

![](00_MN_Final_Project_files/figure-gfm/Gene%20Ontology%20nonbinders,%20normalbinders,%20and%20superbinders-1.png)<!-- -->

``` r
superbinder_promoters_genetype <- ggplot(superbinders_peak_occurrence_df, aes(x = number_of_dbp, fill = gene_type)) +
                                       geom_histogram(bins = 50, position = "dodge")
superbinder_promoters_genetype
```

![](00_MN_Final_Project_files/figure-gfm/Gene%20Ontology%20nonbinders,%20normalbinders,%20and%20superbinders-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/superbinder_promoters_genetype.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Identifying how many protein coding genes are in nonbinders, normalbinders, and superbinders

# Making superbinders data frame, superbinders = 200-400 dbps bound to promoter
superbinders_peak_occurrence_df <- filter(peak_occurrence_df,number_of_dbp %in% c(200:400))
length(superbinders_peak_occurrence_df$gene_type) #11712
```

    ## [1] 11712

``` r
# Promoters that bound no DBPs
nonbinders_peak_occurence_df <- filter(peak_occurrence_df,number_of_dbp %in% (0))
length(nonbinders_peak_occurence_df$gene_type) # 9448
```

    ## [1] 9448

``` r
length(which(nonbinders_peak_occurence_df$gene_type=="protein_coding")) #3402 (36%)
```

    ## [1] 3402

``` r
# Promoters that contain < 100 DBPs
normalbinders_peak_occurrence_df <- filter(peak_occurrence_df,number_of_dbp %in% (1:100))
length(normalbinders_peak_occurrence_df$gene_type) # 12932 (47%)
```

    ## [1] 12932

``` r
length(which(normalbinders_peak_occurrence_df$gene_type=="protein_coding")) #5572 (71%)
```

    ## [1] 5572

``` r
length(superbinders_peak_occurrence_df$gene_type) # 11712 (42%)
```

    ## [1] 11712

``` r
length(which(superbinders_peak_occurrence_df$gene_type=="protein_coding")) #9192 (78%)
```

    ## [1] 9192

``` r
# Promoters that bound a DBP
Allbinders_peak_occurence_df <- filter(peak_occurrence_df,number_of_dbp %in% (1:405))
length(Allbinders_peak_occurence_df$gene_type) # 27366
```

    ## [1] 27366

``` r
length(which(Allbinders_peak_occurence_df$gene_type=="protein_coding")) #16563 (60%)
```

    ## [1] 16563

# Results

# 9448 promoters bound no DBPs (26%), 36% of which were protein coding

# Most promoters that did not bind a DBP are lncRNA promoters

# 27366 promoters bound DBPs (74%), 60% of which were protein coding

# 47% of promoters fall into the “normalbinder” DBP range, 43% of which are protein coding

# 43% of promoters fall into the “superbinder” DBP range, 78% of which are protein coding

# 10% of promters fall into the 100-200 DBP range

# Most superbinder promoters are protein coding

# mRNA and lncRNA promoter and genebody peak relationship

``` r
# Relationship between mRNA and lncRNA promoters 
mRNA_promoter_peaks_vs_lncRNA_promoter_peaks <- ggplot(num_peaks_df, aes(x = peaks_overlapping_mrna_promoters, 
                                      y = peaks_overlapping_lncrna_promoters)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
mRNA_promoter_peaks_vs_lncRNA_promoter_peaks
```

![](00_MN_Final_Project_files/figure-gfm/mRNA%20vs%20lncRNA%20peaks-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mRNA_promoter_peaks_vs_lncRNA_promoter_peaks.pdf")
```

    ## Saving 7 x 5 in image

``` r
# Relationship between mRNA and lncRNA genebodies 
mRNA_genebody_peaks_vs_lncRNA_genebody_peaks <- ggplot(num_peaks_df, aes(x = peaks_overlapping_mrna_genebody, 
                                      y = peaks_overlapping_lncrna_genebody)) +
                                      geom_point(shape = 'circle',
                                      color = 'red')
mRNA_genebody_peaks_vs_lncRNA_genebody_peaks
```

![](00_MN_Final_Project_files/figure-gfm/mRNA%20vs%20lncRNA%20peaks-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mRNA_genebody_peaks_vs_lncRNA_genebody_peaks.pdf")
```

    ## Saving 7 x 5 in image

# Results

# 1) Interesting inflection point where the number of peaks overlapping mRNA and lncRNA promoters becomes exponential. This could indicate a subset of DBPS with permissive binding at promoters

# 2) Linear relationship (mostly, a few outliers) between peaks overlapping mRNA and lncRNA gene bodies.

# Clustering

``` r
# Choosing random promoters from normal promoters (< 100 dbps) and super binding promoters (> 200 dbps)

# CAMP
CAMP_promoter <- promoter_dbps_df %>%
  filter(gene_name == "CAMP")
nrow(CAMP_promoter) # 40 DBPs bound, with 11 unknown DBPs bound containing a ZNF domain (ratio: 28%)
```

    ## [1] 40

``` r
# FAM72A
FAM72A_promoter <- promoter_dbps_df %>%
  filter(gene_name == "FAM72A")
nrow(FAM72A_promoter) # 100 DBPs bound, with 28 unknown DBPs bound containing a ZNF domain (ratio: 28%)
```

    ## [1] 100

``` r
# FOXA3
FOXA3_promoter <- promoter_dbps_df %>%
  filter(gene_name == "FOXA3")
nrow(FOXA3_promoter) # 318 DBPs bound, with 83 unknown DBPs bound containing a ZNF domain (ratio: 26%)
```

    ## [1] 318

``` r
# TBP
TBP_promoter <- promoter_dbps_df %>%
  filter(gene_name == "TBP")
nrow(TBP_promoter) # 343 DBPs bound, with 86 unknown DBPs bound containing a ZNF domain (ratio: 25%)
```

    ## [1] 343

``` r
# GAPDH
GAPDH_promoter <- promoter_dbps_df %>%
  filter(gene_name == "GAPDH")
nrow(GAPDH_promoter) # 360 DBPs bound, with 93 unknown DBPs bound containing a ZNF domain (ratio: 26%)
```

    ## [1] 360

``` r
# ZNF335
ZNF335_promoter <- promoter_dbps_df %>%
  filter(gene_name == "ZNF335")
nrow(ZNF335_promoter) # 371 DBPs bound, with 97 unknown DBPs bound containing a ZNF domain (ratio: 26%)
```

    ## [1] 371

``` r
# HSPE1
HSPE1_promoter <- promoter_dbps_df %>%
  filter(gene_name == "HSPE1")
nrow(HSPE1_promoter) # 399 DBPs bound, with 104 unknown DBPs bound containing a ZNF domain (ratio: 26%)
```

    ## [1] 399

``` r
# Filtering for ZNF proteins on each promoter, removing gene ID column
GAPDH_promoter_ZNF <- GAPDH_promoter %>%  filter(grepl("ZNF*",dbp))
GAPDH_promoter_ZNF <- GAPDH_promoter_ZNF[ -c(2)]
CAMP_promoter_ZNF <- CAMP_promoter %>%  filter(grepl("ZNF*",dbp))
CAMP_promoter_ZNF <- CAMP_promoter_ZNF[ -c(2)]
FAM72A_promoter_ZNF <- FAM72A_promoter %>%  filter(grepl("ZNF*",dbp))
FAM72A_promoter_ZNF <- FAM72A_promoter_ZNF[ -c(2)]
FOXA3_promoter_ZNF <- FOXA3_promoter %>%  filter(grepl("ZNF*",dbp))
FOXA3_promoter_ZNF <- FOXA3_promoter_ZNF[ -c(2)]
TBP_promoter_ZNF <- TBP_promoter %>%  filter(grepl("ZNF*",dbp))
TBP_promoter_ZNF <- TBP_promoter_ZNF[ -c(2)]
ZNF335_promoter_ZNF <- ZNF335_promoter %>%  filter(grepl("ZNF*",dbp))
ZNF335_promoter_ZNF <-ZNF335_promoter_ZNF[ -c(2)]
HSPE1_promoter_ZNF <- HSPE1_promoter %>%  filter(grepl("ZNF*",dbp))
HSPE1_promoter_ZNF <- HSPE1_promoter_ZNF[ -c(2)]


#put all data frames into list
ZNF_promoter_list <- list(CAMP_promoter_ZNF, FAM72A_promoter_ZNF, FOXA3_promoter_ZNF, TBP_promoter_ZNF, GAPDH_promoter_ZNF, ZNF335_promoter_ZNF, HSPE1_promoter_ZNF)

#merge all data frames in list
ZNF_dbps <- ZNF_promoter_list %>% reduce(full_join, by='dbp')
```

# Results

# 1) Regardless of if the promoter is a normal or super binder they bind a similar ratio of these unknown ZNF domain containing DBPs

# 2) It appears there are a subset of six ZNF’s that bind the majority of promoters (caveat: this data is only analyzing 7 different promoters)

# 3) It is likely that these six ZNF’s are ubiquitous DBPs needed for basal transcription to occur

# Metaplots

``` r
# Creating distance matrix of binding profile correlations
metaplot_filtered_matrix <- metaplot_df %>% 
  pivot_wider(names_from = x, values_from = dens) %>%
  column_to_rownames("dbp") %>%
  as.matrix()
mm_scaled <- t(scale(t(metaplot_filtered_matrix)))
metaplot_hclust <- hclust(dist(mm_scaled), method = "complete")

# plotting relationship between binding profiles
plot(metaplot_hclust)
```

![](00_MN_Final_Project_files/figure-gfm/Metaplots-1.png)<!-- -->

``` r
pdf("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/tss_profile_dendrogram.pdf", height = 10, width = 27)
par(cex=0.3)
plot(metaplot_hclust)
dev.off()
```

    ## png 
    ##   2

``` r
# Cutting tree to make clusters

clusters <- cutree(metaplot_hclust, k=30)
table(clusters)
```

    ## clusters
    ##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
    ##  1 93 31 50 12 44 56 79  6 28  7  1  1  1  1  1  1  2  1  1  1  1  4  1  1  1 
    ## 27 28 29 30 
    ##  1  1  1  1

# Results

# 1) The vast majority of dbps clutser relatively close together

``` r
# creating promoters just in case:
lncrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "lncRNA"]
mrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type == "protein_coding"]
```

# Metaplots for each DBP by lncRNA and mRNA promoters

``` r
# lncrna DF
lncrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
for(i in 1:length(filtered_consensus_list)) {
  tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters = lncrna_promoters)
  tmp_df$dbp <- names(filtered_consensus_list)[[i]]
  lncrna_metaplot_df <- bind_rows(lncrna_metaplot_df, tmp_df)
  }

# saving
write_rds(lncrna_metaplot_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/lncRNA_metaplot_df_final.rds")

# mRNAs DF
mrna_metaplot_df <- data.frame(x = integer(), dens = numeric(), dbp = character())
for(i in 1:length(filtered_consensus_list)) {
  tmp_df <- profile_tss(filtered_consensus_list[[i]], lncrna_mrna_promoters = mrna_promoters)
  tmp_df$dbp <- names(filtered_consensus_list)[[i]]
  mrna_metaplot_df <- bind_rows(mrna_metaplot_df, tmp_df)
  }

# Saving 
write_rds(mrna_metaplot_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/mrna_metaplot_df_final.rds")

# Adding the information of gene type
mrna_metaplot_df$gene_type <- "mRNA"
lncrna_metaplot_df$gene_type <- "lncRNA"
combined_metaplot_profile <- bind_rows(mrna_metaplot_df, lncrna_metaplot_df)

# Saving
write_rds(mrna_metaplot_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/metaplot_df_final.rds")
```

``` r
save(filtered_consensus_list, gencode_genes, lncrna_gene_ids, mrna_gene_ids, num_peaks_df, peak_occurrence_df, promoter_dbps_df, metaplot_df, mrna_metaplot_df, lncrna_metaplot_df, mrna_metaplot_df,promoter_peak_occurrence_matrix, lncrna_mrna_promoters, mrna_lncrna_genes, file = "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData")

load("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/results/peak_features.RData", verbose = T)
```

    ## Loading objects:
    ##   filtered_consensus_list
    ##   gencode_genes
    ##   lncrna_gene_ids
    ##   mrna_gene_ids
    ##   num_peaks_df
    ##   peak_occurrence_df
    ##   promoter_dbps_df
    ##   metaplot_df
    ##   mrna_metaplot_df
    ##   lncrna_metaplot_df
    ##   mrna_metaplot_df
    ##   promoter_peak_occurrence_matrix
    ##   lncrna_mrna_promoters
    ##   mrna_lncrna_genes

``` r
# Main data frames for plotting: num_peaks_df, peak_occurrence_df, promoter_dbps_df, metaplot_df, mrna_metaplot_df
```

# Profile plots for lncRNA and mRNA promoters seperated for each DBP

``` r
combined_metaplot_profile <- bind_rows(mrna_metaplot_df, lncrna_metaplot_df)

ggplot(combined_metaplot_profile, 
       aes(x = x, y = dens, color = gene_type )) +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_line(linewidth = 1.5) + 
  facet_wrap(dbp ~ ., scales = "free_y") +
  ggtitle("Promoter Metaplot") + 
  scale_x_continuous(breaks = c(-1000, 0, 1000),
                     labels = c("-1kb", "TSS", "+1kb"),
                     name = "") + 
  ylab("Peak frequency") +
 scale_color_manual(values = c("#424242","#a8404c"))
```

![](00_MN_Final_Project_files/figure-gfm/lncRNA%20and%20mRNA%20profile%20plots%20for%20each%20dbp-1.png)<!-- -->

``` r
# saving
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/mega_meta_plot_lncRNA-mRNA.pdf", width = 49, height = 12)
```

# Results

# Could not get this plot to work for the life of me

# RNA-seq

# Data is from the cell line HEPG2 that is fractionated by cellular compartment (nuclear, cytosolic, total etc)

``` r
# Reading in sample sheet from RNA-seq NF-Core pipeline
samplesheet <- read_rds("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/05_RNAseq/01_differential_expression/results/final_samplesheet.rds")

# Reading in TPM (transcripts per million) values from Salmon for our analyses
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
    ##   promoter_dbps_df
    ##   metaplot_df
    ##   mrna_metaplot_df
    ##   lncrna_metaplot_df
    ##   mrna_metaplot_df
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

# Results

# 1) Bimodal binding distribution at promoters

# 2) Majority of promoters bind &lt;100 dbps, another set of promoters binds 200-400 dbps

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

    ## <environment: 0x24d0b34a8>

``` r
pheatmap::pheatmap(tpm_scaled, show_rownames = FALSE)
```

![](00_MN_Final_Project_files/figure-gfm/TPM%20of%20genes%20in%20each%20fraction-1.png)<!-- -->

``` r
pdf("figures/heatmap_expression.pdf", height =49, width = 12)
```

# Result

# 1) Most RNAs are abundant in the nucleus

# 2) Some RNAs expressed in total that are not in other fractions

# 3) Low abundance of RNAs in the cytosol, as expected, most of it is likely tRNA

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

# Result

# 1) Roughly linear trend when plotting number of dbps bound vs expression levels

# 2) This trend is stronger for protein coding genes as compared to lncRNA

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

# Result

# 1) Trend is nearly identical to total RNA expression vs promoter binding

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

# Results

# 1) Lower abundance when compared to nuclear RNA expression vs promoter binding

# 2) Same linear trend as the first two plots

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

# Results

# 1) More mRNA is expressed than lncRNA

# 2) When looking at just nuclear RNA, the expression levels between mRNA and lncRNA get closer

``` r
superbinders_peak_occurrence_df <- filter(peak_occurrence_df,number_of_dbp %in% c(200:400))
length(superbinders_peak_occurrence_df$gene_type)
```

    ## [1] 11712

``` r
superbinders_tpm_df <- filter(promoter_features_df,number_of_dbp %in% c(200:400))
length(superbinders_tpm_df$gene_name)
```

    ## [1] 11712

``` r
superbinders_binding_density_plot <- ggplot(superbinders_tpm_df, aes(x = number_of_dbp)) +
                            geom_density() 
superbinders_binding_density_plot # Looks correct
```

![](00_MN_Final_Project_files/figure-gfm/High%20binding%20promoters-1.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/superbinders_binding_density_plot.pdf")
```

    ## Saving 7 x 5 in image

``` r
superbinders_vs_expression_total_rna <- ggplot(superbinders_tpm_df, 
                                  aes(y = log2(tpm_homo_sapiens_hepg2 + 0.001), x = number_of_dbp, color = gene_type)) + 
                                  geom_point(data = superbinders_tpm_df %>% filter(tpm_homo_sapiens_hepg2 > 0.001),
                                               shape = 17, alpha = 0.7) +
                                  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
                                  stat_cor() +
                                  geom_smooth(method = "lm") +
                                  scale_x_continuous(expand = c(0,0)) +
                                  scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
                                  ggtitle("Expression vs. promoter binding events") + 
                                  xlab(expression('Number of TFs')) +
                                  ylab(expression(log[2](TPM))) 
superbinders_vs_expression_total_rna
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](00_MN_Final_Project_files/figure-gfm/High%20binding%20promoters-2.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/superbinders_vs_expression_total_rna.pdf")
```

    ## Saving 7 x 5 in image
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
nuclear_expression_vs_superbinders <- ggplot(superbinders_tpm_df, 
                                          aes(y = log2(tpm_homo_sapiens_nuclear_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
                                          geom_point(data = superbinders_tpm_df %>% filter(tpm_homo_sapiens_nuclear_fraction > 0.001),
                                                     shape = 17, alpha = 0.7) +
                                          geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
                                          stat_cor() +
                                          scale_x_continuous(expand = c(0,0)) +
                                          scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
                                          ggtitle("Nuclear Expression vs. promoter binding events") + 
                                          xlab(expression('Number of TFs')) +
                                          ylab(expression(log[2](TPM))) 
nuclear_expression_vs_superbinders
```

![](00_MN_Final_Project_files/figure-gfm/High%20binding%20promoters-3.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/nuclear_expression_vs_superbinders.pdf")
```

    ## Saving 7 x 5 in image

``` r
cytoplasmic_expression_vs_superbinders <- ggplot(superbinders_tpm_df, 
                                            aes(y = log2(tpm_homo_sapiens_cytosolic_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
                                            geom_point(data = superbinders_tpm_df %>% filter(tpm_homo_sapiens_cytosolic_fraction > 0.001),
                                                       shape = 17, alpha = 0.7) +
                                            geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
                                            stat_cor() +
                                            scale_x_continuous(expand = c(0,0)) +
                                            scale_color_manual(values = c("#a8404c", "#424242"), name = "Gene type") + 
                                            ggtitle("Cytoplasmic Expression vs. promoter binding events") + 
                                            xlab(expression('Number of TFs')) +
                                            ylab(expression(log[2](TPM))) 
cytoplasmic_expression_vs_superbinders
```

![](00_MN_Final_Project_files/figure-gfm/High%20binding%20promoters-4.png)<!-- -->

``` r
ggsave("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/06_final_going_global/figures/cytoplasmic_expression_vs_superbinders.pdf")
```

    ## Saving 7 x 5 in image

``` r
superbinders_tpm_no_nuclear_expression_df <- filter(superbinders_tpm_df,tpm_homo_sapiens_nuclear_fraction %in% c(0:1))
length(superbinders_tpm_no_nuclear_expression_df$gene_name)
```

    ## [1] 239

``` r
superbinders_tpm_no_cytoplasmic_expression_df <- filter(superbinders_tpm_df,tpm_homo_sapiens_cytosolic_fraction %in% c(0:1))
length(superbinders_tpm_no_cytoplasmic_expression_df$gene_name)
```

    ## [1] 353

# Results

# 1) Superbinding promoters exhibit a very lienar trend when comparing RNA expression vs number of TFs

# 2) Superbinder expression looks very similar when comparing nuclear vs cytoplasmic expression

# 3) 239 superbinding promoters dont have any nuclear RNA expression

# 4) 353 superbinding promoters dont have any cytoplasmic expression

# 5) However, this only represent roughly 2-3% of the total superbinding promoter population

\`\`\`

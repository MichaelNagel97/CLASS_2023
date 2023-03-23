Class\_exercise\_mina8856
================
Michael Nagel
3/16/2023

# Loading in the pertinent libraries and functions

\#Loading in ChIP peak files of each protein and their replicates

``` r
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/chipseqprofs"
peak_path <- "results/bwa/mergedLibrary/macs/broadPeak_good_files"
broadpeakfilepath <- file.path(basepath,peak_path)

# Importing peaks to get peak list (Grange)
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# Identifying how many peaks are in each file 
peak_numbers <- sapply(peak_list, length) %>% as.data.frame()

# Labeling column as peak_numbers
names(peak_numbers) <- c("peak_numbers")

# Separating DBP name and replicate number
peak_numbers <- peak_numbers %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")

# peak_numbers data frame 
peak_numbers
```

    ##                dbp replicate peak_numbers
    ## 1            BRCA1        R3        44978
    ## 2            BRCA1        R4        44173
    ## 3         H3K36me3        R1       173838
    ## 4         H3K36me3        R2       175056
    ## 5         H3K36me3        R4       148235
    ## 6              JUN        R1        17007
    ## 7              JUN        R2        35109
    ## 8              JUN        R3        40224
    ## 9              JUN        R4        24645
    ## 10             MAX        R3        83146
    ## 11             MAX        R4       113543
    ## 12          POLR2A        R1        28136
    ## 13          POLR2A        R2        29601
    ## 14          POLR2A        R3        16204
    ## 15          POLR2A        R4        15515
    ## 16 POLR2AphosphoS2        R1        64499
    ## 17 POLR2AphosphoS2        R2        99880
    ## 18 POLR2AphosphoS5        R1        91908
    ## 19 POLR2AphosphoS5        R2        89007
    ## 20             TBP        R1        22070
    ## 21             TBP        R2        18302
    ## 22             TBP        R3        51536
    ## 23             TBP        R4        33932

``` r
# Saving peak_numbers data frame
write_csv(peak_numbers, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/peak_numbers_df.csv")
```

# Creating consensus peaks for each protein

``` r
# Generating list of unique DBPs from peak_list
dbps <- unique(sapply(names(peak_list), function(x) {
  unlist(strsplit(x, "_"))[1]
  }))

# Using the consensus_from_reduced function, which finds peaks that overlap between replicates and combines them
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)

#Adding names to dbps
names(consensus_list) <- dbps 

# Making a data frame of the consensus peaks for each DBP
num_consensus_peaks <- sapply(consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(number_consensus_peaks = ".")

# Merging the consensus peaks into the peak_numbers data frame
peak_numbers <- left_join(peak_numbers, num_consensus_peaks)
```

    ## Joining with `by = join_by(dbp)`

``` r
# Printing table
peak_numbers
```

    ##                dbp replicate peak_numbers number_consensus_peaks
    ## 1            BRCA1        R3        44978                  27154
    ## 2            BRCA1        R4        44173                  27154
    ## 3         H3K36me3        R1       173838                  66997
    ## 4         H3K36me3        R2       175056                  66997
    ## 5         H3K36me3        R4       148235                  66997
    ## 6              JUN        R1        17007                   3835
    ## 7              JUN        R2        35109                   3835
    ## 8              JUN        R3        40224                   3835
    ## 9              JUN        R4        24645                   3835
    ## 10             MAX        R3        83146                  69006
    ## 11             MAX        R4       113543                  69006
    ## 12          POLR2A        R1        28136                   4700
    ## 13          POLR2A        R2        29601                   4700
    ## 14          POLR2A        R3        16204                   4700
    ## 15          POLR2A        R4        15515                   4700
    ## 16 POLR2AphosphoS2        R1        64499                  38221
    ## 17 POLR2AphosphoS2        R2        99880                  38221
    ## 18 POLR2AphosphoS5        R1        91908                  49152
    ## 19 POLR2AphosphoS5        R2        89007                  49152
    ## 20             TBP        R1        22070                   8288
    ## 21             TBP        R2        18302                   8288
    ## 22             TBP        R3        51536                   8288
    ## 23             TBP        R4        33932                   8288

``` r
# Saving updates to data frame 
write_csv(peak_numbers, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/peak_numbers_df.csv")



# Exporting consensus peaks to results folder in .bed format

# Setting file path to export
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/mina8856"
consensus_path <- "CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/consensus_peaks/"
exportpath <- file.path(basepath, consensus_path)

# Exporting each DBP consensus peak as .bed file
for(i in 1:length(consensus_list)) {
  rtracklayer::export(consensus_list[[i]], paste0(exportpath, names(consensus_list)[i], "_consensus.bed") )}
```

# Making my consensus peaks compatable with UCSC genome browser

``` r
# Adding column names to .bed consensus files 
# lapply (for loop) across consensus file list to add colnames
# Column names for .broadPeak are: chr, start, end, name, score, strand

consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/consensus_peaks", full.names = T, pattern = ".bed")

peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))

names(peaks) <- dbps

# Removing contigs from peak files
# Make chromosomes a value
canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")

# Using lapply with filter function to cannonical_chr
peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))

# Exporting cleaned up consensus peaks
new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/consensus_peaks/", names(peaks), "_consensus.bed")
for(i in 1:length(peaks)) {
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

# UCSC Header 
headers <- paste0("track type=bed name=", names(peaks))

# creating a path to export after we add header in for loop below
new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/ucsc_consensus_peaks/", names(peaks), ".bed")
for(i in 1:length(peaks)) {
  writeLines(headers[[i]], new_filenames[[i]])
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

# Can now transfer file using file transfer program and look at consensus peaks in UCSC genome browser 

# pSer2 example
knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/pSer2_UCSC.png")
```

<img src="../pSer2_UCSC.png" width="1371" />

# Using human transcription factor (TF) annotations from a cell paper to identify which of my ChIP proteins are TF’s, their DNA binding domain, and their ensemble ID identification

``` r
# Making simplified data frame
# Removing replicates and their peak numbers from "peak_numbers" df
# Adding the total_peak_length

num_peaks_df <- data.frame("dbp" = names(consensus_list),
                           "num_consensus_peaks" = sapply(consensus_list, length),
                           "total_peak_length" = sapply(consensus_list, function(x) sum(width(x))))


# Downloading TF annotations from cell paper
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"
destination_for_url <- "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/TF_annotations.xlsx"
download.file(url, destination_for_url)

#redx1::read_excel to import
suppressWarnings({human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/TF_annotations.xlsx",
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

    ## [1] 5

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
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/num_peaks_df.csv")

# Printing out table
num_peaks_df
```

    ##               dbp num_consensus_peaks total_peak_length      ensembl_id     dbd
    ## 1           BRCA1               27154          21930858 ENSG00000012048 Unknown
    ## 2        H3K36me3               66997         114382144            <NA>    <NA>
    ## 3             JUN                3835           2670049 ENSG00000177606    bZIP
    ## 4             MAX               69006          77629013 ENSG00000125952    bHLH
    ## 5          POLR2A                4700          12715914 ENSG00000181222 Unknown
    ## 6 POLR2AphosphoS2               38221          45021137            <NA>    <NA>
    ## 7 POLR2AphosphoS5               49152          43730816            <NA>    <NA>
    ## 8             TBP                8288           9923640 ENSG00000112592     TBP
    ##     tf
    ## 1   No
    ## 2 <NA>
    ## 3  Yes
    ## 4  Yes
    ## 5   No
    ## 6 <NA>
    ## 7 <NA>
    ## 8  Yes

# Loading in the human genome from gencode

``` r
# Loading in the human genome
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

# Looking at genome gene annotations 
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
#table(gencode_gr$type)
#table(gencode_genes$gene_type)

# Exporting genes file
rtracklayer::export(gencode_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/gencode_genes.gtf")
```

# Creating mRNA and lncRNA annotation files from genes file

``` r
# mRNA genes (called "protein_coding") 
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 

# Exporting annotation file
rtracklayer::export(mrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/mrna_genes.gtf")

# lncRNA genes
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

# Exporting annotation file
rtracklayer::export(lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/lncrna_genes.gtf")

# Combining mRNA and lncRNA annotations files
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
mrna_lncrna_genes_df <- mrna_lncrna_genes %>% as.data.frame()
table(mrna_lncrna_genes$gene_type)
```

    ## 
    ##         lncRNA protein_coding 
    ##          16849          19965

``` r
# Exporting annotation file
rtracklayer::export(mrna_lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/mrna_lncrna_genes.gtf")
```

# Now let’s make a promoter annotation file and gene\_body for future use.

``` r
# Promoter annotations
lncrna_mrna_promoters <- promoters(mrna_lncrna_genes, upstream = 1000, downstream = 1000)

# Checking right size (Should all be 2000)
# width(lncrna_mrna_promoters)
# table(table(lncrna_mrna_promoters))

# Exporting annotation file
rtracklayer::export(lncrna_mrna_promoters, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/lncrna_mrna_promoters.gtf")

# lncRNA gene ID's (for subsetting) 

lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]

# lncRNA gene ID's (for subsetting)

mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]
```

# Importing Consensus Peaks and Annotation Files for ease of access

``` r
# Importing consensus list
consensus_list_files <- list.files("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/consensus_peaks", 
                            pattern = "*.bed",
                            full.names = TRUE)
consensus_list <- lapply(consensus_list_files, rtracklayer::import)
names(consensus_list) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/consensus_peaks/|_consensus.bed","", consensus_list_files)

# Importing annotation files
mrna_lncrna_genes <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/mrna_lncrna_genes.gtf")
lncrna_mrna_promoters <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/gene_annotations/lncrna_mrna_promoters.gtf")

# For subsetting
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]
```

# ChIP Proteins Data Analysis

# Determining how my proteins consensus peaks overlap mRNA and lncRNA genes and promoters

``` r
# Counting promoter overlaps using consensus peaks (columns are promoters, rows are DBPs)
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_list, type = "counts")

# Summing each row for each DBP to identify the amount of promoter overlaps  
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)


# Using gene_id objects (created above) to index into "lncrna" and "mrna" and separate them into two groups


# lncRNA promoter overlaps
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

# mrna promoter overlaps
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

# Printing out table
num_peaks_df
```

    ##               dbp num_consensus_peaks total_peak_length      ensembl_id     dbd
    ## 1           BRCA1               27154          21930858 ENSG00000012048 Unknown
    ## 2        H3K36me3               66997         114382144            <NA>    <NA>
    ## 3             JUN                3835           2670049 ENSG00000177606    bZIP
    ## 4             MAX               69006          77629013 ENSG00000125952    bHLH
    ## 5          POLR2A                4700          12715914 ENSG00000181222 Unknown
    ## 6 POLR2AphosphoS2               38221          45021137            <NA>    <NA>
    ## 7 POLR2AphosphoS5               49152          43730816            <NA>    <NA>
    ## 8             TBP                8288           9923640 ENSG00000112592     TBP
    ##     tf peaks_overlapping_promoters peaks_overlapping_lncrna_promoters
    ## 1   No                       28840                               6776
    ## 2 <NA>                        4358                               3122
    ## 3  Yes                        1400                                462
    ## 4  Yes                       37048                              10342
    ## 5   No                       10204                               2496
    ## 6 <NA>                       16822                               5324
    ## 7 <NA>                       34286                               8826
    ## 8  Yes                       16960                               3512
    ##   peaks_overlapping_mrna_promoters
    ## 1                            22064
    ## 2                             1236
    ## 3                              938
    ## 4                            26706
    ## 5                             7708
    ## 6                            11498
    ## 7                            25460
    ## 8                            13448

``` r
# Exporting updates data frame
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/num_peaks_df.csv")
```

## Results:

# 1) What can you determine from these overlaps?

-   Which proteins preferentially bind to promoter regions (comparing
    num\_consensus\_peaks with peaks\_overlapping\_promoters)

# 2) What is the difference in overlaps between mRNA and lncRNA promoters

-   There are more overlaps at mRNA promoters than lncRNA promoters,
    indicating that mRNA is more frequently transcribed than lncRNA, or
    that the proteins analyzed here are simply more abundant at mRNA
    promoters.
-   However, this can be partially attributed to there being a greater
    number of annotated mRNA genes (mRNA = 19965 annotations vs lncRNA =
    16849 annotations)

<!-- -->

# Identifying amount of protein binding over gene bodies and comparing with promoters

``` r
# Creating genebody variable using count_peaks_per_feature function
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_list, 
                                                type = "counts")

# Extracting the overlaps the same way we did for promoters above
# All gene bodies
num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

# lncRNA gene bodies 
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])

# Printing table
num_peaks_df
```

    ##               dbp num_consensus_peaks total_peak_length      ensembl_id     dbd
    ## 1           BRCA1               27154          21930858 ENSG00000012048 Unknown
    ## 2        H3K36me3               66997         114382144            <NA>    <NA>
    ## 3             JUN                3835           2670049 ENSG00000177606    bZIP
    ## 4             MAX               69006          77629013 ENSG00000125952    bHLH
    ## 5          POLR2A                4700          12715914 ENSG00000181222 Unknown
    ## 6 POLR2AphosphoS2               38221          45021137            <NA>    <NA>
    ## 7 POLR2AphosphoS5               49152          43730816            <NA>    <NA>
    ## 8             TBP                8288           9923640 ENSG00000112592     TBP
    ##     tf peaks_overlapping_promoters peaks_overlapping_lncrna_promoters
    ## 1   No                       28840                               6776
    ## 2 <NA>                        4358                               3122
    ## 3  Yes                        1400                                462
    ## 4  Yes                       37048                              10342
    ## 5   No                       10204                               2496
    ## 6 <NA>                       16822                               5324
    ## 7 <NA>                       34286                               8826
    ## 8  Yes                       16960                               3512
    ##   peaks_overlapping_mrna_promoters peaks_overlapping_genebody
    ## 1                            22064                      53986
    ## 2                             1236                     152886
    ## 3                              938                       6656
    ## 4                            26706                     121446
    ## 5                             7708                      14174
    ## 6                            11498                      84648
    ## 7                            25460                     107476
    ## 8                            13448                      22220
    ##   peaks_overlapping_lncrna_genebody peaks_overlapping_mrna_genebody
    ## 1                             11174                           42812
    ## 2                             10938                          141948
    ## 3                              1684                            4972
    ## 4                             27318                           94128
    ## 5                              2862                           11312
    ## 6                             10082                           74566
    ## 7                             16326                           91150
    ## 8                              4140                           18080

``` r
# Exporting data frame
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/num_peaks_df.csv")
```

## Results:

# 1) Do my proteins have more overlaps with promoters or genebodies?

-   Gene Bodies

# Identifying how many of my proteins bind specific promoters

``` r
# Creating promoter peak occurence variable using count_peaks_per_feature function and type = occurrence
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_list, 
                                               type = "occurrence")

# Let's double check that all lncrna & mrna genes are accounted for:
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

# Great we will use this quite a bit moving forward so let's write it out! 
write.table(promoter_peak_occurence, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Now let's use the 'data.frame()' fucntion. 
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))

# Exporting data frame
write_csv(peak_occurence_df, "/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/peak_occurence_dataframe.csv")
```

## Results: I find the max number of proteins on a promoter to be X

-   RPBS27 and others have 8 proteins (all of them) bound to its
    promoter

# Plotting results (ggplot)

# Relationship between peak number and total DNA covered

``` r
num_peaks_df <- read_csv('/scratch/Shares/rinnclass/CLASS_2023/mina8856/CLASS_2023/CLASSES/05_R_analyses/Class_Exercise/Results/num_peaks_df.csv')
```

    ## Rows: 8 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): dbp, ensembl_id, dbd, tf
    ## dbl (8): num_consensus_peaks, total_peak_length, peaks_overlapping_promoters...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# let's make this simple plot first: number of peaks -vs- total peak length
ggplot(num_peaks_df, aes(x = num_consensus_peaks, 
                 y = total_peak_length)) +
           geom_point(shape = 'circle',
             color = 'red')
```

![](class_exercise_mina8856_files/figure-gfm/ggplot-1.png)<!-- -->

# Coloring plot based on if the protein is a TF or not.

``` r
ggplot(num_peaks_df, aes(x = log2(num_consensus_peaks/1e3), 
                 y = total_peak_length/1e6,
                 color = tf == "Yes")) +
  geom_point()
```

![](class_exercise_mina8856_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

# Making a histogram of the number of peaks for each of my proteins

``` r
ggplot(num_peaks_df, aes(x = num_consensus_peaks, fill = tf)) +
   geom_histogram(bins = 30, position = "dodge")
```

![](class_exercise_mina8856_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# Now I want to facet this by the type of DNA binding domain my protein has.

``` r
ggplot(num_peaks_df, aes(x = num_consensus_peaks, 
                 y = total_peak_length )) +
  facet_wrap(dbd ~ .) +
  geom_point (shape = 'circle',
             color =  'red')
```

![](class_exercise_mina8856_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

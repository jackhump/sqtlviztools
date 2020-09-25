library(dplyr)
library(purrr)
library(data.table)
library(leafviz)
library(stringr)
library(readr)


library(sqtlviztools)

options(echo=TRUE)

# make example files
#rs4733060	8:27302432	PTK2B	chr8:27308412-27309002

# example input files
VCF <- "example/example.vcf"
clusters_table <- "example/clu_52949_cluster_counts.txt"
permutation_res <- "example/example_permutations.all.0.05.bh.txt"
# p values and betas
permutation_full_res <- "example/example_permutations.full.0.05.bh.txt"
  #"clu_52949 8:27302432"

# input files
# sample_table <- "data/ROSMAP/rosmap.id.bam.map.tsv"
# permutation_res <- "data/ROSMAP/permutations.all.ROSMAP.txt.gz.0.05.bh.txt"
# VCF = "data/ROSMAP/all_samples.vcf.gz"
# #clusters_table <- "data/ROSMAP/devsplicing_perind_numers.counts.gz"
# clusters_table <- "data/ROSMAP/ROSMAP_all_samples_cluster_counts_split.Rdata"
# permutation_full_res <- "data/ROSMAP/associations.all.ROSMAP.txt.gz"
# #snp_conversion <- "data/ROSMAP/chrpos_rsid_dbsnp.txt.gz"
# snps_out <- "example/snp_ids_rs_ids.txt"
# #cluster_folder <- "data/ROSMAP/junction_counts/"


# annotation - created by Leafcutter

annotation_code <- "example/PTK2B_gencode_hg19"
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

exons_table <- if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  as.data.frame(fread(cmd = paste("zless",exon_file)) )
} else {
  cat("No exon_file provided.\n")
  NULL
}

# samples

#sampleTable <- read.table(sample_table,header=FALSE, stringsAsFactors = FALSE)

# read in clusters

# if( !file.exists(clusters_table) ){
#
#   # need to read in cluster file for each sample and merge together
#   cluster_files <- list.files(cluster_folder, full.names = TRUE)
#
#   print("reading in clusters")
#   # lovely piece of code this - read in each file and left join.
#   # takes a while but does the job
#   clusters <- cluster_files %>%
#     map(~read_delim(., delim = " ") ) %>%
#     reduce( left_join)
#
#   # save
#   save(clusters, file = "data/ROSMAP/ROSMAP_all_samples_cluster_counts_ratios.Rdata")
#
#   #load("data/ROSMAP/ROSMAP_all_samples_cluster_counts_ratios.Rdata")
#   # split out junction numbers - Towfique gave me the bloody fractions instead
#   clusters_split <- clusters %>%
#     map(~str_split_fixed(. , pattern = "/", n = 2) ) %>%
#     map_df( `[`,1:nrow(clusters) ,1 )
#
#   d <- as.data.frame(clusters_split)
#   # set row names
#   row.names(d) <- d$chrom
#   d$chrom <- NULL
#   # convert to numeric
#   d <- as.data.frame(data.matrix(d, rownames.force = TRUE))
#
#   clusters <- d
#
#   # sort out names
#
#   samples <- names(clusters)
#
#   names(sampleTable) <- c("sample_code", "bam_name")
#
#   # convert sampleTable$V2 to the scheme in the clusters_table
#   sampleTable$cluster_samples <- sampleTable$bam_name %>%
#     gsub("_", ".", ., fixed = TRUE) %>%
#     gsub(".bam", ".final.bam", ., fixed = TRUE) %>%
#     paste0("ROSMAP_STAR_", . )
#
#   # get numeric code for each sample
#   renamed_samples <- sampleTable$sample_code[ match( samples, sampleTable$cluster_samples)]
#
#   # fix clusters with new names
#   names(clusters) <- renamed_samples
#
#
#   save( clusters, file = clusters_table)
#   # find and harmonise sample names
#
#   # write out samples
#   samples_file <- "data/ROSMAP/samples_to_use.txt"
#   writeLines( text=as.character(renamed_samples), con = samples_file)
# }else{
#   message("reading in clusters")
#   load(clusters_table)
# }

clusters <- read.table(clusters_table, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# read in junction x snp results



print("reading in results")

res <- as.data.frame(fread(permutation_res, header =TRUE), stringsAsFactors = FALSE)

genotypes <- unique(res$dummy2)


# genotypes do not have rs numbers - just chr:start values.
# but this matches VCF so maybe not a problem yet
# will need to match in rs numbers at some point, surely


# write out genotypes
#genotypes_file <- "data/ROSMAP/snps_to_use.txt"
#writeLines(text= genotypes, con = genotypes_file)

######
# PREPARE VCF GENOTYPES
######

# use vcftools to filter out just the snps and samples required
# print("filtering VCF")
# vcf_filtered <- "data/ROSMAP/filtered_samples_snps"
# vcf_filtered_full <- "data/ROSMAP/filtered_samples_snps.recode.vcf"
# if( !file.exists(vcf_filtered_full)){
#   cmd <- paste( "vcftools --gzvcf", VCF,  "--snps", genotypes_file, "--keep", samples_file, "--recode --out", vcf_filtered )
#   system(cmd)
# }
#

vcf <- read.table(file =  VCF, header=TRUE,check.names = FALSE)

vcf$REF[vcf$REF == TRUE] <- "T"
vcf$ALT[vcf$ALT == TRUE] <- "T"


###
## FIX SNP IDS

# fix association results
# snps <- read.table(snps_out, header=TRUE)
#
# res$dummy2 <- gsub(":", ".", paste0(snps$RS_id,".", snps$dummy2))
# # fix vcf
# vcf_snp_ids <- gsub(":", ".", paste0(snps$RS_id,".", snps$dummy2))[ match( vcf$ID, snps$dummy2)]
# vcf$ID <- vcf_snp_ids

# get_vcf_meta <- function(vcf){
#   # take a VCF and split out the meta data
#   # to use for querying
#   vcf_meta <-  as.data.frame(str_split_fixed(vcf$ID, "\\.", 3), stringsAsFactors = FALSE)
#   vcf_meta$SNP_pos <- paste(vcf_meta$V2, vcf_meta$V3, sep = ":")
#   vcf_meta <- data.frame( SNP = vcf_meta$V1, SNP_pos = vcf_meta$SNP_pos, REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)
#   return(vcf_meta)
# }
vcf_meta <- sqtlviztools::get_vcf_meta(vcf)


##################
# PREPARE CLUSTERS
##################

# from significant associations
sigClusters <- str_split_fixed(res[,1], ":",4)[,4]
### add Nalls' GWAS SNP clusters and Yang's too!
#sigClusters <- c(sigClusters, gwas_clusters, yang_clusters)

introns <- leafviz::get_intron_meta(row.names(clusters) )
keepClusters <- match(introns$clu,sigClusters)

# remove non-significant (or non-GWAS SNP-associated) clusters
introns <- introns[ !is.na(keepClusters),]
clusters <- clusters[ !is.na(keepClusters),]

# rearrange sample columns in clusters so they match the VCF
samples <- names(vcf)[10:ncol(vcf)]
clusters <- clusters[, samples]

introns_to_plot <- get_intron_meta(row.names(clusters))

# thin down junctions in clusters - remove low contributing junctions

juncProp <- function(cluster){
  cluster$prop <- cluster$meanCount / sum(cluster$meanCount)
  return(cluster)
}

splitClusters <- introns_to_plot %>%
  mutate(
    clu = factor(.$clu, levels = unique(.$clu)),
    meanCount = rowMeans(clusters) ) %>%
  split( .$clu ) %>%
  purrr::map_df( juncProp ) %>%
  mutate( clu = as.character(.$clu))

# thinning out clusters - turn off for now - this needs to be worked on
thinClusters <- FALSE

if( thinClusters == TRUE){
  introns_to_plot <- introns_to_plot[ splitClusters$prop >= 0.01,]
  clusters <- clusters[ splitClusters$prop >= 0.01,]
  introns <- introns[ splitClusters$prop >= 0.01,]
}
####################
# ANNOTATE JUNCTIONS
####################

# functions are now in sqtlviztools
intersects <- sqtlviztools::intersect_introns(introns)

threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")

uniqueClusters <- unique( introns$clu )

# for testing
#save( introns, uniqueClusters, fiveprime_intersect,threeprime_intersect,all.introns_intersect, file = "annotation_test.Rdata")
#load(file="~/Documents/sQTLviz/annotation_test.Rdata")
#annotate_single_cluster(introns, clu = uniqueClusters[1], cluIndex=1)

# for debugging

save.image("debug.RData")

annotatedClusters <- purrr::map_df( seq_along(uniqueClusters),
                                    ~annotate_single_cluster( introns,
                                                              clu = uniqueClusters[.],
                                                              cluIndex = .,
                                                              fiveprime=fiveprime_intersect,
                                                              threeprime=threeprime_intersect,
                                                              bothSS=all.introns_intersect
                                                              )
                                    )

annotatedClusters$gene[ is.na( annotatedClusters$gene) ] <- "."
annotatedClusters$ensemblID[ is.na( annotatedClusters$ensemblID) ] <- "."

#################
# PREPARE RESULTS - MOST SIGNIFICANT SNP x JUNCTION
#################


# get_snp_meta <- function(snps){
#   snp_meta <-  as.data.frame(str_split_fixed(snps, "\\.", 3), stringsAsFactors = FALSE)
#   colnames(snp_meta) <- c("snp_ID", "snp_chr", "snp_pos")
#   # a few snps don't have any IDs - just the coordinates
#   noID <- snp_meta[ grepl("\\.", snp_meta$snp_pos),]
#   noID <- select(noID, snp_pos, snp_ID, snp_chr)
#   names(noID) <- c("snp_ID", "snp_chr","snp_pos")
#   # put back in
#   snp_meta[ grepl("\\.", snp_meta$snp_pos),] <- noID
#
#   snp_meta$snp_pos <- as.numeric(snp_meta$snp_pos)
#   return(snp_meta)
# }

# format of SNPs in res should be rs140285343.chr10.75132726
#snpConversion <- fread(paste("zless",snp_conversion), header=FALSE, data.table=FALSE)
#anames(snpConversion) <- c("coord", "RS_id")
#snps <- res[,2,drop=FALSE]
#snps$RS_id <- snpConversion$V2[ match( snps$dummy2, snpConversion$V1)]
#write.table(snps, file = snps_out, col.names = TRUE, row.names=FALSE, sep = "\t", quote =FALSE)

## the results table should consist of a list of clusters with their most significant SNP
# bind intron_meta, intron, snp, p value, snp_meta
#sigJunctions <- cbind( get_intron_meta( res[,1]), res[, c(1,6,11)], get_snp_meta(res[,6]))

# use associations with q < 0.05 cut-off.
# Bind together metadata with original results
sigJunctions <- cbind( get_intron_meta( res[,1]),
                       res[, c(1,2,3)],
                       get_snp_meta(res[,2]))

# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions)

#names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>%
  dplyr::summarise( chr = first(chr),
                    start = min(start),
                    end = max(end),
                    snp = first(snp_ID),
                    snp_chr = first(snp_chr),
                    pos = first(snp_pos),
                    FDR = first(bpval) ) %>%
  dplyr::arrange(FDR)

####
## PREPARE FOR SHINY
####

code <- "example"
annotation_code <- "gencode_hg19"

resultsByCluster$gene <- annotatedClusters$gene[ match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

# fix coords without "chr"
# if( all( !grepl("chr", sample(resultsByCluster$chr, 100)) ) ){
#   resultsByCluster$chr <- paste0("chr", resultsByCluster$chr)
# }

resultsByCluster$cluster_pos = paste0(resultsByCluster$chr,":", resultsByCluster$start,"-",resultsByCluster$end)

resultsToPlot <- as.data.frame( select( resultsByCluster,
                                        SNP = snp,
                                        SNP_pos,
                                        gene = gene,
                                        cluster_pos,
                                        q = FDR
) )

row.names(resultsToPlot) <- resultsByCluster$clu

resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)

#########
## GET BETAS FOR EACH JUNCTION
#########


# for testing - until I get full assocations table
#junctionTable <- NULL


# get the Betas and per-junction q values
# perm_full <- read_delim( permutation_full_res,
#                          col_names = c("clusterID", "SNP","X","FDR", "Beta"),
#                          col_types = "ccinn",
#                          delim = " "
# )
perm_full <- fread(cmd = paste("zless", permutation_full_res), data.table=FALSE)
names(perm_full) <- c("clusterID", "SNP","X","FDR", "Beta")

junctionsNeeded <- introns_to_plot %>%
    mutate( chr = gsub("chr", "", chr)) %>%
    mutate( cluster = paste( chr, start, end, clu, sep = ":") ) %>%
    pull( cluster )

#perm_clean <- filter( perm_full, clusterID %in% junctionsNeeded )
snps <- read.table(file = "example/snp_ids_rs_ids.txt", stringsAsFactors = FALSE, header=TRUE)

perm_full$RS_id <- snps$RS_id[ match(perm_full$SNP, snps$dummy2)]

#perm_clean$RS_id <- snps$RS_id[ match( perm_clean$SNP, snps$dummy2)]
perm_clean <- perm_full %>%
  filter( clusterID %in% junctionsNeeded ) %>%
  #mutate( RS_id = snps$RS_id[ match( SNP, snps$dummy2)] ) %>%
  filter( !is.na( RS_id)) %>%
  filter( clusterID %in% junctionsNeeded ) %>%
  mutate( SNP_ID = gsub(":", ".", paste0(RS_id,".", SNP)))

# sort out SNP IDs
perm_clean <- select(perm_clean, clusterID, SNP=SNP_ID, Beta, FDR)

perm_clean <- cbind( perm_clean,
                     get_intron_meta(perm_clean$clusterID),
                     get_snp_meta(perm_clean$SNP) )

perm_clean$snp_chr <- paste0("chr", perm_clean$snp_chr)


# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate( clu = row.names(resultsToPlot) ) %>%
  left_join(introns_to_plot, by = "clu" ) %>%
  rename(snp_ID = SNP) %>%
  left_join( perm_clean,
             by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0( chr, ":", start, "-", end)) %>%
  left_join( annotatedClusters,
             by = c("clu" = "clusterID", "coord" )
  ) %>%
  select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate( Beta = signif(Beta, digits = 3),
          q = signif(q, digits = 3)) %>%
  mutate( Beta = ifelse(is.na(Beta), ".", Beta),
          q = ifelse(is.na(q), ".", q))

# why are some of the q values greater than 1?


#########
## SAVE OBJECTS
#########
# too big
rm(perm_full)
#save.image("example_data.Rdata")
print("saving objects")
save( annotatedClusters, # every junction needed
      sigJunctions, # every junction x SNP interaction
      resultsToPlot, #significant clusters and the most significant SNP
      #GWASresults, # associations with SNPs from a PD GWAS
      #YangResults, # associations with Yang's TWAS hit SNPs
      clusters, # junction counts for each sample
      vcf,# the genotypes of each sample
      vcf_meta, # the vcf metadata
      introns_to_plot, # all the intron positions
      #counts,
      #meta,
      exons_table, # the annotation
      junctionTable, # the junctions to display for each cluster
      #pca,
      #intron_summary,
      #cluster_summary,
      #introns_to_plot,
      #cluster_ids,
      #sample_table,
      annotation_code,
      code,
      file = "shiny/sQTL_results.Rdata"
)

#load("data/ROSMAP/all_data.Rdata")

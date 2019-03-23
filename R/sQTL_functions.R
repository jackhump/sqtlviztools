# functions used in preparing sQTLs for browsing

get_snp_meta <- function(snps){
  snp_meta <-  as.data.frame(stringr::str_split_fixed(snps, "\\.", 3), stringsAsFactors = FALSE)
  colnames(snp_meta) <- c("snp_ID", "snp_chr", "snp_pos")
  # a few snps don't have any IDs - just the coordinates
  noID <- snp_meta[ grepl("\\.", snp_meta$snp_pos),]
  noID <- dplyr::select(noID, snp_pos, snp_ID, snp_chr)
  names(noID) <- c("snp_ID", "snp_chr","snp_pos")
  # put back in
  snp_meta[ grepl("\\.", snp_meta$snp_pos),] <- noID

  snp_meta$snp_pos <- as.numeric(snp_meta$snp_pos)
  return(snp_meta)
}


get_vcf_meta <- function(vcf){
  # take a VCF and split out the meta data
  # to use for querying
  vcf_meta <-  as.data.frame(stringr::str_split_fixed(vcf$ID, "\\.", 3), stringsAsFactors = FALSE)
  vcf_meta$SNP_pos <- paste(vcf_meta$V2, vcf_meta$V3, sep = ":")
  vcf_meta <- data.frame( SNP = vcf_meta$V1, SNP_pos = vcf_meta$SNP_pos, REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)
  return(vcf_meta)
}

intersect_introns <- function(introns){
  all.introns <- introns
  all.introns$start <- as.numeric(all.introns$start)
  all.introns$end <- as.numeric(all.introns$end)


  # for each splice site write out a bed file
  all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID = clu)

  all.fiveprime <- data.frame( chr = all.introns$chr,
                               start = all.introns$start,
                               end = all.introns$start + 1,
                               clusterID = all.introns$clu)
  all.threeprime <- data.frame( chr = all.introns$chr,
                                start = all.introns$end,
                                end = all.introns$end + 1,
                                clusterID = all.introns$clu)
  all.file <- "all_junctions.bed"
  all.fiveprime.file <- "all_fiveprime.bed"
  all.threeprime.file <- "all_threeprime.bed"

  write.table( all.junctions, all.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )

  print( "BedTools intersect junctions with list of known splice sites")

  # first match junctions
  all.introns.cmd <- paste0("bedtools intersect -a ", all.file, " -b ", all_introns, " -wa -wb -loj -f 1" )
  all.introns_intersect <- data.table::fread(cmd = all.introns.cmd, header=FALSE)

  # intersect with bedtools to find the annotations of each splice site
  threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ",threeprime_file, " -wa -wb -loj -f 1" )
  threeprime_intersect <- data.table::fread(cmd = threeprime.cmd,header=FALSE)

  fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", fiveprime_file, " -wa -wb -loj -f 1" )
  fiveprime_intersect <- data.table::fread(cmd = fiveprime.cmd,header=FALSE)

  # remove temporary files
  rm.cmd <- paste("rm ", all.file, all.fiveprime.file, all.threeprime.file)
  
  system(rm.cmd)

  return( list(threeprime_intersect, fiveprime_intersect,all.introns_intersect))
}

#' Title
#'
#' @param introns
#' @param clu
#' @param cluIndex
#' @param fiveprime
#' @param threeprime
#' @param bothSS
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples
annotate_single_cluster <- function(introns, clu, cluIndex, fiveprime, threeprime, bothSS){
  #print(clu)
  # for each intron in the cluster, check for coverage of both
  # output a vector of string descriptions
  cluster <- introns[ introns$clu == clu , ]
  cluster$start <- as.integer(cluster$start)
  cluster$end <- as.integer(cluster$end)
  # subset intersects by clusterID (V4)
  # data.table method not working for some reason
  ## fprimeClu <- fiveprime[ V4 == clu,]
  # tprimeClu <- threeprime[ V4 == clu,]
  # bothSSClu <- bothSS[ V4 == clu,]
  fprimeClu <- dplyr::filter( fiveprime, V4 == clu)
  tprimeClu <- dplyr::filter( threeprime, V4 == clu)
  bothSSClu <- dplyr::filter( bothSS, V4 == clu)

  # for each intron in the cluster:
  #   create vector of overlapping splice sites, indexed by the row of the intersect
  # five prime splice sites
  # some clusters only have 1 junction so it's simple
  if( nrow(cluster) == 1 ){
    fprime <- list(fprimeClu)
    tprime <- list(tprimeClu)
    bothSS <- list(bothSSClu)
    cluster_genes <- unique( c(tprime[[1]]$V8, fprime[[1]]$V8, bothSS[[1]]$V8 ) )
    cluster_ensemblIDs <- unique( c(tprime[[1]]$V9, fprime[[1]]$V9, bothSS[[1]]$V9 ) )
  }else{
    fprime <- apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which( names(cluster) == "start" )
      #fprimeClu[ V1 == x[chr] & V2 == as.numeric( x[start] ),]
      dplyr::filter(fprimeClu, V1 == x[chr] & V2 == as.numeric(x [ start]) )
    } )

    # three prime splice sites
    tprime <- apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      end <- which( names(cluster) == "end" )
      #tprimeClu[V1 == x[chr] & V2 == as.numeric( x[end] ),]
      dplyr::filter(tprimeClu, V1 == x[chr] & V2 == as.numeric(x [end]) )
    } )

    # both splice sites
    bothSS <-  apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which(names(cluster) == "start")
      end <- which( names(cluster) == "end" )

      # bothSSClu[
      #   V6 == as.numeric( x[start] ) &
      #     V7 == as.numeric( x[end] ) ,]
      dplyr::filter(bothSSClu, V6 == as.numeric(x[start]) & V7 == as.numeric(x[end]))
    } )
    # find gene and ensemblID by the most represented gene among all the splice sites
    cluster_genes <- names(sort(table(do.call( what = rbind, c(tprime, fprime, bothSS)  )$V8), decreasing = TRUE ))
    cluster_ensemblIDs <- names(sort(table(do.call( what = rbind, c(tprime, fprime, bothSS) )$V9), decreasing = TRUE ))
  }

  cluster_gene <- cluster_genes[ cluster_genes != "." ][1]
  # if no cluster gene found then leave as "."
  if( length(cluster_gene) == 0){
    cluster_gene == "."
  }

  cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]
  if( length( cluster_ensemblID ) == 0 ){
    cluster_ensemblID == "."
  }

  verdict <- c()
  coord <- c()
  gene <- c()
  ensemblID <- c()
  transcripts <- list()

  for( intron in 1:nrow(cluster) ){
    coord[intron] <- paste0(cluster[intron,]$chr,":", cluster[intron,]$start,"-", cluster[intron,]$end )

    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID

    # for each intron create vector of all transcripts that contain both splice sites
    #transcripts[[intron]] <- unique( intersect( tprime[[intron]]$V10, fprime[[intron]]$V10 ) )

    verdict[intron] <- "error"
    if( # if neither are annotated
      all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 == "." )
    ){ verdict[intron] <- "cryptic_unanchored"
    }
    if( # if only one is annotated
      all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 != "." )
    ){ verdict[intron] <- "cryptic_threeprime"
    }
    if(
      all( tprime[[intron]]$V5 != ".") & all( fprime[[intron]]$V5 == "." )
    ){ verdict[intron] <- "cryptic_fiveprime"
    }
    if( # if both splice sites are annotated
      all( tprime[[intron]]$V5 != "." ) & all( fprime[[intron]]$V5 != "." )
    ){
      # test if the splice sites are paired in a known intron
      if( all(bothSS[[intron]]$V5 != ".") ){
        verdict[intron] <- "annotated"
      }else{ # both are annotated but never in the same junction
        verdict[intron] <- "novel annotated pair"
      }
    }

  }
  if( cluIndex %% 500 == 0 ){
    print( paste("processed", cluIndex, "clusters" ))
  }
  #print(clu)

  return(
    data.frame(
      clusterID = clu,
      coord = coord,
      gene = gene,
      ensemblID = ensemblID,
      verdict = verdict,
      stringsAsFactors = FALSE)
  )
}

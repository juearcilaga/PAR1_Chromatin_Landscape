
#------------------------------------------------------------------------------------------------------#
#                             Regions in hg19
#------------------------------------------------------------------------------------------------------#

AKAP17 <- data.frame(Title_plot = "Bait: AKAP17:chrX_1693613_1711602",
                     start = 1693613,
                     end = 1711602)

DHRSX_ZBED1 <-
  data.frame(Title_plot = "Bait: DHRSX;ZBED1: chrX_2417048_2421197",
             start =2417048,
             end = 2421197)

P2RY8_bait<-data.frame(Title_plot="P2RY8-bait chrX_1629232_1662177", 
                    start= 1629232, 
                    end= 1662177)


P2RY8_3<-data.frame(Title_plot="Oe: P2RY8-3 chrX_1594287_1598632", 
                    start= 1594287, 
                    end= 1598632)

P2RY8_6 <-
  data.frame(Title_plot ="Oe: P2RY8-6 chrX_1613326_1613975",
             start = 1613326,
             end = 1613975)



P2RY8_7 <- data.frame(Title_plot =
                        "Oe: P2RY8-7 chrX_1623259_1629231",
                      start = 1623259,
                      end = 1629231)

GTPBP6 <-
  data.frame(Title_plot = "Bait: GTPBP6;LINC00685;PPP2R3B chrX_224862_318353",
             start = 224862,
             end = 318353)

PPP2R3B_1 <- data.frame(Title_plot =
                          "Oe:PPP2R3B-1 chrX_329721_334032",
                        start =329721,
                        end =334032)


PPP2R3B_2<-data.frame(Title_plot = "Oe:PPP2R3B-2 chrX_334033_335968",
             start = 334033,
             end = 335968)
#------------------------------------------------------------------------------------------------------#
#                               Conversion to hg38 using liftOver
#------------------------------------------------------------------------------------------------------#
library(rtracklayer)


RangesofInterest <- data.frame(rbind(
  AKAP17,
  DHRSX_ZBED1,
  P2RY8_bait,
  P2RY8_3,
  P2RY8_6,
  P2RY8_7,
  GTPBP6,
  PPP2R3B_1,
  PPP2R3B_2
))

RangesofInterest$seqnames="chrX"

GRangesofInterest <-
  makeGRangesFromDataFrame(
    RangesofInterest,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqinfo = NULL,
    seqnames.field = c(
      "seqnames",
      "seqname",
      "chromosome",
      "chrom",
      "chr",
      "chromosome_name",
      "seqid"
    ),
    start.field = "start",
    end.field = c("end", "stop"),
    strand.field = "strand",
    starts.in.df.are.0based = FALSE
  )

curl::curl_download(
  "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  destfile = "hg19ToHg38.over.chain.gz")

# gunzip hg19ToHg38.over.chain.gz

ch = import.chain("hg19ToHg38.over.chain")

seqlevelsStyle(GRangesofInterest) = "UCSC"  # necessary
GRangesofInterest38 = liftOver(GRangesofInterest, ch)
GRangesofInterest38 = unlist(GRangesofInterest38)
genome(GRangesofInterest38) = "hg38"

RangesofInterest38<-data.frame(GRangesofInterest38)
#------------------------------------------------------------------------------------------------------#
#                   Getting the chromatin states X chromosome
#------------------------------------------------------------------------------------------------------#


chrlabels<-read.table("labelsChr.txt", header = TRUE,sep="\t")

#Import chromatin states merged
All_data<-readRDS("chromatin_states_merged.rds")

XChr <- rownames(All_data)[grepl("chrX", rownames(All_data))]

#------------------------------------------------------------------------------------------------------#
#                                           Import metadata
#------------------------------------------------------------------------------------------------------#
library(dplyr)

metadata <- read.table("METADATA.txt", sep = "\t")

mycelltypes <- metadata  %>% filter(., DONOR_SEX == "Female")

# Cell types Biola PC-HiC
Biola_PCHiC_cells <-
  unique(mycelltypes$CELL_TYPE)[c(5, 6, 7, 8, 10, 11, 12, 13, 16, 19, 23)]

mycelltypes <-
  mycelltypes %>% filter(., CELL_TYPE %in% Biola_PCHiC_cells)

Replicates <-
  mycelltypes %>% dplyr::count(., CELL_TYPE) %>% dplyr::filter(., n > 1) %>% dplyr::select(CELL_TYPE)

Replicates.names <-
  mycelltypes %>% dplyr::filter(., CELL_TYPE %in% Replicates$CELL_TYPE) %>% dplyr::select(SAMPLE_NAME)

Singletons <-
  mycelltypes %>% dplyr::count(., CELL_TYPE) %>% dplyr::filter(., n == 1) %>% dplyr::select(CELL_TYPE)

Singletons.names <-
  mycelltypes %>% dplyr::filter(., CELL_TYPE %in% Singletons$CELL_TYPE) %>% dplyr::select(SAMPLE_NAME)

#------------------------------------------------------------------------------------------------------#
#                              Conserved chromatin states
#------------------------------------------------------------------------------------------------------#

# 75% of the samples should  be consistent
Sum <- function(x) {
  sum(x) >= round(length(x) * 0.75)
}


FindConserved <- function(chromstate) {
  celltype.chromstate <- celltype == chromstate
  ifelse(as.numeric(dim(data.frame(celltype))[2]) > 1,
         Suma <- apply(celltype.chromstate, 1, Sum),
         Suma <- TRUE)
  
  celltype <- data.frame(celltype, Suma)
  celltype$conserved <- 0
  celltype[celltype$Suma == TRUE, "conserved"] <-  chromstate
  
  print(table(celltype$conserved))
  conserved <-
    data.frame(Conserved = celltype$conserved,
               row.names = rownames(celltype.chromstate))
  colnames(conserved) <- chrlabels$NameState[chromstate]
  return(conserved)
}

#Par1 38 build
PAR1_start<-10001
PAR1_end<-2781479

start <- (PAR1_start / 200) + 1
end <- (PAR1_end / 200) + 1

XPar1 <- XChr[start:end]
Par1States <- All_data[rownames(All_data) %in% XPar1, ]

Conserved = list()

for (cell in unique(mycelltypes$CELL_TYPE)) {
  print(cell)
  celltype <-
    Par1States[, colnames(Par1States) %in%
                 metadata[metadata$CELL_TYPE == cell, "SAMPLE_NAME"]]
  Conservedlist <-
    mclapply(as.numeric(chrlabels$SymbState),
             FindConserved,
             mc.cores = 35)
  Conserved[[as.character(cell)]] <- do.call(cbind, Conservedlist)
}

ConservedDF <- mclapply(Conserved, rowSums, mc.cores = 35)
ConservedDF2 <- data.frame(do.call(cbind, ConservedDF))

#------------------------------------------------------------------------------------------------------#
#                              Import PAR1 Genes annotation
#------------------------------------------------------------------------------------------------------#
library(tidyr)

Par1Genes <- read.delim("PAR1_genes_ENSEMBL", sep = "\t")
colnames(Par1Genes)

Par1Labels <-
  unique(Par1Genes[, c(
    "Gene.name",
    "Gene.stable.ID",
    "Gene.type",
    "Gene.start..bp.",
    "Gene.end..bp.",
    "Strand"
  )])

Par1Labels$Labels <- as.character(Par1Labels$Gene.name)

for (each in 1:length(Par1Labels$Gene.name)) {
  if (Par1Labels$Gene.name[each] == "") {
    Par1Labels$Labels[each] <-
      as.character(Par1Labels$Gene.stable.ID[each])
  }
}


Par1Labels$color <-
  ifelse(
    Par1Labels$Strand == 1,
    Par1Labels$color <-
      "darkorchid2",
    Par1Labels$color <- "darkorchid1"
  )

Par1Labels[Par1Labels$Gene.type != "protein_coding", "color"] <-
  "white"

Par1Labels$pos <- 12
Par1Labels$NameState <- Par1Labels$Gene.type
# Par1Labels$acronym<-"gene"
Par1Labels$acronym <- Par1Labels$Labels
Par1Labels$SymbState <- "gene"
Par1Labels$start <- Par1Labels$Gene.start..bp.
Par1Labels$end <- Par1Labels$Gene.end..bp.
Par1Labels$variable <- "gene"
Par1Labels$Gene.type <- as.character(Par1Labels$Gene.type)
Par1Labels$NameState <- as.character(Par1Labels$NameState)
Par1Labels[Par1Labels$Gene.type != "protein_coding", "NameState"] <-
  "Other"

#------------------------------------------------------------------------------------------------------#
#                             For each of the regions create plots with chrom states
#------------------------------------------------------------------------------------------------------#

library(parallel)
library(MASS)
library(RColorBrewer)
library(scales)
library(ggplot2)



for (i in 1:nrow(RangesofInterest)) {
  # i=3
  Title_plot <- RangesofInterest38$Title_plot[i]
  start1 <- RangesofInterest38$start[i]
  end1 <- RangesofInterest38$end[i]
  
  start <- paste("chrX_", floor((start1 / 200) + 1), sep = "")
  end <- paste("chrX_", ceiling((end1 / 200) + 1), sep = "")
  
  XRegion <-
    ConservedDF2[which(rownames(ConservedDF2) == start):which(rownames(ConservedDF2) == end) ,]
  #------------------------------------------------------------------------------------------------------#
  #                                 Preparing the data for plotting
  #------------------------------------------------------------------------------------------------------#
  
  ensayo <- XRegion
  ensayo$coord <- rownames(ensayo)
  
  ensayo <-
    tidyr::pivot_longer(
      data.frame(ensayo),
      col = 1:10,
      names_to = "CELL_TYPE",
      values_to = "COLOR"
    ) %>% dplyr::select(coord, CELL_TYPE, COLOR)
  
  
  ensayo <- transform(ensayo, ID = as.numeric(factor(CELL_TYPE)))
  ensayo <-
    separate(ensayo, col = "coord", into = c("chr", "region"))
  
  mydata = data.frame(
    variable = ensayo$CELL_TYPE,
    pos = ensayo$ID,
    start = (as.numeric(ensayo$region) * 200) - 199,
    end = as.numeric(ensayo$region) * 200,
    SymbState = ensayo$COLOR
  )
  
  
  celltype_Abb <- data.frame(
    variable = mydata$variable[1:11],
    acronym = c(
      "nCD4",
      "nCD8",
      "Neu",
      "Mon",
      "Ery",
      "MK",
      "M0",
      "M2",
      "M1",
      "nB",
      "EndP"
    )
  )
  
  ### Only Region of interest
  mydata <- left_join(mydata, celltype_Abb)
  chrlabels$SymbState <- as.character(chrlabels$SymbState)
  chrlabels$NameState <- as.character(chrlabels$NameState)
  chrlabels[8, ] <- c("0", "Not conserved")
  mydata$SymbState <- as.character(mydata$SymbState)
  mydata <- left_join(mydata, chrlabels)
  
  chrcolor <- brewer.pal(n = 7, name = 'Dark2')
  chrcolor[c(8, 9, 10)] <- c("#335387", "darkorchid", "black")
  chrcolor[6] <- "grey"
  level_order <-
    c(
      "Other",
      "protein_coding",
      "Not conserved",
      "TransElo",
      "HetChrom",
      "PolyComb",
      "AProm",
      "AeloE",
      "AconE",
      "PoisE"
    )
  
  names(chrcolor) <- rev(level_order)
  
  Fragment_plot <-
    ggplot(mydata,
           aes(
             x = start,
             xend = end,
             y = acronym,
             yend = acronym,
             color = NameState
           ),
           size = 3) +
    geom_segment(size = 5) +
    scale_color_manual(values = chrcolor) +
    theme_bw() +
    xlab("Position in chrX") +
    ggtitle(Title_plot) +
    ylab("Celltype")
  
  print(Fragment_plot)
  
  
  Par1Labels1 <- Par1Labels[Par1Labels$Gene.end..bp. > start1 &
                              Par1Labels$Gene.start..bp. < end1,]
  
  if (nrow(Par1Labels1) > 0) {
    ### The Context of that region of interest inside the gene
    
    mydata2 <- mydata[mydata$start > min(Par1Labels1$start) &
                        mydata$end < max(Par1Labels1$end),]
    
    mydata2 <- rbind(mydata2, Par1Labels1[, colnames(mydata2)])
    unique(mydata2$NameState)
    
    
    Gene_plot <-
      ggplot(
        mydata2,
        aes(
          x = start,
          xend = end,
          y = acronym,
          yend = acronym,
          color = NameState
        ),
        size = 3) +
      geom_segment(size = 5) +
      scale_color_manual(values = chrcolor) +
      theme_bw() +
      xlab("Position in chrX") +
      ggtitle(Title_plot) +
      ylab("Celltype")
    
    print(Gene_plot)
  }
}

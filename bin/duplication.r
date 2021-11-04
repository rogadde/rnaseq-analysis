#!/home/AD/rkgadde/miniconda3/envs/jupyter/bin Rscript
#
# dupRadar shell script
# call dupRadar R package from the shell 
#
# Authors: Rohini Gadde, Holger Klein & Sergi Sayols
#
# https://github.com/ssayols/dupRadar
# http://bioconductor.org/packages/release/bioc/vignettes/dupRadar/inst/doc/dupRadar.html
#
# Usage: ./dupRadar.sh $1 $2 $3 $4
#   This script marks alignments for PCR duplicates then uses those marked
# marked alignments as input to dupRadar, which generates plots about the
# sequence duplication levels.
# 1 - gtf file
# 2 - strandedness (0|1|2)
# 3 - paired-end? (TRUE|FALSE)
# 4 - number of threads

library(dupRadar)

# get name patterns from command line
args <- commandArgs(TRUE)

# usually, same GTF file as used in htseq-count
gtf <- args[1]
# 0(no)|1(yes)|2(reverse)
stranded <- as.integer(args[2])
# is a paired end experiment
paired <- args[3]
# number of threads to be used
threads <- as.integer(args[4])

if(length(args) != 4) { 
  stop (paste0("Usage: ./dupRadar.sh <file.bam> <genes.gtf> ",
               "<0(no)|1(yes)|2(reverse)> <TRUE(PE)|FALSE(SE)> ",
               "<threads>"))
}

if(!file.exists(gtf)) {
  stop(paste("File",gtf,"does NOT exist"))
}

if(is.na(stranded) | !(grepl("0|1|2",stranded))) {
  stop("Stranded has to be 0|1|2")
}

if(is.na(paired) | !(grepl("TRUE|FALSE",paired))) {
  stop("Paired has to be TRUE|FALSE")
}

if(is.na(threads)) {
  stop("Threads has to be an integer number")
}

# end command line parsing

# path based on previous scripts 
files <- list.files(path="./STAR", pattern="sortedByCoord", full.names=TRUE,
  recursive=TRUE)
  
for (alignment in files) {
  # coordinate-sorted bam file to analyze
  bam <- alignment
  outdir <- dirname(bam)
  
  # analyze duprates and create plots
  cat("Processing file ", bam, " with GTF ", gtf, "\n")
  
  # mark duplicates with BamUtil
  bamDuprm <- markDuplicates(dupremover="bamutil",
                bam,
                out=gsub("Aligned.sortedByCoord.out.bam", "_duprm.bam", bam),
                path="/home/AD/rkgadde/miniconda3/pkgs/bamutil-1.0.15-h2e03b76_1/bin",
                rminput=FALSE)
  
  # calculate duplication rate matrix
  dm <- analyzeDuprates(bamDuprm,
                        gtf,
                        stranded,
                        (paired == "TRUE"),
                        threads)
  
  # produce plots
  # duprate vs. expression smooth scatter
  png(file=paste0(outdir,"/",gsub("_duprm.bam","",basename(bamDuprm)),"_dupRadar_drescatter.png"),
      width=1000, height=1000)
  duprateExpDensPlot(dm, main=basename(bamDuprm))
  dev.off()
  
  # expression histogram
  png(file=paste0(outdir,"/",gsub("_duprm.bam","",basename(bamDuprm)),"_dupRadar_ehist.png"),
      width=1000, height=1000)
  expressionHist(dm)
  dev.off()
  
  ## duprate vs. expression boxplot
  png(file=paste0(outdir,"/",gsub("_duprm.bam","",basename(bamDuprm)),"_dupRadar_drebp.png"),
      width=1000, height=1000)
  par(mar=c(10,4,4,2)+.1)
  duprateExpBoxplot(dm, main=basename(bamDuprm))
  dev.off()
}

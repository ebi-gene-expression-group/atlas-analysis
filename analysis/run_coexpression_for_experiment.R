#!/usr/bin/env Rscript
#
# Script to create a gene coexpression matrix for a baseline Expression Atlas experiment.

# Load BiocParallel package.
library(BiocParallel)
#source("kmeanCluster.R")

#############
# Functions #
#############

check_file_exists <- function( filename ) {

	if( !file.exists( filename ) ) {
		stop( paste(
			"Cannot find:",
			filename
		) )
	}
}

# Function that calculates coexpression values between genes and writes a matrix. First line shows default values
kCluster <- function(cD, maxK = 100, matrixFile = NULL, cores = 4, forceReplicates = FALSE, algorithm = "Lloyd") {    
    
    dat <- cD
    message("K-means processing...", appendLF = FALSE) # print a message to the user
    
    # this bit does all possible clusterings, parallelised
    kgene <- do.call("cbind", bplapply(1:nrow(dat), function(ii, maxK) # do all clusterings per gene, and set a default maximum k allowed of 100
        {
            x <- dat[ii,]
            if(sample(1:100, 1) == 1) message(".", appendLF = FALSE)   # print a dot to the screen for each gene evaluated    
        
            sapply(1:min(length(x), maxK), function(k) # iterate from 1 cluster to the maximum length of k or the maxK (default 100) whichever is smaller
                {
                    if(k > length(unique(x))) # return c(NA,NA) for every k> unique(x)
                        return(rep(NA, 2))

                    if(k == length(x))
                        return(c(paste(paste(1:length(x), collapse = ":"), paste(order(x), collapse = "<"), sep = "-"), "0")) # if you have the same number of clusters as data points, then return each data point as a cluster, to avoid giving it to kmeans which sometimes may output an error. The statistic here is 0 because all the data points are equal to the cluster centers
                    
                    clustK <- suppressWarnings(kmeans(x, k, iter.max = 1000, nstart = 100, algorithm = algorithm))            #  run kmeans    

                    if(any(is.na(clustK$centers)))
                        clustK <- suppressWarnings(kmeans(x, k, iter.max = 1000, nstart = 100, algorithm = algorithm)) # another error catching: sometimes kmeans returns an empty cluster, so run it again. The chances are very small for this to happen twice
                    
                                        #store the clustering and orderings as a string     
                    clustering <- paste(match(clustK$cluster, (unique(clustK$cluster))), collapse = ":") # the match bit makes sure that identical clusterings have identical numberings (the output of kmeans gives arbitrary numbers to the clusters but we need them to be similar for comparisons below
                    ordering <- paste(order(clustK$centers[unique(clustK$cluster)]), collapse = "<") # a string encoding the clusters and the order on the centre of those clusters (for comparisons later on)
                    clusterOrdering <- paste(clustering, ordering, sep = "-")

                    stat <- max(sapply(1:k, function(ii)  # iterate over the k clusters
                        {
                            mx <- suppressWarnings(max(abs(split(x, factor(clustK$cluster, levels = 1:k))[[ii]] - clustK$centers[ii])))
                            mx
                        }))
                                        # statistic gets stored as a string too, with marginal loss of precision.
                    return(c(clusterOrdering, stat))                  
                })
        }, maxK = maxK, BPPARAM = MulticoreParam(workers = as.integer(cores), progressbar = TRUE))) 
    message(".done!")
                                        # making sure kgene is a matrix
    kgene <- matrix(kgene, ncol = nrow(dat))
    
                                        # separate out clusterings and stats
    clusterings <- kgene[1:(nrow(kgene) / 2) * 2 - 1,]
    clusterings[is.na(clusterings)] <- FALSE
    stats <- matrix(as.numeric(kgene[1:(nrow(kgene) / 2) * 2,]), nrow = nrow(kgene) / 2)

    stats <- round(stats, 3)

                                        #temporary matrix prototype
    infmat <- matrix(Inf, nrow = nrow(stats), ncol = ncol(stats))    

    # check if a file name is provided in the matrixFile variable. If it isn't null, then prepare the file for writing
    if(!is.null(matrixFile))
        {        
            if(substr(matrixFile, nchar(matrixFile) -2, nchar(matrixFile)) != ".gz")
                {
                    matrixFile <- paste(matrixFile, ".gz", sep = "")

                    message( paste( "Coexpressions matrix will be written to:", matrixFile ) )
                }
            if(file.exists(matrixFile)) {
                
                message( "File with the same name already exists, removing it." )

                file.remove(matrixFile)
            }

            gzOutFile <- gzfile(matrixFile, "w")
        }

        
    
                                        # write one by one the stats for each gene.


    message("Comparing clusterings...", appendLF = FALSE)

    if(!is.null(matrixFile))
        {
            cores <- 1 # here you cannot parallelise, so dropping to 1 core
            writeLines(paste(c("",rownames(dat)), collapse ="\t"),gzOutFile)  # write gene names on first line
        }
    kAM <- do.call("rbind", bplapply(1:ncol(clusterings), function(ii)
        {
            if(sample(1:100, 1) == 1)
                message(".", appendLF = FALSE)       
                                        # fill in temporary matrix with valid statistics for each comparison
            tmat <- infmat; tmat[(clusterings[,ii] == clusterings)] <- stats[clusterings[,ii] == clusterings]
                                        # update temporary matrix with stats from current gene (if bigger)
            tmat <- apply(cbind(stats[,ii], tmat), 1, function(x) pmax(x[-1], x[1]))    
    
                                        # select minimum statistic from each column of temporary matrix and write

            minstats <- do.call(pmin, c(lapply(1:ncol(tmat), function(i)tmat[,i]), list(na.rm = TRUE)))

            if(!is.null(matrixFile))
                writeLines(paste(c(rownames(dat)[ii], minstats), collapse = "\t"), gzOutFile)

            # for clustering, we want to know for the ii'th gene the closest co-expressed gene that has a higher index than ii
            minstats[1:ii] <- Inf # so set all stats for lower indices to Inf
            c(id = ii, pair.id = which.min(minstats), stat = min(minstats)) # and choose the gene with minimum statistic

        },
                                     BPPARAM = MulticoreParam(workers = as.integer(cores), progressbar = TRUE)))

    if(!is.null(matrixFile))
        close(gzOutFile)

    message(".done!")

    kAM
}



###############################
# Script start.


# Get the commandline arguments.
args <- commandArgs( TRUE )

# Stop if we don't have the requisite number of arguments.
# For now this is 1 -- we just want the Atlas experiment directory. We can
# build the other filenames from that.
if( length( args ) == 2 ) {
	
    # Path to the Atlas experiment directory.
	atlasExperimentDirectory <- args[ 1 ]
    
    # Name of file containing expressions matrix
    expressionsFile <- args[ 2 ]
} else {
	# Print a usage message and exit.
	stop( "\nUsage:\n\trun_coexpression_for_experiment.R <Atlas experiment directory> <expressions filename with quartiles>\n\n" )
}

# get the experiment accession
atlasExperimentAccession <- unlist(strsplit(expressionsFile, "[.]"))[1]

message( paste( "Experiment accession is:", atlasExperimentAccession ) )

# combine file path
expressionsFile <- file.path( atlasExperimentDirectory, expressionsFile )

# Check the directory provided exists, die if not.
check_file_exists(expressionsFile)

message( paste( "Reading expression data from file:", expressionsFile ) )

# read file
exp <- read.delim(expressionsFile)
exp[,3:ncol(exp)] <- sapply(exp[3:ncol(exp)], function(x) sub("(^[^,]+[,][^,]+[,])([^,]+)(,.+$)", "\\2",x))  #get the middle value for each gene/tissue
expL <- sapply(exp[,3:ncol(exp)], as.numeric) # make sure the values are numeric
rownames(expL) <- exp[,1]
expL <- log(expL) # get the natural logarithm

                                        
expL[is.na(expL)] <- 0 # turn any NAs to 0

cD <- expL[rowSums(expL) > ncol(expL),] # Filter out non-expressed genes

# preparing the output file name
outFileName <- paste( atlasExperimentAccession, "coexpressions.tsv", sep="-")

# Add full path to experiment directory to output file.
outFileName <- file.path( atlasExperimentDirectory, outFileName )

message( paste( "Coexpression matrix will be written to gzipped file:", outFileName ) )

kClust <- kCluster(cD, cores = 10, matrixFile= outFileName) # run kCluster function to create coexpression matrices, set to use 10 cores. It can be changed to use less or more

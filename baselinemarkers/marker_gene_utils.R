
# Functions modified from MGFR package, initially used for proteomics marker calculation


# Function to test if a given gene is a marker
.isMarker.rnaseq <- function(named.vec, rep.vec) {
    sort.vec <- sort(named.vec, decreasing = TRUE)
    sv.len <- length(sort.vec)
    names.vec <- names(sort.vec)
    rep.fe <- unname(rep.vec[names.vec[1]])
    poss.marker <- (length(unique(names.vec[1:rep.fe])) == 1 && sum((as.numeric(sort.vec)[1:rep.fe]) > 
                                                                        1) == rep.fe)
    if (poss.marker) {
        mean.vec <- c()
        sort.num <- as.numeric(sort.vec)
        mean.vec <- c(mean.vec, mean(sort.num[1:rep.fe]))
        start.p <- rep.fe + 1
        rep.se <- unname(rep.vec[names.vec[start.p]])
        end.pos <- start.p + rep.se - 1
        cp.found <- FALSE
        while (end.pos <= sv.len && !cp.found) {
            len.sa <- length(unique(names.vec[start.p:end.pos]))
            if (sum(unique(names.vec[start.p:end.pos]) %in% names.vec[(end.pos + 1):sv.len]) == 
                0 | end.pos == sv.len) {
                mean.vec <- c(mean.vec, mean(sort.num[start.p:end.pos]))
                cp.found <- TRUE
            } else {
                end.pos <- start.p + sum(rep.vec[unique(names.vec[start.p:end.pos])]) - 1
            }
        }
        #return(c(names.vec[1], (mean.vec[2]/mean.vec[1])))
        return(setNames(c(names.vec[1], mean.vec[2]/mean.vec[1]), c("gene", "score")))
        
        
    }
}

## Function to map gene identifiers (ensembl, refseq or ucsc ids) to gene symbols and entrez
## gene ids using biomaRt.
.get.genes.rnaseq <- function(IDs = "", gene.identifiers = "ensembl") {
    
    # if(!require('biomaRt')) { stop('Please install the R package biomaRt to map the gene
    # identifiers to gene symbols!') } require('biomaRt')
    if (gene.identifiers == "ensembl") {
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))
        ann.df <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol", 
                                                                    "entrezgene_id"), values = IDs, mart = mart)
        return(ann.df)
    }
    if (gene.identifiers == "ucsc") {
        # Convert from UCSC IDs to Gene Name and Description
        ann.df <- getBM(attributes = c("ucsc", "hgnc_symbol", "entrezgene_id"), 
                        filters = "ucsc", values = IDs, mart = mart)  #, uniqueRows=T)
        return(ann.df)
    }
    if (gene.identifiers == "refseq") {
        # Convert RefSEq IDs
        ann.df <- getBM(attributes = c("refseq_mrna", "hgnc_symbol", "entrezgene_id"), 
                        filters = "refseq_mrna", values = IDs, mart = mart)  #, uniqueRows=T)
        return(ann.df)
    }
}

.get.genes.rnaseq2 <- function(IDs = "", gene.identifiers = "ensembl") {
    
    # if(!require('biomaRt')) { stop('Please install the R package biomaRt to map the gene
    # identifiers to gene symbols!') } require('biomaRt')
    if (gene.identifiers == "ensembl") {
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))
        ann.df <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol", 
                                                                    "entrezgene_id", "description"), values = IDs, mart = mart)
        return(ann.df)
    }
    if (gene.identifiers == "ucsc") {
        # Convert from UCSC IDs to Gene Name and Description
        ann.df <- getBM(attributes = c("ucsc", "hgnc_symbol", "entrezgene_id", "description"), 
                        filters = "ucsc", values = IDs, mart = mart)  #, uniqueRows=T)
        return(ann.df)
    }
    if (gene.identifiers == "refseq") {
        # Convert RefSEq IDs
        ann.df <- getBM(attributes = c("refseq_mrna", "hgnc_symbol", "entrezgene_id", "description"), 
                        filters = "refseq_mrna", values = IDs, mart = mart)  #, uniqueRows=T)
        return(ann.df)
    }
}



## Function to get marker genes parameters: data.mat: RNA-Seq gene expression matrix with
## genes corresponding to rows and samples corresponding to columns.  samples2compare: a
## character vector with the sample names to be compared, default is to compare all annotate:
## a boolean value indicating if an annotation (mapping gene identifiers to gene symbols) is
## desired, default is TRUE gene.ids.type: type of the used gene identifiers, the following
## gene identifiers are supported: ensembl, refseq and ucsc gene ids score.cutoff: an integer
## value to filter the resulted markers output: a list with the corresponding markers to each
## sample type

getMarkerGenes.rnaseq <- function(data.mat, class.vec = colnames(data.mat), samples2compare = "all", annotate = FALSE, gene.ids.type = "ensembl", 
                                  score.cutoff = 1) {
    if(length(class.vec) == dim(data.mat)[2]){
        
        colnames(data.mat) <- class.vec
        
    }else{
        stop("class.vec should have the same length as the number of columns of data.mat!!")
    }
    
    if (length(samples2compare) > 1) {
        if (any(!samples2compare %in% colnames(data.mat))) {
            stop("The samples to compare should be included in the data matrix!!")
        } else {
            ii <- which(colnames(data.mat) %in% samples2compare)
            data.mat <- data.mat[, ii]
        }
    }
    ## remove all 0 rows
    x <- data.mat
    data.mat <- x[!apply(x == 0, 1, all), , drop = FALSE]
    markers.list <- list()
    rep.vec <- table(colnames(data.mat))
    cat("Detecting marker genes...\n")
    #res.list <- apply(data.mat, 1, .isMarker.rnaseq, rep.vec = rep.vec)
    #print(res.list)

    res.list <- lapply(rownames(data.mat), function(gene) {
        res <- .isMarker.rnaseq(data.mat[gene, ], rep.vec)
        if (!is.null(res)) {
            names(res) <- c("gene", "score")
            return(res)
        } else {
            return(NULL)
        }
    })
    
    # remove NULLs
    valid.idx <- !sapply(res.list, is.null)
    res.list <- res.list[valid.idx]
    
    # assign correct names after filtering
    names(res.list) <- rownames(data.mat)[valid.idx]
    
    
    
    
    if(length(res.list) > 0) {
        mar.len <- length(unlist(res.list))
        samples.vec <- unname(unlist(res.list))[seq(1, mar.len, 2)]
        scores.vec <- round(as.numeric(unname(unlist(res.list))[seq(2, mar.len, 2)]), digits = 10)
        markers.vec <- names(res.list)
        #markers.vec <- sapply(res.list, function(x) x["gene"])
        #markers.vec <- names(res.list)
        #if (is.null(markers.vec) || any(markers.vec == "")) {
        #    markers.vec <- names(data.mat)[which(!sapply(res.list, is.null))]
        #    }
        
        
        names(scores.vec) <- markers.vec
        u.snames <- unique(colnames(data.mat))
        
        if (annotate) {
            if (!gene.ids.type %in% c("ensembl", "refseq", "ucsc")) {
                stop("To map gene identifiers to gene symbols, only the following gene identifiers are supported: ensembl, refseq and ucsc gene ids!")
            }
            cat("Mapping gene IDs to gene symbols...\n")
            annot.df <- .get.genes.rnaseq(IDs = rownames(data.mat), gene.identifiers = gene.ids.type)
            for (i in seq_along(u.snames)) {
                inds <- which(samples.vec == u.snames[i])
                sort.scores <- sort(scores.vec[inds])
                sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
                genes.inds <- match(names(sort.scores), annot.df[, "ensembl_gene_id"])
                markers.list[[i]] <- paste(names(sort.scores), annot.df[genes.inds, "hgnc_symbol"], 
                                           annot.df[genes.inds, "entrezgene_id"], unname(sort.scores), sep = " : ")
            }
            
        } else {
            
            
            for (i in seq_along(u.snames)) {
                inds <- which(samples.vec == u.snames[i])
                sort.scores <- sort(scores.vec[inds])
                sort.scores <- sort.scores[which(sort.scores <= score.cutoff)]
                markers.list[[i]] <- paste(names(sort.scores), unname(sort.scores), sep = " : ")
            }
        }
        names(markers.list) <- paste(u.snames, "markers", sep = "_")
        cat("Done! \n")
        return(markers.list)
    } else{
        cat("No markers found!!\n This method is designed to be used with few replicates for each sample type!!\n
Please try with fewer replicates!!")
    }
}

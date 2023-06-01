#!/usr/bin/env Rscript

## script to collect devonvolution results from diffrent organism parts accross an atlas
## experiment and create final output. If deconvolution runs was not succesful for 
## an organism part, Cell type column will be empty and sd and proportin will be "0"

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

sdrf = args[1]
accession = args[2]
tissue =  args[3]
sc_reference = args[4]
output = args[5]

if(length(args) != 5) {
   stop("Not correct number of arguments. Please supply five arguments")
}

sdrf = read.csv(sdrf, sep = "\t", row.names = NULL)
# convert "Comment.ENA_RUN." to "Comment..ENA_RUN."
colnames(sdrf) = gsub("[[:punct:]]+", '.', colnames(sdrf))
# remove rows with duplicated run ids in "Comment.ENA_RUN."
sdrf = sdrf[!duplicated(sdrf[ , "Comment.ENA.RUN."]),]
rownames(sdrf) = sdrf$Comment.ENA.RUN.
sdrf$Comment.ENA.RUN. = NULL

#in case whitespace in tissue name was replaced by '_' revert that
tissue_name = gsub('_', ' ', tissue)
# avoid problems if there is '(' or ')' in tissue name
tissue_pattern = gsub("\\(", "\\\\(", tissue)
tissue_pattern = gsub('\\)', '\\\\)', tissue_pattern)

# initialize output file if it does not exist yet
if (!file.exists(output)) {
  write.table(data.frame(ENA_RUN = character(),
                        CL_term = character(),
                        sd = numeric(),
                        proportion = numeric(),
                        organism_part = character(),
                        sc_reference = character()),
                        output, sep = "\t", row.names = FALSE, col.names = TRUE)
}
# check if all three deconvolution results are there
filenames <- paste0('Output/',accession , '/', list.files(paste0("Output/", accession), pattern=paste0(accession,'-', tissue_pattern)))
# check if reference for deconvolution was found... 
if (length(filenames) != 3){ #...if not just append rund ids
    # get the run ids from the runs were we dont have deconvolution results
    runs = rownames(sdrf[sdrf$Characteristics.organism.part. == tissue_name, ])
    MergedDF = data.frame(ENA_RUN =runs,
                          CL_term = character(length(runs)),
                          sd = numeric(length(runs)),
                          proportion = numeric(length(runs)),
                          organism_part = tissue,
                          sc_reference = 'no reference found/deconvolution failed')
    write.table(MergedDF, output, sep = "\t", col.names = FALSE, append = TRUE, row.names = FALSE)
} else { #...if yes, get the average proportions and sds
  # populate file with results
  files <- lapply(filenames, readRDS)
  # calculare pearson corr for each RUN between the three methods
  # to see if the difference between methods is too big
  mean_vector = c()
  for (i in 1:dim(files[[1]])[2]){
    cor_mat = cor(data.frame(files[[1]][,i], files[[2]][,i], files[[3]][,i]))
    # set diagonal to NA
    diag(cor_mat) <- NA
    mean_value <- mean(cor_mat, na.rm = TRUE)
    mean_vector[i] = (mean_value)
  }
  
  print(paste('mean correaltion between the three methods is', (mean(mean_vector))))
  if (mean(mean_vector) < 0.7) {
    print('corelation for methods lower than 0.7')
    # get the run ids from the runs were we dont have deconvolution results
    runs = rownames(sdrf[sdrf$Characteristics.organism.part. == tissue_name, ])
    MergedDF = data.frame(ENA_RUN = runs,
                          CL_term = character(length(runs)),
                          sd = numeric(length(runs)),
                          proportion = numeric(length(runs)),
                          organism_part = tissue,
                          sc_reference = 'no reference found/deconvolution failed')
    write.table(MergedDF, output, sep = "\t", col.names = FALSE, append = TRUE, row.names = FALSE)
  } else {
     # add proportions per run and cell type
     sum_props = Reduce('+',files)
     # get average of the three results
     prop = sum_props/length(files)

     # get standard derivation of the three measurements we have per run and cell type
     vec <- unlist(files, use.names = TRUE)
     DIM <- dim(files[[1]])
     n <- length(files)
     list.mean <- tapply(vec, rep(1:prod(DIM),times = n), mean)
     attr(list.mean, "dim") <- DIM
     list.mean <- as.data.frame(list.mean)

     list.sd <- tapply(vec, rep(1:prod(DIM),times = n), sd)
     attr(list.sd, "dim") <- DIM
     list.sd <- as.data.frame(list.sd)
     colnames(list.sd) = colnames(files[[1]])

     sd_norm = rowMeans(list.sd)/rowMeans(list.mean)
     rownames(list.sd) = rownames(prop)

     list.sd$CL_term = rownames(prop)

     # pivot table
     prop = as.data.frame(prop)
     prop$CL_term = rownames(prop)
     pivoted_mean = pivot_longer(prop, cols = !CL_term, names_to = "ENA_RUN", values_to = "proportion")
     pivoted_sd = pivot_longer(list.sd, cols = !CL_term, names_to = "ENA_RUN", values_to = "sd")

     MergedDF <- merge(pivoted_sd, pivoted_mean, by =c("ENA_RUN", 'CL_term'))
     MergedDF$organism_part = tissue
     MergedDF$sc_reference = sc_reference
  
     write.table(MergedDF, output, sep = "\t", col.names = FALSE, append = TRUE, row.names = FALSE)
     }  
  
}

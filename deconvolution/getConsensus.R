#!/usr/bin/env Rscript

## get consensus of proportions of selected methods
#


suppressMessages(library(dplyr))
#install.packages("Polychrome", repos='http://cran.us.r-project.org')
suppressMessages(library(Polychrome))
suppressMessages(library(scOntoMatch))
suppressMessages(library(plyr))


CiberBarFrazer <- function(ciber, colors, main0, legy){
    ## function to get barplot per sample
    ## from https://github.com/mkrdonovan/gtex_deconvolution
    par(mar = c(2, 4, 2, 0.5))
    
    #rownames(ciber) <- ciber$Input.Sample
    #ciber$Input.Sample <- NULL
    #ciber = ciber[, seq(1, (ncol(ciber) - 3))]
    nsamples =nrow(ciber)
    ciber = as.data.frame(t(ciber) * 100)
    
    ciber$color = colors[ match(colors$celltypes, rownames(ciber)),  "color"]
    #print(ciber)
    barplot(as.matrix(ciber[, seq(1, (ncol(ciber) - 1))]), las = 2, col = ciber$color, border=NA, names.arg = rep("", ncol(ciber) - 1), 
           ylab = "Fraction clusters", main = main0, space=0, cex.main = 2)

    text(nsamples * .05, 90, paste("n = ", nsamples, sep = ""), pos = 4, cex = 1)

    legend(nsamples * .05, legy, gsub("_", " ", colors$celltypes), bty = "n",
           pch = rep(22, nrow(colors)),
           pt.cex = rep(4, nrow(colors)),
           pt.bg = colors$color,
           y.intersp = 1, cex = 1.1
          )
}
suppressMessages({
suppressWarnings({ 
args <- commandArgs(trailingOnly = TRUE)
sampleName = args[1]
tissue = args[2]


tissue = sub("[(]", "[(]", tissue)
tissue = sub("[)]", "[)]", tissue)
tissue = sub("[']", "[']", tissue)
#Get files as list 
#only select DWLS, EpiDISH and FARDEEP proportions for consensus
#filenames <- list.files("Output/", pattern="DWLS|EpiDISH|FARDEEP") %>% lapply(., function(x) { x <- paste0("Output/", x) }) %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) %>% .[lengths(.)!=0] %>% as.character(.)
filenames <- paste0('Output/',list.files("Output/", pattern=paste0(sampleName,'-', tissue)))
print(filenames)
files <- lapply(filenames, readRDS)
names(files) <- lapply(filenames, function(x) gsub("Output/.*res_(.*).rds", "\\1", x))

#Reorder files in case this was without proportions
#files <- lapply(files, function(x) x[order(match(rownames(x), rownames(files[[1]]))),])


#get average proportions per sample and celltype accros the three methods


sum_props = Reduce('+',files)
prop = sum_props/length(files)
#Save consensus proportions
#saveRDS(prop, paste0("Consensus/", sampleName,"/" ,sampleName, "_", tissue, "_consensus.rds"))


#get sd
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

cl_ids = (rownames(prop))

prop = data.frame(t(prop))

#get cell type labels from ids for plot


#read in UBERON ontologys
obo_file = 'deconvolution-analysis/files/cl-basic.obo' #download here: http://purl.obolibrary.org/obo/cl/basic.obo
propagate_relationships = c('is_a', 'part_of', 'relationship: part_of', 'synonym')
# create UBRERON ontology object
ont <- ontologyIndex::get_OBO(obo_file, propagate_relationships =propagate_relationships, extract_tags = 'everything')


cl_terms = getOntologyName(ont = ont, cl_ids)
#ensure right order of colnames
cl_terms = mapvalues(cl_ids, names(cl_terms), unname(cl_terms))

#add sd to columnames 
colnames(prop) = paste0(cl_terms, '(sd=', 100* round(rowMeans(list.sd),3), ')')
#order to have largest proportions first
top = names(prop[order(colSums(prop), decreasing = T)][1])
second_top = names(prop[order(colSums(prop), decreasing = T)][2])
third_top = names(prop[order(colSums(prop), decreasing = T)][3])
color = data.frame(celltypes = colnames(prop),
                   name      =  colnames(prop),
                   color     =  sky.colors(ncol(prop)))

png(filename = paste0("ConsensusPlot/", sampleName,"/",sampleName, "_", tissue, "_consensus.png"), width = 500*3, height = 500*3, res=300)
CiberBarFrazer(prop[order(-prop[top], -prop[second_top], -prop[third_top]),], color, strsplit(sampleName, '_')[[1]][3], 90)
dev.off()



})
})






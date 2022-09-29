#!/usr/bin/Rscript
source("config")


source(file.path(codeDir, "graph.R"))

########################################################################
# deal with town-cell intersections-------------------------------------
########################################################################


load(file.path(dataDir, paste0('intersection_eastern_', productVersion, '.Rda')))

maxPossibleCells <- 70
townCellOverlap <- matrix(0, nTowns, maxPossibleCells)
townCellIds <- matrix(1, nTowns, maxPossibleCells)

for(k in 1:nTowns) {
  tmp <- inter[ , , k]
  ids <- which(tmp > 0)
  nMatch <- length(ids)
  if(nMatch > 0) {
  townCellIds[k, 1:nMatch] <- ids
  townCellOverlap[k, 1:nMatch] <- tmp[ids]
  }
}

tmp <- colSums(townCellOverlap)

townCellOverlap <- townCellOverlap[ , which(tmp > 0)]
townCellIds <- townCellIds[ , which(tmp > 0)]

########################################################################
# read graph for grid --------------------------------------------------
########################################################################


source(file.path(codeDir, 'set_domain.R'))

m1 <- length(easternDomainX)
m2 <- length(easternDomainY)

type <- nbhdStructure
substring(type, 1 ,1) = toupper(substring(type, 1, 1))
fns <- rep("", 2)
fns[1] <- paste('graph', type, '-',  m1, 'x', m2, '.csv', sep='')
fns[2] <- paste('graphCats', type, '-', m1, 'x', m2, '.csv', sep='')

if(!file.exists(file.path(dataDir, fns[1])) || (nbhdStructure != 'bin' && !file.exists(file.path(dataDir, fns[2])))) {
  fns <- graphCreate(m1, m2, type = nbhdStructure, dir = dataDir, fn = fns[1], fn_cats = fns[2])
} 
nbhd <- graphRead(fns[1], fns[2], m1, m2, type = nbhdStructure, dir = dataDir)


########################################################################
# read data ------------------------------------------------------------
########################################################################

easternDataDir <- "eastern"

fn <- file.path(dataDir, paste0(easternVersionID, '-', easternVersion, '.csv'))
data <- read.csv(fn)
data <- data[ , 9:ncol(data)]
## remove the "N.5.2" from the taxa names
names(data) <- sapply(strsplit(names(data), "\\."), `[[`, 1)

cat(paste0("Read ", nrow(data), " rows from ", fn, ", with field names: "))
cat(names(data), sep = ',')
cat("\n")

data[is.na(data)] <- 0
data <- round(data) # there is one count of 3.5

########################################################################
# subset and manipulate taxa ------------------------------------
########################################################################

taxaInfo <- read.csv('taxon_translation.csv')

taxaOther <- taxaInfo[["original"]][taxaInfo[["modelled"]] == "OTHER"]
taxaUse <- unique(taxaInfo[["original"]][!(taxaInfo[["modelled"]] == "OTHER")])


cat("Grouping the following taxa into an 'Other' category:")
cat(taxaOther, sep = ',')
cat("\n")

allTaxa <- c(taxaOther, taxaUse)
zeros <- allTaxa[!(allTaxa %in% names(data))]
data[ , zeros] <- 0

other <- rowSums(data[ , taxaOther])

data <- data[ , taxaUse]
  
data <- data[ , order(names(data))]

data$"OTHER" <- other

taxaUse <- names(data)
nTaxa <- length(taxaUse)

taxa <- data.frame(taxonID = 1:nTaxa, taxonName = taxaUse)


cat("Using the following", nTaxa, "taxa, with", sum(data), "trees, of which", round(100*sum(data$"OTHER")/sum(data), 2), "% are in the 'other' category:")
cat(taxaUse, sep = ",")
cat("\n")

########################################################################
# create data objects for MCMC fitting ---------------------------------
########################################################################


nTrees <- rowSums(data)

town <- rep(seq_len(nTowns), times = nTrees)

tmp <- c(t(as.matrix(data)))

taxon <- rep(rep(1:nTaxa, nTowns), times = tmp)

data <- data.frame(taxon = taxon, town = town)

save(data, townCellOverlap, townCellIds, taxa, nbhd, m1, m2, nTaxa, nTowns, file = file.path(dataDir, paste0('data_eastern_', productVersion, '.Rda')))

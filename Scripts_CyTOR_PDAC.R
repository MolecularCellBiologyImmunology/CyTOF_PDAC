library(CATALYST)
library(flowCore)
library(FlowSOM)

##### Bead Normalization and Debarcoding #####
dir <- #Directory#
files <- list.files(dir, pattern = "fcs$")
ff <- flowCore::read.flowSet(file.path(dir, files))
sce <- prepData(ff)

res <- normCytof(x=sce, beads = "dvs", k=50, plot=T)

sce1 <- res$data
plotScatter(sce, chs = c("Time", "CD45"))
plotScatter(sce1, chs = c("Time", "CD45"))


data("sample_key")
sce2 <- assignPrelim(sce1, sample_key)
sce2 <- estCutoffs(sce2)
plotYields(sce2)

sce2 <- applyCutoffs(sce2)

plotEvents(sce2, which = c(0, "D1"), n = 50)

fs <- sce2fcs(sce2, split_by = "bc_id")

dir <- #Output Directory#
ids <- fsApply(fs, identifier)
for (id in ids) {
  ff <- fs[[id]]
  fn <- sprintf("week1_sample%s.fcs", id)
  fn <- file.path(dir, fn)
  write.FCS(ff, fn)
}


##### Clsutering & Dimensional Reduction #####

list.raw <- list.files(#Folder#, full.names = T, pattern = ".fcs)
panel <- read.table("Panel_CyTOF.txt")
md <- read.table("Metadata_CyTOF.txt")

fs<- read.flowSet(list.raw)
sce <- prepData(fs, panel = panel, md = md)
sce <- cluster(sce, features = "type", xdim = 20, ydim = 20, maxK = 100, verbose = T, seed = 1234)
sce <- runDR(sce, dr = "UMAP", cells = 250, features = "type")
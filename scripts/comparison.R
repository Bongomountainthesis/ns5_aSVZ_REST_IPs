
stringsAsFacctors=FALSE

ns5 <- read.csv(file="data/ns5_nearest_or_overlapping_peak_to_gene_TSS.csv")
asvz <- read.csv(file="data/aSVZ_nearest_or_overlapping_peak_to_gene.TSS.csv")

dim(ns5)
#[1] 1301   20
 dim(asvz)
#[1] 755  20

#cut off all peaks at neg10log10pval >=60

ns5 <- ns5[which(ns5[,"neg10log10pVal"]>=60),]
asvz <- asvz[which(asvz[,"neg10log10pVal"]>=60),]

dim(ns5)
#[1] 1126  20
dim(asvz)
#[1] 477  20

##asvz unique genes

asvz.ids <- unique(asvz[,"EnsemblID"])
ns5.ids <- unique(ns5[,"EnsemblID"])

shared <- intersect(asvz.ids,ns5.ids)

length(shared)
# 373


ns5.unique <- ns5.ids[which(!(ns5.ids %in% asvz.ids))]
ns5.unique.res <- ns5[which(ns5[,"EnsemblID"] %in% ns5.unique),]

dim(ns5.unique.res)
#661

asvz.unique <- asvz.ids[which(!(asvz.ids %in% ns5.ids))]
asvz.unique.res <- asvz[which(asvz[,"EnsemblID"] %in% asvz.unique),]

dim(asvz.unique.res)
#76

shared.ns5 <- ns5[which(ns5[,"EnsemblID"] %in% shared),]
shared.asvz <- asvz[which(asvz[,"EnsemblID"] %in% shared),]

shared.merged <- merge(shared.ns5, shared.asvz, by.x = "EnsemblID", by.y = "EnsemblID",suffixes = c("_NS5","_aSVZ")) 
tidy <- shared.merged[,c(1,19,20,12,13,14,3,4,5,6,7,8,9,10,11,15,16,17,18,22,23,24,25,26,27,28,29,30,34,35,36,37)]

write.csv(tidy, file = "results/peaks_shared_in_asvz_and_ns5.csv")
write.csv(ns5.unique.res,file="results/peaks_unique_to_ns5.csv")
write.csv(asvz.unique.res,file="results/peaks_unique_to_asvz.csv")



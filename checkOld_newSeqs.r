# this script checks the old RData and new RData files in seqData to check if anything has changed
# between the two, and reports on which ones did change.

require(Biostrings)

orgRDat <- file.path("seqData", "orgSequences", dir("seqData/orgSequences/"))
orgRDat <- grep("RData", orgRDat, value=T)

newRDat <- file.path("seqData", dir("seqData"))
newRDat <- grep("RData", newRDat, value=T)

nDat <- length(orgRDat)
allInfo <- vector('list', nDat)
for (iDat in 1:nDat){
	orgE <- new.env()
	load(orgRDat[iDat], env=orgE)
	
	newE <- new.env()
	load(newRDat[iDat], env=newE)
	
	isPM <- all.equal(orgE$pm.mm.dat$pm, newE$pm.mm.dat$pm)
	isMM <- all.equal(orgE$pm.mm.dat$mm, newE$pm.mm.dat$mm)
	
	allInfo[[iDat]] <- list(pm=isPM, mm=isMM)
}

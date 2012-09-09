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
	
	orgPM <- as.character(orgE$pm.mm.dat$pm)
	orgPM <- orgPM[order(names(orgPM))]
	
	orgMM <- as.character(orgE$pm.mm.dat$mm)
	orgMM <- orgMM[order(names(orgMM))]
	
	newPM <- as.character(newE$pm.mm.dat$pm)
	newPM <- newPM[order(names(newPM))]
	
	newMM <- as.character(newE$pm.mm.dat$mm)
	newMM <- newMM[order(names(newMM))]
	
	isPM <- all.equal(orgPM, newPM)
	isMM <- all.equal(orgMM, newMM)
	
	allInfo[[iDat]] <- list(pm=isPM, mm=isMM)
	rm(orgE, newE, orgPM, newPM, orgMM, newMM)
}

names(allInfo) <- basename(orgRDat)
allInfo

evec.lowmem <- read.table("eigenvectors_lowhmem.txt",
   header=FALSE, sep="", stringsAsFactors=FALSE)
evec.highmem <- read.table("eigenvectors_highmem.txt",
   header=FALSE, sep="", stringsAsFactors=FALSE)
evec.plink <- read.table("plink.eigenvec",
   header=FALSE, sep="", stringsAsFactors=FALSE)[, -(1:2)]

r1 <- diag(cor(evec.lowmem, evec.highmem))
stopifnot(all(abs(r1) > 0.99999))

r2 <- diag(cor(evec.lowmem, evec.plink))
stopifnot(all(abs(r1) > 0.99999))



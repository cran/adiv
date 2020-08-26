plot.evoCA <- function(x, xaxis = 1, yaxis = 2, graph = FALSE, ...){

if(ncol(x$co)<max(xaxis, yaxis)) stop(paste(ncol(x$co), " axes only have been saved in evoCA, please change xaxis, yaxis, or run again function evoCA to include more axes"))
 
phyl.phylo <- attributes(x)$phy

tre4 <- as(phyl.phylo, "phylo4")
namesnodes <- nodeLabels(tre4)
pathnodes <- ancestors(tre4, namesnodes, type="ALL")
nodesage <- unlist(lapply(pathnodes, function(x) 
sumEdgeLength(tre4, names(x))))
names(nodesage) <- names(pathnodes)
namestips <- tipLabels(tre4)
pathtips <- ancestors(tre4, namestips, type="ALL")
tipsage <- unlist(lapply(pathtips, function(x) 
sumEdgeLength(tre4, names(x))))
names(tipsage) <- names(pathtips)
ages <- c(tipsage, nodesage)[rownames(x$co)]
ages <- ages - max(ages) 

coordtoutes <- cbind.data.frame(x$co[, c(xaxis, yaxis)], ages)
names(coordtoutes) <- c("X","Y","Z")
coord1 <- coordtoutes[namestips, ]

listchildren <- sapply(namesnodes, function(u) children(tre4, u))
names(listchildren) <- namesnodes

options(rgl.useNULL = !graph)
plot3d(coord1$X, coord1$Y, coord1$Z, xlim=range(coordtoutes$X), ylim=range(coordtoutes$Y), zlim=range(coordtoutes$Z), pch=19, xlab=paste("evoCA Axis", xaxis, sep=" "), ylab=paste("evoCA Axis", yaxis, sep=" "), zlab="Time of evolution", ...)

for(i in 1: length(listchildren)){
for(j in 1: length(listchildren[[i]])){
nom1 <- names(listchildren)[i]
nom2 <- names((listchildren[[i]])[j])
plot3d(coordtoutes[c(nom1, nom2),"X"], coordtoutes[c(nom1, nom2),"Y"], coordtoutes[c(nom1, nom2),"Z"], add=TRUE, type="l")
}
}

}

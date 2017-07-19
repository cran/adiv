dissRicotta <- function(comm, dis)
{
     if(any(dis>1)){
        warning("Phylogenetic dis are not in the range 0-1. They have been normalized by the maximum")
        dis <- dis/max(dis)
    }
if(any(!colnames(comm) %in%rownames(dis))) stop("At least one species in the matrix of comm is missing in the matrix of dis") 
if(any(!colnames(comm) %in%colnames(dis))) stop("At least one species in the matrix of comm is missing in the matrix of dis") 
dis <- dis[colnames(comm), colnames(comm)]
    dataset <- t(comm)
    similarities<-1-as.matrix(dis)
    total <- colSums(dataset)
    rel_abu<- sweep(dataset, 2, total, "/")
    num.plot<-dim(dataset)[2]
    names<-list(colnames(dataset), colnames(dataset))
    dist.matrix<-matrix(0, nrow=num.plot, ncol=num.plot, dimnames=names)
    for (i in 2:num.plot) {
    for (j in 1:(i-1)) {
        mat_folk<-similarities*rel_abu[,i]
        mat_folk2<- similarities*rel_abu[,j]
        xxx<-colSums(mat_folk2)
        yyy<-colSums(mat_folk)
        sub<-abs(xxx-yyy)
        index<-sum(sub)/sum(xxx+yyy)
        dist.matrix[i,j]<-index
    }
    }
    dist.matrix <- dist.matrix + t(dist.matrix)
    semi_matrix<-as.dist(dist.matrix)
    return(semi_matrix)
}

specieseve <- function(comm, method = "full", tol = 1e-8){

if(any(!method%in%c("GiniSimpson", "Simpson", "Shannon", "Heip", "McIntosh", "SmithWilson", "full"))) stop("Your choice for method is not available")
if("full"%in%method) method <- c("GiniSimpson", "Simpson", "Shannon", "Heip", "McIntosh", "SmithWilson")

if(any(comm < (-tol))) stop("Abundance entry in comm must be nonnegative")

comm[comm<tol] <- 0

if(all(rowSums(comm)<tol)) stop("All communities are empty")

if(any(rowSums(comm)<tol)) warning("Empty communities were discarded")

comm <- comm[rowSums(comm)>tol, ]

FUNshannoneq <- function(v){
    if(length(v[v>0])==1) return(0)
    else{
    v <- v[v>0]
    return(-sum(v/sum(v)*log(v/sum(v)))/log(length(v)))
    }
}
FUNSW <- function(v){
    if(length(v[v>0])==1) return(0)
    else{
    v <- v[v>0]
    Evar <- 1-
( 2 / pi *
atan(
sum((log(v)-sum(log(v)/length(v)))^2) /length(v) 
) 
) 
    return(Evar)
    }
}
FUNshannon <- function(v){
    if(length(v[v>0])==1) return(0)
    else{
    v <- v[v>0]
    return(-sum(v/sum(v)*log(v/sum(v))))
    }
}
RES <- matrix(0, nrow(comm), length(method))
rownames(RES) <- rownames(comm)
colnames(RES) <- method
for(i in 1:length(method)){
    if(method[i]=="GiniSimpson")
    RES[,i] <- apply(comm, 1, 
        function(x) (1-sum((x/sum(x))^2)) * 
        length(x[x>0]) / (length(x[x>0])-1))
    else if(method[i]=="Simpson")
    RES[,i] <- apply(comm, 1, function(x) 1/sum((x/sum(x))^2)/length(x[x>0]))
    else if(method[i]=="Shannon")
    RES[,i] <- apply(comm, 1, FUNshannoneq)
    else if(method[i]=="Heip")
    RES[,i] <- apply(comm, 1, function(x) (exp(FUNshannon(x))-1)/(length(x[x>0])-1))
    else if(method[i]=="McIntosh")
    RES[,i] <- apply(comm, 1, function(x) (sum(x)-sqrt(sum(x^2)))/(sum(x)-sum(x)/sqrt(length(x[x>0]))))
    else if(method[i]=="SmithWilson")
    RES[,i] <- apply( comm, 1, FUNSW)
}
return(RES)

}

speciesdiv <- function(comm, method = "full", tol = 1e-8){

if(any(!method%in%c("richness", "GiniSimpson", "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh", "full"))) stop("Your choice for method is not available")
if("full"%in%method) method <- c("richness", "GiniSimpson", "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh")

if(any(comm < (-tol))) stop("Abundance entry in comm must be nonnegative")

comm[comm<tol] <- 0

if(all(rowSums(comm)<tol)) stop("All communities are empty")

if(any(rowSums(comm)<tol)) warning("Empty communities were discarded")

comm <- comm[rowSums(comm)>tol, , drop = FALSE]

FUNshannon <- function(v){
    if(length(v[v>0])==1) return(0)
    else{
    v <- v[v>0]
    return(as.vector(-sum(v/sum(v)*log(v/sum(v)))))
    }
}
RES <- matrix(0, nrow(comm), length(method))
rownames(RES) <- rownames(comm)
colnames(RES) <- method
for(i in 1:length(method)){
    if(method[i]=="richness")
    RES[,i] <- apply(comm, 1, function(x) length(x[x>0]))
    if(method[i]=="GiniSimpson")
    RES[,i] <- apply(comm, 1, function(x) 1-sum((x/sum(x))^2))
    if(method[i]=="Simpson")
    RES[,i] <- apply(comm, 1, function(x) 1/sum((x/sum(x))^2))
    if(method[i]=="Shannon")
    RES[,i] <- apply(comm, 1, FUNshannon)
    if(method[i]=="Margalef")
    RES[,i] <- apply(comm, 1, function(x) (length(x[x>0])-1)/log(sum(x)))
    if(method[i]=="Menhinick")
    RES[,i] <- apply(comm, 1, function(x) length(x[x>0])/sqrt(sum(x)))
    if(method[i]=="McIntosh")
    RES[,i] <- apply(comm, 1, function(x) (sum(x)-sqrt(sum(x^2)))/(sum(x)-sqrt(sum(x))))
}
return(RES)

}


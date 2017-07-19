rare_Rao <-
function(comm, dis = NULL, sim = TRUE, resampling = 999, formula = c("QE", "EDI"))
{
if (!inherits(comm, "matrix") & !inherits(comm, "data.frame")) stop("Non convenient comm")
comm <- t(comm)
if (any(comm<0))
stop("Negative value in comm")
if (any(apply(comm, 2, sum) < 1e-16))
stop("Empty communities")
if(!formula[1]%in%c("QE","EDI")) stop("formula can be either QE or EDI")
formula <- formula[1]
if (!is.null(dis)) {
    if (!inherits(dis, "dist")) 
        stop("Object of class 'dist' expected for distance")
    dis <- as.matrix(dis)
    if (nrow(comm) != nrow(dis)) 
        stop("Non convenient comm")
    dis <- as.dist(dis)
    if(formula=="EDI") dis <- dis^2/2
    if (!is.euclid(sqrt(dis)))
warning("Squared Euclidean or Euclidean property expected for dis")
}
if (is.null(dis)){
    dis <- as.dist((matrix(1, nrow(comm), nrow(comm)) - diag(rep(1, nrow(comm)))))
}
N <- ncol(comm)
rel_abu <- sweep(comm, 2, colSums(comm), "/")
if(sim){
p <- resampling
result <- matrix(0, p, N)
for(j in 1:N){
for(i in 1:p) {
if(j==1) newplot <- as.vector(rel_abu[, sample(1:ncol(rel_abu), j)])
else{
subsample <- rel_abu[, sample(1:ncol(rel_abu), j)]
newplot <- as.vector(apply(subsample, 1, mean))
}
result[i, j] <- (t(newplot) %*% (as.matrix(dis)) %*% newplot)
}
}
aver <- apply(result, 2, mean)
IC_plus <- aver + (1.96*(sd(result)/sqrt(p)))
IC_neg <- aver - (1.96*(sd(result)/sqrt(p)))
rrRao <- data.frame(as.matrix(aver), IC_neg, IC_plus)
names(rrRao) <- c("ExpRao", "LeftIC", "RightIC")
return(rrRao)
}
else{
a <- mean(diag(t(rel_abu)%*%as.matrix(dis)%*%as.matrix(rel_abu)))
f2 <- apply(rel_abu, 1, sum)/sum(rel_abu)
g <- t(f2)%*%as.matrix(dis)%*%f2
b <- g - a
rrARao <- cbind.data.frame(sapply(1:N, function(M) return(a + N * (M-1) * b / M / (N-1))))
names(rrARao) <- c("ExpRao")
return(rrARao)
}
}

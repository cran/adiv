qHdiv <-
function(comm, C, q=2){

    if (!inherits(comm, "data.frame") & !inherits(comm, "matrix")) 
        stop("comm must be a data frame or a matrix")
    if (any(comm < 0))
        stop("non-negative values expected in comm")
    comms <- apply(comm, 1, sum)
    if (any(comms == 0))
        stop("column in comm with zero values only")
    fun <- function(x){
        x <- x/sum(x)
        d <- svd(diag(sqrt(x))%*%C%*%diag(sqrt(x)))$d
        d <- d/sum(d)
        divx <- (sum(d^q))^(1/1-q)
        return(divx)
    }
    div <- apply(comm, 1, fun)
    return(div)

}

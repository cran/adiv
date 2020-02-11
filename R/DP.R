DP<- 
function(mtree, comm, height = NULL, diag = FALSE, upper = FALSE, tol = 0.001){

    if(! (inherits(comm, "matrix") | inherits(comm, "data.frame")))
        stop("comm must be a matrix or a data frame")
    ncom <- nrow(comm)
    if (is.null(colnames(comm))) 
        stop("comm must have names for column")
    if (ncom < 2) 
        stop("At least two rows for comm are required")
    if (is.null(colnames(comm))) 
        stop("comm must have names for column")
    TA <- tecAptree(mtree, tol = tol)
    if (any(!colnames(comm) %in% names(TA$list[[1]]))) 
        stop("comm contains tip names that are not available in mtree")
    if (any(comm < 0)) 
        stop("comm should contain nonnegative values")
    if (any(rowSums(comm) == 0)) 
        stop("empty communities should be discarded in comm")
    sp_names <- names(TA$list[[1]])
    comm <- comm[, sp_names]
    FUN_COM <- function(groups){
        COM <- apply(comm, 1, function(x) tapply(x, groups, sum))
        return(t(COM))
    }
    FUN_BC <- function(tab){
        d <- matrix(0, ncom, ncom)
        funBC <- function(x) {
            p <- tab[x[1], ]
            q <- tab[x[2], ]
            ps <- p[(p + q) > 0]
            qs <- q[(p + q) > 0]
            w <- sum(abs(ps - qs))/sum(ps + qs)
            return(w)
        }
        index <- cbind(col(d)[col(d) < row(d)], 
            row(d)[col(d) < row(d)])
        d <- unlist(apply(index, 1, funBC))
        return(d)
    }
    LISTCOM <- lapply(TA$list, FUN_COM)
    LISTd <- lapply(LISTCOM, FUN_BC)
    d <- LISTd[[1]] * TA$plength[1]
    for(i in 2:length(LISTd)) {
        d <- d + LISTd[[i]] * TA$plength[i]
    }
    if(!is.null(height) && is.numeric(height)&& height >sum(TA$plength))
         d <- d / height
    else 
d <- d / sum(TA$plength)
    attr(d, "Size") <- ncom
    attr(d, "Labels") <- rownames(comm)
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- "DP"
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}

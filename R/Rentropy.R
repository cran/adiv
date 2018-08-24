Rentropy <-
function (comm, dis = NULL, scale = FALSE)
{
    if (!inherits(comm, "data.frame") & !inherits(comm, "matrix")) 
        stop("comm must be a data frame or a matrix")
    comm <- t(comm)
    if (any(comm < 0)) 
        stop("Negative value in comm")
    if (!is.null(dis)) {
        if (!inherits(dis, "dist")) 
            stop("Object of class 'dist' expected for dis")
        dis <- as.matrix(dis)
        if (nrow(comm) != nrow(dis)) 
            stop("Species in comm (columns) must be the same as in dis")
        dis <- as.dist(dis)
    }
    if (is.null(dis)){
        dis <- as.dist((matrix(1, nrow(comm), nrow(comm)) - diag(rep(1, nrow(comm)))))
    }
    div <- as.data.frame(rep(0, ncol(comm)))
    names(div) <- "diversity"
    rownames(div) <- names(comm)
    for (i in 1:ncol(comm)) {
        if (sum(comm[, i]) < 1e-16) 
            div[i, ] <- 0
        else div[i, ] <- (t(sqrt(comm[, i]/sum(comm[, i]))) %*% (as.matrix(dis)) %*% 
            sqrt(comm[, i]/sum(comm[, i])))
    }
    if (scale == TRUE) {
        divmax <- eigen(dis)$value[1]
        div <- div/divmax
    }
    rownames(div) <- colnames(comm)
    return(div)
}

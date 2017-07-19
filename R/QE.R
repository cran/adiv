QE <-
function (comm, dis = NULL, formula = c("QE", "EDI"), scale = FALSE) 
{
    if (!inherits(comm, "data.frame") & !inherits(comm, "matrix")) 
        stop("Non convenient comm")
    comm <- t(comm)
    if (any(comm < 0)) 
        stop("Negative value in comm")
    if(!formula[1]%in%c("QE","EDI")) stop("formula can be either QE or EDI")
    formula <- formula[1]
    if (!is.null(dis)) {
        if (!inherits(dis, "dist")) 
            stop("Object of class 'dist' expected for distance")
        dis <- as.matrix(dis)
        if (nrow(comm) != nrow(dis)) 
            stop("Non convenient comm")
        dis <- as.dist(dis)
        if(formula=="QE") dis <- sqrt(2*dis)
    }
    if (is.null(dis)){
        dis <- as.dist((matrix(1, nrow(comm), nrow(comm)) - diag(rep(1, nrow(comm)))) * sqrt(2))
    }
    div <- as.data.frame(rep(0, ncol(comm)))
    names(div) <- "diversity"
    rownames(div) <- names(comm)
    for (i in 1:ncol(comm)) {
        if (sum(comm[, i]) < 1e-16) 
            div[i, ] <- 0
        else div[i, ] <- (t(comm[, i]) %*% (as.matrix(dis)^2) %*% 
            comm[, i])/2/(sum(comm[, i])^2)
    }
    if (scale == TRUE) {
        divmax <- divcmax(dis)$value
        div <- div/divmax
    }
    return(div)
}

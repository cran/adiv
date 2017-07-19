discomQE <-
function (comm, dis = NULL, structures = NULL, formula = c("QE", "EDI")) 
{
    if (!inherits(comm, "data.frame") & !inherits(comm, "matrix")) 
        stop("comm must be a data frame or a matrix")
    comm <- t(comm)
    if (any(comm < 0)) 
        stop("Negative value in comm")
    if (any(apply(comm, 2, sum) < 1e-16)) 
        stop("Empty comm")
    if(!formula[1]%in%c("QE","EDI")) stop("formula can be either QE or EDI")
    formula <- formula[1]
    if (!is.null(dis)) {
        if (!inherits(dis, "dist")) 
            stop("Object of class 'dist' expected for distance")
        dis <- as.matrix(dis)
        if (nrow(comm) != nrow(dis)) 
            stop("Non convenient comm")
        if(formula=="QE") dis <- sqrt(2*dis)
        if (!is.euclid(as.dist(dis))) 
            stop("Euclidean property is expected for distance")
    }
    if (is.null(dis))
        dis <- (matrix(1, nrow(comm), nrow(comm)) - diag(rep(1, 
            nrow(comm)))) * sqrt(2)
    if (!is.null(structures)) {
        if (!inherits(structures, "data.frame")) 
            stop("Non convenient structures")
        m <- match(apply(structures, 2, function(x) length(x)), 
            ncol(comm), 0)
        if (length(m[m == 1]) != ncol(structures)) 
            stop("Non convenient structures")
        m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)), 
            function(x) is.factor(structures[, x])), TRUE, 0)
        if (length(m[m == 1]) != ncol(structures)) 
            stop("Non convenient structures")
    }
    Structutil <- function(dp2, Np, unit) {
        if (!is.null(unit)) {
            modunit <- model.matrix(~-1 + unit)
            sumcol <- apply(Np, 2, sum)
            Ng <- modunit * sumcol
            lesnoms <- levels(unit)
        }
        else {
            Ng <- as.matrix(Np)
            lesnoms <- colnames(Np)
        }
        sumcol <- apply(Ng, 2, sum)
        Lg <- t(t(Ng)/sumcol)
        colnames(Lg) <- lesnoms
        Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
        rownames(Pg) <- lesnoms
        deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% 
            dp2 %*% x))
        ug <- matrix(1, ncol(Lg), 1)
        dg2 <- t(Lg) %*% dp2 %*% Lg - 1/2 * (deltag %*% t(ug) + 
            ug %*% t(deltag))
        colnames(dg2) <- lesnoms
        rownames(dg2) <- lesnoms
        return(list(dg2 = dg2, Ng = Ng, Pg = Pg))
    }
    Diss <- function(dis, nbhaplotypes, comm, structures) {
        structutil <- list(0)
        structutil[[1]] <- Structutil(dp2 = dis, Np = comm, 
            NULL)
        diss <- list(sqrt(as.dist(structutil[[1]]$dg2)))
        if (!is.null(structures)) {
            for (i in 1:length(structures)) {
                structutil[[i + 1]] <- Structutil(structutil[[1]]$dg2, 
                  structutil[[1]]$Ng, structures[, i])
            }
            diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)), 
                function(x) sqrt(as.dist(structutil[[x + 1]]$dg2))))
        }
        return(diss)
    }
    nbhaplotypes <- sum(comm)
    diss <- Diss(dis^2, nbhaplotypes, comm, structures)
    names(diss) <- c("communities", names(structures))
    if (!is.null(structures)) {
        if(formula=="QE") diss <- lapply(diss, function(x) x^2/2)
        return(diss)
    }
    if(formula=="QE") return(diss$comm^2/2)
    if(formula=="EDI") return(diss$comm)
}

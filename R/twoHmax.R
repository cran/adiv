twoHmax <-
function (C, epsilon = 1e-8, smooth = TRUE, comment = FALSE)
{

    nonnegativedefinite <- function(x){
        eigenval <- eigen(x, only.values = TRUE)$values
        nonzeroeigenval <- eigenval[abs(eigenval)>epsilon]
        return(!any(nonzeroeigenval<0))
    }
    if (!inherits(C, "matrix"))
        stop("Matrix expected")
    if (epsilon <= 0)
        stop("epsilon must be positive")
    if (!nonnegativedefinite(C))
        stop("C must be nonnegative definite")
    C2 <- C^2
    n <- dim(C2)[1]
    relax <- 0
    Z <- diag(1/diag(C))%*%C2%*%diag(1/diag(C))
    D2 <- max(Z) - Z
    x0 <- apply(D2, 1, sum)/sum(D2)
    objective0 <- t(x0) %*% D2 %*% x0
    if (comment == TRUE)
        print("evolution of the objective function:")
    xk <- x0
    repeat {
        repeat {
            maxi.temp <- t(xk) %*% D2 %*% xk
            if (comment == TRUE)
                print(as.character(maxi.temp))
            deltaf <- (-2 * D2 %*% xk)
            sature <- (abs(xk) < epsilon)
            if (relax != 0) {
                sature[relax] <- FALSE
                relax <- 0
            }
            yk <- (-deltaf)
            yk[sature] <- 0
            yk[!(sature)] <- yk[!(sature)] - mean(yk[!(sature)])
            if (max(abs(yk)) < epsilon) {
                break
            }
            alpha.max <- as.vector(min(-xk[yk < 0]/yk[yk < 0]))
            alpha.opt <- as.vector(-(t(xk) %*% D2 %*% yk)/(t(yk) %*%
                D2 %*% yk))
            if ((alpha.opt > alpha.max) | (alpha.opt < 0)) {
                alpha <- alpha.max
            }
            else {
                alpha <- alpha.opt
            }
            if (abs(maxi.temp - t(xk + alpha * yk) %*% D2 %*%
                (xk + alpha * yk)) < epsilon) {
                break
            }
            xk <- xk + alpha * yk
        }
        if (prod(!sature) == 1) {
            if (comment == TRUE)
                print("KT")
            break
        }
        vectD2 <- D2 %*% xk
        u <- 2 * (mean(vectD2[!sature]) - vectD2[sature])
        if (min(u) >= 0) {
            if (comment == TRUE)
                print("KT")
            break
        }
        else {
            if (comment == TRUE)
                print("relaxation")
            satu <- (1:n)[sature]
            relax <- satu[u == min(u)]
            relax <- relax[1]
        }
    }
    if (comment == TRUE)
        print(list(objective.init = objective0, objective.final = maxi.temp))
    result <- as.vector(xk, mode = "numeric")
    result[result < epsilon] <- 0
    if(smooth){
	respos <- (1:n)[result>0]
        ressmooth <- (solve(Z[respos, respos])%*%rep(1, length(respos))/sum(solve(Z[respos, respos])))[, 1]
        if(!any(ressmooth < -epsilon)){
            result <- rep(0, n)
            result[respos] <- ressmooth
            result[result < epsilon] <- 0
        }
    }
    result <- result / diag(C)
    result <- result / sum(result)
    result <- as.data.frame(result)
    names(result) <- "pmax"
    result$pmax[result$pmax < epsilon] <- 0
    restot <- list()
    p <- result[, 1]
    restot$value <- as.vector((t(p)%*%diag(C))^2)/as.vector(t(p)%*%C2%*%p)
    restot$vector <- result
    return(restot)

}

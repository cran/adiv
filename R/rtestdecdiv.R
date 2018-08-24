rtestdecdiv <- function(phyl, vecab, dis = NULL, tol = 1e-08, option = 1:5, formula = c("QE", "EDI"), vranking = c("complexity", "droot"), ties.method = c("average", "first", "last", "random", "max", "min"), statistic = 1:3, optiontest = NULL, nrep = 99) {

    ties.method <- ties.method[1]
    vranking <- vranking[1]
    if(!is.null(optiontest)){
    if(!optiontest%in%c("greater", "less", "two-sided")) stop("Incorrect definition of optiontest")
    if(length(optiontest)!=length(statistic)) stop("Incorrect definition of optiontest")
    }
    mydecdiv <- decdiv(phyl = phyl, 
        comm = vecab, dis = dis, 
        tol = tol, option = option, formula = formula)    
    tre4 <- as(attributes(mydecdiv)$phyl, "phylo4")
    nsp <- length(tipLabels(tre4))
    namesnodes <- nodeLabels(tre4)
    nnodes <- length(namesnodes)
    pathnodes <- ancestors(tre4, namesnodes, type="ALL")
    ties.method <- ties.method[1]
    if(!ties.method%in%c("average", "first", "last", "random", "max", "min")) stop("Incorrect definition of parameter ties.method")
    complexity <- function(phy){
        listno <- lapply(namesnodes, function(x) c(x, names(descendants(phy, x, type="all"))))
        names(listno) <- namesnodes
        lisDD <- lapply(namesnodes, function(x) descendants(phy, x, type="children"))
        nbdes <- unlist(lapply(lisDD, function(x)
        prod(1:length(x))))
        names(nbdes) <- namesnodes
        compl <- unlist(lapply(listno, function(x) prod(nbdes[x[x%in%namesnodes]])))
        names(compl) <- namesnodes
        return(compl)
    }

    droot <- function(phy){
        roottab <- unlist(lapply(pathnodes, function(x) 
            sumEdgeLength(phy, names(x))))
        names(roottab) <- names(pathnodes)
        return(roottab)
    }

    vrank <- rank(get(vranking)(tre4), ties.method = ties.method)
    vrank <- vrank[rownames(mydecdiv)]

    funprincipale <- function(vdecdiv){

    stat1 <- function(v){
        v <- v/sum(v)
        return(max(v))
    }
    stat2 <- function(v){
        v <- v/sum(v)
        fun1 <- function(m){
            return(abs(sum(sort(v)[1:m]) - m/nnodes))
        }
        return(max(unlist(lapply(1:nnodes, fun1))))
    }
    stat3 <- function(v){
        funstat3 <- function(vrank1){
            v <- v/sum(v)
            return(sum(rank(vrank1, ties.method = ties.method)*v)/nnodes)
        }
        return(funstat3(vrank))
    }

    methods <- c("stat1", "stat2", "stat3")[statistic]

    statobs <- unlist(sapply(methods, function(x) get(x)(vdecdiv)))
    return(statobs)

    }

    obs <- funprincipale(mydecdiv[, 1])
    
    funsim <- function(i){
        e <- sample(length(vecab))
        ve <- vecab[e]
        names(ve) <- names(vecab)
        if(!is.null(dis)){
           dise <- as.matrix(dis)[e,e]
           rownames(dise) <- colnames(dise) <- attributes(dis)$Labels
           dise <- as.dist(dise)
        }
        simdecdiv <- decdiv(phyl = phyl, 
            comm = ve, dis = dise, tol = tol, option = option,
            formula = formula)
        simi <- funprincipale(simdecdiv[, 1])
        return(simi)  
    }
    sim <- cbind.data.frame(sapply(1:nrep, funsim))
    
    optiondefault <- c("greater", "two-sided", "two-sided")

    if(length(statistic) == 1)
    {
        if(!is.null(optiontest))
            return(as.randtest(obs = obs, sim = sim[, 1], alter = optiontest, call = "rtest.decdiv"))
        else
        return(as.randtest(obs = obs, sim = sim[, 1], alter = optiondefault[statistic], call = "rtest.decdiv"))
    }
    
    if(!is.null(optiontest))
        return(as.krandtest(obs = obs, sim = t(sim), 
           alter = optiontest, call = "rtest.decdiv"))
    else
        return(as.krandtest(obs = obs, sim = t(sim), 
           alter = optiondefault[statistic], 
            call = "rtest.decdiv"))

}

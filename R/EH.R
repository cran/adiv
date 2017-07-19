EH <-
function(phyl, select = NULL)
  ## phyl can be of class phylo, phylo4, phylo4d, phylog
{
    arg.phyl <- .checkphyloarg(phyl)
    phyl <- arg.phyl$phyl
    phyl.phylo <- arg.phyl$phyl.phylo
    rm(arg.phyl)
    if(!hasEdgeLength(phyl)) stop("Missing branch length")
  
    ## phyl is a phylo4 object, phyl.phylo is a phylo object 
    if(is.null(select)) {
        sum1 <- sumEdgeLength(phyl)
    } 
    else {
        if(is.ultrametric(phyl.phylo)) {
            phyl.D <- cophenetic.phylo(phyl.phylo)/2
            if(length(select)==1)
                sum1 <- max(phyl.D)
            else if(length(select)==2) {
                sum1 <- phyl.D[select[1], select[2]] + max(phyl.D)
            } 
            else {
                fun.EH <- function(i) {
                    min(phyl.D[select[i], select[1:(i - 1)]])
                }
                sum1 <-  phyl.D[select[1], select[2]] + max(phyl.D) + sum(sapply(3:length(select), fun.EH))
            }      
        }
        else {  
            sum1 <- sumEdgeLength(phyl, unique(unlist(ancestors(phyl, select, type="ALL"))))
        }
    }
    return(sum1)
}

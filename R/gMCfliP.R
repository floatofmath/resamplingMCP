##' Perform a graph-based closed test procedure using permutation tests
##'
##' @title graph based multiple comparsions procedure using permutation tests
##' @param permTP flip.object
##' @param G graphMCP
##' @param ... additional arguments to \code{\link{flip::npc}}
##' @return vector with adjustde p-values
##'
##' @importFrom flip npc
##' @importFrom gMCP generateWeights
##' @author Florian Klinglmueller
##'
##' @export
gMCfliP <- function(permTP,G,...){
    weights <- gMCP:::generateWeights(G@m,G@weights)
    n2 <- ncol(weights)
    weight_list <- alply(weights,1,tail,n=n2/2)
    weight_list <- lapply(weight_list,`names<-`,names(permTP))
    subset_list <- alply(weights,1,function(w) names(permTP)[head(w,n=n2/2)>0])
    weight_list <- lapply(1:length(subset_list),function(i) weight_list[[i]][subset_list[[i]]])
    names(weight_list) <- names(subset_list)
    p <- npc(permTP,subsets=subset_list,weights=weight_list,...)
    adjP <- sapply(1:length(permTP),closed_pval,p,m=length(permTP))
    names(adjP) <- names(permTP)
    return(adjP)
}

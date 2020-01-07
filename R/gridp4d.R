gridp4d <- function (p4d, trait = names(tdata(p4d)), center = TRUE, scale = TRUE, 
    tree.ladderize = FALSE, tree.type = "phylogram", tree.ratio = NULL, 
    tree.xlim = NULL, tree.open.angle = 0, tree.open.crown = TRUE, 
    show.tip = TRUE, tip.labels = NULL, tip.col = "black", 
    tip.cex = 1, tip.font = 3, tip.adj = 0, cell.col = topo.colors(100), 
    show.color.scale = TRUE, show.trait = TRUE, trait.labels = NULL, 
    trait.col = "black", trait.cex = 0.7, trait.font = 1, 
    trait.bg.col = "grey90", show.box = FALSE, grid.vertical = FALSE, 
    grid.horizontal = FALSE, grid.col = "grey25", grid.lty = "dashed", 
    ...) 
{
    plot.phylo4d(x = p4d, trait = trait, center = center, scale = scale, 
        plot.type = "gridplot", tree.ladderize = tree.ladderize, 
        tree.type = tree.type, tree.ratio = tree.ratio, tree.xlim = tree.xlim, 
        tree.open.angle = tree.open.angle, tree.open.crown = tree.open.crown, 
        show.tip = show.tip, tip.labels = tip.labels, tip.col = tip.col, 
        tip.cex = tip.cex, tip.font = tip.font, tip.adj = tip.adj, 
        cell.col = cell.col, show.color.scale = show.color.scale, 
        show.trait = show.trait, trait.labels = trait.labels, 
        trait.col = trait.col, trait.cex = trait.cex, trait.font = trait.font, 
        trait.bg.col = trait.bg.col, show.box = show.box, grid.vertical = grid.vertical, 
        grid.horizontal = grid.horizontal, grid.col = grid.col, 
        grid.lty = grid.lty, ...)
}

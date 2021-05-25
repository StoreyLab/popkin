# for determining x-axis limit in plots
phylo_max_edge <- function( tree ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    ## # would validate further with this, but it's too verbose unfortunately :(
    ## ape::checkValidPhylo( tree )
    
    # algorithm is sensitive to edge ordering
    # assumption is that we move from the root up, which is the reverse of postorder
    order_edges <- rev( ape::postorder( tree ) )
    
    # this is what we want, partially, the sum of edges from root
    # initialize to NA to catch issues
    node_edge_from_root <- rep.int( NA, max( tree$edge ) )
    
    # determine root node
    # it is very first parent node (in reverse postorder)
    j_root <- tree$edge[ order_edges[ 1 ], 1 ]
    # only root node is pre-set to zero, or to value of root edge if present
    node_edge_from_root[ j_root ] <- if ( is.null( tree$root.edge ) ) 0 else tree$root.edge
    
    # navigate all edges in reverse postorder!
    for ( e in order_edges ) {
        # get parent and child nodes for this edge
        j_parent <- tree$edge[ e, 1 ]
        j_child <- tree$edge[ e, 2 ]
        
        # get edge length of parent from root
        fst_parent_from_root <- node_edge_from_root[ j_parent ]
        # check
        if ( is.na( fst_parent_from_root ) )
            stop( 'Node index ', j_parent, ' was not assigned edge value from root! (unexpected)' )
        
        # store total edge for child to root (as it may become a parent in later iterations)
        node_edge_from_root[ j_child ] <- fst_parent_from_root + tree$edge.length[ e ]
    }
    
    # just want maximum length
    return( max( node_edge_from_root ) )
}

#' Random Walk with Restart (RWR) methods for functional prediction in a graph
#' @param g An igraph object for networks.
#' @param seedss A character vector of seeds node names (known disease genes or functional genesets).
#' @param restart_prob Numeric (0-1). The restart probability (r). 
#'   - High r (e.g., 0.7): Very local search, stays close to seeds.
#'   - Low r (e.g., 0.1): Global search, diffuses further away.
#'   - Default is 0.7 based on user's image.
#' @param score_name Character. The name of the attribute to store in the graph. Default is "RWR_score".
#' @param max_iter Integer. Max iterations. Default is 100.
#' @param threshold Numeric. Convergence threshold. Default is 1e-6.
#' @importFrom igraph is.igraph V as_adjacency_matrix vertex_attr<-
#' @importFrom Matrix colSums Diagonal Matrix
#' @importFrom methods as
#'
#' @return An igraph object containing RWR score for a specific geneset.
#' @examples
#' \dontrun{
#' data(demo_ppi)
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(dplyr)
#' library(igraph)
#' ppi_g <- igraph::V(demo_ppi)$name
#' x <- enrichGO(gene = ppi_g, 
#'               ont = "BP", 
#'               OrgDb = 'org.Hs.eg.db', 
#'               keyType = "SYMBOL")
#' 
#' go_res <- x@result
#' ir_g <- subset(go_res, ID == "GO:0050727")$geneID |> strsplit("/") |> unlist()
#'
#' ppi <- run_RWR(g = demo_ppi, 
#'                seeds = ir_g, 
#'                score_name = "inflammatory_response")
#' 
#' res_df <- igraph::as_data_frame(ppi, what = "vertices")
#'
#' top_candidates <- res_df %>%
#'   arrange(desc(inflammatory_response)) %>%
#'   filter(is_seeds_inflammatory_response == FALSE) %>%
#'   head(10)
#' print(top_candidates)
#' }
#' @export
run_RWR <- function(g,
                    seeds,
                    restart_prob = 0.7,
                    score_name = "RWR_score",
                    max_iter = 100,
                    threshold = 1e-6){
  stopifnot(is.igraph(g))
  
  all_nodes <- V(g)$name
  n <- length(all_nodes)
  valid_seeds <- seeds[seeds %in% all_nodes]
  
  if (length(valid_seeds) == 0) {
    stop("No valid seeds found in the network. Please check your gene IDs.")
  }
  
  message(sprintf("Running RWR with r=%.2f. %d/%d seeds found in network.", 
                  restart_prob, length(valid_seeds), length(seeds)))
  
  ## initialize the transition matrix
  # normalization for cols
  adj <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  col_sum <- Matrix::colSums(adj)
  col_sum[col_sum == 0] <- 1
  inv_diag <- Matrix::Diagonal(x = 1/col_sum) # faster
  W <- adj %*% inv_diag
  
  ## initialize the probability vector P0
  P0 <- Matrix::Matrix(0, nrow = n, ncol = 1)
  rownames(P0) <- all_nodes
  
  seeds_indices <- which(all_nodes %in% valid_seeds)
  P0[seeds_indices, 1] <- 1 / length(valid_seeds)
  
  P_t <- P0
  
  ## iterations
  for (i in 1:max_iter){
    P_prev <- P_t
    
    # random walk ((1-r) * W * P_t)
    diffusion <- (1 - restart_prob) * (W %*% P_t)
    
    # restart (r * P0)
    restart <- restart_prob * P0
    
    # update
    P_t <- diffusion + restart
    
    # check divergence
    delta <- sum(abs(P_t - P_prev))
    
    if(delta < threshold) {
      message(sprintf("Converged at iteration %d (Delta: %.2e)", i, delta))
      break
    }
  }
  
  ## get results
  # RWR score
  final_scores <- as.numeric(P_t)
  igraph::vertex_attr(g, name = score_name) <- final_scores
  
  # seeds info
  seeds_attr_name <- paste0("is_seeds_", score_name)
  igraph::vertex_attr(g, name = seeds_attr_name) <- V(g)$name %in% valid_seeds
    
  return(g)
}

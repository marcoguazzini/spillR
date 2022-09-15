#' Generate compensated counts dataset
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import tidygraph
#' @import ggplot2
#' @import ggraph
#' @export
#'
#' @param spillover_matrix  Matrix number of markers x number of markers which contains the amount of spillover between markers
#'   
#' @return plot
#'
#' @examples
#' set.seed(23)
#' sm <- spillR::load_spillover()
#' spill_dag(sm)
spill_dag <- function(spillover_matrix){
  
  graph <- as.data.frame(as.table(sm)) %>% 
    setNames(c('from', 'to', 'value')) %>% 
    mutate(percentage = value*100) %>% 
    filter(percentage >= 0.1 & percentage < 100) %>% 
    as_tbl_graph()
  ggraph(graph, layout = "stress") + 
    geom_edge_link(
      aes(
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name),
        color = percentage
      ), 
      arrow = arrow(length = unit(4, 'mm'))) + 
    geom_node_text(aes(label = name))
}

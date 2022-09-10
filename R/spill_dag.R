spill_dag <- function(spillover_matrix){
  library(tidygraph)
  library(ggraph)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
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
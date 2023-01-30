# Library
library(networkD3)
library(dplyr)

# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
  source=c("class_1_blood","class_1_blood", "class_1_blood","class_1_blood", 
           "class_2_blood", "class_2_blood", "class_2_blood", "class_2_blood", 
           "class_3_blood", "class_3_blood", "class_3_blood", "class_3_blood"), 
  target=c("class_1_brain","class_2_brain", "class_3_brain", "none", 
           "class_1_brain","class_2_brain", "class_3_brain", "none", 
           "class_1_brain","class_2_brain", "class_3_brain", "none"), 
  value=c(1167, 679, 1106, 23083,
          236, 341, 743, 17819, 
          307, 360, 2359, 61397)
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 20,
                   sinksRight=FALSE)
p

library(htmlwidgets)
saveWidget(p, file=paste0( getwd(), "/scripts/BRAINcode/blood_and_brain_TNEs_class.html"))

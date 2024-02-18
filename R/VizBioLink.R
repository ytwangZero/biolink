#' @title Visualize the biological link
#' @description Visualize the biological link. It should be noted that the corresponding link 
#'     has been built by BioLink function prior to using it.
#' @usage VizBioLink(eSet = eSet,Mode = "PPI",Layout = "force-directed",Brightness = "dark",Palette = "default2")
#' @param eSet An R6 class object.
#' @param Mode chr. Method to build the biological link between exposure and disease. 
#'    Only PPI is currently supported.
#' @param Layout chr. Visualization layout. Available options include "Force-directed" and "Degree-circle".  
#' @param Brightness chr. Visualization brightness. Available options include "Light" and "Dark".  
#' @param Palette chr. Visualization palette. Available options include "default1", "default2"
#'    and "Journal". The "Journal" option provides several journal preference styles including cell, nature, 
#'    science, lancet, nejm, and jama. 
#' @param Method chr. The method to draw the net plot.users can select "ggraph" or "cytoscape".
#'    
#' @details 
#' @return An R6 class object containing the plot of the biological link.
#' @export
#' @examples eSet <- InitBioLink(1)
#' eSet <- LoadBioLink(eSet = eSet,UseExample = "example#1", FileDirLink = NULL)
#' eSet <- ConvToExpoID(eSet = eSet)
#' eSet <- Biolink(eSet = eSet,Mode = "PPI",ChemCas = "default",ChemInchikey = "default",
#' DiseaseID = "default",MetabolomeID = "default",MetBiospec = "blood",ProteomeID = "default",GenomeID = "default")
#' eSet <- VizBioLink(eSet = eSet,Mode = "PPI",Layout = "force-directed",Brightness = "dark",Palette = "default2")
#' @author Mingliang Fang, Bin Wang (corresponding author), Yuting Wang, Ting Wu
VizBioLink <- function(eSet,
                       Mode = "PPI",
                       Method = "ggraph",#cytoscape
                       Layout = "all",
                       Brightness = "light",
                       Palette = 'cell'){ 
  
  tictoc::tic()
  
  lubridate::now() %>% 
    stringr::str_replace_all(":",".") %>% 
    stringr::str_replace_all("-",".") -> NowTime 
  
  # Path = stringr::str_c(getwd(),"/output_",eSet$PID, "/Biolink")
  # if(!file.exists(Path)) {dir.create(Path)}
  
  # Path2 = stringr::str_c(getwd(),"/output_",eSet$PID, "/Biolink/Legend")
  # if(!file.exists(Path2)) {dir.create(Path2)}
  message("It will take several minutes, please wait patiently.")
  switch(Method,
         "cytoscape" = {
           RCy3::cytoscapePing() 
           RCy3::closeSession(F)
           
           #save plot=----------------------
           func_save_ppi <- function(y){
             ddpcr::quiet(
               switch(Layout,
                      "force-directed" = {
                        RCy3::layoutNetwork(Layout)
                        RCy3::saveSession(stringr::str_c(Path, "/PPI_",Method,"_", Layout,"_",Brightness,"_",Palette))
                        RCy3::exportImage(stringr::str_c(Path, "/PPI_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                    overwriteFile = T,
                                    type = "PNG",
                                    zoom=300) 
                        RCy3::exportImage(stringr::str_c(Path, "/PPI_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                          overwriteFile = T,
                                          type = "SVG",
                                          zoom=300) 
                      },
                      "degree-circle" = {
                        RCy3::layoutNetwork(Layout)
                        RCy3::saveSession(stringr::str_c(Path, "/PPI_",Method,"_", Layout,"_",Brightness,"_",Palette))
                        RCy3::exportImage(stringr::str_c(Path, "/PPI_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                    overwriteFile = T,
                                    type = "PNG",
                                    zoom=300)
                        RCy3::exportImage(stringr::str_c(Path, "/PPI_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                          overwriteFile = T,
                                          type = "SVG",
                                          zoom=300)
                      })
             )
           }
           func_save_go <- function(z){
             ddpcr::quiet(
               switch(Layout,
                      "force-directed" = {
                        RCy3::layoutNetwork(Layout)
                        RCy3::saveSession(stringr::str_c(Path, "/GO_",Method,"_", Layout,"_",Brightness,"_",Palette))
                        RCy3::exportImage(stringr::str_c(Path, "/GO_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                    overwriteFile = T,
                                    type = "PNG",
                                    zoom=300) 
                        RCy3::exportImage(stringr::str_c(Path, "/GO_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                          overwriteFile = T,
                                          type = "SVG",
                                          zoom=300) 
                      },
                      "degree-circle" = {
                        RCy3::layoutNetwork(Layout)
                        RCy3::saveSession(stringr::str_c(Path, "/GO_",Method,"_", Layout,"_",Brightness,"_",Palette))
                        RCy3::exportImage(stringr::str_c(Path, "/GO_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                    overwriteFile = T,
                                    type = "PNG",
                                    zoom=300) 
                        RCy3::exportImage(stringr::str_c(Path, "/GO_",Method,"_", Layout,"_",Brightness,"_",Palette), 
                                          overwriteFile = T,
                                          type = "SVG",
                                          zoom=300)
                      })
             )
           }
           
           # color----------------------
           Node_Color <- NULL
           NodeLabel_Color <- NULL
           NodeBorder_Color <- NULL
           Node_Shape <- NULL
           Edge_Color <- NULL
           
           # node color----------------------
           c("#E85472","#64B5F6","#d27eef","#ffa3d8","#d27eef","#d27eef") -> Node_Color[["light_default1"]]
           c("#D81C38","#00B64C","#2196F3","#ee88fd","#2196F3","#2196F3") -> Node_Color[["light_default2"]]
           c("#ff8033","#64acbf","#1b6393","#edab63","#1b6393","#1b6393") -> Node_Color[["light_cell"]]
           c("#ff0000","#4DBBD5","#5c50e6","#F39B7F","#5c50e6","#5c50e6") -> Node_Color[["light_nature"]]
           c("#EE0000","#20B2AA","#5000B8","#B85798","#5000B8","#5000B8") -> Node_Color[["light_science"]]
           c("#ED0000","#0099B4","#ffb366","#925E9F","#ffb366","#ffb366") -> Node_Color[["light_lancet"]]
           c("#EE4C97","#E18727","#5A9F77","#4A91C5","#5A9F77","#5A9F77") -> Node_Color[["light_nejm"]]
           c("#ff4d00","#FFA54F","#79af97","#00a1d5","#79af97","#79af97") -> Node_Color[["light_jama"]]
           
           c("#00441B","#67001F","#08306B","#A63603","#08306B","#08306B") -> Node_Color[["dark_default1"]]
           c("#600707","#0C1B46","#33691E","#7B1FA2","#33691E","#33691E") -> Node_Color[["dark_default2"]]
           c("#A51526","#323695","#698B69","#8B4500","#698B69","#698B69") -> Node_Color[["dark_cell"]]
           c("#A51526","#3C5488","#8B7500","#00A087","#8B7500","#8B7500") -> Node_Color[["dark_nature"]]
           c("#bb0021","#191970","#008280","#CD6600","#008280","#008280") -> Node_Color[["dark_science"]]
           c("#AD002A","#00468B","#ADB6B6","#800080","#ADB6B6","#ADB6B6") -> Node_Color[["dark_lancet"]]
           c("#A52A2A","#CC5500","#6F99AD","#2A52BE","#6F99AD","#6F99AD") -> Node_Color[["dark_nejm"]]
           c("#b24745","#cc5500","#6a6599","#374e55","#6a6599","#6a6599") -> Node_Color[["dark_jama"]]
           
           # node label color----------------------
           c("#E85472","#64B5F6","#d27eef","#ffa3d8","#d27eef","#d27eef") -> NodeLabel_Color[["light_default1"]]
           c("#D81C38","#00B64C","#2196F3","#ee88fd","#2196F3","#2196F3") -> NodeLabel_Color[["light_default2"]]
           c("#ff8033","#64acbf","#1b6393","#edab63","#1b6393","#1b6393") -> NodeLabel_Color[["light_cell"]]
           c("#dc0000","#4DBBD5","#5c50e6","#F39B7F","#5c50e6","#5c50e6") -> NodeLabel_Color[["light_nature"]]
           c("#EE0000","#20B2AA","#5000B8","#B85798","#5000B8","#5000B8") -> NodeLabel_Color[["light_science"]]
           c("#ED0000","#0099B4","#d2691e","#925E9F","#d2691e","#d2691e") -> NodeLabel_Color[["light_lancet"]]
           c("#EE4C97","#E18727","#5A9F77","#4A91C5","#5A9F77","#5A9F77") -> NodeLabel_Color[["light_nejm"]]
           c("#ff4d00","#FFA54F","#2e8b57","#00a1d5","#2e8b57","#2e8b57") -> NodeLabel_Color[["light_jama"]]
           
           c("#00441B","#67001F","#08306B","#A63603","#08306B","#08306B") -> NodeLabel_Color[["dark_default1"]]
           c("#600707","#0C1B46","#33691E","#7B1FA2","#33691E","#33691E") -> NodeLabel_Color[["dark_default2"]]
           c("#A51526","#323695","#004d99","#8B4500","#004d99","#004d99") -> NodeLabel_Color[["dark_cell"]]
           c("#A51526","#3C5488","#8B7500","#00A087","#8B7500","#8B7500") -> NodeLabel_Color[["dark_nature"]]
           c("#bb0021","#191970","#008280","#CD6600","#008280","#008280") -> NodeLabel_Color[["dark_science"]]
           c("#AD002A","#00468B","#404040","#800080","#404040","#404040") -> NodeLabel_Color[["dark_lancet"]]
           c("#A52A2A","#CC5500","#50788C","#2A52BE","#50788C","#50788C") -> NodeLabel_Color[["dark_nejm"]]
           c("#b24745","#cc5500","#483D8B","#374e55","#483D8B","#483D8B") -> NodeLabel_Color[["dark_jama"]]
           
           # node border color----------------------
           c("#f1c9cf","#BFD0F6","#C6C4E1","#F8BBD0","#C6C4E1","#C6C4E1") -> NodeBorder_Color[["light_default1"]]
           c("#D81C38","#00B64C","#2196F3","#ee88fd","#2196F3","#2196F3") -> NodeBorder_Color[["light_default2"]]
           c("#cd853f","#5686bf","#5686bf","#ff8033","#5686bf","#5686bf") -> NodeBorder_Color[["light_cell"]]
           c("#ff6347","#87cefa","#6495ed","#ff8033","#6495ed","#6495ed") -> NodeBorder_Color[["light_nature"]]
           c("#F08080","#33E6CC","#9932CC","#EE82EE","#9932CC","#9932CC") -> NodeBorder_Color[["light_science"]]
           c("#ED0000","#0099B4","#ff8c69","#925E9F","#ff8c69","#ff8c69") -> NodeBorder_Color[["light_lancet"]]
           c("#FF80BF","#FF9900","#16982B","#6495ED","#16982B","#16982B") -> NodeBorder_Color[["light_nejm"]]
           c("#cd5c5c","#e9967a","#2e8b57","#87ceeb","#2e8b57","#2e8b57") -> NodeBorder_Color[["light_jama"]]
           
           c("#00441B","#67001F","#08306B","#A63603","#08306B","#08306B") -> NodeBorder_Color[["dark_default1"]]
           c("#600707","#0C1B46","#33691E","#7B1FA2","#33691E","#33691E") -> NodeBorder_Color[["dark_default2"]]
           c("#b22222","#2a52be","#698B69","#d2691e","#698B69","#698B69") -> NodeBorder_Color[["dark_cell"]]
           c("#A51526","#3C5488","#8B7500","#00A087","#8B7500","#8B7500") -> NodeBorder_Color[["dark_nature"]]
           c("#E32636","#0000CD","#48D1CC","#FF8C00","#48D1CC","#48D1CC") -> NodeBorder_Color[["dark_science"]]
           c("#E32636","#007FFF","#ADB6B6","#9932CC","#ADB6B6","#ADB6B6") -> NodeBorder_Color[["dark_lancet"]]
           c("#E60000","#FF7300","#4798B3","#4169E1","#4798B3","#4798B3") -> NodeBorder_Color[["dark_nejm"]]
           c("#CD5C5C","#F28500","#8674A1","#006374","#8674A1","#8674A1") -> NodeBorder_Color[["dark_jama"]]
           
           # edge color----------------------
           c("#f1c9cf","#C6C4E1","#BFD0F6","#F8BBD0") -> Edge_Color[["light_default1"]]
           c("#EFA5B0","#2196F3","#98CB9A","#C097C6") -> Edge_Color[["light_default2"]]
           c("#f5cead","#1b6393","#64acbf","#edab63") -> Edge_Color[["light_cell"]]
           c("#EA7D7C","#5c50e6","#4DBBD5","#F39B7F") -> Edge_Color[["light_nature"]]
           c("#ECC1C7","#5000B8","#20B2AA","#B85798") -> Edge_Color[["light_science"]]
           c("#f9abab","#ff8c69","#0099B4","#925E9F") -> Edge_Color[["light_lancet"]]
           c("#F8BDD7","#5A9F77","#E18727","#4A91C5") -> Edge_Color[["light_nejm"]]
           c("#e1bcbc","#79af97","#FFA54F","#00a1d5") -> Edge_Color[["light_jama"]]
           c("#956171","#42546f","#799c87","#cd9880") -> Edge_Color[["dark_default1"]]
           c("#fd8a8a","#81C784","#1C7CD1","#C097C6") -> Edge_Color[["dark_default2"]]
           c("#ff8099","#4575B4","#323695","#8B4500") -> Edge_Color[["dark_cell"]]
           c("#EA7D7C","#8B7500","#6B7CA2","#00A087") -> Edge_Color[["dark_nature"]]
           c("#ECC1C7","#008280","#191970","#CD6600") -> Edge_Color[["dark_science"]]
           c("#E0ABB5","#ADB6B6","#00468B","#800080") -> Edge_Color[["dark_lancet"]]
           c("#E5B8B5","#ADB6B6","#CC5500","#2A52BE") -> Edge_Color[["dark_nejm"]]
           c("#D29898","#6a6599","#cc5500","#374e55") -> Edge_Color[["dark_jama"]]
           # node shape----------------------
           c("TRIANGLE", "DIAMOND", "ELLIPSE", "ROUND_RECTANGLE", "ELLIPSE", "ELLIPSE") -> Node_Shape
           
           # plot----------------------
           switch (Mode,
                   "PPI" = {
                     RCy3::createNetworkFromDataFrames(eSet$Biolink$Nodes_PPI,
                                                 eSet$Biolink$Edges_PPI, 
                                                 title = "PPI",
                                                 collection = "biolinker")
                     #default style--------------------------------------------------------
                     style = "Style1"
                     defaults <- list(NODE_SHAPE="diamond",NODE_SIZE=20,
                                      EDGE_TRANSPARENCY=225,NODE_LABEL_POSITION="S,W,c,0.00,0.00")
                     nodeLabels <- RCy3::mapVisualProperty('node label','label','p')
                     nodeShape <- RCy3::mapVisualProperty('Node Shape','group',"d",
                                                    c("exposure","disease","ppi", "chem_by_met","prot_by_metabolome", "prot_by_proteome"),
                                                    c("TRIANGLE", "DIAMOND", "ELLIPSE", "ROUND_RECTANGLE", "ELLIPSE", "ELLIPSE")
                                                    #Node_Shape[1:length(Node_Group)]
                     )
                     nodeSize <- RCy3::mapVisualProperty('Node Size','group',"d",
                                                   c("exposure","disease","ppi", "chem_by_met","prot_by_metabolome", "prot_by_proteome"),
                                                   c(40, 40, 25, 35, 40, 40))
                     nodeLableFontSize <- RCy3::mapVisualProperty('Node Label Font Size','group',"d",
                                                            c("exposure","disease","ppi", "chem_by_met","prot_by_metabolome", "prot_by_proteome"),
                                                            c(30, 30, 25, 30, 30, 30))
                     nodeWidth <- RCy3::mapVisualProperty('Node Width','group',"d",
                                                    c("exposure","disease","ppi", "chem_by_met","prot_by_metabolome", "prot_by_proteome"),
                                                    c(40, 40, 25, 35, 40, 40))
                     arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','edge_type','d',
                                                      c(1,2,3,4),
                                                      c("Arrow","Arrow","Arrow","Arrow"))   
                     arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','edge_type','d',
                                                      c(1,2,3,4), 
                                                      c("Arrow","Arrow","Arrow","Arrow"))   
                     edgeWidth <- RCy3::mapVisualProperty('edge width','edge_type','d',
                                                    c(1,2,3,4), 
                                                    c(2,  #uniprotKB
                                                      1,  #protein_protein
                                                      2,  #chem_protein
                                                      2)  #metabolite_protein
                     )  
                     RCy3::createVisualStyle(style, 
                                       defaults, 
                                       list(nodeLabels,
                                            #nodeFills,
                                            nodeShape,
                                            nodeSize,
                                            nodeLableFontSize,
                                            #nodeLableColor,
                                            nodeWidth,
                                            arrowShapes,
                                            edgeWidth
                                       ))
                     RCy3::setVisualStyle(style)
                     #plot---------------------
                     RCy3::setNodeColorMapping(mapping.type = 'd',
                                         table.column = 'group',
                                         table.column.values = c("exposure","disease","ppi", "chem_by_met","protein_by_metabolome", "protein_by_proteome"),
                                         Node_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                         style.name = style)
                     RCy3::setNodeLabelColorMapping(mapping.type = 'd',
                                              table.column = 'group',
                                              table.column.values = c("exposure","disease","ppi", "metabolite","protein_by_metabolome", "protein_by_proteome"),
                                              NodeLabel_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                              style.name = style)
                     RCy3::setNodeBorderWidthDefault(new.width = 4,style.name = style)
                     RCy3::setNodeBorderColorMapping(mapping.type = 'd',
                                               table.column = 'group',
                                               table.column.values = c("exposure","disease","ppi", "metabolite","protein_by_metabolome", "protein_by_proteome"),
                                               NodeBorder_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                               style.name = style)
                     RCy3::setEdgeColorMapping('edge_type', c(1,2,3,4), 
                                         Edge_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                         style.name = style)
                     RCy3::setEdgeTargetArrowColorMapping('edge_type', c(1,2,3,4), 
                                                    Edge_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                                    style.name = style)
                     
                     # node group----------------------
                     eSet$Biolink$Nodes_PPI %>% 
                       dplyr::distinct(group) %>% 
                       as.matrix() %>% as.vector() -> Node_Group
                     # match group and color legend -----------------------------------------
                     Legend_Shape <- NULL
                     24 -> Legend_Shape[["exposure"]]
                     23 -> Legend_Shape[["disease"]]
                     21 -> Legend_Shape[["ppi"]]
                     22 -> Legend_Shape[["metabolite"]]
                     21 -> Legend_Shape[["protein_by_metabolome"]]
                     21 -> Legend_Shape[["protein_by_proteome"]]
                     
                     Legend_Color <- NULL
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][1] -> Legend_Color[["exposure"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][2] -> Legend_Color[["disease"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][3] -> Legend_Color[["ppi"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][4] -> Legend_Color[["metabolite"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][5] -> Legend_Color[["protein_by_metabolome"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][6] -> Legend_Color[["protein_by_proteome"]]
                     
                     # create legend -----------------------------------------
                     ggplot2::ggplot(data = eSet$Biolink$Nodes_PPI %>% 
                              dplyr::mutate(Group = group),
                            aes(x = id,y = label,
                                #size = MinusLog10P,
                                fill = Group))+
                       ggplot2::geom_point(aes(shape = Group),
                                  size = 5)+
                       ggplot2::scale_shape_manual(values = Legend_Shape[Node_Group])+
                       ggplot2::scale_fill_manual(values = Legend_Color[Node_Group]) -> Plot_Legend
                     ggplot2::ggsave(stringr::str_c(Path,"/Legend/",Mode,"_",Brightness,"_",Palette,"_PlotLegend.png"),
                                     cowplot::get_legend(Plot_Legend),
                                     #limitsize = FALSE,
                                     width = 2,
                                     height = 2
                     )
                     
                     # add legend into cytoscape -----------------------------------------------
                     RCy3::deleteAnnotation(names = "legend")
                     
                     RCy3::addAnnotationImage(url = str_c(Path,"/Legend/",Mode,"_",Brightness,"_",Palette,"_PlotLegend.png"),
                                        network = "PPI",
                                        x.pos = 950,
                                        y.pos = -300,
                                        name = "legend",
                                        height = 400,
                                        width = 300,
                                        canvas = "background",
                                        borderColor = "#FFFFFF")
                     
                     # save---------------------
                     func_save_ppi(y)
                   },
                   
                   "GO" = {
                     RCy3::createNetworkFromDataFrames(eSet$Biolink$Nodes_GO,
                                                 eSet$Biolink$Edges_GO, 
                                                 title = "GO",
                                                 collection = "biolinker")
                     
                     #default style--------------------------------------------------------
                     style = "Style1"
                     defaults <- list(NODE_SHAPE="diamond",NODE_SIZE=20,
                                      EDGE_TRANSPARENCY=225,NODE_LABEL_POSITION="S,W,c,0.00,0.00")
                     nodeLabels <- RCy3::mapVisualProperty('node label','label','p')
                     nodeShape <- RCy3::mapVisualProperty('Node Shape','group',"d",
                                                    c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                                    c("TRIANGLE", "DIAMOND", "ELLIPSE", "ROUND_RECTANGLE", "ELLIPSE", "ELLIPSE"))
                     nodeSize <- RCy3::mapVisualProperty('Node Size','group',"d",
                                                   c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                                   c(40, 40, 25, 35, 30, 40))
                     nodeLableFontSize <- RCy3::mapVisualProperty('Node Label Font Size','group',"d",
                                                            c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                                            c(30, 30, 25, 30, 30, 30))
                     nodeWidth <- RCy3::mapVisualProperty('Node Width','group',"d",
                                                    c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                                    c(40, 40, 25, 35, 30, 40))
                     arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','edge_type','d',
                                                      c(1,2,3,4),
                                                      c("Arrow","Arrow","Arrow","Arrow"))   
                     arrowShapes <- RCy3::mapVisualProperty('Edge Target Arrow Shape','edge_type','d',
                                                      c(1,2,3,4), 
                                                      c("Arrow","Arrow","Arrow","Arrow"))   
                     edgeWidth <- RCy3::mapVisualProperty('edge width','edge_type','d',
                                                    c(1,2,3,4), 
                                                    c(2,  #uniprotKB
                                                      1,  #protein_protein
                                                      2,  #chem_protein
                                                      2)  #metabolite_protein
                     )  
                     RCy3::createVisualStyle(style, 
                                       defaults, 
                                       list(nodeLabels,
                                            #nodeFills,
                                            nodeShape,
                                            nodeSize,
                                            nodeLableFontSize,
                                            #nodeLableColor,
                                            nodeWidth,
                                            arrowShapes,
                                            edgeWidth
                                       ))
                     RCy3::setVisualStyle(style)
                     #plot---------------------
                     RCy3::setNodeColorMapping(mapping.type = 'd',
                                         table.column = 'group',
                                         table.column.values = c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                         Node_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                         style.name = style)
                     RCy3::setNodeLabelColorMapping(mapping.type = 'd',
                                              table.column = 'group',
                                              table.column.values = c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                              NodeLabel_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                              style.name = style)
                     RCy3::setNodeBorderWidthDefault(new.width = 8,style.name = style)
                     RCy3::setNodeBorderColorMapping(mapping.type = 'd',
                                               table.column = 'group',
                                               table.column.values = c("exposure","disease","go", "chem_by_met","go_by_metabolome", "go_by_proteome"),
                                               NodeBorder_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                               style.name = style)
                     RCy3::setEdgeColorMapping('edge_type', c(1,2,3,4), 
                                         Edge_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                         style.name = style)
                     RCy3::setEdgeTargetArrowColorMapping('edge_type', c(1,2,3,4), 
                                                    Edge_Color[[stringr::str_c(Brightness,"_",Palette)]],
                                                    style.name = style)
                     
                     # node group----------------------
                     eSet$Biolink$Nodes_GO %>% 
                       dplyr::distinct(group) %>% 
                       as.matrix() %>% as.vector() -> Node_Group
                     # match group and color legend -----------------------------------------
                     Legend_Shape <- NULL
                     24 -> Legend_Shape[["exposure"]]
                     23 -> Legend_Shape[["disease"]]
                     21 -> Legend_Shape[["go"]]
                     22 -> Legend_Shape[["chem_by_met"]]
                     21 -> Legend_Shape[["go_by_metabolome"]]
                     21 -> Legend_Shape[["go_by_proteome"]]
                     
                     Legend_Color <- NULL
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][1] -> Legend_Color[["exposure"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][2] -> Legend_Color[["disease"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][3] -> Legend_Color[["go"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][4] -> Legend_Color[["chem_by_met"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][5] -> Legend_Color[["go_by_metabolome"]]
                     Node_Color[[stringr::str_c(Brightness,"_",Palette)]][6] -> Legend_Color[["go_by_proteome"]]
                     
                     # create legend -----------------------------------------
                     ggplot2::ggplot(data = eSet$Biolink$Nodes_GO %>% 
                              dplyr::mutate(Group = group),
                            aes(x = id,y = label,
                                #size = MinusLog10P,
                                fill = Group))+
                       ggplot2::geom_point(aes(shape = Group),
                                  size = 5)+
                       ggplot2::scale_shape_manual(values = Legend_Shape[Node_Group])+
                       ggplot2::scale_fill_manual(values = Legend_Color[Node_Group]) -> Plot_Legend
                     ggplot2::ggsave(stringr::str_c(Path,"/Legend/",Mode,"_",Brightness,"_",Palette,"_PlotLegend.png"),
                                     cowplot::get_legend(Plot_Legend),
                                     #limitsize = FALSE,
                                     width = 1.5,
                                     height = 2
                     )
                     
                     # add legend into cytoscape -----------------------------------------------
                     RCy3::deleteAnnotation(names = "legend")
                     
                     RCy3::addAnnotationImage(url = str_c(Path,"/Legend/",Mode,"_",Brightness,"_",Palette,"_PlotLegend.png"),
                                        network = "GO",
                                        x.pos = 450,
                                        y.pos = -100,
                                        name = "legend",
                                        height = 350,
                                        width = 200,
                                        canvas = "background",
                                        borderColor = "#FFFFFF")
                     
                     # save---------------------
                     func_save_go(z)
                     
                   }
           )
         },
         "ggraph" = {
           switch(Mode,
                  'PPI' = {
                    eSet$Biolink$Nodes_PPI %>% 
                      dplyr::distinct(group) %>% 
                      as.matrix() %>% as.vector() -> GroupName
                    
                    eSet$Biolink$Edges_PPI %>% 
                      dplyr::select(source,target,edge_type) %>% 
                      setNames(c('from','to',"value")) -> connect_data
                    
                    rbind(
                      data.frame(from = "origin",
                                 to = GroupName),
                      data.frame(to = unique(c(as.character(connect_data$to),as.character(connect_data$to)))) %>% 
                        dplyr::left_join(eSet$Biolink$Nodes_PPI %>% 
                                           dplyr::select(id,group),
                                         by = c("to" = "id")) %>% 
                        dplyr::distinct(to,.keep_all = TRUE) %>% 
                        dplyr::arrange(group) %>% 
                        setNames(c("to","from")) %>% 
                        dplyr::select(from,to)
                    ) -> Edges
                    
                    data.frame(
                      name = unique(c(as.character(Edges$from),as.character(Edges$to)))
                    ) -> vertices_data
                    
                    vertices_data %>% 
                      dplyr::left_join(eSet$Biolink$Nodes_PPI %>% 
                                         dplyr::select(id),
                                       by = c("name" = "id")) -> vertices_data
                    
                    vertices_data %>%
                      dplyr::mutate(Group = Edges$from[ match( vertices_data$name, Edges$to ) ],
                             id = NA) -> vertices_data
                    Leaves <- which(is.na(match(vertices_data$name,Edges$from)))
                    NLeaves <- length(Leaves)
                    vertices_data$id[Leaves] <- seq(1:NLeaves)
                    vertices_data$angle <- 90 - 360 * vertices_data$id / NLeaves
                    vertices_data$hjust <- ifelse(vertices_data$angle < -90,1,0)
                    vertices_data$angle <- ifelse(vertices_data$angle < -90,vertices_data$angle+180,
                                                  vertices_data$angle)
                    
                    c(as.character(connect_data$from),as.character(connect_data$to)) %>% 
                      as_tibble() %>% 
                      dplyr::group_by(value) %>% 
                      dplyr::summarize(n=n()) %>% 
                      # dplyr::left_join(eSet$Expo$Voca %>% 
                      #                    dplyr::select(SerialNo_Raw,
                      #                                  GroupName),by = c("value" = "SerialNo_Raw")) %>% 
                      dplyr::left_join(eSet$Biolink$Nodes_PPI,
                                       by = c("value" = "id")) %>% 
                      dplyr::distinct(value,.keep_all = TRUE) -> Vertices

                    colnames(Vertices) <- c("name", "n","GroupName","label")
                    
                    mygraph <- igraph::graph_from_data_frame( connect_data, vertices = Vertices, directed = FALSE )
                    
                    com <- igraph::walktrap.community(mygraph)
                    
                    Vertices1 <- Vertices %>% 
                      dplyr::mutate( group = com$membership) %>%
                      dplyr::mutate(group=as.numeric(factor(group,
                                                     levels=sort(summary (as.factor(group)),index.return=TRUE,decreasing = T)$ix,
                                                     order=TRUE)))%>%
                      dplyr::filter( group<10) %>%
                      dplyr::arrange(group,desc(n)) %>%
                      dplyr::mutate(name=factor(name, name))
                    
                    Connect_data <- connect_data %>%
                      dplyr::filter(from %in% Vertices1$name) %>%
                      dplyr::filter(to %in% Vertices1$name) %>%
                      dplyr::left_join(Vertices1,by=c('from'='name')) 
                    
                    Vertices1 %>% 
                      dplyr::distinct(label) %>% 
                      dplyr::arrange(label) %>% 
                      as.matrix() %>% as.vector() -> GroupName
                    
                    Vertices1$label <- factor(Vertices1$label,
                                              levels = GroupName)
                    
                    mygraph <- igraph::graph_from_data_frame( Connect_data, vertices = Vertices1, directed = TRUE )
                    #mycolor <- wesanderson::wes_palette("Darjeeling1", max(Vertices1$group), type = "continuous")
                    
                    # color -------------------------------------------------------------------
                    #Darjeeling1 BottleRocket2
                    # match group and color legend -----------------------------------------
                    Node_Shape <- data.frame(
                      exposure = 24,
                      disease = 23,
                      metabolite = 22,
                      gene = 25,
                      protein_by_metabolome = 21,
                      ppi = 21,
                      protein_by_proteome = 21,
                      protein_by_protome_metabolome = 21,
                      protein_by_genome = 21,
                      protein_by_protome_genome = 21,
                      protein_by_protome_metabolome_genome = 21
                    )
                    
                    Node_Shape %>% 
                      dplyr::select(all_of(GroupName)) %>% 
                      as.matrix() %>% 
                      as.vector() -> nodeshape
                    # # Legend_Color <- NULL
                    # # Node_Color[[stringr::str_c(Brightness,"_",Palette)]][1] -> Legend_Color[["exposure"]]
                    # # Node_Color[[stringr::str_c(Brightness,"_",Palette)]][2] -> Legend_Color[["disease"]]
                    # # Node_Color[[stringr::str_c(Brightness,"_",Palette)]][3] -> Legend_Color[["ppi"]]
                    # # Node_Color[[stringr::str_c(Brightness,"_",Palette)]][4] -> Legend_Color[["metabolite"]]
                    # # Node_Color[[stringr::str_c(Brightness,"_",Palette)]][5] -> Legend_Color[["protein_by_metabolome"]]
                    # # Node_Color[[stringr::str_c(Brightness,"_",Palette)]][6] -> Legend_Color[["protein_by_proteome"]]
                    # # 
                    # mycolor <- wesanderson::wes_palette("Darjeeling1", 
                    #                                     length(Vertices1 %>% dplyr::distinct(label) %>% 
                    #                                              as.matrix() %>% as.vector()),
                    #                                     type = "continuous")
                    # mycolor <- sample(mycolor, length(mycolor))
                    # 
                    mycolor <- c("#FF0000","turquoise1","#FBA3F3","#B7EF97",
                                 "#F2AD00","#F98400","#FDF7BB") %>% 
                      head(length(GroupName))
                    #scales::show_col(mycolor)
                    # labelsize ---------------------------------------------------------------
                    if(nrow(Vertices1) <= 100){
                      LabelSize = 2.5
                    }else if(nrow(Vertices1) > 100 & nrow(Vertices1) <= 300){
                      LabelSize = 2
                    }else if(nrow(Vertices1) >300 & nrow(Vertices1) <= 450){
                      LabelSize = 1.5
                    }else if(nrow(Vertices1) >450 & nrow(Vertices1) <= 600){
                      LabelSize = 1
                    }else{
                      LabelSize = 0.7
                    }
                    if(nrow(Vertices1) <= 100){
                      NodeSize = 4
                    }else if(nrow(Vertices1) > 100 & nrow(Vertices1) <= 300){
                      NodeSize = 3.5
                    }else if(nrow(Vertices1) >300 & nrow(Vertices1) <= 450){
                      NodeSize = 3
                    }else if(nrow(Vertices1) >450 & nrow(Vertices1) <= 600){
                      NodeSize = 2.5
                    }else{
                      NodeSize = 2
                    }
                    # Plot --------------------------------------------------------------------
                    switch(Layout,
                           "force-directed" = {
                             set.seed(123)
                             ddpcr::quiet(
                               ggraph::ggraph(mygraph,layout='fr') + 
                                 ggraph::geom_edge_link(edge_colour="grey50",
                                                        arrow = arrow(length = unit(2, 'mm'),
                                                                      #type = 'closed',
                                                                      ends = "last"),
                                                        edge_width=0.2) +
                                 ggraph::geom_node_point(ggplot2::aes(shape=label,
                                                              fill=label), 
                                                 size = NodeSize,
                                                 color = 'grey20',
                                                 #size=absCoef,color='grey20',
                                                 alpha=0.9) +
                                 ggplot2::scale_size_continuous(range=c(0.5,10)) +
                                 ggplot2::scale_shape_manual(values = nodeshape)+
                                 ggplot2::scale_fill_manual(values=mycolor) +
                                 ggraph::geom_node_text(ggplot2::aes(label = GroupName), 
                                                nudge_x = 0.22,
                                                nudge_y = -0.22,
                                                size=LabelSize, 
                                                color="grey30") +
                                 ggplot2::expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+
                                 ggplot2::theme_bw() +
                                 ggplot2::labs(fill = "GroupName",
                                      shape = "GroupName")+
                                 ggplot2::theme(panel.grid = element_blank(),
                                       panel.border = element_blank(),
                                       axis.line = element_blank(),
                                       axis.ticks =element_blank(),
                                       axis.text =element_blank(),
                                       axis.title = element_blank()) -> networkPlot
                             )
                           },
                           "degree-circle" = {
                             #构造数据标签的放置角度
                             number_of_bar<-nrow(Vertices1)
                             Vertices1$id<-seq(1, nrow(Vertices1))
                             angle<- 360 * (Vertices1$id-0.5) /number_of_bar 
                             Vertices1$hjust<-ifelse(angle>180, 1, 0)
                             Vertices1$angle<-ifelse(angle>180, 90-angle+180, 90-angle)
                             # 重新构造 graph 图结构的数据
                             mygraph <- igraph::graph_from_data_frame( Connect_data, vertices = Vertices1, directed = TRUE )
                             mycolor <- wesanderson::wes_palette("Darjeeling1", 
                                                                 length(Vertices1 %>% dplyr::distinct(label) %>% 
                                                                          as.matrix() %>% as.vector()),
                                                                 type = "continuous")
                             
                             ggraph::ggraph(mygraph,layout = 'linear', circular = TRUE) +
                               ggraph::geom_edge_arc(aes(edge_colour=label), 
                                             arrow = arrow(length = unit(2, 'mm'),
                                                           ends = "last"),
                                             edge_alpha=0.8, edge_width=0.5,
                                             show.legend = FALSE) +
                               ggraph::geom_node_point(aes(fill=label, 
                                                   shape=label),
                                               color = 'grey20',
                                               size = NodeSize) +
                               ggplot2::scale_size_continuous(range=c(0.5,10)) +
                               ggraph::geom_node_text(aes(x = x*1.05, y=y*1.05, label=GroupName, angle=angle,hjust=hjust,
                                                  color=label),
                                              size=LabelSize,
                                              show.legend = FALSE) +
                               ggplot2::scale_shape_manual(values = nodeshape)+
                               ggplot2::scale_fill_manual(values=mycolor) +
                               ggraph::scale_edge_color_manual(values=mycolor) +
                               ggplot2::labs(fill = "GroupName",
                                    shape = "GroupName")+
                               ggplot2::expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))+
                               ggplot2::coord_fixed()+
                               ggplot2::theme_bw()+
                               ggplot2::theme(panel.grid = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_blank(),
                                     axis.ticks =element_blank(),
                                     axis.text =element_blank(),
                                     axis.title = element_blank()) -> networkPlot
                           })

                    # save --------------------------------------------------------------------
                    Path = str_c(getwd(), "/output_", eSet$PID)
                    ggplot2::ggsave(stringr::str_c(Path, "/PPI_",Method,"_",Layout,"_",Brightness,"_",Palette,'.png'),
                            plot = networkPlot,
                            width =  20, 
                            height = 20,
                            units = "cm")
                    # ggplot2::ggsave(stringr::str_c(Path, "/PPI_",Method,"_",Layout,"_",Brightness,"_",Palette,'.svg'),
                    #                 plot = networkPlot,
                    #                 width =  20, 
                    #                 height = 20,
                    #                 units = "cm")
                    
                  },
                  'GO' = {
                    eSet$Biolink$Nodes_GO %>% 
                      dplyr::distinct(group) %>% 
                      as.matrix() %>% as.vector() -> GroupName
                    
                    eSet$Biolink$Edges_GO %>% 
                      dplyr::select(source,target,edge_type) %>% 
                      setNames(c('from','to',"value")) -> connect_data
                    
                    rbind(
                      data.frame(from = "origin",
                                 to = GroupName),
                      data.frame(to = unique(c(as.character(connect_data$to),as.character(connect_data$to)))) %>% 
                        dplyr::left_join(eSet$Biolink$Nodes_GO %>% 
                                           dplyr::select(id,group),
                                         by = c("to" = "id")) %>% 
                        dplyr::distinct(to,.keep_all = TRUE) %>% 
                        dplyr::arrange(group) %>% 
                        setNames(c("to","from")) %>% 
                        dplyr::select(from,to)
                    ) -> Edges
                    
                    data.frame(
                      name = unique(c(as.character(Edges$from),as.character(Edges$to)))
                    ) -> vertices_data
                    
                    vertices_data %>% 
                      dplyr::left_join(eSet$Biolink$Nodes_GO %>% 
                                         dplyr::select(id),
                                       by = c("name" = "id")) -> vertices_data
                    
                    vertices_data %>%
                      dplyr::mutate(Group = Edges$from[ match( vertices_data$name, Edges$to ) ],
                             id = NA) -> vertices_data
                    Leaves <- which(is.na(match(vertices_data$name,Edges$from)))
                    NLeaves <- length(Leaves)
                    vertices_data$id[Leaves] <- seq(1:NLeaves)
                    vertices_data$angle <- 90 - 360 * vertices_data$id / NLeaves
                    vertices_data$hjust <- ifelse(vertices_data$angle < -90,1,0)
                    vertices_data$angle <- ifelse(vertices_data$angle < -90,vertices_data$angle+180,
                                                  vertices_data$angle)
                    
                    c(as.character(connect_data$from),as.character(connect_data$to)) %>% 
                      as_tibble() %>% 
                      dplyr::group_by(value) %>% 
                      dplyr::summarize(n=n()) %>% 
                      # dplyr::left_join(eSet$Expo$Voca %>% 
                      #                    dplyr::select(SerialNo_Raw,
                      #                                  GroupName),by = c("value" = "SerialNo_Raw")) %>% 
                      dplyr::left_join(eSet$Biolink$Nodes_GO,
                                       by = c("value" = "id")) %>% 
                      dplyr::distinct(value,.keep_all = TRUE) -> Vertices
                    
                    colnames(Vertices) <- c("name", "n","GroupName","label")
                    
                    mygraph <- igraph::graph_from_data_frame( connect_data, vertices = Vertices, directed = FALSE )
                    
                    com <- igraph::walktrap.community(mygraph)
                    
                    Vertices1 <- Vertices %>% 
                      dplyr::mutate( group = com$membership) %>%
                      dplyr::mutate(group=as.numeric(factor(group,
                                                     levels=sort(summary (as.factor(group)),index.return=TRUE,decreasing = T)$ix,
                                                     order=TRUE)))%>%
                      dplyr::filter( group<10) %>%
                      dplyr::arrange(group,desc(n)) %>%
                      dplyr::mutate(name=factor(name, name))
                    
                    Connect_data <- connect_data %>%
                      dplyr::filter(from %in% Vertices1$name) %>%
                      dplyr::filter(to %in% Vertices1$name) %>%
                      dplyr::left_join(Vertices1,by=c('from'='name')) 
                    
                    Vertices1 %>% 
                      dplyr::distinct(label) %>% 
                      dplyr::arrange(label) %>% 
                      as.matrix() %>% as.vector() -> GroupName
                    
                    Vertices1$label <- factor(Vertices1$label,
                                              levels = GroupName)
                    
                    mygraph <- igraph::graph_from_data_frame( Connect_data, vertices = Vertices1, directed = TRUE )
                    #mycolor <- wesanderson::wes_palette("Darjeeling1", max(Vertices1$group), type = "continuous")
                    
                    # color -------------------------------------------------------------------
                    Node_Shape <- data.frame(
                      exposure = 24,
                      disease = 23,
                      go_by_proteome = 21,
                      go_by_metabolome = 21,
                      go = 21,
                      go_by_metabolome_proteome = 21
                    )
                    
                    Node_Shape %>% 
                      dplyr::select(all_of(GroupName)) %>% 
                      as.matrix() %>% 
                      as.vector() -> nodeshape
                    mycolor <- c("#FF0000","turquoise1","#FBA3F3","#B7EF97",
                                 "#F2AD00","#F98400","#FDF7BB") %>% 
                      head(length(GroupName))
                    #Darjeeling1
                    # mycolor <- wesanderson::wes_palette("Darjeeling1", 
                    #                                     length(Vertices1 %>% dplyr::distinct(label) %>% 
                    #                                              as.matrix() %>% as.vector()),
                    #                                     type = "discrete")
                    # mycolor <- sample(mycolor, length(mycolor))
                    #scales::show_col(mycolor)
                    # labelsize ---------------------------------------------------------------
                    if(nrow(Vertices1) <= 100){
                      LabelSize = 3
                    }else if(nrow(Vertices1) > 100 & nrow(Vertices1) <= 300){
                      LabelSize = 2.5
                    }else if(nrow(Vertices1) >300 & nrow(Vertices1) <= 450){
                      LabelSize = 2
                    }else if(nrow(Vertices1) >450 & nrow(Vertices1) <= 600){
                      LabelSize = 1.5
                    }else{
                      LabelSize = 1
                    }
                    if(nrow(Vertices1) <= 100){
                      NodeSize = 4
                    }else if(nrow(Vertices1) > 100 & nrow(Vertices1) <= 300){
                      NodeSize = 3.5
                    }else if(nrow(Vertices1) >300 & nrow(Vertices1) <= 450){
                      NodeSize = 3
                    }else if(nrow(Vertices1) >450 & nrow(Vertices1) <= 600){
                      NodeSize = 2.5
                    }else{
                      NodeSize = 2
                    }
                    # Plot --------------------------------------------------------------------
                    switch(Layout,
                           "force-directed" = {
                             set.seed(123)
                             ddpcr::quiet(
                               ggraph::ggraph(mygraph,layout='fr') + 
                                 ggraph::geom_edge_link(edge_colour="grey50",
                                                        arrow = arrow(length = unit(2, 'mm'),
                                                                      #type = 'closed',
                                                                      ends = "last"),
                                                        edge_width=0.4) +
                                 ggraph::geom_node_point(ggplot2::aes(shape=label,
                                                              fill=label), 
                                                 size = NodeSize,
                                                 color = 'grey20',
                                                 #size=absCoef,color='grey20',
                                                 alpha=0.9) +
                                 ggplot2::scale_size_continuous(range=c(0.5,10)) +
                                 ggplot2::scale_shape_manual(values = nodeshape)+
                                 ggplot2::scale_fill_manual(values=mycolor) +
                                 ggraph::geom_node_text(ggplot2::aes(label = GroupName), 
                                                nudge_x = 0.22,
                                                nudge_y = -0.12,
                                                size=LabelSize, 
                                                color="grey30") +
                                 ggplot2::expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+
                                 ggplot2::theme_bw() +
                                 ggplot2::labs(fill = "GroupName",
                                      shape = "GroupName")+
                                 ggplot2::theme(panel.grid = element_blank(),
                                       panel.border = element_blank(),
                                       axis.line = element_blank(),
                                       axis.ticks =element_blank(),
                                       axis.text =element_blank(),
                                       axis.title = element_blank()) -> networkPlot
                             )
                           },
                           "degree-circle" = {
                             #构造数据标签的放置角度
                             number_of_bar<-nrow(Vertices1)
                             Vertices1$id<-seq(1, nrow(Vertices1))
                             angle<- 360 * (Vertices1$id-0.5) /number_of_bar 
                             Vertices1$hjust<-ifelse(angle>180, 1, 0)
                             Vertices1$angle<-ifelse(angle>180, 90-angle+180, 90-angle)
                             # 重新构造 graph 图结构的数据
                             mygraph <- igraph::graph_from_data_frame( Connect_data, vertices = Vertices1, directed = TRUE )
                             mycolor <- wesanderson::wes_palette("Darjeeling1", 
                                                                 length(Vertices1 %>% dplyr::distinct(label) %>% 
                                                                          as.matrix() %>% as.vector()),
                                                                 type = "continuous")
                             
                             ggraph(mygraph,layout = 'linear', circular = TRUE) +
                               ggraph::geom_edge_arc(aes(edge_colour=label), 
                                             arrow = arrow(length = unit(2, 'mm'),
                                                           ends = "last"),
                                             edge_alpha=0.8, edge_width=0.5,
                                             show.legend = FALSE) +
                               ggraph::geom_node_point(aes(fill=label, 
                                                   shape=label),
                                               color = 'grey20',
                                               size = NodeSize) +
                               ggplot2::scale_size_continuous(range=c(0.5,10)) +
                               ggraph::geom_node_text(aes(x = x*1.05, y=y*1.05, label=GroupName, angle=angle,hjust=hjust,
                                                  color=label),
                                              size=LabelSize,
                                              show.legend = FALSE) +
                               ggplot2::scale_shape_manual(values = nodeshape)+
                               ggplot2::scale_fill_manual(values=mycolor) +
                               ggraph::scale_edge_color_manual(values=mycolor) +
                               ggplot2::labs(fill = "GroupName",
                                    shape = "GroupName")+
                               ggplot2::expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))+
                               ggplot2::coord_fixed()+
                               ggplot2::theme_bw()+
                               ggplot2::theme(panel.grid = element_blank(),
                                     panel.border = element_blank(),
                                     axis.line = element_blank(),
                                     axis.ticks =element_blank(),
                                     axis.text =element_blank(),
                                     axis.title = element_blank()) -> networkPlot
                           })
                    
                    # save --------------------------------------------------------------------
                    # ggplot2::ggsave(stringr::str_c(Path, "/GO_",Method,"_",Layout,"_",Brightness,"_",Palette,'.png'),
                    #        plot = networkPlot,
                    #        width =  20, 
                    #        height = 20,
                    #        units = "cm")
                    # ggplot2::ggsave(stringr::str_c(Path, "/GO_",Method,"_",Layout,"_",Brightness,"_",Palette,'.svg'),
                    #                 plot = networkPlot,
                    #                 width =  20, 
                    #                 height = 20,
                    #                 units = "cm")
                    
                  })
         })
  
  #save R package data-------------------
  networkPlot -> eSet$Plot[[stringr::str_c(Mode,"_",Layout,"_",Brightness,"_",Palette)]]
  
  #save R command log-------------------
  eSet$RCommandLog <- eSet$AddCommand(stringr::str_c("\n eSet <- VizBioLink(eSet = eSet, \n",
                                                     " Mode = '", Mode, "', \n",
                                                     " Layout = '", Layout, "', \n",
                                                     " Brightness = '", Brightness, "', \n",
                                                     " Palette = '", Palette,"') \n")
                                      )
  
  eSet$RCommandLog %>% 
    as_tibble() %>% 
    purrr::set_names("R commands") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/rcommands log.txt"))
  
  #print message and save log
  message("Complete VizBiolink in ", Mode, " mode. ", NowTime, "\n")
  eSet$ExcecutionLog  <- eSet$AddLog(stringr::str_c("Complete VizBiolink in ", Mode, " mode. ", NowTime))
  
  eSet$ExcecutionLog %>% 
    as_tibble() %>% 
    purrr::set_names("running log") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/running log.txt"))  
  
  # eSet %>% 
  #   save(file = stringr::str_c(eSet$FileDirOut,"/eSet.Rdata"))
  
  tictoc::toc()
  
  return(eSet)
}


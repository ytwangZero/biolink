
#' @title End the module analysis
#' @description End the module analysis
#' @usage Exit(eSet = eSet)
#' @param eSet An R6 class object.
#'
#' @details 
#' @return
#' @export
#' @examples eSet <- InitBioLink(1)
#' eSet <- LoadBioLink(eSet = eSet,UseExample = "example#1", FileDirLink = NULL)
#' eSet <- ConvToExpoID(eSet = eSet)
#' eSet <- Biolink(eSet = eSet,Mode = "PPI",ChemCas = "default",ChemInchikey = "default",
#' DiseaseID = "default",MetabolomeID = "default",MetBiospec = "blood",ProteomeID = "default",GenomeID = "default")
#' eSet <- VizBioLink(eSet = eSet,Mode = "PPI",Layout = "force-directed",Brightness = "dark",Palette = "default2")
#' Exit(eSet = eSet)
#' @author Mingliang Fang, Bin Wang (corresponding author), Yuting Wang, Ting Wu
Exit <- function(eSet = eSet){
  
  #if(file.exists(eSet$FileDirIn)){
  #  dir_delete(eSet$FileDirIn)
  #}
  
  if(file.exists(eSet$FileDirOut)){
    fs::file_delete(eSet$FileDirOut)
  }
  
  #if(file.exists(eSet$FolderOut)){
  #  fs::file_delete(eSet$FolderOut)
  #}
}
  
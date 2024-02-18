
#' @title Initialize BioLink module
#' @param seednum set a random seed number.
#'
#' @description Initialize BioLink module analysis. It can generate an R6 class object 
#'     integrating all the analysis information.
#' @usage InitBioLink(123)
#' @param 
#' @details BioLink module is designed to find the biological relationships between 
#'     exposure factors and health outcome. This module adopts the most frequently-used and 
#'     authoritative databases, e.g., T3DB, CTD, ToxCast, StringDB, STITCH, KEGG, and GO.
#' @return An R6 class object. 
#' @export
#' @examples eSet <- InitBioLink(123)
#' @author Mingliang Fang, Bin Wang (corresponding author), Yuting Wang, Ting Wu
InitBioLink <- function(seednum){

tictoc::tic()
lubridate::now() %>% 
  stringr::str_replace_all(":",".") %>% 
  stringr::str_replace_all("-",".") -> NowTime 

  #create R6 class eSet ----------------------------------------------------
    eSet <- R6::R6Class(
      "eSet",
      public = list(
        PID = NULL,
        EpiDesign = NULL,
        DataStr = NULL,
        FileDirIn = NULL,
        FileDirOut = NULL,
        ExcecutionLog = NULL,
        RCommandLog = NULL,
        Expo = list(Data = NULL),
        AddLog = function(x){
          self$ExcecutionLog = c(self$ExcecutionLog, 
                                 stringr::str_c(x))},
        
        AddCommand = function(x){
          self$RCommandLog = c(self$RCommandLog, 
                                 stringr::str_c(x))}
        ),
      
      lock_class = FALSE,
      lock_objects = FALSE
      )

    eSet <- eSet$new()

    # generate default PID
    set.seed(seednum)
    pid_char = sample(c("A","B","C","D","E","F","G","H","J","K","L","M","N","P","Q","R","S","T","U","V","W","X","Y","Z"),
                      6, replace = TRUE) %>% 
      stringr::str_c(collapse = "")
    
    eSet$PID <- stringr::str_c(lubridate::now() %>% 
                                 stringr::str_remove_all(" ") %>% 
                                 stringr::str_remove_all("-") %>%
                                 stringr::str_remove_all(":") %>%
                                 stringr::str_sub(9,-1),
                                 collapse = "") %>% 
      stringr::str_c(pid_char)
    

    # initialize the folder ---------------------------------------------------
    stringr::str_c(getwd(), "/output_", eSet$PID) -> eSet$FileDirOut
    dir.create(eSet$FileDirOut)

    stringr::str_c(getwd(), "/database") %>% 
      dir.create()
    message("Download data files from our GitHub and move them to the newly generated folder named database.")

  # save R command log ----
  eSet$RCommandLog <- eSet$AddCommand(stringr::str_c("\n eSet <- InitExpoData(PID = 'Any ID your like', \n FileDirIn = 'Any input file directory your like', \n FileDirOut = 'Any output file directory your like')\n"))
  
  eSet$RCommandLog %>% 
    as_tibble() %>% 
    purrr::set_names("R commands") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/rcommands log.txt"))
  
  #save log 
  message("Complete initializing the BioLink module.", NowTime, "\n")
  
  eSet$ExcecutionLog  <- eSet$AddLog(stringr::str_c("Complete initializing the BioLink module.", 
                                                    NowTime))
  ddpcr::quiet(eSet$ExcecutionLog %>% 
          as_tibble() %>% 
          purrr::set_names("running log") %>% 
          data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/running log.txt"))
        )
  
  # save eSet ----
  eSet %>% 
    save(file = str_c(eSet$FileDirOut,"/eSet.Rdata"))
  
  tictoc::toc()  
  
  return(eSet)
}

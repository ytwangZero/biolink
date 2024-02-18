#' @title Load data file for BioLink module
#' @description Load data file for BioLink module. 
#' @usage LoadBioLink(PID, UseExample = "default", FileDirIn=NULL)
#'
#' @param eSet An R6 class object.
#' @param UseExample chr. Method of uploading data. If "default",user should upload their own data files,
#'                   or use "example#1" provided by this module.  
#' @param FileDirLink  chr. the name of the data file provided by user. Remember move it into the working 
#'                     directory(e.g. "eg_biolink.xlsx").
#'    
#' @details 
#' @return An R6 class object containing the input data.
#' @export
#' @examples eSet = InitBioLink()
#'           eSet = LoadBioLink(PID = res$PID, UseExample = "example#1")
#' @author Mingliang Fang, Bin Wang (corresponding author), Yuting Wang, Ting Wu
LoadBioLink <- function(eSet,
                        UseExample = "default",
                        FileDirLink
                        ){
  tictoc::tic()
  
  lubridate::now() %>% 
    stringr::str_replace_all(":",".") %>% 
    stringr::str_replace_all("-",".") -> NowTime 

  #find the detailed FileDirExpo
  switch(UseExample,
         "default" = {
           readxl::read_xlsx(Expodata) -> eSet$Expo$Data
         },
         "example#1" = {
           "database/eg_biolink.xlsx" -> FileDirExpo
           readxl::read_xlsx(FileDirExpo) -> eSet$Expo$Data
         }
  )
  
   #check the name of Exposome file-------------------------------------
   all(eSet$Expo$Data %>% 
         names() == c("SerialNo", 
                "FullName",
                "GroupName", 
                "DiseaseID", 
                "ExposureID", 
                "MetabolomeID", 
                "ProteomeID",
                "GenomeID")) -> flag1
  
  if(!flag1){
    message("Error: Please use the templete file of Biolink!")
    print(Flag_Error)
    }
     
  #check GroupName-----------------------------------------------------------------
  eSet$Expo$Data %>% 
    dplyr::filter(str_detect(SerialNo, "X")) %>% 
    dplyr::filter(is.na(GroupName)) %>% 
    nrow() == 0 -> flag2
    
  if(!flag2){
    message("Error: Please fill the GroupName of all X-variables")
    print(Flag_Error)
    }
  
  #check GroupName-----------------------------------------
  c("disease", "exposure", "metabolome", "proteome","genome") -> groupname_fixed
    
  eSet$Expo$Data %>% 
    dplyr::distinct(GroupName) %>% 
      .$GroupName -> Var_GroupName
  
    if(all(Var_GroupName %in% groupname_fixed)){
      flag3 = T
      }else{
        message("Error: The GroupNames are not 'disease', exposure', metabolome', 'proteome'")
        flag3 = F
        print(Flag_Error)
        }
    

  #summary the flags
  if(all(c(flag1,flag2,flag3))){

    ddpcr::quiet(
      eSet$Expo$Data %>%
        vroom::vroom_write(stringr::str_c(eSet$FileDirOut, "/ExpoData.csv"),
                           delim = ",")
    )

     }else{
       message("Error: Fail to load Exposome file! ", NowTime, "\n")
       eSet$ExcecutionLog  <- eSet$AddLog(stringr::str_c("Error: Fail to load Exposome file!", NowTime))
       eSet$Expo$Data <- NULL
       
       print(Flag_Error)
     }


  #save R command and running logs
  eSet$RCommandLog <- eSet$AddCommand(stringr::str_c("\n eSet <- LoaBioLink(eSet = eSet,  \n", 
                                                     "UseExample = '", UseExample, "',  \n",
                                                     "Expodata = Expodata'", "') \n")
                                      )
  
  eSet$RCommandLog %>% 
    as_tibble() %>% 
    purrr::set_names("R commands") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/rcommands log.txt"))
  
  
  message("Complete the data loading! ", NowTime, "\n")
  eSet$ExcecutionLog  <- eSet$AddLog(stringr::str_c("Complete the data loading! ", NowTime))
  
  eSet$ExcecutionLog %>% 
    as_tibble() %>% 
    purrr::set_names("running log") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/running log.txt")) 
  
  eSet %>% 
    save(file = str_c(eSet$FileDirOut,"/eSet.Rdata"))
  
  tictoc::toc()
  
  return(eSet)
  }




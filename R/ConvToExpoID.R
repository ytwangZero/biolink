
#' @title Convert different IDs to the unified Biolink ID
#' @description Convert the IDs of exposure, chemicals, metabolites, genes, or proteins to the unified Biolink ID, 
#' @param eSet An R6 class object.
#'   
#' @usage ConvToExpoID(eSet = eSet)
#' @details 
#' @return An R6 class object containing the converted ID information.
#' @export
#' @examples eSet <- InitBioLink(1)
#' eSet <- LoadBioLink(eSet = eSet,UseExample = "example#1", FileDirLink = NULL)
#' eSet <- ConvToExpoID(eSet = eSet)
#' @author Mingliang Fang, Bin Wang (corresponding author), Yuting Wang, Ting Wu
ConvToExpoID <- function(eSet){
  
  tictoc::tic()
  
  lubridate::now() %>% 
    stringr::str_replace_all(":",".") %>% 
    stringr::str_replace_all("-",".") -> NowTime 
  
  #check the file existance
  file.exists(str_c(eSet$FileDirOut, "/ExpoData.csv")) -> flag0
  
  if(!flag0){
    
    #save R command log
    eSet$RCommandLog <- eSet$AddCommand(str_c("\n eSet <- ConvToExpoID(eSet = eSet) \n")
    )
    
    eSet$RCommandLog %>% 
      as_tibble() %>% 
      purrr::set_names("R commands") %>% 
      data.table::fwrite(str_c(eSet$FileDirOut,"/rcommands log.txt"))
    
    #save running log
    message("Error: Biolink data file is missing! Please add them, or the following functions cann't run. ", NowTime, "\n")
    eSet$EXEecutionLog  <- eSet$AddLog(str_c("Error: Biolink data or vocabulary file is missing! ", NowTime))
    
    eSet$EXEecutionLog %>% 
      as_tibble() %>% 
      purrr::set_names("running log") %>% 
      data.table::fwrite(str_c(eSet$FileDirOut,"/running log.txt")) 
    
    eSet %>% 
      save(file = str_c(eSet$FileDirOut,"/eSet.Rdata"))
    
    print(Flag_Error)
    
    tictoc::toc()
    
    return(eSet)
    
  }else{
    
    #load database
    ddpcr::quiet(
    vroom::vroom("database/DB_Disease.csv",
                   show_col_types = F) %>% 
      dplyr::select(-disease) -> DB_Disease
     )
    
    ddpcr::quiet(
      vroom::vroom( "database/DB_Protein_Unfold_EX.csv",
                   show_col_types = F) -> DB_Protein
      )
    
    ddpcr::quiet(
      vroom::vroom("database/DB_Exposure_Unfold_EX.csv",
                   show_col_types = F) -> DB_Exposure
      )
    
    ddpcr::quiet(
      vroom::vroom("database/DB_Metabolite_Unfold.csv",
                   show_col_types = F)  %>% 
        distinct(EXM, Entry) -> DB_Metabolite
    )
    
    ddpcr::quiet(
      vroom::vroom("database/DB_Gene.csv",
                   show_col_types = F) -> DB_Gene
    )

    #convert disease ID
    eSet$Expo$Data %>% 
      dplyr::filter(GroupName %in% c("disease")) %>% 
      dplyr::left_join(DB_Disease %>% 
                         na.omit(), 
                       by = c("DiseaseID" = "disease.id")) %>% 
      dplyr::rename(EX = EXD) -> df_disease

    #convert exposure ID
    eSet$Expo$Data %>% #by cas
      dplyr::filter(!ExposureID %in% c(NA)) %>% 
      dplyr::filter(ExposureID %in% unique(DB_Exposure$cas.rn)) -> temp1
    
    if(length(temp1) > 0){
      temp1 %>% 
        dplyr::left_join(DB_Exposure %>% 
                           dplyr::select(EXE, cas.rn) %>% 
                           dplyr::distinct(EXE, cas.rn) %>% 
                           dplyr::filter(cas.rn %in% eSet$Expo$Data$ExposureID) %>% 
                           na.omit(), 
                         by = c("ExposureID" = "cas.rn")) %>% 
        dplyr::rename(EX = EXE) -> DB_exposure_cas
      }else{
        DB_exposure_cas <- NULL
    }

    eSet$Expo$Data %>% #by inchikey
      dplyr::filter(!ExposureID %in% c(NA)) %>% 
      dplyr::filter(ExposureID %in% unique(DB_Exposure$inchikey)) -> temp2
    
    if(length(temp2) > 0){
      temp2 %>% 
        dplyr::left_join(DB_Exposure %>% 
                           dplyr::select(EXE, inchikey) %>% 
                           dplyr::distinct(EXE, inchikey) %>% 
                           dplyr::filter(inchikey %in% eSet$Expo$Data$ExposureID) %>% 
                           na.omit(), 
                         by = c("ExposureID" = "inchikey")) %>% 
        dplyr::rename(EX = EXE) -> DB_exposure_inchikey
      }else{
        DB_exposure_inchikey <- NULL
    }

    #convert metabolite ID
    eSet$Expo$Data %>% #by cas
      dplyr::filter(!MetabolomeID %in% c(NA)) %>% 
      dplyr::filter(MetabolomeID %in% unique(DB_Metabolite$Entry)) -> temp3
    
    if(nrow(temp3) > 0){
      temp3 %>% 
        dplyr::left_join(DB_Metabolite %>% 
                           dplyr::filter(Entry %in% eSet$Expo$Data$MetabolomeID) %>% 
                           na.omit(), 
                         by = c("MetabolomeID" = "Entry")) %>% 
        dplyr::rename(EX = EXM) -> df_metabolite
    }else{
      df_metabolite <- NULL
    }

    
    #convert protein ID
    eSet$Expo$Data %>% #by uniprot
      dplyr::filter(!ProteomeID %in% c(NA)) %>% 
      dplyr::filter(ProteomeID %in% unique(DB_Protein$Uniprot_KB)) -> temp4
    
    if(nrow(temp4) > 0){
      temp4 %>% 
        dplyr::left_join(DB_Protein %>% 
                           dplyr::select(EXP, Uniprot_KB) %>% 
                           dplyr::distinct(EXP, Uniprot_KB) %>% 
                           na.omit(), 
                         by = c("ProteomeID" = "Uniprot_KB")) %>% 
        dplyr::rename(EX = EXP) -> df_protein_uniprot
    }else{
      df_protein_uniprot <- NULL
    }
    
    eSet$Expo$Data %>% #by ensp
      dplyr::filter(!ProteomeID %in% c(NA)) %>% 
      dplyr::filter(ProteomeID %in% unique(DB_Protein$ENSP)) -> temp5
    
    if(length(temp5) > 0){
      temp5 %>% 
        dplyr::left_join(DB_Protein %>% 
                           dplyr::select(EXP, ENSP) %>% 
                           dplyr::distinct(EXP, ENSP) %>% 
                           na.omit(), 
                         by = c("ProteomeID" = "ENSP")) %>% 
        dplyr::rename(EX = EXP) -> df_protein_ensp
    }else{
      df_protein_ensp <- NULL
    }

    ##convert gene ID
    as.character(DB_Gene$GeneID) -> DB_Gene$GeneID
    eSet$Expo$Data %>% 
      dplyr::filter(GroupName %in% c("genome")) %>% 
      dplyr::left_join(DB_Gene %>% 
                         na.omit(), 
                       by = c("GenomeID" = "GeneID")) %>% 
      dplyr::rename(EX = EXG) %>% 
      select(-c(Gene_Names, database)) -> df_gene

    rbind(df_disease %>% dplyr::select(-PreferredName),
          DB_exposure_cas,
          DB_exposure_inchikey,
          df_metabolite,
          df_protein_uniprot,
          df_protein_ensp,
          df_gene) -> eSet$Expo$Data 
    
    #save R exposome data
    # ddpcr::quiet(
    #   eSet$Expo$Data %>% 
    #     vroom::vroom_write(str_c(eSet$FileDirOut, "/ExpoData.csv"),
    #                        delim = ",")
    # )
    
    #save R command and running logs
    eSet$RCommandLog <- eSet$AddCommand(str_c("\n eSet <- ConvToExpoID(eSet = eSet) \n")
    )
    
    eSet$RCommandLog %>% 
      as_tibble() %>% 
      purrr::set_names("R commands") %>% 
      data.table::fwrite(str_c(eSet$FileDirOut,"/rcommands log.txt"))
    
    message("Complete the ID concersion! ", NowTime, "\n")
    eSet$EXEecutionLog  <- eSet$AddLog(str_c("Complete the ID concersion! ", NowTime))
    
    eSet$EXEecutionLog %>% 
      as_tibble() %>% 
      purrr::set_names("running log") %>% 
      data.table::fwrite(str_c(eSet$FileDirOut,"/running log.txt")) 
    
    eSet %>% 
      save(file = str_c(eSet$FileDirOut,"/eSet.Rdata"))
    
    tictoc::toc()
    
    return(eSet)
  }
}




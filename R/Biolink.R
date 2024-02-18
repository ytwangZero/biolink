#' @title Build the biological link 
#' @description Build the biological link between the exposures and diseases
#' @usage BioLink(eSet, Mode="PPI",ChemCas="default", ChemInchikey="default",
#'    DiseaseID="default", MetabolomeID="default", MetBiospec="blood", ProteomeID="default",GenomeID = "default")
#' @param eSet An R6 class object.
#' @param Mode chr. Method to build the biological link between exposure and disease. 
#'             Only PPI is currently supported.
#' @param ChemCas chr. CAS Registry Number of chemicals. Default means using the values in the input data file. 
#'    Users can also copy the part of them by clicking "Available vars". It should be noted that 
#'    there is fixed format for the entering characters separated with comma and without space, 
#'    e.g., "7440-43-9,333-41-5,20461-54-5". 
#' @param ChemInchikey chr. InChiKey serial number of chemicals. Default means using the values in the input data file. 
#'    It should be noted that there is fixed format for the entering characters separated with comma and without space, 
#'    e.g., "WABPQHHGFIMREM-UHFFFAOYSA-N,BTAGRXWGMYTPBY-UHFFFAOYSA-N,IAKOZHOLGAGEJT-UHFFFAOYSA-N". 
#' @param DiseaseID chr. ID of the concerned diseases. Both IDs from OMIM (e.g., OMIM:220100) and 
#'    MESH (e.g., MESH:C536409) are accepted. It should be noted that there is fixed format 
#'    for the entering characters separated with comma and without space, 
#'    e.g., "OMIM:244200,MESH:C536409,OMIM:181500". 
#' @param MetabolomeID chr. KEGG entry number of metabolites. Default means using the values 
#'    in the input data file. Users can also copy the part of them by clicking "Available vars". 
#'    It should be noted that there is fixed format for the entering characters separated with comma 
#'    and without space, e.g., "C00022,C00117,C00794". 
#' @param MetBiospec chr. Biological sample matrix for the metabolome analysis. 
#'    Options include "Blood" and "Urine".
#' @param ProteomeID chr. Protein ID. Both IDs of Ensembl and UniProt are accepted. 
#'    Default means using the values in the input data file. Users can also copy the part of them 
#'    by clicking "Available vars". It should be noted that there is fixed format for the entering characters
#'    separated with comma and without space, e.g., "Q9Y3X0,Q8N5I3,ENSP00000000233,ENSP00000000412". 
#' @param GenomeID chr. Gene ID. NCBI Gene identifier is accepted. users also can get it from UniProt or CTD databases.
#'    
#' @details 
#' @return An R6 class object containing the edges and nodes of the biological link.
#' @export
#' @examples eSet <- InitBioLink(1)
#' eSet <- LoadBioLink(eSet = eSet,UseExample = "example#1", FileDirLink = NULL)
#' eSet <- ConvToExpoID(eSet = eSet)
#' eSet <- Biolink(eSet = eSet,Mode = "PPI",ChemCas = "default",ChemInchikey = "default",
#' DiseaseID = "default",MetabolomeID = "default",MetBiospec = "blood",ProteomeID = "default",GenomeID = "default")
#' @author Mingliang Fang, Bin Wang (corresponding author), Yuting Wang, Ting Wu
Biolink <- function(eSet,
                    Mode,
                    ChemCas = "default",
                    ChemInchikey = "default",
                    DiseaseID = "default",
                    MetabolomeID = "default",
                    MetBiospec = "blood",
                    ProteomeID = "default",
                    GenomeID = "default"
                    ){
  
  tictoc::tic()
  
  lubridate::now() %>% 
    stringr::str_replace_all(":",".") %>% 
    stringr::str_replace_all("-",".") -> NowTime 
  
  # Path = stringr::str_c(getwd(),"/output_",eSet$PID, "/Biolink")
  # if(!file.exists(Path)) {dir.create(Path)}
  message("It will take 5 to 10 minutes, please wait patiently.")
  #define the variable
    if(ChemCas != "default"){
      ChemCas -> eSet$Expo$Data$ChemCas}

    if(ChemInchikey != "default"){
      ChemInchikey -> eSet$Expo$Data$ChemInchikey}
 
    if(ChemCas != "default" & ChemInchikey == "default"){
       eSet$Expo$Data$ChemCas  -> eSet$Expo$Data$ExposureID
      }else if (ChemCas == "default" & ChemInchikey != "default") {
       eSet$Expo$Data$ChemInchikey  -> eSet$Expo$Data$ExposureID
      }else if (ChemCas != "default" & ChemInchikey != "default") {
     c(eSet$Expo$Data$ChemCas,
     eSet$Expo$Data$ChemInchikey) -> eSet$Expo$Data$ExposureID
      }

    if(DiseaseID != "default"){
       DiseaseID -> eSet$Expo$Data$DiseaseID
    }

    if(MetabolomeID != "default"){
       MetabolomeID -> eSet$Expo$Data$MetabolomeID
    }
    
    if(ProteomeID != "default"){
       ProteomeID -> eSet$Expo$Data$ProteomeID
    }
    if(GenomeID != "default"){
       GenomeID -> eSet$Expo$Data$GenomeID
    }

  #Build Biolink model
  switch (Mode, 
    "PPI" = { 
  # load databases -----------------------------------------------
  # ddpcr::quiet(
  #   vroom::vroom(str_c(getwd(), "/data/DB_Disease.csv"),
  #                  show_col_types = F) %>% 
  #     dplyr::select(-disease) -> DB_DiseaseDB
  #       ) 
  
  # ddpcr::quiet(
  #   vroom::vroom(str_c(getwd(), "/data/DB_Protein_Unfold.csv"),
  #                show_col_types = F) %>% 
  #     dplyr::distinct(EXP, ENSP, Uniprot_KB) -> DB_Protein 
  # )
  
  # ddpcr::quiet(
  #   vroom::vroom(str_c(getwd(), "/data/DB_Exposure_Unfold.csv"),
  #                #col_select = c("EXE", "inchikey", "cas.rn"),
  #                show_col_types = F)  %>% 
  #     distinct(EXE, inchikey, cas.rn) -> DB_Exposure
  # )
  
  ddpcr::quiet(
    vroom::vroom("database/DB_Metabolite_Unfold.csv",
                 col_select = c("EXM", "Entry"),
                 show_col_types = F)  %>% 
      distinct(EXM, Entry) -> DB_Metabolite)
  
  ddpcr::quiet(
    vroom::vroom("database/DB_Enzyme_Unfold.csv", #to be improved
                     show_col_types = F) -> DB_Enzyme)


  # ddpcr::quiet(
  #   vroom::vroom("data/DB_Exposure_Protein.csv", #to be improved
  #                           #n_max = 5000000,
  #                    show_col_types = F) -> DB_Exposure_Protein
  #   )
  
  ddpcr::quiet(vroom::vroom("database/DB_Exposure_Protein_EX.csv", 
                   show_col_types = F) -> DB_Exposure_Protein_EX)

  ddpcr::quiet(vroom::vroom("database/DB_Protein_Protein_EX.csv",
                     show_col_types = F) -> DB_Protein_Protein_EX)   
  
  ddpcr::quiet(vroom::vroom("database/DB_Protein_Disease_EX.csv",
                     show_col_types = F) -> DB_Protein_Disease_EX) 
  
  ddpcr::quiet(vroom::vroom("database/DB_Gene_Disease_EX.csv",
                     show_col_types = F) -> DB_Gene_Disease_EX) 

  ddpcr::quiet(vroom::vroom("database/DB_Gene_Protein_EX.csv",
                     show_col_types = F) -> DB_Gene_Protein_EX) 
  
  ddpcr::quiet(vroom::vroom("database/DB_Exposure_Gene_EX.csv",
                     show_col_types = F) -> DB_Exposure_Gene_EX)

  #step1: filter EXD of uniprotKB by EX
  DB_Protein_Disease_EX %>% 
    dplyr::filter(EXD %in% eSet$Expo$Data$EX) %>% 
  dplyr::rename(source = EXP,
     target = EXD,
     database = Database
       ) %>% 
  dplyr::mutate(interaction = "involvement",
  source.class = "protein",
  target.class = "disease"
     ) %>% 
  dplyr::select(source,
                  target,
                  interaction,
                  source.class,
                  target.class,
                  database
    ) %>% 
    dplyr::distinct(source,
             target,
             .keep_all = T) %>% 
    na.omit()  ->  edges_protein_disease  

  #step2: filter edges by matching "EXP.to" of DB_Protein_Protein to "source" of edges_protein_disease
  DB_Protein_Protein_EX %>% 
    dplyr::filter(EXP.to %in% unique(edges_protein_disease$source)) %>% 
    dplyr::rename(source = EXP.from,
           target = EXP.to,
           interaction = remarks
    ) %>% 
    dplyr::mutate(source.class = "protein",
    target.class = "protein"
    ) %>% 
     dplyr::select(source,
                  target,
                  interaction,
                  source.class,
                  target.class,
                  database
    ) %>% 
    dplyr::distinct(source,
             target,
             .keep_all = T
             ) %>% 
    na.omit() -> edges_protein_protein
  
  #step3: filter "EXE" of DB_Exposure_Protein_EX by EX
  DB_Exposure_Protein_EX %>% 
    dplyr::filter(EXE %in%  eSet$Expo$Data$EX
    )  %>% #filter based on the target chemical
    dplyr::rename(source = EXE,
           target = EXP,
           interaction = remarks
    ) %>% 
    na.omit() %>% 
    dplyr::mutate(source.class = "chemical",
    target.class = "protein"
    ) %>% 
    dplyr::select(source,
                  target,
                  interaction,
                  source.class,
                  target.class,
                  database
    ) %>% 
    dplyr::distinct(source,
             target,
             .keep_all = T) -> edges_exposure_protein

  #step4: filter edges by matching "EXP" of DB_Exposure_Protein_EX with  either "source" or "target" of edges_protein_disease
  edges_exposure_protein %>% 
    dplyr::filter(target  %in%  unique(c(edges_protein_protein$source, #by ppi
                                         edges_protein_protein$target,
                                         edges_protein_disease$source))
    )  %>% 
    na.omit() -> edges_exposure_protein_1
  
  #step5: re-filter "source" of edges_protein_protein edges by matching "target" of edges_exposure_protein_1
  edges_protein_protein %>%  
    dplyr::filter(source %in% unique(edges_exposure_protein_1$target) &
    target %in% unique(edges_protein_disease$source)) -> edges_protein_protein_1

  #step6: re-filter "source" of edges_protein_disease by matching "source" or "target" of edges_protein_protein_1
  edges_protein_disease %>% 
    dplyr::filter(source %in% unique(edges_exposure_protein_1$target) |
    source %in% unique(edges_protein_protein_1$target)) -> edges_protein_disease_1

  #step7: re-filter "target" of edges_protein_protein_1 by matching "source" of edges_protein_disease_1
  edges_protein_protein_1 %>% 
    dplyr::filter(target %in% unique(edges_protein_disease_1$source)) -> edges_protein_protein_2

  #step8: re-filter "target" of edges_exposure_protein_1 by matching "source" and "target" of edges_protein_protein_2
  edges_exposure_protein_1 %>% 
    dplyr::filter(target %in% unique(c(edges_protein_disease_1$source, edges_protein_protein_2$source))
    ) -> edges_exposure_protein_2

  #step9: filter EXG, EXD of CTD by EX
  DB_Gene_Disease_EX %>% 
    dplyr::filter(EXD %in% eSet$Expo$Data$EX) %>% 
    dplyr::filter(EXG %in% eSet$Expo$Data$EX) %>%
    dplyr::rename(source = EXG,
     target = EXD,
     interaction = DirectEvidence) %>% 
    dplyr::mutate(source.class = "gene",
     target.class = "disease") %>% 
    dplyr::select(source,
                  target,
                  interaction,
                  source.class,
                  target.class,
                  database
    ) %>% 
    dplyr::distinct(source,target,.keep_all = T) %>% 
    na.omit()  ->  edges_gene_disease_1  

  #step10: filter "EXE" of DB_Exposure_Gene_EX by EX
  DB_Exposure_Gene_EX %>% 
    dplyr::filter(EXE %in%  eSet$Expo$Data$EX
    )  %>% #filter based on the target chemical
    dplyr::rename(source = EXE,
           target = EXG,
           interaction = Interaction
    ) %>% 
    na.omit() %>% 
    dplyr::mutate(source.class = "chemical",
    target.class = "gene"
    ) %>% 
    dplyr::select(source,
                  target,
                  interaction,
                  source.class,
                  target.class,
                  database
    ) %>% 
    dplyr::distinct(source,
             target,
             .keep_all = T) -> edges_exposure_gene_1
  
  #step11: re-filter "source" of edges_gene_disease_1 by matching "target" of edges_exposure_gene_1
  edges_gene_disease_1 %>% 
    dplyr::filter(source %in% unique(edges_exposure_gene_1$target)) -> edges_gene_disease_2
  #step12: re-filter "target" of edges_exposure_gene_1 by matching "source" of edges_gene_disease_2
  edges_exposure_gene_1 %>% 
    dplyr::filter(target %in% unique(edges_gene_disease_2$source)) -> edges_exposure_gene_2
 

  #combine all edges
  rbind(edges_exposure_protein_2,
        edges_protein_protein_2,
        edges_protein_disease_1,
        edges_gene_disease_2,
        edges_exposure_gene_2)  %>% 
    dplyr::mutate(edge_type = case_when(database == "uniprotKB" ~ 3,
                                 database == "stringdb" ~ 2,
                                 T ~ 1)
    )  %>% 
    dplyr::distinct(source, target, .keep_all = T)  %>% 
    dplyr::filter(source != target) -> edges_all_1

  #find edges_metabolome to the proteins of edges_all
  if(nrow(edges_all_1) > 0){ #skip this line if debug line-by-line
    DB_Enzyme %>% 
      dplyr::filter(!uniprot %in% c(NA)) %>% 
      dplyr::mutate(Metabolites = str_extract_all(Product,"(\\C\\d{5})")) %>%
      tidyr::unnest(Metabolites) %>% 
      dplyr::filter(Metabolites %in% na.omit(eSet$Expo$Data$MetabolomeID)
      ) %>% 
      dplyr::filter(EXP %in% unique(c(edges_all_1$source, edges_all_1$target)) # screen the target uniprot in PPI
              ) -> edges_metabolome

     if(nrow(edges_metabolome) > 0){
       edges_metabolome %>% 
         dplyr::distinct(EXP, Metabolites,Reaction
         ) %>% 
         left_join(eSet$Expo$Data %>% 
         dplyr::distinct(MetabolomeID, EX) %>% 
         dplyr::filter(!MetabolomeID %in% c(NA)),
         by = c("Metabolites" = "MetabolomeID")
         ) %>% 
         dplyr::rename(source = EX,
                target = EXP,
                interaction = Reaction
                )  %>% 
         dplyr::mutate(
                source.class = "metabolite",
                target.class = "protein",
                database = "kegg",
                edge_type = 4
                )  %>% 
         dplyr::select(source,
                       target,
                       interaction,
                       source.class,
                       target.class,
                       database,
                       edge_type)  %>% 
         dplyr::distinct(source, target,
                  .keep_all = T) -> edges_metabolome
     
         if(length(edges_metabolome) > 0){
           edges_metabolome %>% 
             dplyr::distinct(target) %>% 
             na.omit() %>% 
             .$target -> nodes.metabolome.sig
         }
    }else{
      nodes.metabolome.sig <- NULL
    }

    #find edges_gene to the proteins of edges_all
    DB_Gene_Protein_EX %>% 
    filter(EXG %in% edges_exposure_gene_2$target) %>% 
    filter(EXP %in% unique(c(edges_protein_protein_2$source,
     edges_protein_protein_2$target, edges_protein_disease_1))) %>% 
         dplyr::rename(source = EXG,
                target = EXP
                )  %>% 
         dplyr::mutate(
                source.class = "gene",
                target.class = "protein",
                edge_type = 5
                )  %>% 
         dplyr::select(source,
                       target,
                       interaction,
                       source.class,
                       target.class,
                       database,
                       edge_type)  %>% 
         dplyr::distinct(source, target,
                  .keep_all = T) -> edges_gene

    if(length(edges_gene) > 0){
           edges_gene %>% 
             dplyr::distinct(target) %>% 
             na.omit() %>% 
             .$target -> nodes.genome.sig
    }else{
      nodes.genome.sig <- NULL
    }


  #find edges_proteome to the proteins of edges_all
    edges_all_1 %>% 
      dplyr::filter(edge_type == 2) %>% 
      dplyr::select(source, target) %>% 
      tidyr::pivot_longer(source:target,
                   values_to = "EX") %>% 
      dplyr::distinct(EX) %>% 
      dplyr::filter(EX %in% eSet$Expo$Data$EX) %>% 
      na.omit() %>% 
      .$EX -> nodes.proteome.sig
   
     #add edges_metabolome to edges_all
     rbind(edges_all_1, edges_metabolome, edges_gene) %>% 
      dplyr::mutate(source.label = source,
       target.label = target) -> edges_all_2

     #generate nodes
     edges_all_2 %>% 
       tidyr::pivot_longer(c("source","target"),
                           values_to = "node") %>% 
       dplyr::mutate(group = case_when(name == "source" ~ source.class,
                                name == "target" ~ target.class)
                                ) %>% 
       dplyr::distinct(node) %>% 
       dplyr::left_join(vroom::vroom("database/DB_Alias_new.csv", delim = ",", show_col_types = F),
       by = c("node" = "ID")) %>% 
       dplyr::rename(id = node,
       label = Alias) %>% 
       dplyr::mutate(group = case_when(str_detect(id, "EX:E") ~ "exposure",
                                       str_detect(id, "EX:D")  ~ "disease",
                                       str_detect(id, "EXM:M") ~ "metabolite",
                                       str_detect(id, "EX:G") ~ "gene",
                                       id  %in%  nodes.metabolome.sig ~ "protein_by_metabolome",
                                       id  %in%  nodes.genome.sig ~ "protein_by_genome",
                                       id  %in%  nodes.proteome.sig ~ "protein_by_proteome",
                                       id  %in%  intersect(nodes.metabolome.sig, 
                                           nodes.proteome.sig) ~ "protein_by_protome_metabolome",
                                       id  %in%  intersect(nodes.genome.sig, 
                                           nodes.proteome.sig) ~ "protein_by_protome_genome",
                                       id  %in%  intersect(intersect(nodes.genome.sig, 
                                           nodes.metabolome.sig), nodes.proteome.sig) ~ "protein_by_protome_metabolome_genome",
                                       T ~ "ppi")
       ) -> nodes_all
     # 
     # edges_all_2 %>% 
     #   dplyr::arrange(edge_type) %>% 
     #   fwrite(str_c(Path, "/", Mode, "_edges.csv"))
     # 
     # nodes_all %>% 
     #   dplyr::arrange(group) %>% 
     #   fwrite(str_c(Path, "/", Mode, "_nodes.csv"))
     
     edges_all_2 %>% 
       distinct(source, target, 
                .keep_all = T)  %>% 
       dplyr::filter(!source == target) -> eSet$Biolink$Edges_PPI
     nodes_all -> eSet$Biolink$Nodes_PPI
     
     
     #save data for R package --------------------------------------------------------------------------------
     edges_all_2 -> eSet$Edges[[Mode]]
     nodes_all -> eSet$Nodes[[Mode]] 

nodes_all %>% 
  filter(group == "gene")

     
    }else{
      message("No matched PPI pathway! ")
    }
  },

  "GO" = {
  #read database
  vroom::vroom("data/DB_Exposure_GO_EX.csv",
  show_col_types = F) -> DB_Exposure_GO_EX

  vroom::vroom("data/DB_GO_Disease_EX.csv",
  show_col_types = F) -> DB_GO_Disease_EX

  vroom::vroom("data/DB_GO.csv",
  show_col_types = F) -> DB_GO

  #filter EXE of DB_Exposure_GO_EX by EX (ExposureID)
  ddpcr::quiet(
    DB_Exposure_GO_EX %>% 
      dplyr::rename(database = source) %>% 
      dplyr::filter(EXE  %in%  na.omit(eSet$Expo$Data$EX)) %>%
      dplyr::rename(source = EXE,
           target = go.id) %>% 
      dplyr::mutate(source.class = "chemical",
           target.class = "GO",
           interaction = "association",
             edge_type = 1) %>% 
      dplyr::select(source,
                    target,
                    interaction,
                    source.class,
                    target.class,
                    database,
                    edge_type) %>% 
      dplyr::distinct(source, target, .keep_all = T) -> edges_chemical.go_ctd
  )
  
  #filter DB_GO_Disease_EX by DiseaseID
  ddpcr::quiet(
  DB_GO_Disease_EX %>% 
    dplyr::rename(database = source) %>% 
    dplyr::filter(EXD  %in%  na.omit(eSet$Expo$Data$EX)) %>% 
    dplyr::mutate(interaction = "association",
           edge_type = 2) %>% 
    dplyr::rename(source = go.id,
                  target = EXD) %>% 
    dplyr::mutate(source.class = "GO",
                  target.class = "disease",
                  database = "ctd") %>% 
    dplyr::select(source,
                  target,
                  interaction,
                  source.class,
                  target.class,
                  database,
                  edge_type) %>% 
    dplyr::distinct(source, 
             target, 
             .keep_all = T) -> edges_go.disease_ctd
  )

  #cbind edges
  intersect(edges_chemical.go_ctd$target,
   edges_go.disease_ctd$source) -> go.intersect 
  
  rbind(edges_chemical.go_ctd %>% 
          dplyr::filter(target %in% go.intersect),
        edges_go.disease_ctd %>% 
          dplyr::filter(source %in% go.intersect)
          ) -> edges_all_1
  
  edges_all_2 <- NULL
  
  #enrich proteome using stringdb.enrichment (enrich->go)
  if(length(eSet$Expo$Data$ProteomeID) > 0){
    ddpcr::quiet(
      STRINGdb::STRINGdb$new(version = "11.0b",
                           species = 9606, #9606 for Human, 10090 for mouse, NCBI taxonomy identiers
                           score_threshold = 200,
                           input_directory= "")$get_enrichment -> eSet$Biolink$Enrich
      )

    eSet$Biolink$Enrich(eSet$Expo$Data$ProteomeID)$term %>%
      as_tibble() %>%
      dplyr::filter(str_detect(value, "GO:.*")) %>%
      .$value -> nodes.proteome.sig

    if(length(nodes.proteome.sig) > 0){
      tibble(nodes.proteome.sig) %>%
        purrr::set_names("id") %>%
        dplyr::left_join(DB_GO %>% 
                           dplyr::filter(!str_detect(Synonym, "GO")) %>% 
                           dplyr::distinct(go_id, .keep_all = T) %>% 
                           dplyr::select(go_id, Synonym),
                         by = c("id" = "go_id")) %>% 
        dplyr::rename(label = Synonym) %>% 
        dplyr::mutate(group = "proteome") -> nodes_proteome
      
      
      if(length(nodes_proteome) > 0){
      edges_all_1 %>% 
        dplyr::filter(source %in% c(nodes_proteome$id) |
                        target %in% c(nodes_proteome$id)) -> edges_all_proteome
      
      edges_all_1 %>% 
        dplyr::filter(source %in% c(nodes_proteome$id)) %>% 
        distinct(source) %>% 
        .$source -> nodes_proteome_included
      }
    }
  }else{
      nodes.proteome.sig <- NULL
    }

  
  #enrich metabolome using metaboanalyst (pathway->go.name->go.id)
  if(length(eSet$Expo$Data$MetabolomeID %>% 
  as_tibble() %>% 
  na.omit() %>% 
  .$value) > 0){
    switch (MetBiospec,
            "blood" = {vroom::vroom("data/DB_Enrich_Blood_Metabolite_50sets.csv", 
                                    show_col_types = F)  %>% 
                distinct(pathway, kegg.id,
                         .keep_all = T) -> df_background
                         },
            "urine" = {vroom::vroom("data/DB_Enrich_Urine_Metabolite_45sets.csv", 
                                    show_col_types = F) %>% 
                distinct(pathway, kegg.id,
                         .keep_all = T) -> df_background
                         }
            )
    #metabolites enrichment analysis
    # pp = "hsa00030"
    FuncEnrichMet <- function(pp){
      intersect(eSet$Expo$Data$MetabolomeID, df_background$kegg.id) -> mm.incl
      
      df_background %>% 
        dplyr::filter(pathway %in% pp) -> Path.target
      
      df_background %>% 
        dplyr::filter(!pathway %in% pp) -> Path.untarget
      
      aa = Path.target %>% 
        dplyr::filter(kegg.id %in% mm.incl) %>% 
        nrow()
      bb = nrow(Path.target) - aa
      cc = Path.untarget %>% 
        dplyr::filter(kegg.id %in% mm.incl) %>% 
        nrow()
      dd = nrow(Path.untarget) - cc
      
      if(aa/nrow(Path.target) > cc/nrow(Path.untarget)){
        tibble(sig = c(aa, cc),
               unsig = c(bb, dd)) %>% 
          fisher.test(alternative = "two.sided") %>% 
          .$p.value -> p.value
      }else{
        p.value = 1
      }
      
      Path.target %>% 
        dplyr::mutate(P.value = p.value) %>% 
        dplyr::distinct(go_id, pathway,
                        .keep_all = T) -> temp
        
      return(temp)
    }
    
    furrr::future_map_dfr(unique(df_background$pathway), 
                          FuncEnrichMet) -> res.enrich
    
    res.enrich %>% 
     dplyr::filter(P.value < 0.05) %>% 
      dplyr::select(go_id, Synonym) %>%
      purrr::set_names("id", "label") %>%
      dplyr::mutate(group = "metabolome") %>% 
      dplyr::select(id, group, label) -> nodes_metabolome
    
    if(length(nodes_metabolome) > 0){
      edges_all_1 %>% 
        dplyr::filter(source %in% c(nodes_metabolome$id) |
                        target %in% c(nodes_metabolome$id)) -> edges_all_metabolome
      
      edges_all_1 %>% 
        dplyr::filter(source %in% c(nodes_metabolome$id)) %>% 
        distinct(source) %>% 
        .$source -> nodes_metabolome_included
      }
  }else{
    res.enrich <- NULL
  }

  rbind(edges_all_proteome,
   edges_all_metabolome) -> edges_all_2

  if(length(edges_all_2) > 0){  
  #generate nodes
        edges_all_2 %>% 
           tidyr::pivot_longer(c("source","target"),
                        values_to = "node") %>% 
                        dplyr::left_join(vroom::vroom("database/DB_Alias_new.csv", delim = ",",  show_col_types = FALSE),
                        by = c("node" = "ID")) %>% 
           dplyr::rename(label = Alias) %>% 
           dplyr::select(node,
                         label,
                         database) %>% 
           dplyr::distinct(node, label, database) %>% 
           dplyr::rename(id = node,
                  group = database) %>% 
           dplyr::mutate(group = case_when(str_detect(id, "EX:E") ~ "exposure",
                                           str_detect(id, "EX:D") ~ "disease",
                                           (id  %in%  nodes_proteome_included) | !(id  %in%  nodes_metabolome_included) ~ "go_by_proteome",
                                           !(id  %in%  nodes_proteome_included) | (id  %in%  nodes_metabolome_included) ~ "go_by_metabolome",
                                           (id  %in%  nodes_proteome_included) | (id  %in%  nodes_metabolome_included) ~ "go_by_metabolome_proteome",
                                           T ~ "go")
           ) -> nodes_all
         
        edges_all_2 %>% 
          distinct(source, target, 
                   .keep_all = T)  %>% 
          dplyr::filter(!source == target) -> eSet$Biolink$Edges_GO

        nodes_all -> eSet$Biolink$Nodes_GO
         
         #save data --------------------------------------------------------------------------------
         # edges_all_2 %>% 
         #  dplyr::arrange(edge_type) %>% 
         #  data.table::fwrite(str_c(Path, "/", Mode, "_edges.csv"))
         # 
         # #save data --------------------------------------------------------------------------------
         # nodes_all %>% 
         #  dplyr::arrange(group) %>% 
         #  data.table::fwrite(str_c(Path, "/", Mode, "_nodes.csv"))
         # 
         #save data for R package --------------------------------------------------------------------------------
         edges_all_2 -> eSet$Edges[[Mode]]
         nodes_all -> eSet$Nodes[[Mode]] 
        
        }else{
           message("No matched GO pathway! ")
         }
      }
  )

  #save R command log
  eSet$RCommandLog <- eSet$AddCommand(stringr::str_c("\n eSet <- Biolink(eSet = eSet, \n",
                                                     " Mode = '", Mode, "', \n",
                                                     " ChemCas = '", ChemCas, "', \n",
                                                     " ChemInchikey = '", ChemInchikey, "', \n",
                                                     " DiseaseID = '", DiseaseID, "', \n",
                                                     " MetabolomeID = '", MetabolomeID, "', \n",
                                                     " MetBiospec = '", MetBiospec, "', \n",
                                                     " ProteomeID = '", ProteomeID, "') \n",
                                                     " GenomeID = '", GenomeID, "') \n")
                                      )
  
  eSet$RCommandLog %>% 
    as_tibble() %>% 
    purrr::set_names("R commands") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/rcommands log.txt"))
  
  #print message and save log
  message("Complete Biolink modeling in ", Mode," mode. ", NowTime, "\n")
  eSet$EXEecutionLog  <- eSet$AddLog(stringr::str_c("Complete Biolink modeling in ", Mode," mode. ", NowTime))
  
  eSet$EXEecutionLog %>% 
    as_tibble() %>% 
    purrr::set_names("running log") %>% 
    data.table::fwrite(stringr::str_c(eSet$FileDirOut,"/running log.txt"))  
  
  eSet %>% 
    save(file = str_c(eSet$FileDirOut,"/eSet.Rdata"))
  
  tictoc::toc()
  
  return(eSet)
}


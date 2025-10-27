#------------------------------------------------------------------------------#
#-------------------Functions to compare two datababases ----------------------#
#------------------------------------------------------------------------------#
#
# AUTHORS:
# Paula Pappalardo - Smithsonian Environmental Research Center
# email: paulapappalardo@gmail.com
#
# Please see README file for background and citation instructions. 

addScinameLevel <- function(mydf){
  # Finds the taxonomic resolution for a scientific name
  # Args:
  #  mydf: dataframe with taxonomic columns kingdom, phylum, class, order, family, genus, and species 
  #  TODO: consider activating subfamily if you are dealing with BOLD output
  # Ouput:
  #   same dataframe with an additional column named sciname_level
  mydf_ed <- mydf %>% 
    mutate(sciname_level = case_when(sciname %in% na.omit(species) ~ "species",
                                     sciname %in% na.omit(genus) ~ "genus",
                                     sciname %in% na.omit(family) ~ "family",
                                     #sciname %in% na.omit(subfamily) ~ "subfamily", # if using
                                     sciname %in% na.omit(order) ~ "order",
                                     sciname %in% na.omit(class) ~ "class",
                                     sciname %in% na.omit(phylum) ~ "phylum",
                                     sciname %in% na.omit(kingdom) ~ "kingdom")) %>% 
    relocate(sciname_level, .after = sciname)
  # return modified dataframe
  return(mydf_ed)
}

labelFinalMatch <- function(mydf, label){
  # condenses taxonomy to one tax string and adds reference db label
  # Args:
  #   mydf: blast final dataset (after quality filters, tags, and identical matches were solved)
  #   label: name to represent the origin database (e.g., "mlml" or "midori")
  # Output: the function will remove higher taxonomy columns and add a taxonomy string column
  # Usage: for dplyr use - labelFinalMatch(., label = "mlml")
  # 
  # combine all taxonomic levels into one string
  mydf_ed <- mydf %>% 
    mutate(tax_string = paste(kingdom, phylum, class, order, family, genus, species, sep = ";")) %>% 
    select(-kingdom, -phylum, -class, -order, -family, -genus, -species)
  # rename columns to indicate database source
  newnames <- paste(names(mydf_ed), label, sep = "_")
  names(mydf_ed) <- newnames
  return(mydf_ed)
}

compareBOLDandMIDORI <- function(mydf_bold, mydf_midori, identity_th){
  # Compare matches to bold and midori for same sequence
  # Args
  # mydf_bold: dataframe with blast results from BOLD
  # mydf_midori: dataframe with blast results from MIDORI
  # identity_th: which identity threshold you want to assume "good matches"
  # both dataframes need columns for similarity metrics and taxonomic levels 
  # prepare data to be compared
  # midori
  mydf_midori_ed <- mydf_midori %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "midori") %>% 
    dplyr::rename(query_seqid = query_seqid_midori)
  # bold
  mydf_bold_ed <- mydf_bold %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "bold") %>% 
    dplyr::rename(query_seqid = query_seqid_bold)
  
  # merge and compare matches
  mydf_ed <- mydf_midori_ed %>% 
    full_join(mydf_bold_ed, by = "query_seqid") %>%
    # first check for agreement on sciname and sciname level
    mutate(sciname_agree = ifelse(sciname_bold == sciname_midori, TRUE, FALSE),
           sciname_level_agree = ifelse(sciname_level_bold == sciname_level_midori, TRUE, FALSE)) %>% 
    # then identify the best identity for later
    mutate(best_db_identity  = case_when(is.na(percent_identity_bold) & !is.na(percent_identity_midori) ~ "keep midori - bold NA",
                                         is.na(percent_identity_midori) & !is.na(percent_identity_bold) ~ "keep bold - midori NA",
                                         percent_identity_midori > percent_identity_bold ~ "keep midori",
                                         percent_identity_midori <= percent_identity_bold ~ "keep bold",
                                         percent_identity_midori == percent_identity_bold ~ "identical",
                                         T ~ "CHECK")) %>% 
    # now try to come up with the outcomes based on the info we have:
    mutate(match_outcome = case_when(sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_midori >= identity_th & percent_identity_bold >= identity_th ~ "tricky and relevant - separate function",
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_midori <= identity_th & percent_identity_bold <= identity_th ~ "tricky and irrelevant - separate function",
                                     percent_identity_bold > identity_th & percent_identity_midori <= identity_th ~ "keep bold - bold higher",
                                     percent_identity_midori > identity_th & percent_identity_bold <= identity_th ~ "keep midori - midori higher",
                                     is.na(percent_identity_bold) & !is.na(percent_identity_midori) ~ "keep midori - bold NA",
                                     is.na(percent_identity_midori) & !is.na(percent_identity_bold)  ~ "keep bold - midori NA",
                                     is.na(tax_string_midori) & is.na(tax_string_bold) ~ "no match from either",
                                     sciname_agree == TRUE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == FALSE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == TRUE & sciname_level_agree == FALSE ~ best_db_identity,# few cases, e.g. Oomycota
                                     T ~ "CHECK"))
  # create ordered factor of the taxonomic level
  higher_tax <- factor(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                       levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), ordered = T)
  
  # separate "Tricky" cases where we want to compare taxonomic resolution achieved by each database
  mydf_ed_1 <- mydf_ed %>% 
    filter(grepl("tricky", match_outcome)) %>% 
    rowwise() %>% 
    mutate(match_outcome = ifelse(which(higher_tax == sciname_level_bold) > which(higher_tax == sciname_level_midori), "keep bold", "keep midori")) %>% 
    ungroup()
  
  mydf_ed_2 <- mydf_ed %>% 
    filter(!grepl("tricky", match_outcome))
  
  mydf_ed_final <- bind_rows(mydf_ed_1, mydf_ed_2)
  
  return(mydf_ed_final)
  
}

pickFinalTax_BoldVsMidori <- function(mydf){
  # Pick best final id and keep final percent identity
  # Args
  #   mydf: classified dataframe after running compareDatabases(), compareTaxMatches(), and classifyAssignments()
  # Output
  #   The function will keep only the higher taxonomy of the best hit, and will add column database_source to track where the match came from
  # Usage: pickFinalTax(.) 
  #
  mydf_ed <- mydf %>% 
    dplyr::mutate(tax_string = case_when(grepl("keep midori", match_outcome) ~ tax_string_midori,
                                         grepl("keep bold", match_outcome) ~ tax_string_bold,
                                         T ~ "CHECK"),
                  sciname = case_when(grepl("keep midori", match_outcome) ~ sciname_midori,
                                      grepl("keep bold", match_outcome) ~ sciname_bold,
                                      T ~ "CHECK"),
                  sciname_level = case_when(grepl("keep midori", match_outcome) ~ sciname_level_midori,
                                            grepl("keep bold", match_outcome) ~ sciname_level_bold,
                                            T ~ "CHECK"),
                  taxonomy_source = case_when(grepl("keep midori", match_outcome) ~ taxonomy_source_midori,
                                              grepl("keep bold", match_outcome) ~ taxonomy_source_bold,
                                              T ~ "CHECK"),
                  database_source = case_when(grepl("keep midori", match_outcome) ~ "MIDORI2",
                                              grepl("keep bold", match_outcome)~ "BOLD",
                                              T ~ "CHECK"),
                  match_name = case_when(grepl("keep midori", match_outcome) ~ result_seqid_midori,
                                         grepl("keep bold", match_outcome)~ result_seqid_bold,
                                         T ~ "CHECK"),
                  bitscore = case_when(grepl("keep midori", match_outcome) ~ bitscore_midori,
                                       grepl("keep bold", match_outcome) ~ bitscore_bold,
                                       T ~ NA),
                  percent_coverage = case_when(grepl("keep midori", match_outcome) ~ percent_coverage_midori,
                                               grepl("keep bold", match_outcome) ~ percent_coverage_bold,
                                               T ~ NA),
                  percent_identity = case_when(grepl("keep midori", match_outcome) ~ percent_identity_midori,
                                               grepl("keep bold", match_outcome) ~ percent_identity_bold,
                                               T ~ NA)) %>%
    dplyr::select(query_seqid, database_source, match_name, bitscore, percent_coverage, percent_identity, sciname, sciname_level, tax_string, taxonomy_source) 
  # return modified dataframe that includes only the final chosen hit from each database
  return(mydf_ed)
}

compareMLMLandLavrador<- function(mydf_mlml, mydf_lavrador, identity_th){
  # Compare matches between MLML and Lavrador for each unknown sequence
  # Args
  # mydf_mlml: dataframe with blast results from BOLD
  # mydf_labrador: dataframe with blast results from MIDORI
  # identity_th: which identity threshold you want to assume "good matches"
  # both dataframes need columns for similarity metrics and taxonomic levels 
  # prepare data to be compared
  # midori
  mydf_mlml_ed <- mydf_mlml %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "mlml") %>% 
    dplyr::rename(query_seqid = query_seqid_mlml)
  # bold
  mydf_lavrador_ed <- mydf_lavrador %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "lavrador") %>% 
    dplyr::rename(query_seqid = query_seqid_lavrador)
  
  # merge and compare matches
  mydf_ed <- mydf_mlml_ed %>% 
    full_join(mydf_lavrador_ed, by = "query_seqid") %>%
    # first check for agreement on sciname and sciname level
    mutate(sciname_agree = ifelse(sciname_mlml == sciname_lavrador, TRUE, FALSE),
           sciname_level_agree = ifelse(sciname_level_mlml == sciname_level_lavrador, TRUE, FALSE)) %>% 
    # then identify the best identity for later
    mutate(best_db_identity  = case_when(is.na(percent_identity_mlml) & !is.na(percent_identity_lavrador) ~ "keep lavrador - mlml NA",
                                         is.na(percent_identity_lavrador) & !is.na(percent_identity_mlml) ~ "keep mlml - lavrador NA",
                                         percent_identity_lavrador > percent_identity_mlml ~ "keep lavrador",
                                         percent_identity_lavrador <= percent_identity_mlml ~ "keep mlml",
                                         percent_identity_lavrador == percent_identity_mlml ~ "identical",
                                         T ~ "CHECK")) %>% 
    # now try to come up with the outcomes based on the info we have:
    mutate(match_outcome = case_when(sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_lavrador >= identity_th & percent_identity_mlml >= identity_th ~ "tricky and relevant - separate function",
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_lavrador <= identity_th & percent_identity_mlml <= identity_th ~ "tricky and irrelevant - separate function",
                                     percent_identity_mlml > identity_th & percent_identity_lavrador <= identity_th ~ "keep mlml - mlml higher",
                                     percent_identity_lavrador > identity_th & percent_identity_mlml <= identity_th ~ "keep lavrador - lavrador higher",
                                     is.na(percent_identity_mlml) & !is.na(percent_identity_lavrador) ~ "keep lavrador - mlml NA",
                                     is.na(percent_identity_lavrador) & !is.na(percent_identity_mlml)  ~ "keep mlml - lavrador NA",
                                     is.na(tax_string_lavrador) & is.na(tax_string_mlml) ~ "no match from either",
                                     sciname_agree == TRUE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == FALSE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == TRUE & sciname_level_agree == FALSE ~ best_db_identity,# few cases, e.g. Oomycota
                                     T ~ "CHECK"))
  # create ordered factor of the taxonomic level
  higher_tax <- factor(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                       levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), ordered = T)
  
  # separate "Tricky" cases where we want to compare taxonomic resolution achieved by each database
  mydf_ed_1 <- mydf_ed %>% 
    filter(grepl("tricky", match_outcome)) %>% 
    rowwise() %>% 
    mutate(match_outcome = ifelse(which(higher_tax == sciname_level_mlml) > which(higher_tax == sciname_level_lavrador), "keep mlml", "keep lavrador")) %>% 
    ungroup()
  
  mydf_ed_2 <- mydf_ed %>% 
    filter(!grepl("tricky", match_outcome))
  
  mydf_ed_final <- bind_rows(mydf_ed_1, mydf_ed_2)
  
  return(mydf_ed_final)
  
}

pickFinalTax_MLMLvsLavrador <- function(mydf){
  # Pick best final id and keep final percent identity
  # Args
  #   mydf: classified dataframe after running compareDatabases(), compareTaxMatches(), and classifyAssignments()
  # Output
  #   The function will keep only the higher taxonomy of the best hit, and will add column database_source to track where the match came from
  # Usage: pickFinalTax(.) 
  #
  mydf_ed <- mydf %>% 
    mutate(tax_string = case_when(grepl("keep mlml", match_outcome) ~ tax_string_mlml,
                                  grepl("keep lavrador", match_outcome) ~ tax_string_lavrador,
                                  T ~ "CHECK"),
           sciname = case_when(grepl("keep mlml", match_outcome) ~ sciname_mlml,
                               grepl("keep lavrador", match_outcome) ~ sciname_lavrador,
                               T ~ "CHECK"),
           sciname_level = case_when(grepl("keep mlml", match_outcome) ~ sciname_level_mlml,
                                     grepl("keep lavrador", match_outcome) ~ sciname_level_lavrador,
                                     T ~ "CHECK"),
           taxonomy_source = case_when(grepl("keep mlml", match_outcome) ~ taxonomy_source_mlml,
                                       grepl("keep lavrador", match_outcome) ~ taxonomy_source_lavrador,
                                       T ~ "CHECK"),
           database_source = case_when(grepl("keep mlml", match_outcome) ~ "mlml",
                                       grepl("keep lavrador", match_outcome)~ "lavrador",
                                       T ~ "CHECK"),
           match_name = case_when(grepl("keep mlml", match_outcome) ~ result_seqid_mlml,
                                  grepl("keep lavrador", match_outcome)~ result_seqid_lavrador,
                                  T ~ "CHECK"),
           bitscore = case_when(grepl("keep mlml", match_outcome) ~ bitscore_mlml,
                                grepl("keep lavrador", match_outcome) ~ bitscore_lavrador,
                                T ~ NA),
           percent_coverage = case_when(grepl("keep mlml", match_outcome) ~ percent_coverage_mlml,
                                        grepl("keep lavrador", match_outcome) ~ percent_coverage_lavrador,
                                        T ~ NA),
           percent_identity = case_when(grepl("keep mlml", match_outcome) ~ percent_identity_mlml,
                                        grepl("keep lavrador", match_outcome) ~ percent_identity_lavrador,
                                        T ~ NA)) %>%
    select(query_seqid, database_source, match_name, bitscore, percent_coverage, percent_identity, sciname, sciname_level, tax_string, taxonomy_source) 
  # return modified dataframe that includes only the final chosen hit from each database
  return(mydf_ed)
}

compareGlobalVsCurated <- function(mydf_global, mydf_curated, identity_th){
  # Compare matches from global (midori, bold) to curated (mlml, lavrador)
  # Args
  # mydf_global: dataframe from comparison of global databases (e.g., bold vs midori)
  # mydf_curated: dataframe from comparison of curated databases (e.g., mlml vs Lavrador)
  # identity_th: which identity threshold you want to assume "good matches"
  #
  # Join dataframes and compare matches
  #
  mydf_ed <- mydf_global %>%
    full_join(mydf_curated, by = "query_seqid") %>% 
    # first check for agreement on sciname and sciname level
    mutate(sciname_agree = ifelse(sciname.x == sciname.y, TRUE, FALSE),
           sciname_level_agree = ifelse(sciname_level.x == sciname_level.y, TRUE, FALSE)) %>% 
    # then identify the best identity for later
    mutate(best_db_identity  = case_when(!is.na(percent_identity.x) & is.na(percent_identity.y) ~ "keep x - y NA",
                                         is.na(percent_identity.x) & !is.na(percent_identity.y) ~ "keep y - x NA",
                                         percent_identity.x > percent_identity.y ~ "keep x",
                                         percent_identity.x < percent_identity.y ~ "keep y",
                                         percent_identity.x == percent_identity.y ~ "keep y - identical",
                                         T ~ "CHECK")) %>% 
    # now try to come up with the outcomes based on the info we have:
    mutate(match_outcome = case_when(sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       (percent_identity.x >= identity_th) & (percent_identity.y >= identity_th) ~ "tricky and relevant - separate function",
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       (percent_identity.x <= identity_th) & (percent_identity.y <= identity_th) ~ best_db_identity,
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity.x > identity_th & percent_identity.y <= identity_th ~ "keep x - x higher",
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity.y > identity_th & percent_identity.x <= identity_th ~ "keep y - y higher",
                                     !is.na(percent_identity.x) & is.na(percent_identity.y) ~ "keep x - y NA",
                                     is.na(percent_identity.x) & !is.na(percent_identity.y) ~ "keep y - x NA",
                                     is.na(tax_string.x) & is.na(tax_string.y) ~ "no match from either",
                                     sciname_agree == TRUE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == FALSE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == TRUE & sciname_level_agree == FALSE ~ best_db_identity,# few cases, e.g. Oomycota
                                     T ~ "CHECK"))
  # create ordered factor of the taxonomic level
  higher_tax <- factor(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                       levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), ordered = T)
  
  # separate "Tricky" cases where we want to compare taxonomic resolution achieved by each database
  mydf_ed_1 <- mydf_ed %>% 
    filter(grepl("tricky", match_outcome)) %>% 
    rowwise() %>% 
    mutate(match_outcome = ifelse(which(higher_tax == sciname_level.x) > which(higher_tax == sciname_level.y), "keep x", "keep y")) %>% 
    ungroup()
  
  mydf_ed_2 <- mydf_ed %>% 
    filter(!grepl("tricky", match_outcome))
  
  mydf_ed_final <- bind_rows(mydf_ed_1, mydf_ed_2)
  
  return(mydf_ed_final)
  
}

pickFinalTax_GlobalVsCurated <- function(mydf){
  # Pick best final id and keep final percent identity
  # Args
  #   mydf: final dataframe after running compareGlobalVsCurated()
  # Output
  #   The function will keep only the higher taxonomy of the best hit, and
  # it also adds column database_source to track where the match came from
  # Usage: pickFinalTax(.) 
  #
  mydf_ed <- mydf %>% 
    mutate(tax_string = case_when(grepl("keep x", match_outcome, fixed = T) ~ tax_string.x,
                                  grepl("keep y", match_outcome, fixed = T) ~ tax_string.y,
                                  T ~ "CHECK"),
           sciname = case_when(grepl("keep x", match_outcome) ~ sciname.x,
                               grepl("keep y", match_outcome) ~ sciname.y,
                               T ~ "CHECK"),
           sciname_level = case_when(grepl("keep x", match_outcome) ~ sciname_level.x,
                                     grepl("keep y", match_outcome) ~ sciname_level.y,
                                     T ~ "CHECK"),
           taxonomy_source = case_when(grepl("keep x", match_outcome) ~ taxonomy_source.x,
                                       grepl("keep y", match_outcome) ~ taxonomy_source.y,
                                       T ~ "CHECK"),
           database_source = case_when(grepl("keep x", match_outcome) ~ database_source.x,
                                       grepl("keep y", match_outcome) ~ database_source.y,
                                       T ~ "CHECK"),
           match_name = case_when(grepl("keep x", match_outcome) ~ match_name.x,
                                  grepl("keep y", match_outcome) ~ match_name.y,
                                  T ~ "CHECK"),
           bitscore = case_when(grepl("keep x", match_outcome) ~ bitscore.x,
                                grepl("keep y", match_outcome) ~ bitscore.y,
                                T ~ NA),
           percent_coverage = case_when(grepl("keep x", match_outcome) ~ percent_coverage.x,
                                        grepl("keep y", match_outcome) ~ percent_coverage.y,
                                        T ~ NA),
           percent_identity = case_when(grepl("keep x", match_outcome) ~ percent_identity.x,
                                        grepl("keep y", match_outcome) ~ percent_identity.y,
                                        T ~ NA)) %>%
    select(query_seqid, database_source, match_name, bitscore, percent_coverage, percent_identity, sciname, sciname_level, tax_string, taxonomy_source) 
  # return modified dataframe that includes only the final chosen hit from each database
  return(mydf_ed)
}
#------------------------------------------------------------------------------#
#-------------------Functions to compare three databases-----------------------#
#------------------------------------------------------------------------------#
#
# AUTHORS:
# Emma Palmer - Smithsonian Environmental Research Center
# email: palmerem@si.edu
#
# Paula Pappalardo - Smithsonian Environmental Research Center
# email: paulapappalardo@gmail.com
#
# CITATION: Please see README file for background and citation instructions.
#
# Notes: requires addScinameLevel() and labelFinalMatch() from compare-two-databases.R
#------------------------------------------------------------------------------#

addScinameLevel <- function(mydf){
  # Finds the taxonomic resolution for a scientific name
  # Args:
  #  mydf: dataframe with taxonomic columns kingdom, phylum, class, order, family, genus, and species 
  #  TODO: consider activating subfamily if you are dealing with BOLD output
  # Ouput:
  #   same dataframe with an additional column named sciname_level
  if(!"sciname" %in% names(mydf)) {
    mydf <- mydf %>% 
      rowwise() %>% 
      mutate(sciname = dplyr::last(na.omit(c(kingdom, phylum, class, order, family, genus, species)))) %>% 
      ungroup()
  } 
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

compare3Databases <- function(mydf_local, mydf_global1, mydf_global2, identity_th){
  # Compares similarity measures between local and global databases
  # Args
  #   mydf_local: merged dataframe with blast results from the local database (most trusted database)
  #   mydf_global1: merged dataframe with blast results from the 1st global database (second most trusted database)
  #   mydf_global2: merged dataframe with blast results from the 2nd global database (least trusted database)
  #   identity_th: percent identity threshold we define as a "good" match - can vary with project & PI & genetic marker
  # Output
  #   the function will add columns useful for later comparison of matches to the local and global databases
  #   function can also determine the databases with the best coverage and/or length (later functions do not use these, therefore they are commented out)
  #
  mydf_local_ed <- mydf_local %>% 
    addScinameLevel() %>% 
    unite("tax_string", c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = TRUE) %>% 
    rename_all(~paste(.x, "local", sep = "_")) %>% 
    dplyr::rename(query_seqid = query_seqid_local)
  
  mydf_global1_ed <- mydf_global1 %>% 
    addScinameLevel() %>% 
    unite("tax_string", c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = TRUE) %>% 
    rename_all(~paste(.x, "global1", sep = "_")) %>% 
    dplyr::rename(query_seqid = query_seqid_global1)
  
  mydf_global2_ed <- mydf_global2 %>% 
    addScinameLevel() %>%
    unite("tax_string", c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = TRUE) %>% 
    rename_all(~paste(.x, "global2", sep = "_")) %>% 
    dplyr::rename(query_seqid = query_seqid_global2)
  
  mydf_ed <- mydf_local_ed %>%
    full_join(mydf_global1_ed, by = "query_seqid") %>%
    full_join(mydf_global2_ed, by = "query_seqid") %>%
    mutate(sciname_agree = case_when(sciname_local == sciname_global1 & sciname_global1 == sciname_global2 ~ "all_agree", 
                                     sciname_local == sciname_global1 & sciname_global1 != sciname_global2 ~ "local&global1 agree",
                                     sciname_local != sciname_global1 & sciname_local == sciname_global2 ~ "local&global2 agree",
                                     sciname_local != sciname_global1 & sciname_global1 == sciname_global2 ~ "global1&global2 agree",
                                     sciname_local == sciname_global1 & is.na(sciname_global2) ~ "local&global1 agree",
                                     sciname_local == sciname_global2 & is.na(sciname_global1) ~ "local&global2 agree",
                                     sciname_global1 == sciname_global2 & is.na(sciname_local) ~ "global1&global2 agree",
                                     T ~ "no_agreement"),
           sciname_level_agree = case_when(sciname_level_local == sciname_level_global1 & sciname_level_global1 == sciname_level_global2 ~ "all_agree", 
                                           sciname_level_local == sciname_level_global1 & sciname_level_global1 != sciname_level_global2 ~ "local&global1 agree",
                                           sciname_level_local != sciname_level_global1 & sciname_level_local == sciname_level_global2 ~ "local&global2 agree",
                                           sciname_level_local != sciname_level_global1 & sciname_level_global1 == sciname_level_global2 ~ "global1&global2 agree",
                                           sciname_level_local == sciname_level_global1 & is.na(sciname_level_global2) ~ "local&global1 agree",
                                           sciname_level_local == sciname_level_global2 & is.na(sciname_level_global1) ~ "local&global2 agree",
                                           sciname_level_global1 == sciname_level_global2 & is.na(sciname_level_local) ~ "global1&global2 agree",
                                           T ~ "no_agreement")) %>%
    mutate(best_db_identity  = case_when(!is.na(percent_identity_local) & is.na(percent_identity_global1) & is.na(percent_identity_global2) ~ "local only",
                                         is.na(percent_identity_local) & !is.na(percent_identity_global1) & is.na(percent_identity_global2) ~ "global1 only",
                                         is.na(percent_identity_local) & is.na(percent_identity_global1) & !is.na(percent_identity_global2) ~ "global2 only",
                                         percent_identity_local >= percent_identity_global1 & percent_identity_local >= percent_identity_global2 ~ "local",
                                         percent_identity_local < percent_identity_global1 & percent_identity_global1 >= percent_identity_global2 ~ "global1",
                                         percent_identity_local < percent_identity_global2 & percent_identity_global1 < percent_identity_global2 ~ "global2",
                                         percent_identity_local >= percent_identity_global1 & is.na(percent_identity_global2) ~ "local",
                                         percent_identity_local >= percent_identity_global2 & is.na(percent_identity_global1) ~ "local",
                                         percent_identity_global1 > percent_identity_local & is.na(percent_identity_global2) ~ "global1",
                                         percent_identity_global1 >= percent_identity_global2 & is.na(percent_identity_local) ~ "global1",
                                         percent_identity_global2 > percent_identity_local & is.na(percent_identity_global1) ~ "global2",
                                         percent_identity_global2 > percent_identity_global1 & is.na(percent_identity_local) ~ "global2",
                                         T ~ "CHECK"),
           match_outcome = case_when(sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local >= identity_th & 
                                       percent_identity_global1 >= identity_th & percent_identity_global2 >= identity_th ~ "tricky and relevant - all - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local >= identity_th & 
                                       percent_identity_global1 >= identity_th & percent_identity_global2 < identity_th ~ "tricky and relevant - local&global1 - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local >= identity_th & 
                                       percent_identity_global1 < identity_th & percent_identity_global2 >= identity_th ~ "tricky and relevant - local&global2 - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local < identity_th & 
                                       percent_identity_global1 >= identity_th & percent_identity_global2 >= identity_th ~ "tricky and relevant - global1&global2 - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local >= identity_th & 
                                       percent_identity_global1 >= identity_th & is.na(percent_identity_global2) ~ "tricky and relevant - local&global1 - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local >= identity_th & 
                                       is.na(percent_identity_global1) & percent_identity_global2 >= identity_th ~ "tricky and relevant - local&global2 - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & is.na(percent_identity_local) & 
                                       percent_identity_global1 >= identity_th & percent_identity_global2 >= identity_th ~ "tricky and relevant - global1&global2 - separate function",
                                     sciname_agree == "no_agreement" & sciname_level_agree == "no_agreement" & percent_identity_local < identity_th & 
                                       percent_identity_global1 < identity_th & percent_identity_global2 < identity_th ~ "tricky and irrelevant - separate function",
                                     is.na(percent_identity_local) & is.na(percent_identity_global1) & is.na(percent_identity_global2) ~ "no match from any",
                                     T ~ best_db_identity)) %>% # includes when sciname and level do not agree or when only one agrees (scinames can agree at name but not at level in a few cases like Oomycota)
    relocate(c("best_db_identity", "match_outcome"), .after = query_seqid) %>%  #"best_db_coverage", "best_db_length", 
    mutate(across(starts_with("sciname_level"), ~ factor(.x, levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), ordered = T)))
  
  # separate "Tricky" cases where we want to compare taxonomic resolution achieved by each database
  mydf_ed_1 <- mydf_ed %>% 
    filter(grepl("tricky", match_outcome)) %>% 
    rowwise() %>% 
    mutate(match_outcome = case_when(sciname_level_local == sciname_level_global1 & sciname_level_local == sciname_level_global2 ~ "keep local - synonym?",
                                     sciname_level_local >= sciname_level_global1 & sciname_level_local >= sciname_level_global2 ~ "keep local",
                                     sciname_level_local < sciname_level_global1 & sciname_level_global1 > sciname_level_global2 ~ "keep global1",
                                     sciname_level_local < sciname_level_global2 & sciname_level_global1 < sciname_level_global2 ~ "keep global2",
                                     sciname_level_local < sciname_level_global1 &  sciname_level_global1 == sciname_level_global2 ~ "keep global1 - synonym?",
                                     sciname_level_local > sciname_level_global1 &  is.na(sciname_level_global2) ~ "keep local",
                                     sciname_level_local < sciname_level_global1 &  is.na(sciname_level_global2) ~ "keep global1",
                                     sciname_level_local > sciname_level_global2 &  is.na(sciname_level_global1)  ~ "keep local",
                                     sciname_level_local < sciname_level_global2 &  is.na(sciname_level_global1)  ~ "keep global2",
                                     sciname_level_global1 > sciname_level_global2 &  is.na(sciname_level_local)  ~ "keep global1",
                                     sciname_level_global1 < sciname_level_global2 &  is.na(sciname_level_local)  ~ "keep global2",
                                     T ~ match_outcome)) %>% 
    ungroup()
  
  mydf_ed_2 <- mydf_ed %>% 
    filter(!grepl("tricky", match_outcome))
  
  mydf_ed_final <- bind_rows(mydf_ed_1, mydf_ed_2) %>% 
    mutate(match_outcome = ifelse(str_detect(match_outcome, "keep"), match_outcome, paste("keep", match_outcome)))
  
  return(mydf_ed_final)
}

pick3FinalTax <- function(mydf){
  # Pick best final id and keep final percent identity
  # Args
  #   mydf: classified dataframe after running compare3Databases() 
  # Output
  #   The function will keep only the higher taxonomy of the best hit, and will add column database_source to track where the match came from
  # Usage: pick3FinalTax(.) 
  # Notes:
  #   Currently there is no function (work in progress) for dealing with assignments to equal, but varying taxonomic levels (any possible values are likely synonyms)
  mydf_ed <- mydf %>% 
    mutate(final_tax = case_when(str_detect(match_outcome, "synonym") ~ NA_character_,
                                 str_detect(match_outcome, "local") ~ tax_string_local,
                                 str_detect(match_outcome, "global1") ~ tax_string_global1,
                                 str_detect(match_outcome, "global2") ~ tax_string_global2,
                                 T ~ "CHECK"),
           taxonomy_source = case_when(str_detect(match_outcome, "synonym") ~ "attention_required",
                                       str_detect(match_outcome, "local") ~ taxonomy_source_local,
                                       str_detect(match_outcome, "global1") ~ taxonomy_source_global1,
                                       str_detect(match_outcome, "global2") ~ taxonomy_source_global2,
                                       T ~ "CHECK"),
           database_source = case_when(str_detect(match_outcome, "synonym") ~ "attention_required",
                                       str_detect(match_outcome, "local") ~ local_name,
                                       str_detect(match_outcome, "global1") ~ global1_name,
                                       str_detect(match_outcome, "global2") ~ global2_name,
                                       T ~ "CHECK"),
           match_name = case_when(str_detect(match_outcome, "synonym") ~ "attention_required",
                                  str_detect(match_outcome, "local") ~ result_seqid_local,
                                  str_detect(match_outcome, "global1") ~ result_seqid_global1,
                                  str_detect(match_outcome, "global2") ~ result_seqid_global2,
                                  T ~ "CHECK"),
           percent_identity_final = case_when(str_detect(match_outcome, "synonym") ~ NA,
                                              str_detect(match_outcome, "local") ~ percent_identity_local,
                                              str_detect(match_outcome, "global1") ~ percent_identity_global1,
                                              str_detect(match_outcome, "global2") ~ percent_identity_global2,
                                              T ~ NA)) %>%
    select(query_seqid, database_source, match_name, percent_identity_final, final_tax, taxonomy_source) %>% 
    unpackTaxString()
  
  if(any(mydf_ed$taxonomy_source == "attention_required")) {
    
    warning("Some taxonomic assignments differ at the lowest taxonomic level. Please manually check for synonyms or another lowest common taxonomic level before proceeding.")
    
  }
  # return modified dataframe that includes only the final chosen hit from each database
  return(mydf_ed)
  
}


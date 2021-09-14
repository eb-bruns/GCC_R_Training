################################################################################
### Author: Emily Beckman  ###  Date: 09/14/2021

### DESCRIPTION:
  ## Takes SurveyMonkey Excel export and format for analysis. Specially-designed
  #   to work with a "conservation action questionnaire" that has multiple
  #   matrix questions for multiple species and regions. SurveyMonkey exports
  #   such questions in a very wide (horizontal) format, and this script puts
  #   the data in a more vertical format and aligns results with respondent
  #   metadata.

### INPUTS:
  ## SurveyMonkey Excel (xlsx) export from survey with multiple matrices,
  #   species, and regions.

### OUTPUTS:
  ## Summary CSV for each matrix question (Question1Results.csv,
  #   Question2Results.csv, etc.)
  ## One CSV with results from all matrix questions (AllResults.csv)
  ## Notes from open-ended question about additional target species
  #   (AdditionalSpeciesNotes.csv)


################################################################################


## SET UP YOUR WORKSPACE

# set working directory
main_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/Magnolia_ConservationActionQuestionnaire"
survey_export_file <- "Magnolia_Conservation_Actions_Questionnaire.xlsx"

# load libraries
library("tidyverse")
library("readxl")


## READ IN DATA

# read in SurveyMonkey export
raw_data <- read_excel(file.path(main_dir,survey_export_file))


## CREATE RESPONDENT INFO DATAFRAME

# list of names for the respondent info columns (can use X1, X2, X3, etc. for
#   cols you don't want)
respondent_col_names <- c("respondent_id","collector_id","start_date",
  "end_date","ip_address","X1","X2","X3","X4","dont_make_public","name",
  "institution","X5","X6","X7","X8","X9","country","email_address","X10",
  "institution_category","institution_category_other")

# select and rename respondent info columns
respondent_info <- raw_data[2:nrow(raw_data),1:length(respondent_col_names)]
for(i in 1:ncol(respondent_info)){
  names(respondent_info)[i] <- respondent_col_names[i] }

# remove respondent info columns with "X" in the name
t <- grepl("^X",names(respondent_info))
if(length(unique(t))>1){
  respondent_info <- respondent_info[, -grep("^X", names(respondent_info))] }
str(respondent_info)


## ORGANIZE DATA BY QUESTIONN TYPE AND ADD REGION COLUMN

## change the following as needed, based on your quesionnaire:

# list of regions (in order!)
regions <- c("Mexico & Central America","Caribbean","South America","East Asia",
  "South and Southeast Asia")

# unique part of each main question asked
select_sp <- "Select all species from"
Q1_match <- "Select all conservation activities your institution participates in"
Q2_match <- "Select what you see as the most urgent conservation activity"
Q3_match <- "Select what you see as the most significant threat"

# create dataframes for collecting results, one for each main question
  # QUESTION1
Q1_results <- data.frame(
  "species" = as.character(NA),
  "Collect and distribute germplasm" = as.character(NA),
  "Conservation horticulture" = as.character(NA),
  "Cryopreservation and/or micropropagation" = as.character(NA),
  "Habitat restoration" = as.character(NA),
  "Implement protection policies or regulations" = as.character(NA),
  "Occurrence surveys or population monitoring" = as.character(NA),
  "Pollen and/or seed banking" = as.character(NA),
  "Population reinforcement or introduction" = as.character(NA),
  "Protect and/or manage habitat" = as.character(NA),
  "Public awareness or education" = as.character(NA),
  "Research: Climate Change" = as.character(NA),
  "Research: Genetics" = as.character(NA),
  "Research: Pests & Pathogens" = as.character(NA),
  "Research: Taxonomy" = as.character(NA),
    # metadata you want to add
  "region" = as.character(NA),
  "respondent_id" = as.double(NA))
head(Q1_results)
  # QUESTION2
Q2_results <- data.frame(
  "species" = as.character(NA),
  "Collect and distribute germplasm" = as.character(NA),
  "Conservation horticulture" = as.character(NA),
  "Cryopreservation and/or micropropagation" = as.character(NA),
  "Habitat restoration" = as.character(NA),
  "Implement protection policies or regulations" = as.character(NA),
  "Occurrence surveys or population monitoring" = as.character(NA),
  "Pollen and/or seed banking" = as.character(NA),
  "Population reinforcement or introduction" = as.character(NA),
  "Protect and/or manage habitat" = as.character(NA),
  "Public awareness or education" = as.character(NA),
  "Research: Climate Change" = as.character(NA),
  "Research: Genetics" = as.character(NA),
  "Research: Pests & Pathogens" = as.character(NA),
  "Research: Taxonomy" = as.character(NA),
  "Unknown" = as.character(NA),
  "None (no conservation actions currently needed)" = as.character(NA),
    # metadata you want to add
  "region" = as.character(NA),
  "respondent_id" = as.double(NA))
head(Q2_results)
  # QUESTION3
Q3_results <- data.frame(
  "species" = as.character(NA),
  "Agriculture, silviculture, and/or ranching" = as.character(NA),
  "Climate change" = as.character(NA),
  "Development, mining, and/or roads" = as.character(NA),
  "Disturbance regime modification" = as.character(NA),
  "Inbreeding or introgression" = as.character(NA),
  "Invasive species competition" = as.character(NA),
  "Pests or pathogens" = as.character(NA),
  "Tourism or recreation" = as.character(NA),
  "Wild harvesting" = as.character(NA),
  "Unknown" = as.character(NA),
    # metadata you want to add
  "region" = as.character(NA),
  "respondent_id" = as.double(NA))
head(Q3_results)

# create list of results dataframes to cycle through
results <- list(Q1_results,Q2_results,Q3_results)

# create list of column numbers to cycle through
start <- list(
  grep(Q1_match,colnames(raw_data)),
  grep(Q2_match,colnames(raw_data)),
  grep(Q3_match,colnames(raw_data)),
  ncol(raw_data)+1)

num_sp <- grep(select_sp,colnames(raw_data))


### LOOP THROUGH AUTOMATICALLY ###

# FOR EACH QUESTIONN
for(q in 1:length(results)){
  # get the number of categories for the question
  num_fields <- ncol(results[[q]])-3 # -3 removes metadata cols (not categories)
  # results dataframe for the question
  df_out <- results[[q]]

  # FOR EACH REGION
  for(r in 1:length(regions)){
    # get the raw data for the current question and region only
    data <- raw_data[,start[[q]][r]:(start[[q]][r]+((start[[1]][r]-num_sp[[r]]-1)*num_fields)-1)]
    # use the first row to assign column names
    colnames(data) = data[1, ]
    # leave just species name in first row (cut off rest)
    data <- as.data.frame(lapply(data,function(x) gsub("( -).*", "", x)),
      stringsAsFactors=F)
    # cycle through each species and add the responses to a dataframe
    start_col <- 1
    end_col <- num_fields

    # FOR EACH SPECIES
    for(s in 1:(start[[1]][r]-num_sp[[r]]-1)){ #number of species
      add <- data[2:nrow(data),start_col:end_col]
      sp_add <- data.frame("species" = rep(data[1,start_col],
        times=nrow(data)-1))
      add <- cbind(sp_add,add)
      # add respondent info and region
      add$region <- regions[[r]]
      add$respondent_id <- respondent_info$respondent_id
      # add species to running dataframe for the question
      df_out <- rbind(df_out,setNames(add,names(df_out)))
      # set up index for next species
      start_col <- start_col + num_fields
      end_col <- end_col + num_fields
    }
    # remove rows with no data
    df_out <- df_out[rowSums(is.na(df_out[,2:ncol(df_out)]))<ncol(df_out)-3,]
  }
  results[[q]] <- df_out
}


## SAVE ALL RESULTS IN ONE FILE:

# join all questions
join_all <- results[[1]]
join_all <- join_all[,c(ncol(join_all),(ncol(join_all)-1),1,2:(ncol(Q1_results)-3))]
for(i in 1:(length(results)-1)){
  join_all <- left_join(join_all,results[[i+1]],
    by=c("respondent_id","region","species"))
}
join_all <- right_join(respondent_info,join_all)
write.csv(join_all,file.path(main_dir,"AllResults.csv"),row.names = F)


## SAVE RESULTS IN SEPARATE FILES, BY QUESTIONN:

# join all respondent info to each dataframe and write files
for(i in 1:length(results)){
  results[[i]] <- left_join(results[[i]],respondent_info)
  write.csv(results[[i]],file.path(main_dir,paste0("Question",i,"Results.csv")),
    row.names = F)
}


## GET NOTES ABOUT ADDITIONAL SPECIES

# dataframe with notes about additional target species suggested, one per region
  # phrase to match
Q4_match <- "You may also provide names of other species you would add to the above list"
  # find locations of phrase
t <- grep(Q4_match,colnames(raw_data))
  # create dataframe
add_species <- data.frame(
  "Additional species and reason" = as.character(NA),
  "region" = as.character(NA),
  "respondent_id" = as.double(NA))
for(i in 1:length(t)){
  data <- raw_data[2:nrow(raw_data),t[i]]
  colnames(data) <- "Additional.species.and.reason"
  data$region <- regions[[i]]
  data$respondent_id <- respondent_info$respondent_id
  add_species <- rbind(add_species,data)
}
  # remove rows with no data
add_species <- add_species[!is.na(add_species[,1]),]
  # join to rest of respondent info
add_sp_join <- left_join(add_species,respondent_info)
head(add_sp_join)
  # write file
write.csv(add_sp_join,file.path(main_dir,"AdditionalSpeciesNotes.csv"),
  row.names = F)

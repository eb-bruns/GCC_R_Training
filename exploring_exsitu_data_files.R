


################################################################################
# Load libraries
################################################################################

rm(list=ls())
my.packages <- c('plyr', 'tidyverse', 'textclean', 'spatialEco', 'maps',
                 'measurements', 'CoordinateCleaner', 'rnaturalearth', 'raster')
select <- dplyr::select
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
# when running above script you may get an error that naniar cant be loaded.
# Load separately:
  # install.packages("naniar")
  # library(nanier)
rm(my.packages)


################################################################################
# Set working directory
################################################################################

# set manually:
main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/Global Tree Conservation Program/4. GTCP_Projects/Gap Analyses/ALL ACCESSION-LEVEL DATA"


################################################################################
# Load functions
################################################################################

# function to read in ex situ files from different folders/years and stack
read.exsitu.csv <- function(path,submission_year){
  # create list of paths to ex situ accessions CSV files in folder
  file_list <- list.files(path=path,pattern=".csv",full.names=TRUE)
  # read in each csv in path list to create list of dataframes
  file_dfs <- lapply(file_list,read.csv,header=TRUE,fileEncoding="LATIN1",
                     strip.white=TRUE,colClasses="character",na.strings=c("","NA"))
  print(paste0("Number of files: ",length(file_dfs)))
  #sapply(file_dfs, nrow) # can look at number of rows in each csv
  for(file in seq_along(file_dfs)){
    df <- file_dfs[[file]]
    # add file name as column, to record home institution for each record
    df$filename <- rep(file_list[file],nrow(df))
    # remove file path portion
    df$filename <- mgsub(df$filename,c(paste0(path,"/"),".csv"),"")
    # add year of submission
    df$submission_year <- submission_year
    # remove extra blank columns
    t <- grepl("^X",names(df))
    if(length(unique(t))>1){
      #print(df$filename[1])
      df <- df[, -grep("^X", names(df))]
    }
    # add accession number if there isn't one
    if("acc_num" %in% names(df) & nrow(df[which(is.na(df$acc_num)),]) > 0){
      df[which(is.na(df$acc_num)),]$acc_num <- paste0("added",
                                                      sprintf("%04d", 1:nrow(df[which(is.na(df$acc_num)),])))
      #print(nrow(df))
    } else if ("acc_no" %in% names(df) & nrow(df[which(is.na(df$acc_no)),]) > 0){
      df[which(is.na(df$acc_no)),]$acc_no <- paste0("added",
                                                    sprintf("%04d", 1:nrow(df[which(is.na(df$acc_no)),])))
      #print(nrow(df))
    } else if (!("acc_num" %in% names(df)) & !("acc_no" %in% names(df))){
      df$acc_num <- paste0("added", sprintf("%04d", 1:nrow(df)))
      #print(nrow(df))
    } else {
      #print(paste("NO ACC NUM EDITS:",df$filename[1]))
    }
    # replace old df with new df
    file_dfs[[file]] <- df
    #print(head(file_dfs[[file]],n=2))
  }
  # stack all datasets using rbind.fill, which keeps non-matching columns
  #   and fills with NA; 'Reduce' iterates through and merges with previous
  # this may take a few minutes if you have lots of data
  all_data <- Reduce(rbind.fill, file_dfs)
  print(paste0("Number of rows: ",nrow(all_data)))
  print(paste0("Number of columns: ",ncol(all_data)))
  return(all_data)
}

# remove duplicates in network datasets (e.g., PCN) if institution submitted
#  their own data separately
remove.network.dups <- function(df,rm_inst_names,file_name){
  df <- df[!(df$filename == file_name & df$inst_short %in% rm_inst_names),]
  return(df)
}


################################################################################
################################################################################
# Read in and stack all accessions data
################################################################################

## read in data from multiple surveys and stack, or just read in from one folder
##    this function also adds columns for 1) the file name [equivalent to the
##    "inst_short" institution nickname] 2) a sumbission year, 3) an accession
##    number if one isn't given
## Warnings are usually ok here, but you can look at the file causing the
#     warning to see if there is an obvious formatting issue
raw_2017 <- read.exsitu.csv(file.path(main_dir,
  "exsitu_standard_column_names","2017"), "2017")
  # if genus is blank in 2017 data, make it "Quercus"
  raw_2017[which(is.na(raw_2017$genus)),]$genus <- "Quercus"
raw_2018 <- read.exsitu.csv(file.path(main_dir,
  "exsitu_standard_column_names","2018"), "2018")
raw_2019 <- read.exsitu.csv(file.path(main_dir,
  "exsitu_standard_column_names","2019"), "2019")
raw_2020 <- read.exsitu.csv(file.path(main_dir,
  "exsitu_standard_column_names","2020"), "2020")
raw_2021 <- read.exsitu.csv(file.path(main_dir,
  "exsitu_standard_column_names","2021"), "2021")
raw_2022 <- read.exsitu.csv(file.path(main_dir,
  "exsitu_standard_column_names","2022"), "2022")

# stack all data
to_stack <- list(raw_2022,raw_2021,raw_2020,raw_2019,raw_2018,raw_2017)
all_data_raw <- Reduce(rbind.fill,to_stack)

# create new version before big changes, so can easily go back to original
all_data <- all_data_raw

# check out column names
sort(colnames(all_data))
# merge similar columns (you may not need to do this if no schema
#   mistakes were made when manually editing column names).
  all_data <- tidyr::unite(all_data,"inst_short",
    c("inst_short","ï..inst_short","inst_short2","inst_short"),
    sep=";",remove=T,na.rm=T)
  all_data <- tidyr::unite(all_data,"taxon_full_name",
    c("taxon_full_name","ï..taxon_full_name","full_sp_name","ï..sp_full_name","sp_full_name","taxon_name"),
    sep=";",remove=T,na.rm=T)
  all_data <- tidyr::unite(all_data,"infra_name",
    c("infra_name","intra_name","specific_name","specific","Infraspecific.Epithet"),
    sep=";",remove=T,na.rm=T)
  all_data <- tidyr::unite(all_data,"infra_rank",
    c("infra_rank","intra_rank","specific_rank","specific2","Infraspecific.Rank"),
    sep=";",remove=T,na.rm=T)
  all_data <- tidyr::unite(all_data,"cultivar",
    c("cultivar","Cultivar.Epithet"),
    sep=";",remove=T,na.rm=T)
  all_data <- tidyr::unite(all_data,"hybrid",
    c("hybrid","Hybrid","hybrid.1"),
    sep=";",remove=T,na.rm=T)
# check again
sort(colnames(all_data))

# fill in inst_short column with filename if none provided
all_data$inst_short[all_data$inst_short==""] <-
  all_data$filename[all_data$inst_short==""]
# remove rows with no inst_short
all_data <- all_data[which(all_data$inst_short!=""),]

### CHECK ALL INSTITUTIONS ARE HERE ###
sort(unique(all_data$inst_short)) #336

# remove leading, trailing, and middle (e.g., double space) whitespace,
#   to prevent future errors
all_data <- as.data.frame(lapply(all_data, function(x) str_squish(x)),
  stringsAsFactors=F)
# replace "" cells with NA in whole dataset
all_data[all_data == ""] <- NA


################################################################################
# Compile taxon names
################################################################################

# preserve original taxon name
all_data$taxon_full_name_orig <- all_data$taxon_full_name

# fill genus column if not already filled
all_data <- all_data %>% separate("taxon_full_name","genus_temp",sep=" ",remove=F)
all_data[which(is.na(all_data$genus)),]$genus <- all_data[which(is.na(all_data$genus)),]$genus_temp
# standardize capitalization
all_data$genus <- str_to_title(all_data$genus)

### MAKE SURE NO GENUS MISSPELLINGS OR ABBREVIATIONS ###
sort(unique(all_data$genus))
all_data$genus <- mgsub(all_data$genus,
  c("^Q$","^Q\\.$","Querces","Quercues","Querus","Qurercus","Querucs","Cyclobalanopsis"),
  "Quercus",fixed=F)

# remove rows not in target genus/genera
all_data_preQ <- all_data
target_genera <- c("Quercus")
all_data <- all_data %>% filter(genus %in% target_genera)
nrow(all_data_preQ) #260356
nrow(all_data) #79886

# create concatenated taxon_full_name column
all_data <- tidyr::unite(all_data, "taxon_full_name_concat",
  c(genus,hybrid,species,infra_rank,infra_name,cultivar),
  sep=" ", remove=F, na.rm=T)

# when blank, fill taxon_full_name column with concatenated full name
all_data[is.na(all_data$taxon_full_name),]$taxon_full_name <-
  all_data[is.na(all_data$taxon_full_name),]$taxon_full_name_concat
unique(all_data$taxon_full_name)
all_data$taxon_full_name <- str_squish(all_data$taxon_full_name)

sort(unique(all_data$inst_short)) #4326


################################################################################
# Look at which instituions have Quercus, and which years the data are from
################################################################################

summary <- unique(data.frame(inst_short = all_data$inst_short,
  genera = all_data$genus,
  submission_year = all_data$submission_year))
summary <- summary %>%
  arrange(genera, submission_year) %>%
  #group_by(inst_short) %>%
  #mutate(#genera = paste(genera,collapse = ", "),
  #       submission_year = paste(submission_year,collapse = ", ")) %>%
  #ungroup() %>%
  arrange(inst_short) %>%
  distinct(inst_short,submission_year)#.keep_all=T)
as.data.frame(summary)

# since there are so many, lets see which do NOT have Quercus data
#all_inst <- unique(all_data_preQ$inst_short)
#Q_inst <- unique(all_data$inst_short)
#sort(all_inst[!all_inst %in% Q_inst])

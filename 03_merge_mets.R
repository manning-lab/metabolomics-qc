library(dplyr)
library(tidyr)

##Read data over
cn <- read.csv("results/no_batch_adj/plasma/QCd_data/inv_norm/DoD_c18_plasma_QCd_inv_norm.csv")
cp <- read.csv("results/no_batch_adj/plasma/QCd_data/inv_norm/DoD_c8_plasma_QCd_inv_norm.csv")
hn <- read.csv("results/no_batch_adj/plasma/QCd_data/inv_norm/DoD_hilic_neg_plasma_QCd_inv_norm.csv")
hp <- read.csv("results/no_batch_adj/plasma/QCd_data/inv_norm/DoD_hilic_pos_plasma_QCd_inv_norm.csv")
output_fnm <-  "results/no_batch_adj/plasma/merged_data/inv_norm_merged_QCd_plasma_knowns.csv"
met_info <- fread("results/no_batch_adj/plasma/met_info.csv",na.strings="")
#cn <- read.csv("results/no_batch_adj/muscle/QCd_data/inv_norm/DoD_c18_muscle_QCd_inv_norm.csv")
#cp <- read.csv("results/no_batch_adj/muscle/QCd_data/inv_norm/DoD_c8_muscle_QCd_inv_norm.csv")
#hn <- read.csv("results/no_batch_adj/muscle/QCd_data/inv_norm/DoD_hilic_neg_muscle_QCd_inv_norm.csv")
#hp <- read.csv("results/no_batch_adj/muscle/QCd_data/inv_norm/DoD_hilic_pos_muscle_QCd_inv_norm.csv")
#output_fnm <-  "results/no_batch_adj/muscle/merged_data/inv_norm_merged_QCd_muscle_knowns.csv"
#met_info <- fread("results/no_batch_adj/muscle/met_info.csv",na.strings="")
dir.create(dirname(output_fnm))

names(cn)[1] <- "Compound_Id"
names(cp)[1] <- "Compound_Id"
names(hn)[1] <- "Compound_Id"
names(hp)[1] <- "Compound_Id"

# Deal with muscle "_repeat" samples in muscle. Remove duplicates (preferring the samples with _repeat), then rename the ids to remove "_repeat".
fix_repeats <- function(df) {
  id_no_repeat <- sub("_repeat","",colnames(df))
  cols2keep <- grepl("_repeat",colnames(df))   |   !( duplicated(id_no_repeat) | duplicated(id_no_repeat,fromLast=T) )
  df <- df[,cols2keep]
  colnames(df) <- sub("_repeat","",colnames(df))
  df
}

cn <- cn %>%
  mutate(Compound_Id = paste(Compound_Id, "cn", sep = "_")) %>% fix_repeats

cp <- cp %>%
  mutate(Compound_Id = paste(Compound_Id, "cp", sep = "_")) %>% fix_repeats

hn <- hn %>%
  mutate(Compound_Id = paste(Compound_Id, "hn", sep = "_")) %>% fix_repeats

hp <- hp %>%
  mutate(Compound_Id = paste(Compound_Id, "hp", sep = "_")) %>% fix_repeats


## Merging the four methods
# Assuming df1 and df2 are your data frames
merged_df <- bind_rows(cn, cp, hn, hp)

#Selecting only relevant columns for merging
info <- met_info[, c("Compound_Id","Name", "HMDB_Id")]

## Merging with metab info
data <- merge(info, merged_df, by = "Compound_Id" )

#Write the files
#write.csv(data, "results/merged_data/inv_norm_merged_QCd_plasma.csv")

## Selecting only known metabolites
selected_rows <- data[!is.na(HMDB_Id) & HMDB_Id!="NA",]

## Remove columns and keeping only HMDB_Id column
df <- subset(selected_rows, select = -c(Compound_Id,Name))

#Remove the controls
df <- df[!grepl("PREF",HMDB_Idl)]
# Removing metabolite duplicates
df <- distinct(df, HMDB_Id, .keep_all = TRUE)

# Transposing data samples as rows, metbaolites as columns
transposed_df <- pivot_longer(df, cols = -HMDB_Id, names_to = "sample_id", values_to = "Value")
transposed_df <- pivot_wider(transposed_df, names_from = HMDB_Id, values_from = Value)

fwrite(transposed_df, output_fnm)

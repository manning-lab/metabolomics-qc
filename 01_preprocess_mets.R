invisible(sapply(c("data.table","dplyr","tidyverse"), library, character.only=T))

options(stringsAsFactors=FALSE)
res_path     <- "results/plasma/"
data_path    <- "data/plasma/"
#res_path     <- "results/muscle/"
#data_path    <- "data/muscle/"
dir.create(file.path(paste0(res_path,"/formatted_data/")), recursive=T)


# Functions ---------------------------------------------------------------
`%!in%` <- Negate(`%in%`)

rmDtRows <- function(DT, del.idxs) {
  keep.idxs <- setdiff(DT[, .I], del.idxs); # row indexes to keep
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]); # this is the subsetted table
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}

keep_first <- function(dt, i, IDs) {
  IDs <- IDs[grep("HMDB0.*", IDs)]
  for (ID in IDs){
    inds <- which(dt[[i]][,hmdb_id_cols[i],with=F]==ID); dt[[i]][inds,]
    dt[[i]] <- rmDtRows(dt[[i]], inds[2]) # Just keep the first entry
  }
  return(dt[[i]])
}

selectMetInfoCols  <- function(dts, col_inds) {
  lapply(seq_along(dts), function(i) {
    if(is.na(col_inds[i])) { rep(NA,times=nrow(dts[[i]])-data_start_rows[i]+1) }
    else                   { dts[[i]][data_start_rows[i]:nrow(dts[[i]]), col_inds[i], with=F] }
  })
}

renameMethods <- function(method_vec) {
  ifelse(method_vec=="C8-pos", "cp",
  ifelse(method_vec=="C18-neg", "cn",
  ifelse(method_vec=="HIL-neg" | method_vec=="HILIC-neg", "hn",
  ifelse(method_vec=="HIL-pos" | method_vec=="HILIC-pos", "hp", NA))))
}

addMethodSuffix <- function(compound_id_vec, method_vec) {
  mapply(compound_id_vec, method_vec, FUN=function(id,method) {
    if(is.na(id)) return(NA)
    if(is.na(method)) return(id)
    
    if(method=="C8-pos") return(paste0(id,"_cp"))
    if(method=="C18-neg") return(paste0(id,"_cn"))
    if(method=="HIL-pos") return(paste0(id,"_hp"))
    if(method=="HILIC-pos") return(paste0(id,"_hp"))
    if(method=="HILIC-neg") return(paste0(id,"_hn"))
    if(method=="HIL-neg") return(paste0(id,"_hn"))
  })
}

# Set up ------------------------------------------------------------------

# Read in data
data_paths <- list.files(data_path, pattern = ".csv", full.names = T)
dts <- lapply(data_paths, fread)
names(dts) <- gsub("\\..*", "", basename(data_paths))

# Data row start
data_start_rows <- unlist(lapply(dts, function(x) which(x[,1]!="")[1]+1))
data_start_cols <- unlist(lapply(seq_along(dts), function(x) which(as.numeric(dts[[x]][data_start_rows[x],])%%1==0 &
                                                                   as.numeric(dts[[x]][data_start_rows[x],]) > 2)[1]))

# Sample-related stuff all starts @ data_start_COL
sample_id_rows <- unlist(lapply(dts, function(x) which(x[,1]!="")[1]))
extr_date_rows <- rep(c(3), length(dts))

# Metabolite-related stuff all starts @ data_start_ROW
method_cols      <- unlist(lapply(dts, function(x) which(x == "Method"|x == "method", arr.ind = TRUE)[2]))
compound_id_cols <- unlist(lapply(dts, function(x) which(grepl("compound", x, ignore.case = T), arr.ind = TRUE)[1]))
mrm_cols <- unlist(lapply(dts, function(x) which(grepl("mrm", x, ignore.case = T), arr.ind = TRUE)[1]))
mz_cols  <- unlist(lapply(dts, function(x) which(grepl("mz", x, ignore.case = T), arr.ind = TRUE)[1]))
rt_cols  <- unlist(lapply(dts, function(x) which(grepl("rt", x, ignore.case = T), arr.ind = TRUE)[1]))
hmdb_id_cols  <- unlist(lapply(dts, function(x) which(grepl("hmdb", x, ignore.case = T), arr.ind = TRUE)[1]))
met_name_cols <- unlist(lapply(dts, function(x) which(grepl("metabolite", x, ignore.case = T), arr.ind = TRUE)[1]))

# Clean up data ------------------------------------------------------------------

for(dti in 1:length(dts)) {
  hmdbs_col <- unlist(dts[[dti]][,hmdb_id_cols[dti],with=F])
  hmdbs_col <- sub("\\*","",hmdbs_col)
  for(r in 1:nrow(dts[[dti]]))
    set(dts[[dti]], i=r, j=as.integer(hmdb_id_cols[dti]), hmdbs_col[r])
}

# Duplicated metabolites
ID_list <- lapply(seq_along(dts), function(i) {unname(unlist(unique(dts[[i]][,hmdb_id_cols[i], with=F])))})
dts <- lapply(seq_along(dts), function(i) {keep_first(dts, i, ID_list[[i]])})

dts_clean <- list()
sample_idss <- lapply(seq_along(dts), function(i) unlist(dts[[i]][ sample_id_rows[i], data_start_cols[i]:ncol(dts[[i]]) ]))

# Use compound id w/ method suffix, but for amines must use HMDB id.
met_idss <- lapply(seq_along(dts), function(i) {
  tmp <- if(is.na(compound_id_cols[i])) {unlist( dts[[i]][data_start_rows[i]:nrow(dts[[i]]), hmdb_id_cols[i]    , with=F] )}
  else                    {unlist( dts[[i]][data_start_rows[i]:nrow(dts[[i]]), compound_id_cols[i], with=F] )}
  tmp[duplicated(tmp)] <- paste0(tmp[duplicated(tmp)], "_dup") # If an ID is repeated
  return(tmp)
})

for(i in seq_along(dts)) {
  dts_clean[[i]] <- dts[[i]][data_start_rows[i]:nrow(dts[[i]]),
                               data_start_cols[i]:ncol(dts[[i]])]
  
  # Now data is numbers only, convert to numeric to save memory
  # (but not integer, b/c some counts are > MAX_INT and also amine measurements are floats)
  for(col in 1:ncol(dts_clean[[i]])) set(dts_clean[[i]], j=col, value=as.numeric(dts_clean[[i]][[col]]))
  
  # tranform, and naming
  dts_clean[[i]] <- t(dts_clean[[i]])
  
  colnames(dts_clean[[i]]) <- met_idss[[i]]
  dts_clean[[i]] <- data.table(sample_id=sample_idss[[i]], dts_clean[[i]])
  dts_clean[[i]] <- dts_clean[[i]][!is.na(sample_id)]
}

names(dts) <- gsub("\\..*", "", basename(data_paths))
names(dts_clean) <- names(dts)


# Write cleaned data ------------------------------------------------------

for(i in seq_along(dts_clean)) { write.table(dts_clean[[i]], paste0(res_path, "/formatted_data/", names(dts_clean)[i],"_formatted.txt"), row.names=F) }

# Get metabolite and sample info ------------------------------------------------------------------

extr_datess <- lapply(seq_along(dts), function(i) unlist(dts[[i]][ extr_date_rows[i], data_start_cols[i]:ncol(dts[[i]]) ]))
sample_idss <- lapply(seq_along(dts), function(i) unlist(dts[[i]][ sample_id_rows[i], data_start_cols[i]:ncol(dts[[i]]) ]))
unique_ids <- unique(unlist(sample_idss))

# sample id & cohort cols
sample_info <- data.table(matrix( "", nrow=length(unique_ids), ncol=length(dts)+1 ))
colnames(sample_info) <- c("sample_id", paste0(names(dts),"_batch"))
sample_info$sample_id <- unique_ids

for(i in seq_along(extr_datess)) {
  for(j in seq_along(extr_datess[[i]])) {
    sample_info[sample_id==sample_idss[[i]][j], i+1] <- extr_datess[[i]][j]
}}

sample_info[, c(2:(length(dts)+1) ):= lapply(.SD, as.factor), .SDcols=2:(length(dts)+1)] # Convert batch columns to factors

sample_info[, is_control := !grepl("MCD",sample_id)]
sample_info <- sample_info[!is.na(sample_info$sample_id)]

met_info <- data.table(
  Compound_Id = unlist(selectMetInfoCols(dts, compound_id_cols)),
  HMDB_Id     = unlist(selectMetInfoCols(dts, hmdb_id_cols)),
  Name        = unlist(selectMetInfoCols(dts, met_name_cols)),
  MRM         = unlist(selectMetInfoCols(dts, mrm_cols)),
  MZ          = as.numeric(unlist(selectMetInfoCols(dts,mz_cols))),
  RT          = as.numeric(unlist(selectMetInfoCols(dts,rt_cols))),
  Method      = unlist(selectMetInfoCols(dts, method_cols))
)

met_info[met_info==""] <- NA
met_info[grepl("tandard",HMDB_Id), HMDB_Id:=NA]
met_info[met_info=="n/a"] <- NA
met_info <- met_info[!is.na(Compound_Id)]
met_info[nchar(HMDB_Id)!=11, HMDB_Id := sub("HMDB","HMDB00", HMDB_Id)]

# Merging knowns across cohorts where the matches exactly (both HMDB and method). 
met_info$Compound_Id <- addMethodSuffix(met_info$Compound_Id,met_info$Method )
met_info$HMDB_Id <- addMethodSuffix(met_info$HMDB_Id,met_info$Method )

# Change bulky method names in the Method column.
met_info$Method  <- renameMethods(met_info$Method)

# Write sample and metabolite info files -------------------------------------------------------------

write.csv(sample_info %>% as.data.frame(), paste0(res_path, "/sample_info.csv"),  row.names=F)
fwrite(met_info, paste0(res_path,"met_info.csv"))

# Messages ---------------------------------------------------------

print(paste0("Info written to:", res_path, ""))
print(paste0("Data written to:", res_path, "formatted_data/"))

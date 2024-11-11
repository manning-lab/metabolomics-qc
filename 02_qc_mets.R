invisible(sapply(c("GGally","matrixStats","ggplot2","data.table","tibble","tidyverse","janitor","gridExtra","swamp"), library, character.only=T))

options(stringsAsFactors=FALSE)
res_path     <- "results/no_batch_adj/plasma/"
#res_path     <- "results/no_batch_adj/muscle/"
dir.create(file.path(paste0(res_path,"/QCd_data/")), recursive=T)

# Functions ---------------------------------------------------------------
stat.mode <- function(x) { names(sort(-table(x)))[1] } # Statistical mode. Tabulate then sort then pick the first item.

t2 <- function(x) {
  rn <- rownames(x); cn <- colnames(x)
  x <- t(x)
  rownames(x) <- cn; colnames(x) <- rn
  return(x)
}

#Limited to categorical variables for now
PCAPlot <- function(df, color_var, n_PC, title="") {
  pca <- prcomp(df, center=T, scale.=T)
  
  var_explained <- scales::percent(summary(pca)$importance[2,], 0.1)
  for(i in 1:n_PC) colnames(pca$x)[i] <- paste0("PC",i," (",var_explained[i],")")
  
  ggpairs(data = data.frame(pca$x[,1:n_PC], check.names=F),
          lower = list(continuous = wrap("points", alpha=0.5, size=0.5, pch=1)),
          diag  = list(continuous = wrap("densityDiag", alpha=0.5)),
          #upper= list(continuous = wrap("cor"),
          axisLabels = "none",
          mapping = aes(color = color_var),
          legend  = c(2,1) ) + labs(color="legend_title") + ggtitle(title)
}
rm0VarRows <- function(df) {
  rows2keep <- rowVars(df, na.rm=T) > 0
  rows2keep[is.na(rows2keep)] <- FALSE
  print(paste(nrow(df)-sum(rows2keep),"/",nrow(df),"rows will be removed for having 0 variance"))
  df[rows2keep,]
}

rmHighMissingness <- function(df, thresh, verbosity=1) {
  rows2keep <- rowSums(is.na(df)) < thresh*ncol(df)
  if(verbosity>0) print(paste( sum(!rows2keep),"/",nrow(df),"rows will be removed for >",scales::percent(thresh),"missingness" ))
  if(verbosity>2 & sum(!rows2keep>0)) print(paste( "Rows removed:",names(df)[!rows2keep] ))
  df <- df[rows2keep,]
  
  cols2keep <- colSums(is.na(df)) < thresh*nrow(df)
  if(verbosity>0) print(paste( sum(!cols2keep),"/",ncol(df),"cols will be removed for >",scales::percent(thresh),"missingness" ))
  if(verbosity>1 & sum(!cols2keep>0)) print(paste( "Cols removed:",names(df)[!cols2keep] ))
  df <- df[,cols2keep]
}
winsorizeBySd <- function(df, n_sds) {
  t2(apply(df,1, function(r) {
    u <- mean(r) + sd(r)*n_sds
    l <- mean(r) - sd(r)*n_sds
    r <- sapply(r, function(x) {
      if(x > u) {x <- u}
      if(x < l) {x <- l}
      else      {x}})}))
}
batchAdj <- function(df, batch_factor) { # Median scaling. batch_factor should have same cols as df.
  if(length(levels(batch_factor)) < 2) { print("Need at least 2 batch levels to adjust for batch, returning df untouched."); return(df) }
  swamp::quickadjust.ref(df, batches=batch_factor, refbatch=stat.mode(batch_factor))$adjusted.data
}
zScoreRows <- function(df) {
  t2(apply(df,1, function(r) {
    m <- mean(r)
    s <- sd(r)
    r <- sapply(r, function(x) (x-m)/s)
}))}
plotSkewKurt <- function(df, title) {
  n <- ncol(df)
  ms <- rowMeans(df)
  sds <- rowSds(df)
  skews <- sapply(1:nrow(df), function(r) (sum((df[r,]-ms[r])^3)/sds[r]^3)/n   )
  kurts <- sapply(1:nrow(df), function(r) (sum((df[r,]-ms[r])^4)/sds[r]^4)/n-3 )
  sk <- data.frame(skews=skews, kurts=kurts)
  
  ggplot(sk, aes(x=skews, y=kurts)) +
    geom_point(size=0.5) + ggtitle(title) +
    geom_vline(xintercept = -0.5, linetype="dotted", color="blue", size=1) +
    geom_vline(xintercept =  0.5, linetype="dotted", color="blue", size=1) + 
    geom_hline(yintercept = -2.0, linetype="dotted", color="blue", size=1) +
    geom_hline(yintercept =  2.0, linetype="dotted", color="blue", size=1)
}

# Read in data ------------------------------------------------------------
data_files <- list.files(path=paste0(res_path, "/formatted_data/"), pattern=".txt", full.names = T)
dict_files <- list.files(path=res_path, pattern="info.csv", full.names = T)
sample_info <- fread(dict_files[2])

# Format ------------------------------------------------------------------

dfs <- lapply(data_files, function(filename) {
  df <- as.matrix(fread(filename))
  rownames(df) <- df[,"sample_id"]; df <- df[,-1]
  mode(df) <- "numeric"
  df <- t2(df)
  non_control <- colnames(df) %in% sample_info$sample_id[!sample_info$is_control]
  if(sum(non_control)>0) df <- df[,non_control] # only non-control samples
  return(df)
})


batch_info_cols <- grep("batch", colnames(sample_info), value = T)
batch_list <- lapply(1:length((batch_info_cols)), function(i) {
  df <- dfs[[i]]
  batch <- sample_info[sample_id %in% colnames(df), sample_id, eval(batch_info_cols[i])]
  batch <- batch[complete.cases(batch),]
  batch <- batch[match(colnames(df),sample_id),] # Reorder batch to match the order of samples
  batch <- factor(unlist(batch[,1])) # Need only the dates column
})

df_labels <- gsub("_batch", "", batch_info_cols)
names(dfs) <- df_labels
names(batch_list) <- names(dfs)

# QC -----------------------------------------------------------------------

tmp <- mapply(dfs,names(dfs),batch_list, SIMPLIFY=F, FUN=function(df,nm,batch_factor) {
  print(nm)
  df <- rm0VarRows(df) # 1
  df <- rmHighMissingness(df, 0.25) # 2
  df <- t2(apply(df,1, function(r) { r[is.na(r)] <- min(r,na.rm=T)/2; r })) # 3
  
  df_l2 <- log2(df+1) # 5a
  df_l2_z <- zScoreRows(log2(df+1)) # 5b
  df_inv_norm <- t2(apply(df,1, function(r) qnorm( (rank(r)-0.5)/length(r) ) )) # 5c
  df_ln <- log(df+1) # 5d
  df_ln_z <- zScoreRows(log(df+1))

  dfs <- list(default=df, l2=df_l2, l2_z=df_l2_z, inv_norm=df_inv_norm, ln=df_ln, ln_z=df_ln_z)

  dfs <- lapply(dfs, winsorizeBySd, 5) # 4
  dfs <- lapply(dfs, batchAdj, batch_factor) # Comment out for non-adjusted QC'd data for PCA plots
})


# Transform ---------------------------------------------------------------

# tmp is a list of QC transforms (l2, l2_z, etc.) per cohort. Rearrange to List of cohorts per transform. 
dfss <- list(default  = lapply(tmp, function(cohort) cohort[["default" ]]),
             l2       = lapply(tmp, function(cohort) cohort[[   "l2"   ]]),
             l2_z     = lapply(tmp, function(cohort) cohort[[  "l2_z"  ]]),
             inv_norm = lapply(tmp, function(cohort) cohort[["inv_norm"]]),
             ln       = lapply(tmp, function(cohort) cohort[[   "ln"   ]]),
             ln_z     = lapply(tmp, function(cohort) cohort[[  "ln_z"  ]]))

# Figures -----------------------------------------------------------------
outfile <- paste0(res_path, "/QCd_data/QC_figures.pdf")
pdf(outfile)
for(i in seq_along(dfss[1])) {
  plot_data_column = function (col) {
    dfss[[col]][[i]] %>% as.data.frame() %>% mutate(mean=rowMeans(., na.rm = TRUE)) %>% 
      ggplot(aes(x=mean)) + 
      geom_histogram(color="black", fill="white", bins = 30) +
      ggtitle(paste(names(dfss[[col]])[i],names(dfss)[col],"distribution")) +
      xlab(names(dfss)[col])
  }
  myplots <- lapply(seq_along(dfss), plot_data_column)
  do.call("grid.arrange", c(myplots, ncol=2))
}

for(i in seq_along(dfss[1])) {
    plot_data_column = function (col) {
      plotSkewKurt(dfss[[col]][[i]], 
                   title=paste(names(dfss[[col]])[i],names(dfss)[col]))
    }
    myplots <- lapply(seq_along(dfss), plot_data_column)
    do.call("grid.arrange", c(myplots, ncol=2))
  }

for(i in seq_along(dfss)) {
  for(j in seq_along(dfss[[i]])) {
    df <- as.data.frame(dfss[[i]][[j]])
    df_t <- t(df) %>% as.data.frame
    data <- rownames_to_column(df_t, var = "Sample_Id") %>% 
      mutate(group = ifelse(str_detect(Sample_Id,'QC'),1,0))
    
    plot_data_column = function (data, column) {
      ggplot(data, aes(x = 1:nrow(data), y = get(column), col=group))+
          geom_point() + ylab("sample") + xlab(column) + 
          ggtitle(paste(names(dfss[[i]])[j],names(dfss)[i], ":", column)) + 
         theme_bw() + theme(legend.position="none")
    }
    myplots <- lapply(colnames(data)[sample(c(2:length(colnames(data))),4)], plot_data_column, data = data)
    do.call("grid.arrange", c(myplots, ncol=2))
  }}
garbage <- dev.off()

# Write files -------------------------------------------------------------

for(type in names(dfss)) {dir.create(file.path(paste0(res_path, "/QCd_data/",type)), showWarnings=F)}

for(i in seq_along(dfss)) {
  for(j in seq_along(dfss[[i]])) {
    write.csv(dfss[[i]][[j]], paste0(res_path, "/QCd_data/", names(dfss)[i],"/", names(dfss[[i]])[j],"_QCd_",names(dfss)[i],".csv"))
}}

# Message -----------------------------------------------------------------
print(paste("QC written to:", paste0(res_path, "/QCd_data/")))

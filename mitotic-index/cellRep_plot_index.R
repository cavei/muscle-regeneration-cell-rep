computeMitoticIndex <- function(exprs, migenes) {
  # migenes <- c("Cdkn3","Ilf2","Kdelr2","Rfc4","Top2a","Mcm3","Kpna2","Cks2","Cdk1")
  missing <- setdiff(migenes, row.names(exprs))
  if (length(missing)>0)
    warning(paste0(paste(missing, collapse=", "), " gene was not found in "))
  usable <- intersect(migenes, row.names(exprs))
  colMeans(exprs[usable, ])
}

splitNames <- function(names) {
  lapply(strsplit(names, "_"), function(x) {c(time=x[2], type=x[3])})
}

extract_days <- function(names) {
  unname(unlist(lapply(strsplit(names, "_"), function(x) x[2])))
}

extract_conditions <- function(names) {
  unname(unlist(lapply(strsplit(names, "_"), function(x) x[3])))
}

create_df_barplot <- function(profile, cell, force_cats=NULL){
  if (is.null(profile))
    return(NULL)
  
  day_replicate_cell <- extract_days(names(profile))
  
  if (!is.null(force_cats)){
    fake_days <- as.numeric(day_replicate_cell)
    rep <- cut(fake_days, breaks=force_cats)
    levels(rep) <- names(force_cats[-length(force_cats)])
    day_replicate_cell <- as.character(rep)
  }
  
  condition_replicate_cell <- extract_conditions(names(profile))
  
  dd <- tapply(seq_along(day_replicate_cell), day_replicate_cell, function(idx) {
    if (is.null(force_cats)){
      day=as.numeric(unique(day_replicate_cell[idx]))
    } else {
      day=unique(day_replicate_cell[idx])
    }
  
    c(MI=mean(profile[idx]),
      sd=sd(profile[idx]),
      condition=gsub(paste0(cell, "-"), "", unique(condition_replicate_cell[idx])),
      day=day)
  })
  
  df <- data.frame(do.call(rbind, dd), stringsAsFactors = F)
  df$MI <- as.numeric(df$MI)
  df$sd   <- as.numeric(df$sd)
  df
}

create_full_df_for_anova <- function(profile, cell, force_cats=NULL){
  if (is.null(profile))
    return(NULL)
  
  day_replicate_cell <- extract_days(names(profile))
  
  if (!is.null(force_cats)){
    fake_days <- as.numeric(day_replicate_cell)
    rep <- cut(fake_days, breaks=force_cats)
    levels(rep) <- names(force_cats[-length(force_cats)])
    day_replicate_cell <- as.character(rep)
  }
  
  condition_replicate_cell <- extract_conditions(names(profile))
  full_df <- data.frame(MI=profile,
                        condition=gsub("-", "_", condition_replicate_cell),
                        day=day_replicate_cell,
                        cnt = paste0(
                          gsub("-", "_", condition_replicate_cell),
                          "_",
                          day=day_replicate_cell))
  if (!is.null(force_cats)){
    full_df$day <- factor(day_replicate_cell, levels=names(force_cats[-length(force_cats)]))
  }
  full_df
}


anova_multi_test <- function(df, formula, contrasts=NULL){
  formula = as.formula(formula)
  a <- aov(formula = as.formula(formula), data = df)
  var <- all.vars(formula)[-1]
  pvalues <- TukeyHSD(a)[[var]]
  if (!is.null(contrasts)){
    pvalues <- pvalues[contrasts, , drop=F]
  }
  pvalues
  }

add_fake_day <- function(days, condition) {
  lapply(days, function(d) {
    data.frame(MI=NA, sd=NA, condition=condition, day=d, stringsAsFactors = F)
  })
}

create_mitotic_index_barplot_data <- function(cell, wt, ko, removeCondition=NULL, force_cats=NULL) {
  library(ggplot2)
  profile_wt <- wt[[cell]]
  profile_ko <- ko[[cell]]
  
  df_wt <- create_df_barplot(profile_wt, cell, force_cats = force_cats)
  df_ko <- create_df_barplot(profile_ko, cell, force_cats = force_cats)
  
  full_df_wt <- create_full_df_for_anova(profile_wt, cell, force_cats = force_cats)
  full_df_ko <- create_full_df_for_anova(profile_ko, cell, force_cats = force_cats)
  full_df <- rbind(full_df_wt, full_df_ko)
  contrasts <- sapply(levels(full_df$day), function(x) {
    paste(
      c(paste(c(cell, "ko", x), collapse = "_"),
        paste(c(cell, "wt", x), collapse = "_")),collapse="-")
  })
  anova_turkey_test <- anova_multi_test(full_df, "MI~cnt", contrasts = contrasts)
  pvalues <- anova_turkey_test[, "p adj"]
  names(pvalues) <- names(contrasts)
  sig_cmp <- names(which(pvalues<=0.05))
  
  days <- unique(df_wt$day, df_ko$day)
  
  add_wt <- do.call(rbind, add_fake_day(setdiff(days, df_wt$day), "wt"))
  add_ko <- do.call(rbind, add_fake_day(setdiff(days, df_ko$day), "ko"))
  
  df <- rbind(df_wt, df_ko, add_ko, add_wt)
  
  if(!is.null(removeCondition)) {
    df <- df[!df$condition==removeCondition, ]
  }
  if (is.null(force_cats)) {
    levels <- as.character(sort(as.numeric(unique(df$day))))
  } else {
    levels <- names(force_cats[-c(length(force_cats))])
  }
  
  df$day <- factor(df$day, levels = levels)
  fcondition <- factor(df$condition)
  df$condition <- relevel(fcondition, "wt")
  df$sig <- ""
  df$sig[df$day %in% sig_cmp] <- "**"
  df
}

create_mitotic_index_ggpubr_data <- function(cell, wt, ko, removeCondition=NULL, force_cats=NULL,
                                             condition_order=NULL) {
  library(ggplot2)
  profile_wt <- wt[[cell]]
  profile_ko <- ko[[cell]]
  
  full_df_wt <- create_full_df_for_anova(profile_wt, cell, force_cats = force_cats)
  full_df_ko <- create_full_df_for_anova(profile_ko, cell, force_cats = force_cats)
  full_df <- rbind(full_df_wt, full_df_ko)
  
  levels_wt <- unique(full_df_wt$condition)
  levels_ko <- unique(full_df_ko$condition)
  
  full_df$condition <- factor(full_df$condition, levels=c(levels_wt, levels_ko))
  
  if(!is.null(condition_order))
    full_df$condition <- factor(full_df$condition, levels=condition_order)
  
  
  full_df
}

plot_barplot <- function(df, palette) {
  ggplot(df, aes(x=day,y=MI,fill=condition)) +
    geom_bar(position="dodge",stat="identity",color="black") +
    scale_fill_manual(values=palette) +
    geom_errorbar(aes(ymin=MI-sd,ymax=MI+sd), width=0.2, position=position_dodge(0.9)) + 
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black"))
}

#### Deprecated function

mergeProfiles <- function(x, y) {
  d <- unique(c(names(x), names(y)))
  profile <- rep(NA, length(d))
  names(profile) <- d
  profileX <- profile
  profileY <- profile
  profileX[names(x)] <- x
  profileY[names(y)] <- y
  rbind(profileX, profileY)
}

plot_and_test <- function(cell, wt, ko) {
  profile_wt <- wt[[cell]]
  profile_ko <- ko[[cell]]
  day_replicate_wt <- extract_days(names(wt[[cell]]))
  day_replicate_ko <- extract_days(names(ko[[cell]]))
  days <- unique(unique(day_replicate_wt, day_replicate_ko))
  w <- lapply(days, function(day){
    day_sel_wt <- day == day_replicate_wt
    day_sel_ko <- day == day_replicate_ko
    if (any(day_sel_wt) & any(day_sel_ko))
      wilcox.test(profile_wt[day_sel_wt], profile_ko[day_sel_ko])
    
  })
  lapply(w, function(x) x$p.value)
}

compute_sd <- function(cell, profile) {
  profile_cell <- profile[[cell]]
  day_replicate_cell <- extract_days(names(profile_cell))
  tapply(seq_along(day_replicate_cell), day_replicate_cell, function(idx) {
    c(mean=mean(profile_cell[idx]), sd=sd(profile_cell[idx]))
  })
}

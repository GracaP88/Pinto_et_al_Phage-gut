#' @title Compare Distances (Beta-diversity)
#' @description Calculated between and within group distances.
#' @details Provided with \code{\link{phyloseq-class}} object and distance matrix, 
#' this wrapper with calculate between and within group distances.  
#' @param x \code{\link{phyloseq-class}} object. 
#' @param dist.matrix User provided distance matrix. Preferable of class 'dist'.
#' @param method Default= "mean". Either mean or median.  
#' @param group Groups to compare distances.
#' @param type Default= "between". Whether to compare distances 
#'             'within' or 'between' else `all` comparison are returned groups.
#' @param na.rm Default=TRUE. Will return 'NA' if there are NA values. Sanity check. 
#' @examples 
#' library(tibble)
#' library(microbiomeutilities)
#' data("zackular2014")
#' ps <- zackular2014
#' ps <- microbiome::transform(ps, "compositional")
#' dist.method= "bray"
#' ps.dist <- distance(ps, "bray")
#' beta_df <- compare_distances(ps, 
#'                              dist.matrix=ps.dist,
#'                              method="median",
#'                              group="DiseaseState",
#'                              na.rm=TRUE,
#'                              type="all")
#' head(beta_df)
#' 
#' @keywords Analysis  
#' @export
compare_distances <- function(x,
                              dist.matrix,
                              group=NULL,
                              method = "mean",
                              na.rm =TRUE,
                              type="between"){
  
  S1 <- S2 <-  comparison <- group_var_1 <- group_var_2 <- ps.dist <- value <- NULL
  if(is.null(group) | is.na(group) | isFALSE(any(group %in% colnames(meta(x))))){  
    
    stop(" Please specify correct 'group' argument")
  }
  
  #x <- ps
  ps.meta <- get_metadf(x)
  
  ps.meta$group_var <- ps.meta[,group]
  
  dist_melt <- suppressMessages(reshape2::melt(
    as.matrix(dist.matrix),
    varnames = c("S1", "S2")) %>% 
      mutate_if(is.factor, as.character) %>%
      left_join(ps.meta, by = c("S1" = "Var2")) %>% 
      filter(S1 != S2))
  
  if (type == "within") {
    
    dist_df <- suppressMessages(calculate_within_dist(dist_melt=dist_melt, 
                                                      ps.meta=ps.meta, 
                                                      method=method, 
                                                      na.rm=na.rm))
    #return(dist_df)
    
  } else if (type == "between") {
    
    dist_df <- suppressMessages(calculate_between_dist(dist_melt=dist_melt, 
                                                       ps.meta=ps.meta, 
                                                       method=method, 
                                                       na.rm=na.rm))
    #return(dist_df)
    
  } else if (type=="all") {
    
    dist_df1 <- suppressMessages(calculate_within_dist(dist_melt=dist_melt, 
                                                       ps.meta=ps.meta, 
                                                       method=method, 
                                                       na.rm=na.rm))
    
    dist_df2 <- suppressMessages(calculate_between_dist(dist_melt=dist_melt, 
                                                        ps.meta=ps.meta, 
                                                        method=method, 
                                                        na.rm=na.rm))
    
    colnames(dist_df1)[3] <- paste0( method,".distance")
    colnames(dist_df2)[3] <- paste0( method,".distance")
    
    dist_df <- bind_rows(dist_df1,dist_df2)
    
    #return(dist_all)
    
  }
  return(dist_df)
}

#' @keywords utilities 
get_metadf<- function(x){
  df <- meta(x) %>% 
    rownames_to_column("Var2")
}


#' @keywords utilities 
calculate_within_dist <- function(dist_melt, 
                                  ps.meta, 
                                  method, 
                                  na.rm){
  
  #group_var_1 <- group_var_2 <- S2 <- comparison <- NULL
  
  if (method=="mean"){
    
    dist_within <- dist_melt %>% 
      left_join(ps.meta,
                by = c("S2" = "Var2"), 
                suffix = c("_1", "_2")) %>% 
      filter(group_var_1 == group_var_2) %>% 
      mutate(comparison= paste0(group_var_1, " vs ", group_var_2)) %>% 
      group_by(S2, comparison) %>% 
      summarise(within.mean.dist=mean(value, na.rm = na.rm)) 
    
    return(dist_within)
    
  } else if (method=="median") {
    
    dist_within <- dist_melt %>% 
      left_join(ps.meta,
                by = c("S2" = "Var2"), 
                suffix = c("_1", "_2")) %>% 
      filter(group_var_1 == group_var_2) %>% 
      mutate(comparison= paste0(group_var_1, " vs ", group_var_2)) %>% 
      group_by(S2, comparison) %>%
      summarise(within.median.dist=median(value, na.rm = na.rm)) 
    
    return(dist_within)
  } 
  
  
  
}



#' @keywords utilities 
calculate_between_dist <- function(dist_melt, 
                                   ps.meta, 
                                   method, 
                                   na.rm){
  
  #group_var_1 <- group_var_2 <- S2 <- comparison <- value <- NULL
  
  if (method=="mean") {
    dist_between <- dist_melt %>% 
      left_join(ps.meta,
                by = c("S2" = "Var2"), 
                suffix = c("_1", "_2")) %>% 
      filter(group_var_1 != group_var_2) %>% 
      mutate(comparison= paste0(group_var_1, " vs ", group_var_2)) %>% 
      group_by(S2, comparison) %>% 
      summarise(between.mean.dist=mean(value, na.rm = na.rm))
    
    return(dist_between)
    
  } else if (method=="median") {
    
    dist_between <- dist_melt %>% 
      left_join(ps.meta,
                by = c("S2" = "Var2"), 
                suffix = c("_1", "_2")) %>% 
      filter(group_var_1 != group_var_2) %>% 
      mutate(comparison= paste0(group_var_1, " vs ", group_var_2)) %>% 
      group_by(S2, comparison) %>% 
      summarise(between.median.dist=median(value, na.rm = na.rm))
    
    return(dist_between)
  }
  
  
  
}




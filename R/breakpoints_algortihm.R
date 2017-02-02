# Breakpoints algorithm ---------------------------------------------------

library(magrittr)
library(dplyr)
library(strucchange)
library(tsoutliers)

#' Title
#'
#' @param d the data frame containing the target variable
#' @param x The column in d which represents x variable in the time-series 
#' @param target_var The column in d that will be the y variable in the time-series 
#' @param freq Frequency of the time-series - set to 12 by default for monthly data
#' @param min_seg The minimum number of of points in a segment in case floor(n * h) < 1 - defaults to 4
#' @param h minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment.
#' @param debug_info Whether to print the error message from the outliers package
#' @param log_file File location for the outliers package to throw output
#' @param cval A value for the t-test statistics for accepting outliers in the tsoutliers package
#'
#' @return a list of objects - a breakpoint object and an outlier object
#' @export
#'
#' @examples
breakpoints_algorithm <- function(d, 
                                  x = "Date",  
                                  target_var,  
                                  freq = 12, 
                                  min_seg = 4, 
                                  h = 0.15,
                                  use_tsoutliers = F, # not recommended for performance reasons - very slow
                                  debug_info = F,
                                  log_file = NULL,
                                  cval = 10,
                                  nstdev = 3)
{
  outlier_obj <- NA
  out_liers_idxs <- NA
  breakpoints_on_outliers <- NA
  bp_idxs <- NA
  y_var <- d[[target_var]]
  x_var <- d[[x]]
  
  
  
  if(use_tsoutliers) {
    # if there are any breakpoints run an outlier detection algorithm 
    tryCatch(
      outlier_obj <- tsoutliers::tso(ts_obj, logfile = log_file, cval = cval),
      error = function(c) {
        if(debug_info) {
          print(c)
        }
        outlier_obj <- NA
      }
    )
    
    if(!is.na(outlier_obj) && nrow(outlier_obj$outliers) > 0) {
      ao_outliers <- outlier_obj$outliers %>% filter(type == "AO")
      if(nrow(ao_outliers) > 0) {
        filter_cond <- !(1:nrow(d) %in% ao_outliers$ind)
        filtered_d <- d[filter_cond,]
        out_liers_idxs <- ao_outliers$ind
      }
    }
    
  } else {
    med <- median(y_var)
    dev <- sd(y_var)
    out_liers <- abs(y_var - med) > nstdev * dev
    out_liers_idxs <- which(out_liers)
    filtered_d <- d %>% filter_(lazyeval::interp(~ abs(var - med) < nstdev * dev, var = as.name(target_var) ))
  }
  
  if(nrow(filtered_d) == 0) stop("The time-series contains zero points after filtering for outliers.")
  
  ts_obj <- ts(filtered_d[[target_var]])
  
  # This part determines the actual minimal segment size.
  # By default this will be floor(15 % of the number of points in the data)
  
  n <- nrow(filtered_d) # number of points in dataset
  k <- 1       # number of regressors
  if (is.null(h)) 
    h <- k + 1 # if h is set to null make it one larger than the number of breakpoints
  if (h < 1) 
    h <- floor(n * h) # set minimal number of points in a segment size to h * n 
  if (h <= k) {
    if(!min_seg > 1) {
      return(NA)
      warning("Minimum segment size must be greater than the number of regressors (1).")
    }
    h <- min_seg
  }
  if (h > floor(n/2)) {
    return(NA)
    warning("Minimum segment size must be smaller than half the number of observations.")
  }
  
  bp <- breakpoints(ts_obj ~ 1, h = h)
  bp_idxs <- bp$breakpoints
  if(any(!is.na(out_liers_idxs)) & !any(is.na(bp_idxs))) {
    for(i in out_liers_idxs) {
      shifted_idxs <- which(bp_idxs >= i)
      bp_idxs[shifted_idxs] <- bp_idxs[shifted_idxs] + 1
    }
  }
  
  return(list(bp_obj = bp, # object of class breakpoint
              bp_idxs = bp_idxs, # the indexes of the breakpoints from 1 to n time points
              outlier_obj = outlier_obj, # An outlier object if the tsoutliers package is used, otherwise NA
              out_liers_idxs = out_liers_idxs)) # The list of indexes indicating outliers
}



#' Returns a factor vector of the segments to which a point in the time series belongs
#'
#' @param d # the data frame from which the time series was computed
#' @param obj # an object of class breakpoint
#' @param bp_idxs # The indexes of the breakpoints
#' @param out_liers_idxs # The indexes of outliers if there were any in the data
#' @param n # the number of points in the time series
#' @param labels # labels to be used for the segments
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
breakfactors_with_outliers <- function (d,
                                        obj, 
                                        bp_idxs, 
                                        n = NULL, 
                                        out_liers_idxs = NA, 
                                        labels = NULL, ...) 
{
  if(!is.null(n)) n = nrow(d)
  if(all(is.na(out_liers_idxs))) {
    break_factors <- breakfactor(obj)
  } else {
    nbreaks <- length(bp_idxs)
    fac <- rep(1:(nbreaks + 1), c(bp_idxs[1], diff(c(bp_idxs, n))))
    fac[out_liers_idxs] <- NA
    if (is.null(labels)) 
      labels <- paste("segment", 1:(nbreaks + 1), sep = "")
    break_factors <- factor(fac, labels = labels, ...)
  }
  return(break_factors)
}

# TODO: generalize for arbitrary column instead of break_factors column
#' Title
#'
#' @param d # The data frame over which models will be computed - must contain a break_factors column
#' @param f 
#'
#' @return
#' @export
#'
#' @examples


compute_segment_models <- function(d, f) {
  if(!"break_factors" %in% colnames(d)) stop("Break factor column doesn't exist")
  # Break up d by break_factors, then fit the specified model to each piece and
  # return a list
  d %<>% filter(!is.na(break_factors))
  models <- dlply(d, "break_factors", function(df) lm(f, data = df))
  models
}



#' Computes the characterstics of each linear segment in a time series with breakpoints - 
#'
#' @param models - a list of linear models for each segment
#' @param angle_in_deg - whether to display the coefficient of the slop in angles or points per period
#'
#' @return
#' @export
#'
#' @examples
extract_segment_coefs <- function(models, angle_in_deg = F) {
  # get the prediction at the last point of each segment
  last_points <- ldply(models, function(x) {last_point <- last_elem(x$fitted.values); return(last_point) })
  colnames(last_points)[2] <- "last_point"
  
  # get the coefficients of each segment
  segment_models <- ldply(models, coef)
  segment_models %<>% mutate(slope = x_seg_idx)
  if(angle_in_deg) {
    segment_models %<>% mutate(angles = atan(x_seg_idx) * 180 / pi)
  } else {
    segment_models %<>% mutate(angles = x_seg_idx)
  }
  suppressMessages({segment_models %<>% left_join(last_points)})
  segment_models
}




#' Title
#'
#' @param segment_models - a dataframe contianing the information for the coefficients of the linear models in each segment of the time series
#'
#' @return - a data frame with absolute and relative differences in the 
#' @export
#'
#' @examples
compute_breakpoint_changes <- function(segment_models) {
  
  if(nrow(segment_models) < 2) return(NA) # if only one segment is present
  
  breakpoints_shifts_abs <- numeric(0)
  breakpoints_shifts_rel <- numeric(0)
  breakpoints_angle_bef <- numeric(0)
  breakpoints_angle_aft <- numeric(0)
  
  for(i in 1:(nrow(segment_models) -1)) {
    breakpoints_shifts_abs %<>% append(., segment_models$`(Intercept)`[i+1] - segment_models$last_point[i])
    breakpoints_shifts_rel %<>% append(., (segment_models$`(Intercept)`[i+1] - segment_models$last_point[i]) / segment_models$last_point[i] * 100)
    breakpoints_angle_bef  %<>% append(., segment_models$angles[i])
    breakpoints_angle_aft  %<>% append(., segment_models$angles[i+1])
  }
  
  angle_diff_abs <- diff(segment_models$angles)
  angle_diff_rel <- diff(segment_models$angles) / abs(segment_models$angles[1 : (nrow(segment_models) - 1)]) * 100
  
  breakpoint_changes <- data.frame(breakpoints_shifts_abs, 
                                   breakpoints_shifts_rel, 
                                   angle_diff_abs, 
                                   angle_diff_rel, 
                                   breakpoints_angle_bef,
                                   breakpoints_angle_aft)
  
  breakpoint_changes
}


#' Formats a label to print the shift in intercept over every breakpoint
#'
#' @param breakpoint_changes - a dataframe containing the information abouth the changes in the linear segments, before and after
#' @param digits - number of digits for displaying
#'
#' @return
#' @export
#'
#' @examples
shift_diff_format <- function(breakpoint_changes, digits = 2) {
  x <- breakpoint_changes
  qchange_desc <- character(0)
  for(i in 1:nrow(x)) {
    qchange_desc %<>% c(paste("Diff =", round(x$breakpoints_shifts_abs[i], digits = digits)))
  }
  return(qchange_desc)
}


#' Returns the confidence intervals of the breakpoints in case there were outliers
#'
#' @param bp_obj 
#' @param out_liers_idxs 
#'
#' @return
#' @export
#'
#' @examples
compute_bp_ci_with_outliers <- function(bp_obj, out_liers_idxs) {
  bp_ci <- confint(bp_obj)$confint
  
  if(any(!is.na(out_liers_idxs))) {
    for(i in out_liers_idxs) {
      shifted_idxs <- bp_ci >= i
      bp_ci[shifted_idxs] <- bp_ci[shifted_idxs] + 1
    }
  }
  bp_ci
}



#' Title
#'
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
compute_bp_obj <- function(d, target_var) {
  
  res <- breakpoints_algorithm2(d)
  if(all(is.na(res)) | all(is.na(res$bp_idxs))) return(NA)
  
  bp_obj          <- res$bp_obj
  bp_idxs         <- res$bp_idxs
  out_liers_idxs  <- res$out_liers_idxs
  res$bp_ci       <- compute_bp_ci_with_outliers(bp_obj, out_liers_idxs)
  
  # enrich d with breakpoints(is.na(res) | all(is.na(res$bp_idxs))
  d$break_factors <- breakfactors_with_outliers(d,
                                                bp_obj,
                                                bp_idxs,
                                                out_liers_idxs,
                                                n = n)
  
  # compute models on each segment
  d %<>% group_by(break_factors) %>% mutate(x_seg_idx = 1:n()) %>% ungroup()
  f <- paste(target_var, "~ x_seg_idx")
  models <- compute_segment_models(d, f)
  
  # compute coeffiecient segments on each chunk
  segment_models <- extract_segment_coefs(models)
  
  # compute the changes occuring at the breakpoints
  breakpoint_changes <- compute_breakpoint_changes(segment_models)
  
  res$models <- models
  res$segment_models <- segment_models
  res$breakpoint_changes <- breakpoint_changes
  res$d <- d
  
  res
}






#' Generates a ggplot2 plot of a time series with breakpoints
#'
#' @param bp_out_list - a list containing an object of class breakpoint, the indexes of the breakpoints & the indexes of outliers locations
#' @param d - the dataframe used to compute the breakpoints object
#' @param titl - title for the plot - defaults to empty title
#' @param x - the x variable in the data vrame - defaults to Date
#' @param y - the y variable in the data vrame - defaults to Value
#' @param xlabel - the actual label to display on the x axis
#' @param ylabel - the actual label to display on the x axis
#' @param y_scale_factor - A scaling factor for extending the y axis - defaults to 5%
#' @param coef_digits - the number of digits for the segment coefficients to be displayed
#' @param models - a list of linear models for each segment
#'
#' @return
#' @export
#'
#' @examples
bp_ggplot2 <- function(bp_out_list,
                       # d,
                       x, 
                       y, 
                       titl = "", 
                       xlabel = NULL, 
                       ylabel = NULL,
                       y_scale_factor = 0.05,
                       coef_digits = 3,
                       models = NULL) {
  
  bp_obj          <- bp_out_list$bp_obj
  bp_idxs         <- bp_out_list$bp_idxs
  out_liers_idxs  <- bp_out_list$out_liers_idxs
  bp_ci           <- bp_out_list$bp_ci # compute_bp_ci_with_outliers(bp_obj, out_liers_idxs)
  models          <- bp_out_list$models
  segment_models  <- bp_out_list$segment_models
  breakpoint_changes <- bp_out_list$breakpoint_changes
  d <- bp_out_list$d
  if(is.null(xlabel)) xlabel <- x
  if(is.null(ylabel)) ylabel <- y
  
  ts_x <- d[[x]]
  myplot = function(df, x_string, y_string) {
    ggplot(df, aes_string(x = x_string, y = y_string))
  }
  # print(titl)
  n <- nrow(d)
  
  # if(!"break_factors" %in% colnames(d)) {
  #   d$break_factors <- breakfactors_with_outliers(d,
  #                                                 bp_obj,
  #                                                 bp_idxs,
  #                                                 out_liers_idxs,
  #                                                 n = n)
  # }
  # 
  # if(is.null(models)) {
  #   d %<>% group_by(break_factors) %>% mutate(x_seg_idx = 1:n()) %>% ungroup()
  #   f <- paste(y, "~ x_seg_idx")
  #   models <- compute_segment_models(d, f)
  # }
  # 
  # segment_models <- extract_segment_coefs(models)
  # add a label
  segment_models %<>% mutate(lab = paste("Angle = ", round(x_seg_idx, coef_digits), " per period", sep = ""))
  
  # breakpoint_changes <- compute_breakpoint_changes(segment_models)
  # add poisition locations and labels
  breakpoint_changes$x_pos_bp <- d$Date[bp_idxs]
  breakpoint_changes$y_pos_bp <- rep(max(d[[y]]) * 1.05 , times = length(bp_idxs))
  breakpoint_changes$lab <- shift_diff_format(breakpoint_changes)
  
  
  # This creates the initial plot with linear segments
  p <- myplot(d, x_string = x, y_string = y) + 
    geom_point(shape = 1) + 
    geom_line() + 
    scale_x_date(breaks =  d[[x]][seq(n, 1, -3)], 
                 labels = scales::date_format("%y %b")) +
    geom_smooth(aes(group = break_factors), 
                method ="lm", 
                data = subset(d, !is.na(break_factors))) +
    geom_point(color = 'red', # the outliers have na in break_factors
               size = 3, 
               data = subset(d, is.na(break_factors)))
  
  
  # This part formats the positioning of the breakpoints vertical dashed lines and confidence intervals
  ggb <- ggplot_build(p)
  nbp <- nrow(bp_ci)
  yrange <- ggb$panel$ranges[[1]]$y.range
  at  <- yrange
  at  <- diff(at)/1.08 * 0.02 + at[1]
  at2  <- at * 0.97
  if(length(at) < nbp) at <- rep(at, length.out = nbp)
  bp_ci[,1][bp_ci[,1] <= 0] <- 1
  bp_ci[,3][bp_ci[,3] > n] <- n
  segments <- data.frame(x1 = ts_x[bp_ci[,1]], x2 = ts_x[bp_ci[,3]], y1 = at, y2 = at)
  
  p <- p + geom_segment(aes(x = x1, 
                            y = y1, 
                            xend = x2, 
                            yend = y2), 
                        data = segments, 
                        colour = 'red')
  
  p <- p + geom_vline(xintercept = ggb$data[[1]]$x[bp_ci[,2]], 
                      linetype = 'longdash')
  
  
  
  # This part formats the positioning of the segment coeficients labels
  # Create a dataframe for printing the segment coefficients
  segment_coefs <- subset(d, !is.na(break_factors)) %>% group_by(break_factors) %>% summarise(x2 = floor(n()/2), 
                                                                                              med_x = Date[x2], 
                                                                                              med_y = median(Value_Share_Shampoo), 
                                                                                              max_y = max(Value_Share_Shampoo),
                                                                                              min_y = min(Value_Share_Shampoo))
  med <- subset(d, !is.na(break_factors)) %>% .[[y]] %>% median()
  segment_coefs %<>% mutate(med_y = if_else(med_y < med, max_y + 0.1 * max_y, min_y - 0.1 * min_y))
  if(nrow(segment_coefs) <= 3) {
    segment_coefs$med_y <- at2
  }
  segment_coefs$lab <- segment_models$lab
  
  
  # Add segment coefficients labels
  p <- p + ggplot2::annotate("text", 
                             x = segment_coefs$med_x, 
                             y = segment_coefs$med_y, 
                             label = segment_coefs$lab,
                             fontface = 'bold',
                             color = 'Black',
                             size = 3)
  
  # Add breakpoints labels
  p <- p + ggplot2::annotate("text", 
                             x = breakpoint_changes$x_pos_bp, 
                             y = breakpoint_changes$y_pos_bp, 
                             label = breakpoint_changes$lab,
                             fontface = 'bold',
                             color = 'Black')
  # Extend the y axis by a certain factor
  p <- p + scale_y_continuous(limits = c(yrange[1]- yrange[1] * y_scale_factor, 
                                         yrange[2] * (1 + y_scale_factor)))
  
  # Set title of plot
  p <- p + ggtitle(titl) + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    theme(axis.text.y = element_text(face = "bold", size = 13))
  
  return(p)
}





qualitatitve_change <- function(breakpoint_changes, digits = 2, stable_tresh = 0.3) {
  x <- breakpoint_changes
  qchange_desc <- character(0)
  for(i in 1:nrow(x)) {
    # Enter the region of shifts
    if(abs(x$breakpoints_shifts_rel[i]) > 1) {
      # stable to stable
      if(abs(x$breakpoints_angle_bef[i]) < stable_tresh && abs(x$breakpoints_angle_aft[i]) < stable_tresh) {
        # qchange_desc %<>% c(paste("Stable to stable, <br> I diff.: ", round(x$breakpoints_shifts_abs[i],digits = digits)))
        qchange_desc %<>% c(paste("Stable to <br> stable"))
        
      } else if(abs(x$breakpoints_angle_bef[i]) < stable_tresh && abs(x$breakpoints_angle_aft[i]) > stable_tresh) { # stable to trend
        qchange_desc %<>% c(paste("Stable to <br> trend"))
        # qchange_desc %<>% c(paste("Stable to trend, <br> I diff.: ", round(x$breakpoints_shifts_abs[i],digits = digits), 
        # "<br> Ang diff: ", round(x$angle_diff_abs[i], digits = digits)))
      } else if(abs(x$breakpoints_angle_bef[i]) > stable_tresh && abs(x$breakpoints_angle_aft[i]) < stable_tresh) { # Trend to stable
        qchange_desc %<>% c(paste("Trend to <br> stable"))
        # qchange_desc %<>% c(paste("Trend to stable, <br> I diff. ", round(x$breakpoints_shifts_abs[i],digits = digits),
        # "<br> Ang diff: ", round(x$angle_diff_abs[i], digits = digits)))
      } else if(abs(x$breakpoints_angle_bef[i]) > stable_tresh && abs(x$breakpoints_angle_aft[i]) > stable_tresh) { # Trend to Trend
        # qchange_desc %<>% c(paste("Trend to trend, <br> I diff. ", round(x$breakpoints_shifts_abs[i],digits = digits), 
        # "<br> Ang diff: ", round(x$angle_diff_abs[i], digits = digits)))
        qchange_desc %<>% c(paste("Trend to <br> trend"))
      } 
    } else { # enter the region of acceleration/deceleration
      
    }
  }
  if(length(qchange_desc) < nrow(x)){
    empty_str <- rep("", times = nrow(x) - length(qchange_desc))
    qchange_desc %<>% c(empty_str)
  }
  return(qchange_desc)
  
}

e1_test <- function(p4, mid_line, ucl, lcl) {
  test_fail <- F
  abs_diff_to_ucl <- abs(p4 - ucl)
  abs_diff_to_mid <- abs(p4 - mid_line)
  abs_diff_to_lcl <- abs(p4 - lcl)
  if(abs_diff_to_mid[4] < min(abs_diff_to_lcl[4], abs_diff_to_ucl[4])) return(test_fail)
  dp4 <- diff(p4)
  test_fail <- all(dp4 > 0) | all(dp4 < 0) | all(p4 > mid_line) | all(p4 < mid_line)
  return(test_fail)
}

e1_test_part1 <- function(p4) {
  test_fail <- F
  dp4 <- diff(p4)
  test_fail <- all(dp4 > 0) | all(dp4 < 0)
  return(test_fail)
}

e2_test <- function(p3, mid_line, ucl, lcl) {
  test_fail <- F
  abs_diff_to_ucl <- abs(p3 - ucl)
  abs_diff_to_mid <- abs(p3 - mid_line)
  abs_diff_to_lcl <- abs(p3 - lcl)
  if(all(p3 > mid_line)) { # if all above midline test distance to ucl
    if(sum(abs_diff_to_ucl > abs_diff_to_mid) >= 2) test_fail = T
  } else if (all(p3 < mid_line)) {
    if(sum(abs_diff_to_lcl > abs_diff_to_mid) >= 2) test_fail = T
  }
  return(test_fail)
}


e2_test_v1 <- function(p3, mid_line, ucl, lcl) {
  test_fail <- F
  # if points 1 and 3 are on the same side of the midline but the mid point is on the oppoosite side
  if(sign(p3[1] - mid_line[1]) == sign(p3[3] - mid_line[3]) && sign(p3[1] - mid_line[1]) != sign(p3[2] - mid_line[2])) return(test_fail)
  
  up_lim <- 0.8 * (ucl - mid_line)
  lw_lim <- 0.8 * (mid_line - lcl)
  abs_diff_to_ucl <- abs(p3 - ucl)
  abs_diff_to_mid <- abs(p3 - mid_line)
  abs_diff_to_lcl <- abs(p3 - lcl)
  # if the last point is very close to the midpoint then things got back to a reasonable value and the test fails
  if(abs_diff_to_mid[3] < min(abs(up_lim[3]), abs(lw_lim[3]))) return(test_fail)
  # else if in addition to the last point one of the other two points are running away from mid_line
  if(abs_diff_to_mid[2] > min(abs(up_lim[2]), abs(lw_lim[2])) |
     abs_diff_to_mid[1] > min(abs(up_lim[1]), abs(lw_lim[1]))) test_fail = T
  return(test_fail)
}


early_warning_tests <- function(obj, r2cond = 0.9) {
  
  early_warning <- list()
  
  
  
  if(!is.null(obj$models)) {
    last_segment_model <- last_elem(obj$models)[[1]]
    r2 <- summary(last_segment_model)$r.squared
    if(r2 > r2cond) {
      obj$early_warning <- early_warning
      return(obj) # no early warning
    }
    
    df_points <- stats::predict(last_segment_model, interval = "confidence") %>% as.data.frame()
    df_points$points <- last_segment_model$model[[1]]
    
    p4 <- tail(df_points, 4)
    e1_test_res <- e1_test(p4 = p4$points, 
                           mid_line = p4$fit,
                           ucl = p4$upr,
                           lcl = p4$lwr)
    
    if(e1_test_res) {
      early_warning$test_type <- "e1"
      early_warning$test_res <- T
      early_warning$early_warning_idx <- 4
      obj$early_warning <- early_warning
      return(obj)
    }
    
    
    p3 <- tail(df_points, 3)
    e2_test_res <- e2_test_v1(p3 = p3$points, 
                              mid_line = p3$fit,
                              ucl = p3$upr,
                              lcl = p3$lwr)
    if(e2_test_res) {
      early_warning$test_type <- "e2"
      early_warning$test_res <- T
      early_warning$early_warning_idx <- 3
      obj$early_warning <- early_warning
      return(obj)
    }
    
    
    
    obj$early_warning <- early_warning
    return(obj)
  }
}


transform_bp_df <- function(obj) {
  d <- obj$bp_obj$d
  print(unique(d$id))
  d$breakpoint_indicator <- FALSE
  d$breakpoint_indicator[obj$bp_obj$bp_idxs] <- T
  d$outlier_indicator <- FALSE
  if(length(obj$bp_obj$out_liers_idxs)) {
    d$outlier_indicator[obj$bp_obj$out_liers_idxs] <- T
  }
  
  d$ci_indicator <- FALSE
  for(ci_idx in 1:nrow(obj$bp_obj$bp_ci)) {
    ci_row <- obj$bp_obj$bp_ci[ci_idx,]
    
    ci_idxs <- seq(from = ci_row[1], to = ifelse(ci_row[3] < nrow(d),ci_row[3], nrow(d)))
    d$ci_indicator[ci_idxs] <- TRUE
  }
  
  d %<>% left_join(obj$bp_obj$segment_models[, c("break_factors", "(Intercept)", "slope")], by = "break_factors")
  d$breakpoint_shifts <- NA
  d$breakpoint_shifts[d$breakpoint_indicator] <- obj$bp_obj$breakpoint_changes$breakpoints_shifts_abs
  
  d$early_warning_indicator <- FALSE
  if(length(obj$bp_obj$early_warning)) {
    ew_idx <- nrow(d) - obj$bp_obj$early_warning$early_warning_idx + 1
    d$early_warning_indicator[ew_idx] <- TRUE
  }
  
  return(d)
  
}



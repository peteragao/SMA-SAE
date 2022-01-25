library(dplyr)
library(tidyr)

get_direct <- function(formula, by, sample_des, CI = .90) {
  res <- svyby(formula, by, design = sample_des, svymean, na.rm = T)
  out_dat <- data.frame(
    region = as.vector(res[[all.vars(by)[1]]]),
    est = as.vector(model.matrix(formula, res)[, 2]),
    var = res$se ^ 2
  ) %>% mutate(
    median = est,
    lower = est + qnorm((1-CI)/2) * res$se,
    upper = est + qnorm(1 - (1-CI)/2) * res$se,
    method = "HT"
  )
  return(out_dat)
}

get_greg <- function(working_fit, formula, by,
                     pop_dat, sample_des, CI = .90) {
  pop_unit_ests <- as.vector(predict(working_fit, pop_dat,
                                     type = "response")) 
  area_ests <-
    aggregate(pop_unit_ests, 
              list(region = as.vector(pop_dat[[all.vars(by)[1]]])),
              mean)
  colnames(area_ests)[2] <- "working_est"
  sample_des$variables$res <- 
    sample_des$variables[[all.vars(formula)[1]]] -
    as.vector(predict(working_fit, sample_des$variables, type = "response")) 
  sample_des$variables$region <- 
    as.vector(sample_des$variables[[all.vars(by)[1]]])
  res_ht <- svyby(~res, ~region, sample_des, svymean)
  out_dat <- left_join(area_ests, res_ht, by = "region")
  out_dat$est = out_dat$working_est + out_dat$res
  out_dat$median <- out_dat$est
  out_dat$var = out_dat$se ^ 2
  out_dat$method = "GREG"
  out_dat$lower = out_dat$est + qnorm((1-CI)/2) * out_dat$se
  out_dat$upper = out_dat$est + qnorm(1 - (1-CI)/2) * out_dat$se
  out_dat <- dplyr::select(out_dat, region, median, est, var, lower, upper, method)
  return(out_dat)
}

get_greg <- function(working_fit, formula, by, 
                     pop_dat, sample_des, pop = NULL, CI = .90) {
  pop_unit_ests <- as.vector(predict(working_fit, pop_dat,
                                     type = "response")) 
  if (is.null(pop)) {
    area_ests <-
      aggregate(pop_unit_ests, 
                list(region = as.vector(pop_dat[[all.vars(by)[1]]])),
                mean)
  } else {
    pop_dat$pop <- pop_dat[[all.vars(pop)[1]]]
    wt_area_ests <-
      aggregate(pop_unit_ests * pop_dat$pop, 
                list(region = as.vector(pop_dat[[all.vars(by)[1]]])),
                sum)
    area_pops <-
      aggregate(pop_dat$pop, 
                list(region = as.vector(pop_dat[[all.vars(by)[1]]])),
                sum)
    area_ests <- wt_area_ests
    area_ests[, 2] <- wt_area_ests[, 2] / area_pops[, 2]
  }
  colnames(area_ests)[2] <- "working_est"

  sample_des$variables$res <- 
    sample_des$variables[[all.vars(formula)[1]]] -
    as.vector(predict(working_fit, sample_des$variables, type = "response")) 
  sample_des$variables$region <- 
    as.vector(sample_des$variables[[all.vars(by)[1]]])
  res_ht <- svyby(~res, ~region, sample_des, svymean)
  out_dat <- left_join(area_ests, res_ht, by = "region")
  out_dat$est = out_dat$working_est + out_dat$res
  out_dat$median <- out_dat$est
  out_dat$var = out_dat$se ^ 2
  out_dat$method = "GREG"
  out_dat$lower = out_dat$est + qnorm((1-CI)/2) * out_dat$se
  out_dat$upper = out_dat$est + qnorm(1 - (1-CI)/2) * out_dat$se
  out_dat <- dplyr::select(out_dat, region, est, median, var, lower, upper, method)
  return(out_dat)
}
get_iid_sdir <- function(direct_est,
                         pc_u = 5,
                         pc_alpha = 0.01, 
                         CI = .90) {
  hyperpc_iid_int <- list(prec = list(prior = "pc.prec",
                                      param = c(pc_u , pc_alpha)))
  sd_dat <- direct_est %>%
    mutate(est = ifelse(est != 0 & est != 1 & var > 1e-5, est, NA),
           prec = 1 / var)
  sd_fit <-
    INLA::inla(est ~ f(region, model = "iid", 
                             hyper = hyperpc_iid_int),
               family = "gaussian", data = sd_dat, 
               scale = sd_dat$prec,
               control.family = 
                 list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  sd_fit_sample <- inla.posterior.sample(n = 1000, sd_fit,
                                         list(Predictor = 1:nrow(sd_dat)))
  sd_est_mat <- do.call(cbind, lapply(sd_fit_sample,
                                      function(x) x$latent[1:nrow(sd_dat)]))
  out_dat <- data.frame(region = direct_est$region,
                        est = rowMeans(sd_est_mat),
                        median = apply(sd_est_mat, 1,
                                       function(x) median(x, na.rm = T)),
                        var = apply(sd_est_mat, 1, var),
                        lower = apply(sd_est_mat, 1,
                                      function(x) quantile(x, (1-CI)/2)),
                        upper = apply(sd_est_mat, 1,
                                      function(x) quantile(x, 1-(1-CI)/2)),
                        method = paste0("iidS", direct_est$method[1]))
  return(out_dat)
}
get_iid_sdir_logit <- function(direct_est,
                               pc_u = 5,
                               pc_alpha = 0.01, 
                               CI = .90) {
  hyperpc_iid_int <- list(prec = list(prior = "pc.prec",
                                      param = c(pc_u , pc_alpha)))
  sd_dat <- direct_est %>%
    mutate(est = ifelse(est != 0 & est != 1 & var > 1e-5, est, NA)) %>%
    mutate(logit_est = SUMMER::logit(est),
           logit_prec = (est^2*(1-est)^2) / var)
  sd_fit <-
    INLA::inla(logit_est ~ f(region, model = "iid", 
                             hyper = hyperpc_iid_int),
               family = "gaussian", data = sd_dat, 
               scale = sd_dat$logit_prec,
               control.family = 
                 list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  sd_fit_sample <- inla.posterior.sample(n = 1000, sd_fit,
                                         list(Predictor = 1:nrow(sd_dat)))
  sd_est_mat <- do.call(cbind, lapply(sd_fit_sample,
                                      function(x) x$latent[1:nrow(sd_dat)]))
  sd_est_mat <- SUMMER::expit(sd_est_mat)
  out_dat <- data.frame(region = direct_est$region,
                        est = rowMeans(sd_est_mat),
                        median = apply(sd_est_mat, 1,
                                       function(x) median(x, na.rm = T)),
                        var = apply(sd_est_mat, 1, var),
                        lower = apply(sd_est_mat, 1,
                                      function(x) quantile(x, (1-CI)/2)),
                        upper = apply(sd_est_mat, 1,
                                      function(x) quantile(x, 1-(1-CI)/2)),
                        method = paste0("iidSLogit", direct_est$method[1]))
  return(out_dat)
}
get_bym2_sdir <- function(direct_est, 
                                adj_mat,
                                pc_u = 5,
                                pc_alpha = 0.01, 
                                pc_u_phi = 0.5,
                                pc_alpha_phi = 2/3,
                                CI = .90) {
  hyperpc_bym_int <- list(
    prec = list(prior = "pc.prec", param = c(pc_u , pc_alpha)),  
    phi = list(prior = 'pc', param = c(pc_u_phi , pc_alpha_phi))
  )
  sd_dat <- direct_est %>%
    mutate(est = ifelse(est != 0 & est != 1 & var > 1e-5, est, NA)) %>%
    mutate(prec = 1 / var,
           region = match(region, rownames(adj_mat)))
  sd_fit <-
    INLA::inla(est ~ f(region, model = "bym2", 
                       graph = adj_mat, 
                       hyper = hyperpc_bym_int, 
                       scale.model = TRUE),
               family = "gaussian", data = sd_dat, 
               scale = sd_dat$prec,
               control.family = 
                 list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  
  sd_fit_sample <-
    inla.posterior.sample(n = 1000, sd_fit,
                          list(region = 1:nrow(adj_mat), "(Intercept)" = 1))
  sd_est_mat <-
    do.call(cbind, lapply(sd_fit_sample,
                          function(x) x$latent[1:nrow(adj_mat)] +
                            x$latent[nrow(adj_mat) + 1]))
  out_dat <- data.frame(region = rownames(adj_mat),
                        est = rowMeans(sd_est_mat),
                        median = apply(sd_est_mat, 1,
                                       function(x) median(x, na.rm = T)),
                        var = apply(sd_est_mat, 1, var),
                        lower = apply(sd_est_mat, 1,
                                      function(x) quantile(x, (1-CI)/2)),
                        upper = apply(sd_est_mat, 1,
                                      function(x) quantile(x, 1-(1-CI)/2)),
                        method = paste0("bymS", direct_est$method[1]))
  return(out_dat)
}
get_bym2_sdir_logit <- function(direct_est, 
                                adj_mat,
                                pc_u = 5,
                                pc_alpha = 0.01, 
                                pc_u_phi = 0.5,
                                pc_alpha_phi = 2/3,
                                CI = .90) {
  hyperpc_bym_int <- list(
    prec = list(prior = "pc.prec", param = c(pc_u , pc_alpha)),  
    phi = list(prior = 'pc', param = c(pc_u_phi , pc_alpha_phi))
  )
  sd_dat <- direct_est %>%
    mutate(est = ifelse(est != 0 & est != 1 &var > 1e-5, est, NA)) %>%
    mutate(logit_est = SUMMER::logit(est),
           logit_prec = (est^2*(1-est)^2) / var,
           region = match(region, rownames(adj_mat)))
  sd_fit <-
    INLA::inla(logit_est ~ f(region, model = "bym2", 
                             graph = adj_mat, 
                             hyper = hyperpc_bym_int,
                             scale.model = TRUE),
               family = "gaussian", data = sd_dat, 
               scale = sd_dat$logit_prec,
               control.family = 
                 list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  
  sd_fit_sample <-
    inla.posterior.sample(n = 1000, sd_fit,
                          list(region = 1:nrow(adj_mat), "(Intercept)" = 1))
  sd_est_mat <-
    do.call(cbind, lapply(sd_fit_sample,
                          function(x) x$latent[1:nrow(adj_mat)] +
                            x$latent[nrow(adj_mat) + 1]))
  
  sd_est_mat <- SUMMER::expit(sd_est_mat)
  out_dat <- data.frame(region = rownames(adj_mat),
                        est = rowMeans(sd_est_mat),
                        median = apply(sd_est_mat, 1,
                                       function(x) median(x, na.rm = T)),
                        var = apply(sd_est_mat, 1, var),
                        lower = apply(sd_est_mat, 1,
                                      function(x) quantile(x, (1-CI)/2)),
                        upper = apply(sd_est_mat, 1,
                                      function(x) quantile(x, 1-(1-CI)/2)),
                        method = paste0("bymSLogit", direct_est$method[1]))
  return(out_dat)
}
#### GAUSSIAN ULM ####
get_iid_gau_ulm <- function(formula, by, 
                            pop_dat, sample_dat, 
                            h = NULL,
                            pc_u = 5,
                            pc_alpha = 0.01, 
                            CI = .90) {
  hyperpc_iid_int <- list(prec = list(prior = "pc.prec",
                                      param = c(pc_u , pc_alpha)))
  
  pop_dat$region <-  as.vector(pop_dat[[all.vars(by)[1]]])
  sample_dat$region <-  as.vector(sample_dat[[all.vars(by)[1]]])
  
  reg_vec <- unique(pop_dat$region)
  pop_dat$reg_id <- match(pop_dat$region, reg_vec)
  sample_dat$reg_id <- match(sample_dat$region, reg_vec)
  inla_frm <- 
    update(formula, ~ . + f(reg_id, model = "iid", hyper = hyperpc_iid_int))
  ulm_fit <-
    INLA::inla(inla_frm,
               family = "gaussian", data = sample_dat,
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  ulm_fit_sample <- inla.posterior.sample(n = 1000, ulm_fit)
  pop_X <- model.matrix(formula, pop_dat)
  fe_idx <- grep(colnames(pop_X)[1], rownames(ulm_fit_sample[[1]]$latent))
  fe_idx <- fe_idx:(fe_idx + ncol(pop_X) - 1)
  re_idx <- grep("reg_id", x = rownames(ulm_fit_sample[[1]]$latent))
  summary_sample <- function(x) {
    pop_unit_ests <-
      x$latent[re_idx][pop_dat$reg_id] + pop_X %*% x$latent[fe_idx] +
      rnorm(nrow(pop_X), sd = sqrt(1 / x$hyperpar[1]))
    if (!is.null(h)) {
      pop_unit_ests <- h(pop_unit_ests)
    }
    area_ests <- aggregate(pop_unit_ests, list(region = pop_dat$reg_id), mean)
    return(area_ests[match(1:length(re_idx), area_ests[, 1]), 2])
  }
  ulm_est_mat <- do.call(cbind, lapply(ulm_fit_sample, summary_sample))
  out_dat <-
    data.frame(region = reg_vec,
               est = rowMeans(ulm_est_mat),
               median = apply(ulm_est_mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(ulm_est_mat, 1, var),
               lower = apply(ulm_est_mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(ulm_est_mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = "iidGauULM")
  return(out_dat)
}


get_bym_gau_ulm <- function(formula, by, 
                            adj_mat,
                            pop_dat, sample_dat, 
                            h = NULL,
                            pc_u = 5,
                            pc_alpha = 0.01, 
                            pc_u_phi = 0.5,
                            pc_alpha_phi = 2/3,
                            CI = .90) {
  hyperpc_bym_int <- list(
    prec = list(prior = "pc.prec", param = c(pc_u , pc_alpha)),  
    phi = list(prior = 'pc', param = c(pc_u_phi , pc_alpha_phi))
  )
  
  pop_dat$region <-  as.vector(pop_dat[[all.vars(by)[1]]])
  sample_dat$region <-  as.vector(sample_dat[[all.vars(by)[1]]])

  pop_dat$reg_id <- match(pop_dat$region, rownames(adj_mat))
  sample_dat$reg_id <- match(sample_dat$region, rownames(adj_mat))
  inla_frm <- 
    update(formula, ~ . + f(reg_id, model = "bym2", 
                            graph = adj_mat, 
                            hyper = hyperpc_bym_int))
  ulm_fit <-
    INLA::inla(inla_frm,
               family = "gaussian", data = sample_dat,
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  ulm_fit_sample <- inla.posterior.sample(n = 1000, ulm_fit)
  pop_X <- model.matrix(formula, pop_dat)
  fe_idx <- grep(colnames(pop_X)[1], rownames(ulm_fit_sample[[1]]$latent))
  fe_idx <- fe_idx:(fe_idx + ncol(pop_X) - 1)
  re_idx <- grep("reg_id", x = rownames(ulm_fit_sample[[1]]$latent))
  re_idx <- re_idx[1:(length(re_idx) / 2)]
  summary_sample <- function(x) {
    pop_unit_ests <-
      x$latent[re_idx][pop_dat$reg_id] + pop_X %*% x$latent[fe_idx] +
      rnorm(nrow(pop_X), sd = sqrt(1 / x$hyperpar[1]))
    if (!is.null(h)) {
      pop_unit_ests <- h(pop_unit_ests)
    }
    area_ests <- aggregate(pop_unit_ests, list(region = pop_dat$reg_id), mean)
    return(area_ests[match(1:length(re_idx), area_ests[, 1]), 2])
  }
  ulm_est_mat <- do.call(cbind, lapply(ulm_fit_sample, summary_sample))
  out_dat <-
    data.frame(region = rownames(adj_mat),
               est = rowMeans(ulm_est_mat),
               median = apply(ulm_est_mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(ulm_est_mat, 1, var),
               lower = apply(ulm_est_mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(ulm_est_mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = "bymGauULM")
  return(out_dat)
}
#### BINOMIAL ULM ####
get_iid_bin_ulm_fit <- function(formula, by, 
                                sample_dat,
                                pc_u = 5,
                                pc_alpha = 0.01) {
  hyperpc_iid_int <- list(prec = list(prior = "pc.prec",
                                      param = c(pc_u , pc_alpha)))
  
  sample_dat$region <-  as.vector(sample_dat[[all.vars(by)[1]]])
  
  reg_vec <- unique(sample_dat$region)
  sample_dat$reg_id <- match(sample_dat$region, reg_vec)
  inla_frm <- 
    update(formula, ~ . + f(reg_id, model = "iid", hyper = hyperpc_iid_int))
  ulm_fit <-
    INLA::inla(inla_frm,
               family = "binomial", data = sample_dat,
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  ulm_fit_sample <- inla.posterior.sample(n = 1000, ulm_fit)
  return(list(sample_dat = sample_dat,
              ulm_fit = ulm_fit,
              ulm_fit_sample = ulm_fit_sample))
}

get_iid_bin_ulm_est <- function(formula, by, 
                                sample_dat,
                                pop_dat,
                                ulm_fit,
                                ulm_fit_sample,
                                pop = NULL, CI = .90) {
  pop_dat$region <-  as.vector(pop_dat[[all.vars(by)[1]]])
  reg_vec <- unique(sample_dat$region)
  pop_dat$reg_id <- match(pop_dat$region, reg_vec)
  cov_frm <- update(formula, NULL ~ .)
  pop_X <- model.matrix(cov_frm, pop_dat)
  fe_idx <- grep(colnames(pop_X)[1], rownames(ulm_fit_sample[[1]]$latent))
  fe_idx <- fe_idx:(fe_idx + ncol(pop_X) - 1)
  re_idx <- grep("reg_id", x = rownames(ulm_fit_sample[[1]]$latent))
  summary_sample <- function(x) {
    pop_unit_ests <- x$latent[re_idx][pop_dat$reg_id] 
    pop_unit_ests[is.na(pop_unit_ests)] <-
      rnorm(sum(is.na(pop_unit_ests)), sd = sqrt(1 / x$hyperpar[1]))
    pop_unit_ests <- pop_unit_ests + pop_X %*% x$latent[fe_idx]
    pop_unit_ests <- SUMMER::expit(pop_unit_ests)
    
    if (is.null(pop)) {
      area_ests <- aggregate(pop_unit_ests, list(region = pop_dat$region), mean)
    } else {
      pop_dat$pop <- pop_dat[[all.vars(pop)[1]]]
      wt_area_ests <-
        aggregate(pop_unit_ests * pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_pops <-
        aggregate(pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_ests <- wt_area_ests
      area_ests[, 2] <- wt_area_ests[, 2] / area_pops[, 2]
    }
    
    area_ests <- area_ests[match(unique(pop_dat$region), area_ests[, 1]), ]
    return(area_ests[, 2, drop=F])
  }
  ulm_est_mat <- do.call(cbind, lapply(ulm_fit_sample, summary_sample))
  out_dat <-
    data.frame(region = unique(pop_dat$region),
               est = rowMeans(ulm_est_mat),
               median = apply(ulm_est_mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(ulm_est_mat, 1, var),
               lower = apply(ulm_est_mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(ulm_est_mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = "iidBinULM")
  return(out_dat)
}

get_iid_bin_ulm <- function(formula, by,
                            sample_dat,
                            pop_dat, 
                            pop = NULL,
                            pc_u = 5,
                            pc_alpha = 0.01, 
                            CI = .90) {
  fit_obj <- get_iid_bin_ulm_fit(formula, by,
                                 sample_dat)
  out_dat <-
    get_iid_bin_ulm_est(formula, by, 
                        fit_obj$sample_dat, 
                        pop_dat,
                        fit_obj$ulm_fit, 
                        fit_obj$ulm_fit_sample,
                        pop, CI)
  
  return(out_dat)
}


get_bym_bin_ulm_fit <- function(formula, by,
                                adj_mat,
                                sample_dat,
                                pc_u = 5,
                                pc_alpha = 0.01,
                                pc_u_phi = 0.5,
                                pc_alpha_phi = 2/3) {
  hyperpc_bym_int <- list(
    prec = list(prior = "pc.prec", param = c(pc_u , pc_alpha)),  
    phi = list(prior = 'pc', param = c(pc_u_phi , pc_alpha_phi))
  )

  sample_dat$region <-  as.vector(sample_dat[[all.vars(by)[1]]])
  sample_dat$reg_id <- match(sample_dat$region, rownames(adj_mat))
  inla_frm <- 
    update(formula, ~ . + f(reg_id, model = "bym2", 
                            graph = adj_mat, 
                            hyper = hyperpc_bym_int))
  ulm_fit <-
    INLA::inla(inla_frm,
               family = "binomial", data = sample_dat,
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  ulm_fit_sample <- inla.posterior.sample(n = 1000, ulm_fit)
  return(list(sample_dat = sample_dat,
              ulm_fit = ulm_fit,
              ulm_fit_sample = ulm_fit_sample))
}
get_bym_bin_ulm_est <- function(formula, by,
                                adj_mat,
                                sample_dat,
                                pop_dat,
                                ulm_fit,
                                ulm_fit_sample,
                                pop = NULL, CI = .90) {
  pop_dat$region <-  as.vector(pop_dat[[all.vars(by)[1]]])
  pop_dat$reg_id <- match(pop_dat$region, rownames(adj_mat))
  
  cov_frm <- update(formula, NULL ~ .)
  pop_X <- model.matrix(cov_frm, pop_dat)
  fe_idx <- grep(colnames(pop_X)[1], rownames(ulm_fit_sample[[1]]$latent))
  fe_idx <- fe_idx:(fe_idx + ncol(pop_X) - 1)
  re_idx <- grep("reg_id", x = rownames(ulm_fit_sample[[1]]$latent))
  re_idx <- re_idx[1:(length(re_idx) / 2)]
  summary_sample <- function(x) {
    pop_unit_ests <- x$latent[re_idx][pop_dat$reg_id] 
    pop_unit_ests[is.na(pop_unit_ests)] <-
      rnorm(sum(is.na(pop_unit_ests)), sd = sqrt(1 / x$hyperpar[1]))
    pop_unit_ests <- pop_unit_ests + pop_X %*% x$latent[fe_idx]
    pop_unit_ests <- SUMMER::expit(pop_unit_ests)
    
    if (is.null(pop)) {
      area_ests <- aggregate(pop_unit_ests, list(region = pop_dat$region), mean)
    } else {
      pop_dat$pop <- pop_dat[[all.vars(pop)[1]]]
      wt_area_ests <-
        aggregate(pop_unit_ests * pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_pops <-
        aggregate(pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_ests <- wt_area_ests
      area_ests[, 2] <- wt_area_ests[, 2] / area_pops[, 2]
    }
    area_ests <- area_ests[match(unique(pop_dat$region), area_ests[, 1]), ]
    return(area_ests[, 2, drop=F])
  }
  ulm_est_mat <- do.call(cbind, lapply(ulm_fit_sample, summary_sample))
  out_dat <-
    data.frame(region = unique(pop_dat$region),
               est = rowMeans(ulm_est_mat),
               median = apply(ulm_est_mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(ulm_est_mat, 1, var),
               lower = apply(ulm_est_mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(ulm_est_mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = "bymBinULM")
  return(out_dat)
  
}

get_bym_bin_ulm <- function(formula, by,
                            adj_mat,
                            sample_dat,
                            pop_dat, 
                            pop = NULL,
                            pc_u = 5,
                            pc_alpha = 0.01, 
                            CI = .90) {
  fit_obj <- get_bym_bin_ulm_fit(formula, by,
                                 adj_mat,
                                 sample_dat)
  out_dat <-
    get_bym_bin_ulm_est(formula, by, 
                        adj_mat,
                        fit_obj$sample_dat, 
                        pop_dat,
                        fit_obj$ulm_fit, 
                        fit_obj$ulm_fit_sample,
                        pop, CI)
  
  return(out_dat)
}

#### BETABINOMIAL ULM ####
get_iid_bbin_ulm_fit <- function(formula, by, 
                                sample_dat,
                                pc_u = 5,
                                pc_alpha = 0.01) {
  hyperpc_iid_int <- list(prec = list(prior = "pc.prec",
                                      param = c(pc_u , pc_alpha)))
  
  sample_dat$region <-  as.vector(sample_dat[[all.vars(by)[1]]])
  
  reg_vec <- unique(sample_dat$region)
  sample_dat$reg_id <- match(sample_dat$region, reg_vec)
  inla_frm <- 
    update(formula, ~ . + f(reg_id, model = "iid", hyper = hyperpc_iid_int))
  ulm_fit <-
    INLA::inla(inla_frm,
               family = "betabinomial", data = sample_dat,
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  ulm_fit_sample <- inla.posterior.sample(n = 1000, ulm_fit)
  return(list(sample_dat = sample_dat,
              ulm_fit = ulm_fit,
              ulm_fit_sample = ulm_fit_sample))
}

get_iid_bbin_ulm_est <- function(formula, by, 
                                sample_dat,
                                pop_dat,
                                ulm_fit,
                                ulm_fit_sample,
                                pop = NULL, CI = .90) {
  pop_dat$region <-  as.vector(pop_dat[[all.vars(by)[1]]])
  reg_vec <- unique(sample_dat$region)
  pop_dat$reg_id <- match(pop_dat$region, reg_vec)
  cov_frm <- update(formula, NULL ~ .)
  pop_X <- model.matrix(cov_frm, pop_dat)
  fe_idx <- grep(colnames(pop_X)[1], rownames(ulm_fit_sample[[1]]$latent))
  fe_idx <- fe_idx:(fe_idx + ncol(pop_X) - 1)
  re_idx <- grep("reg_id", x = rownames(ulm_fit_sample[[1]]$latent))
  summary_sample <- function(x) {
    pop_unit_ests <- x$latent[re_idx][pop_dat$reg_id] 
    pop_unit_ests[is.na(pop_unit_ests)] <-
      rnorm(sum(is.na(pop_unit_ests)), sd = sqrt(1 / x$hyperpar[1]))
    pop_unit_ests <- pop_unit_ests + pop_X %*% x$latent[fe_idx]
    pop_unit_ests <- SUMMER::expit(pop_unit_ests)
    
    if (is.null(pop)) {
      area_ests <- aggregate(pop_unit_ests, list(region = pop_dat$region), mean)
    } else {
      pop_dat$pop <- pop_dat[[all.vars(pop)[1]]]
      wt_area_ests <-
        aggregate(pop_unit_ests * pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_pops <-
        aggregate(pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_ests <- wt_area_ests
      area_ests[, 2] <- wt_area_ests[, 2] / area_pops[, 2]
    }
    
    area_ests <- area_ests[match(unique(pop_dat$region), area_ests[, 1]), ]
    return(area_ests[, 2, drop=F])
  }
  ulm_est_mat <- do.call(cbind, lapply(ulm_fit_sample, summary_sample))
  out_dat <-
    data.frame(region = unique(pop_dat$region),
               est = rowMeans(ulm_est_mat),
               median = apply(ulm_est_mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(ulm_est_mat, 1, var),
               lower = apply(ulm_est_mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(ulm_est_mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = "iidBBinULM")
  return(out_dat)
}

get_iid_bbin_ulm <- function(formula, by,
                            sample_dat,
                            pop_dat, 
                            pop = NULL,
                            pc_u = 5,
                            pc_alpha = 0.01, 
                            CI = .90) {
  fit_obj <- get_iid_bbin_ulm_fit(formula, by,
                                 sample_dat)
  out_dat <-
    get_iid_bbin_ulm_est(formula, by, 
                        fit_obj$sample_dat, 
                        pop_dat,
                        fit_obj$ulm_fit, 
                        fit_obj$ulm_fit_sample,
                        pop, CI)
  
  return(out_dat)
}


get_bym_bbin_ulm_fit <- function(formula, by,
                                adj_mat,
                                sample_dat,
                                pc_u = 5,
                                pc_alpha = 0.01,
                                pc_u_phi = 0.5,
                                pc_alpha_phi = 2/3) {
  hyperpc_bym_int <- list(
    prec = list(prior = "pc.prec", param = c(pc_u , pc_alpha)),  
    phi = list(prior = 'pc', param = c(pc_u_phi , pc_alpha_phi))
  )
  
  sample_dat$region <-  as.vector(sample_dat[[all.vars(by)[1]]])
  sample_dat$reg_id <- match(sample_dat$region, rownames(adj_mat))
  inla_frm <- 
    update(formula, ~ . + f(reg_id, model = "bym2", 
                            graph = adj_mat, 
                            hyper = hyperpc_bym_int))
  ulm_fit <-
    INLA::inla(inla_frm,
               family = "betabinomial", data = sample_dat,
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  ulm_fit_sample <- inla.posterior.sample(n = 1000, ulm_fit)
  return(list(sample_dat = sample_dat,
              ulm_fit = ulm_fit,
              ulm_fit_sample = ulm_fit_sample))
}
get_bym_bbin_ulm_est <- function(formula, by,
                                adj_mat,
                                sample_dat,
                                pop_dat,
                                ulm_fit,
                                ulm_fit_sample,
                                pop = NULL, CI = .90) {
  pop_dat$region <-  as.vector(pop_dat[[all.vars(by)[1]]])
  pop_dat$reg_id <- match(pop_dat$region, rownames(adj_mat))
  
  cov_frm <- update(formula, NULL ~ .)
  pop_X <- model.matrix(cov_frm, pop_dat)
  fe_idx <- grep(colnames(pop_X)[1], rownames(ulm_fit_sample[[1]]$latent))
  fe_idx <- fe_idx:(fe_idx + ncol(pop_X) - 1)
  re_idx <- grep("reg_id", x = rownames(ulm_fit_sample[[1]]$latent))
  re_idx <- re_idx[1:(length(re_idx) / 2)]
  summary_sample <- function(x) {
    pop_unit_ests <- x$latent[re_idx][pop_dat$reg_id] 
    pop_unit_ests[is.na(pop_unit_ests)] <-
      rnorm(sum(is.na(pop_unit_ests)), sd = sqrt(1 / x$hyperpar[1]))
    pop_unit_ests <- pop_unit_ests + pop_X %*% x$latent[fe_idx]
    pop_unit_ests <- SUMMER::expit(pop_unit_ests)
    
    if (is.null(pop)) {
      area_ests <- aggregate(pop_unit_ests, list(region = pop_dat$region), mean)
    } else {
      pop_dat$pop <- pop_dat[[all.vars(pop)[1]]]
      wt_area_ests <-
        aggregate(pop_unit_ests * pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_pops <-
        aggregate(pop_dat$pop, 
                  list(region = pop_dat$region),
                  sum)
      area_ests <- wt_area_ests
      area_ests[, 2] <- wt_area_ests[, 2] / area_pops[, 2]
    }
    area_ests <- area_ests[match(unique(pop_dat$region), area_ests[, 1]), ]
    return(area_ests[, 2, drop=F])
  }
  ulm_est_mat <- do.call(cbind, lapply(ulm_fit_sample, summary_sample))
  out_dat <-
    data.frame(region = unique(pop_dat$region),
               est = rowMeans(ulm_est_mat),
               median = apply(ulm_est_mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(ulm_est_mat, 1, var),
               lower = apply(ulm_est_mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(ulm_est_mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = "bymBBinULM")
  return(out_dat)
  
}

get_bym_bbin_ulm <- function(formula, by,
                            adj_mat,
                            sample_dat,
                            pop_dat, 
                            pop = NULL,
                            pc_u = 5,
                            pc_alpha = 0.01, 
                            CI = .90) {
  fit_obj <- get_bym_bbin_ulm_fit(formula, by,
                                 adj_mat,
                                 sample_dat)
  out_dat <-
    get_bym_bbin_ulm_est(formula, by, 
                        adj_mat,
                        fit_obj$sample_dat, 
                        pop_dat,
                        fit_obj$ulm_fit, 
                        fit_obj$ulm_fit_sample,
                        pop, CI)
  
  return(out_dat)
}

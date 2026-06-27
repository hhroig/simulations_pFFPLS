library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)
library(plotly)
library(writexl)


beta_num_to_text <- function(beta_num_in) {
  if (beta_num_in == 1) {
    beta_txt = "symm" # symmetrical function
  }else if (beta_num_in == 2) {
    beta_txt = "exp" # single exponential top right corner
  }else if (beta_num_in == 3) {
    beta_txt = "saddle" # a horse saddle
  }else if (beta_num_in == 4) {
    beta_txt = "dbl_exp" # a double exponential top right and bottom left
  }
  
  return(beta_txt)
  
}

compare_methods_fun <- function(input_folder, 
                                zoom_r2_lower = 0, 
                                do_rough_r2 = TRUE,
                                theta = 30,   # Angle for viewing (rotation beta surface)
                                phi = 30    # Angle for viewing (tilt beta surface)
){
  
  out_folder <- paste0(input_folder, "results_plots/")
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  
  # show_col(hue_pal()(6))
  
  color_codes <- c(
    "pFFPLS" = hue_pal()(6)[1],
    "pFFR_I" = hue_pal()(6)[2],
    "pFFR_RS" = hue_pal()(6)[3],
    "FFPLS_OB" = hue_pal()(6)[4],
    "FFPLS" = hue_pal()(6)[5]
  )
  
  ## Final Models Files  ---------------------------------------------------
  
  # IMSE for Beta and validation Y:
  all_final_res <- data.frame()
  
  final_res_files <- list.files(path = input_folder, pattern = "final_models")
  
  for (ind_file in final_res_files) {
    
    all_final_res <- rbind(
      all_final_res,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # CVEs:
  all_cves <- data.frame()
  
  all_cves_files <- list.files(path = input_folder, pattern = "cves_rep")
  
  for (ind_file in all_cves_files) {
    
    all_cves <- rbind(
      all_cves,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # Computation time
  all_computation_times <- data.frame()
  
  all_computation_times_files <- list.files(path = input_folder, pattern = "computation_times")
  
  for (ind_file in all_computation_times_files) {
    
    all_computation_times <- rbind(
      all_computation_times,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  # Best number bases for FFPLS:
  
  all_best_num_bases_FFPLS <- data.frame()
  
  all_best_num_bases_FFPLS_files <- list.files(path = input_folder, pattern = "best_num_bases_FFPLS")
  
  for (ind_file in all_best_num_bases_FFPLS_files) {
    
    all_best_num_bases_FFPLS <- rbind(
      all_best_num_bases_FFPLS,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  # Best Lambdas files:
  
  all_best_lambdas <- data.frame()
  
  all_best_lambdas_files <- list.files(path = input_folder, pattern = "best_lambdas")
  
  for (ind_file in all_best_lambdas_files) {
    
    all_best_lambdas <- rbind(
      all_best_lambdas,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # All R^2 files:
  
  all_best_r2 <- data.frame()
  
  all_best_r2_files <- list.files(path = input_folder, pattern = "R2s_rep")
  
  for (ind_file in all_best_r2_files) {
    
    all_best_r2 <- rbind(
      all_best_r2,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # Betas files:
  
  all_betas <- data.frame()
  
  all_betas_files <- list.files(path = input_folder, pattern = "betas_")
  
  for (ind_file in all_betas_files) {
    
    all_betas <- rbind(
      all_betas,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  all_final_res <- all_final_res %>%
    mutate(nComp = as.factor(nComp))
  
  all_cves <- all_cves %>%
    mutate(nComp = as.factor(nComp))
  
  all_best_lambdas <- all_best_lambdas %>%
    mutate(nComp = as.factor(nComp)) 
  
  all_computation_times <- all_computation_times %>%
    mutate(nComp = as.factor(nComp)) 
  
  all_betas <- all_betas %>%
    mutate(nComp = as.factor(nComp))
  
  if (nrow(all_best_num_bases_FFPLS) > 0) {
    all_best_num_bases_FFPLS <- all_best_num_bases_FFPLS %>%
      mutate(nComp = as.factor(nComp)) 
  }
  
  
  
  # Limits:
  
  # cve_limits <- range(all_cves$CVE)
  # imse_limits <- range(all_final_res$imse )
  # mean_imse_Y_val_limits <- range(all_final_res$mean_imse_Y_val )
  
  
  # IMSE + MSE --------------------------------------------------------------
  
  
  out_folder_IMSE_CVEs <- paste0(out_folder, "IMSEs_CVEs_Excel/")
  
  if (!dir.exists(out_folder_IMSE_CVEs)) {
    dir.create(out_folder_IMSE_CVEs)
  }
  
  list_of_beta_paths = list()
  
  for (uniq_beta in unique(all_cves$beta.num)) {
    
    list_of_beta_paths[[uniq_beta]] <- paste0(out_folder, "IMSEs_CVEs_", 
                                              beta_num_to_text(uniq_beta), "/")
    
    
    if (!dir.exists(list_of_beta_paths[[uniq_beta]])) {
      dir.create(list_of_beta_paths[[uniq_beta]])
    }
  }
  
  # Excel tables
  
  all_cves_summ <- all_cves %>%
    group_by(method, nComp, beta.num) %>%
    summarize(
      mean_CVE = mean(CVE, na.rm = TRUE),
      median_CVE = median(CVE, na.rm = TRUE),
      sd_CVE = sd(CVE, na.rm = TRUE),
      iqr_CVE = IQR(CVE, na.rm = TRUE)
    )
  write_xlsx(all_cves_summ, paste0(out_folder_IMSE_CVEs, "summary_training_cves.xlsx"))
  
  
  all_imse_beta_summ <- all_final_res %>%
    group_by(method, nComp, beta.num) %>%
    summarize(
      mean_IMSE = mean(imse, na.rm = TRUE),
      median_IMSE = median(imse, na.rm = TRUE),
      sd_IMSE = sd(imse, na.rm = TRUE),
      iqr_IMSE = IQR(imse, na.rm = TRUE)
    )
  write_xlsx(all_imse_beta_summ, paste0(out_folder_IMSE_CVEs, "summary_imse_beta.xlsx"))
  
  
  all_val_imse_summ <- all_final_res %>%
    group_by(method, nComp, beta.num) %>%
    summarize(
      mean_mean_imse_Y_val  = mean(mean_imse_Y_val , na.rm = TRUE),
      median_mean_imse_Y_val  = median(mean_imse_Y_val , na.rm = TRUE),
      sd_mean_imse_Y_val  = sd(mean_imse_Y_val , na.rm = TRUE),
      iqr_mean_imse_Y_val  = IQR(mean_imse_Y_val , na.rm = TRUE)
    )
  write_xlsx(all_val_imse_summ, paste0(out_folder_IMSE_CVEs, "summary_imse_val_y.xlsx"))
  
  
  
  # Plots
  
  for (beta_num in unique(all_final_res$beta.num)) {
    
    p_cve <- ggplot(all_cves %>% filter(beta.num == beta_num),
                    aes(x = nComp, y = CVE, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: CVE(Y)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="none", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_imse <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                     aes(x = nComp, y = imse, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "IMSE(" * beta *")" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20))+
      labs(fill = "")
    
    p_imse_val <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                         aes(x = nComp, y = mean_imse_Y_val, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "Test: IMSE(Y)" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    # get legend
    leg <- get_legend(p_imse_val)
    # remove from last plot
    p_imse_val <- p_imse_val + theme(legend.position="none")
    
    p_both <- grid.arrange(p_cve, p_imse_val, nrow = 1, bottom = leg)
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_both,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_cveY_val_imseY_beta_", beta_num_to_text(beta_num),".png")  ),
           width = 16, height = 8 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_both,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_cveY_val_imseY_beta_", beta_num_to_text(beta_num),".pdf")  ),
           width = 16, height = 8 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_imse,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("imse_beta_", beta_num_to_text(beta_num),".png")  ),
           width = 12, height = 8 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_imse,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("imse_beta_", beta_num_to_text(beta_num),".pdf")  ),
           width = 12, height = 8 )
    
    
    
    # Log-scale
    p_cve_log <- ggplot(all_cves %>% filter( beta.num == beta_num),
                        aes(x = nComp, y = log(CVE, base = 10), fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: log{ CVE(Y) }") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="none", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_imse_log <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                         aes(x = nComp, y = log(imse, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab(expression( "log { IMSE(" * beta *") }" ) ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20))+
      labs(fill = "")
    
    p_imse_val_log <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                             aes(x = nComp, y = log(mean_imse_Y_val, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab(expression( "Test: log { IMSE(Y) }" ) ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20))+
      labs(fill = "")
    
    # get legend
    leg <- get_legend(p_imse_val_log)
    # remove from last plot
    p_imse_val_log <- p_imse_val_log + theme(legend.position="none")
    
    p_both_log <- grid.arrange(p_cve_log, p_imse_val_log, nrow = 1, bottom = leg)
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_both_log,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_cveY_val_imseY_log_beta_", beta_num_to_text(beta_num),
                                    ".pdf")  ),
           width = 16, height = 8 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_both_log,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_cveY_val_imseY_log_beta_", beta_num_to_text(beta_num),
                                    ".png")  ),
           width = 16, height = 8 )
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_imse_log,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("log_imse_beta_", beta_num_to_text(beta_num),
                                    ".pdf")  ),
           width = 12, height = 6 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_imse_log,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("log_imse_beta_", beta_num_to_text(beta_num),
                                    ".png")  ),
           width = 12, height = 6 )
    
    
    
  } # loop "beta.num"
  
  
  ### Disaggregated (original scale) -----
  
  for (beta_num in unique(all_final_res$beta.num)) {
    
    p_cve <- ggplot(all_cves %>% 
                      filter(beta.num == beta_num,
                             str_detect(method, "PLS" )),
                    aes(x = nComp, y = CVE, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: CVE(Y)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_cve,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_PLS_cveY_", beta_num_to_text(beta_num),".png")  ),
           width = 8, height = 6 )
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_cve,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_PLS_cveY_", beta_num_to_text(beta_num),".pdf")  ),
           width = 8, height =6 )
    
    
    
    p_imse_val <- ggplot(all_final_res %>% filter( beta.num == beta_num,
                                                   str_detect(method, "PLS")  ),
                         aes(x = nComp, y = mean_imse_Y_val, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "Test: IMSE(Y)" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_imse_val,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("test_PLS_imseY_", beta_num_to_text(beta_num),".png")  ),
           width = 8, height = 6 )
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_imse_val,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("test_PLS_imseY_", beta_num_to_text(beta_num),".pdf")  ),
           width = 8, height = 6 )
    
    
    
    for (opt_comp in unique(all_final_res$nComp)) {
      
      # IMSE on Beta
      
      p_imse <- ggplot(all_final_res %>% 
                         filter(beta.num == beta_num, 
                                nComp == opt_comp),
                       aes(x = method, y = imse, fill = method)) +
        geom_boxplot(position=position_dodge(0.8)) +
        ylab( expression( "IMSE(" * beta *")" )  ) +
        xlab("") +
        scale_fill_manual(values = color_codes) +
        theme_bw() +
        theme(legend.position="none", text = element_text(size = 20))+
        labs(fill = "", subtitle =paste0("PLS methods using ", opt_comp ," components")) 
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("imse_beta_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".png")  ),
             width = 8, height = 6 )
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("imse_beta_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".pdf")   ),
             width = 8, height = 6 )
      
      
      # IMSE on test set Y
      
      p_imse_val <- ggplot(all_final_res %>% filter( beta.num == beta_num, 
                                                     nComp == opt_comp),
                           aes(x = method, y = mean_imse_Y_val, fill = method)) +
        geom_boxplot(position=position_dodge(0.8)) +
        ylab( expression( "Test: IMSE(Y)" )  ) +
        xlab("") +
        scale_fill_manual(values = color_codes) +
        theme_bw() +
        theme(legend.position="none", text = element_text(size = 20)) +
        labs(fill = "", subtitle =paste0("PLS methods using ", opt_comp ," components"))
      
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse_val,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("test_imseY_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".png")  ),
             width = 8, height = 6 )
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse_val,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("test_imseY_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".pdf")   ),
             width = 8, height = 6 )
      
      
      
    } # end beta loop on original scale
  } # end beta loop on original scale
  
  
  
  ### Disaggregated (log) -----
  
  for (beta_num in unique(all_final_res$beta.num)) {
    
    
    p_cve_log <- ggplot(all_cves %>% 
                          filter(beta.num == beta_num,
                                 str_detect(method, "PLS" )),
                        aes(x = nComp, y = CVE, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      scale_y_log10() +  # Log-10 scale transformation
      ylab("Training: log{CVE(Y)}") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_cve_log,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_log_cveY_", beta_num_to_text(beta_num),".png")  ),
           width = 8, height = 6 )
    
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_cve_log,
           filename = paste0(list_of_beta_paths[[beta_num]],
                             paste0("train_log_cveY_", beta_num_to_text(beta_num),".pdf")  ),
           width = 8, height = 6 )
    
    
    for (opt_comp in unique(all_final_res$nComp)) {
      
      # IMSE on Beta
      
      p_imse_log <- ggplot(all_final_res %>% 
                             filter(beta.num == beta_num, 
                                    nComp == opt_comp),
                           aes(x = method, y = imse, fill = method)) +
        geom_boxplot(position=position_dodge(0.8)) +
        ylab( expression( " log{ IMSE(" * beta *") }" )  ) +
        xlab("") +
        scale_y_log10() +  # Log-10 scale transformation
        scale_fill_manual(values = color_codes) +
        theme_bw() +
        theme(legend.position="none", text = element_text(size = 20))+
        labs(fill = "", subtitle =paste0("PLS methods using ", opt_comp ," components")) 
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse_log,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("log_imse_beta_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".png")  ),
             width = 8, height = 6 )
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse_log,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("log_imse_beta_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".pdf")   ),
             width = 8, height = 6 )
      
      
      # IMSE on test set Y
      
      p_imse_val_log <- ggplot(all_final_res %>% filter( beta.num == beta_num, 
                                                         nComp == opt_comp),
                               aes(x = method, y = mean_imse_Y_val, fill = method)) +
        geom_boxplot(position=position_dodge(0.8)) +
        ylab( expression( "Test: log { IMSE(Y) }" )  ) +
        scale_y_log10() +  # Log-10 scale transformation
        xlab("") +
        scale_fill_manual(values = color_codes) +
        theme_bw() +
        theme(legend.position="none", text = element_text(size = 20)) +
        labs(fill = "", subtitle =paste0("PLS methods using ", opt_comp ," components"))
      
      
      if (!is.null(dev.list())) dev.off()
      
      
      ggsave(p_imse_val_log,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("log_test_imseY_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".png")  ),
             width = 8, height = 6 )
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_imse_val_log,
             filename = paste0(list_of_beta_paths[[beta_num]],
                               paste0("log_test_imseY_", 
                                      beta_num_to_text(beta_num), 
                                      "_ncomp", 
                                      opt_comp,
                                      ".pdf")   ),
             width = 8, height = 6 )
      
      
      
    } # end loop on n components
    
    
    
    
    
  } # end beta loop on original and log scales
  
  
  # Computation time --------------------------------------------------------------
  
  out_folder_computation_times <- paste0(out_folder, "computation_times/")
  
  if (!dir.exists(out_folder_computation_times)) {
    dir.create(out_folder_computation_times)
  }
  
  for (beta_num in unique(all_computation_times$beta.num)) {
    
    p_elapsed <- ggplot(all_computation_times %>% filter(beta.num == beta_num),
                        aes(x = nComp, y = elapsed_time, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Cross-Validation Time") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_elapsed_log <- ggplot(all_computation_times %>% filter(beta.num == beta_num),
                            aes(x = nComp, y = log(elapsed_time, base = 10), fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Cross-Validation Time (log-scale)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_elapsed,
           filename = paste0(out_folder_computation_times,
                             paste0("computaion_times_beta", beta_num,".png")  ),
           width = 8, height = 6 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_elapsed,
           filename = paste0(out_folder_computation_times,
                             paste0("computaion_times_beta", beta_num,".pdf")  ),
           width = 8, height = 6 )
    
    
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_elapsed_log,
           filename = paste0(out_folder_computation_times,
                             paste0("log_computaion_times_beta", beta_num,".png")  ),
           width = 8, height = 6 )
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_elapsed_log,
           filename = paste0(out_folder_computation_times,
                             paste0("log_computaion_times_beta", beta_num,".pdf")  ),
           width = 8, height = 6 )
    
    
  } # loop "beta.num"
  
  # Save as Excel file to report a single repetition:
  
  rownames(all_computation_times) <- NULL
  
  write_xlsx(all_computation_times, paste0(out_folder_computation_times, "all_computation_times.xlsx"))
  
  
  
  # Penalties ---------------------------------------------------------------
  
  df <- all_best_lambdas %>% 
    pivot_longer(cols = starts_with("lambda"),
                 names_to = "target",
                 values_to = "penalty") 
  
  out_folder_penalties <- paste0(out_folder, "best_penalties/")
  
  if (!dir.exists(out_folder_penalties)) {
    dir.create(out_folder_penalties)
  }
  
  number_of_pen_methods <- length(unique(df$method))
  
  for (beta_num in unique(df$beta.num)) {
    
    p_lamb <- ggplot(df %>% filter(beta.num == beta_num),
                     aes(x = nComp, y = penalty)) +
      facet_wrap(~method, ncol = number_of_pen_methods) +
      geom_boxplot() +
      ylab(expression(lambda)) +
      xlab("# of components") +
      theme_bw() +
      theme(text = element_text(size = 20))
    
    
    if (!is.null(dev.list())) dev.off()
    
    
    ggsave(p_lamb,
           filename = paste0(out_folder_penalties, "2d_best_lambdas_",  beta_num, ".png"),
           width = 8*number_of_pen_methods, height = 8 )
    
    
    if (!is.null(dev.list())) dev.off()
    
    ggsave(p_lamb,
           filename = paste0(out_folder_penalties, "2d_best_lambdas_",  beta_num, ".pdf"),
           width = 8*number_of_pen_methods, height = 8 )
  }
  
  
  
  
  # R^2 ------------------------------------------
  
  
  out_folder_R2 <- paste0(out_folder, "R2/")
  
  if (!dir.exists(out_folder_R2)) {
    dir.create(out_folder_R2)
  }
  
  
  
  summ_all_r2 <- all_best_r2 %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_OB",
                                              "pFFPLS",
                                              "pFFR_I",
                                              "pFFR_RS",
                                              "True Beta") )) %>%
    group_by(method, beta.num, nComp, q) %>%
    summarise(Training = mean(R2_train), Test = mean(R2_val))
  
  summ_all_r2_long <- summ_all_r2 %>%
    pivot_longer(cols = c(Training, Test),
                 names_to = "partition",
                 values_to = "r2")
  
  
  for ( n_comps_loop in unique(summ_all_r2$nComp) ){
    
    for (beta_num in unique(summ_all_r2$beta.num)) {
      
      if (!is.null(dev.list())) dev.off()
      
      p_r2_both <- ggplot(summ_all_r2_long %>% 
                            filter(beta.num == beta_num, nComp == n_comps_loop),
                          aes(x = q, y = r2, color = method)) +
        facet_wrap(~partition) +
        geom_smooth(se = F, linewidth = 1, alpha = 0.8)  +
        ylab("R^2") +
        xlab("q") +
        scale_color_manual(values = color_codes)+
        ylim(0, 1) +
        theme_bw()  +
        theme(legend.position="bottom", text = element_text(size = 20)) +
        labs(color = "")
      
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_r2_both,
             filename = paste0(out_folder_R2,
                               paste0("R2_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
             width = 12, height = 6 )
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_r2_both,
             filename = paste0(out_folder_R2,
                               paste0("R2_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
             width = 12, height = 6 )
      
      # Limit lower ylim for more details:
      
      p_r2_both_zoom <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                               aes(x = q, y = r2, color = method)) +
        facet_wrap(~partition) +
        geom_smooth(se = F, linewidth = 1, alpha = 0.8)  +
        ylab("R^2") +
        xlab("q") +
        scale_color_manual(values = color_codes)+
        ylim(zoom_r2_lower, 1) +
        theme_bw()  +
        theme(legend.position="bottom", text = element_text(size = 20)) +
        labs(color = "")
      
      if (!is.null(dev.list())) dev.off()
      
      
      ggsave(p_r2_both_zoom,
             filename = paste0(out_folder_R2,
                               paste0("R2_zoom_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
             width = 12, height = 6 )
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_r2_both_zoom,
             filename = paste0(out_folder_R2,
                               paste0("R2_zoom_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
             width = 12, height = 6 )
      
      
      
      
      if (do_rough_r2) {
        p_r2_both_rough <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                                  aes(x = q, y = r2, color = method)) +
          facet_wrap(~partition) +
          geom_line(linewidth = 1, alpha = 0.8)  +
          ylab("R^2") +
          xlab("q") +
          scale_color_manual(values = color_codes)+
          ylim(0, 1) +
          theme_bw()  +
          theme(legend.position="bottom", text = element_text(size = 20)) +
          labs(color = "")
        
        if (!is.null(dev.list())) dev.off()
        
        
        ggsave(p_r2_both_rough,
               filename = paste0(out_folder_R2,
                                 paste0("R2_rough_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
               width = 12, height = 6 )
        
        if (!is.null(dev.list())) dev.off()
        
        ggsave(p_r2_both_rough,
               filename = paste0(out_folder_R2,
                                 paste0("R2_rough_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
               width = 12, height = 6 )
        
        p_r2_both_rough_zoom <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                                       aes(x = q, y = r2, color = method)) +
          facet_wrap(~partition) +
          geom_line(linewidth = 1, alpha = 0.8)  +
          ylab("R^2") +
          xlab("q") +
          scale_color_manual(values = color_codes)+
          ylim(zoom_r2_lower, 1) +
          theme_bw()  +
          theme(legend.position="bottom", text = element_text(size = 20)) +
          labs(color = "")
        
        
        if (!is.null(dev.list())) dev.off()
        
        ggsave(p_r2_both_rough_zoom,
               filename = paste0(out_folder_R2,
                                 paste0("R2_rough_zoom_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
               width = 12, height = 6 )
        
        if (!is.null(dev.list())) dev.off()
        
        ggsave(p_r2_both_rough_zoom,
               filename = paste0(out_folder_R2,
                                 paste0("R2_rough_zoom_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
               width = 12, height = 6 )
        
      } # if do_rough_r2
      
      
    } # loop R^2
    
  } # loop ncomponents
  
  
  
  # Repeat each train/test graph separately
  
  for ( n_comps_loop in unique(summ_all_r2$nComp) ){
    
    for (beta_num in unique(summ_all_r2$beta.num)) {
      
      for (partition_out in unique(summ_all_r2_long$partition)) {
        
        
        p_r2_partition <- ggplot(summ_all_r2_long %>% 
                                   filter(beta.num == beta_num, 
                                          nComp == n_comps_loop,
                                          partition == partition_out),
                                 aes(x = q, y = r2, color = method)) +
          geom_smooth(se = F, linewidth = 1, alpha = 0.8)  +
          ylab("R^2") +
          xlab("q") +
          scale_color_manual(values = color_codes)+
          ylim(0, 1) +
          theme_bw()  +
          theme(legend.position="bottom", text = element_text(size = 20)) +
          labs(color = "")
        
        if (!is.null(dev.list())) dev.off()
        
        
        ggsave(p_r2_partition,
               filename = paste0(out_folder_R2,
                                 paste0("R2_", partition_out, 
                                        "_beta_", beta_num_to_text(beta_num), 
                                        "_nComp", n_comps_loop,".png")  ),
               width = 8, height = 6 )
        
        if (!is.null(dev.list())) dev.off()
        
        ggsave(p_r2_partition,
               filename = paste0(out_folder_R2,
                                 paste0("R2_", partition_out, 
                                        "_beta_", beta_num_to_text(beta_num),
                                        "_nComp", n_comps_loop,".pdf")  ),
               width = 8, height = 6 )
        
        # Limit lower ylim for more details:
        
        p_r2_partition_zoom <- ggplot(summ_all_r2_long %>% 
                                        filter(beta.num == beta_num, 
                                               nComp == n_comps_loop,
                                               partition == partition_out),
                                      aes(x = q, y = r2, color = method)) +
          geom_smooth(se = F, linewidth = 1, alpha = 0.8)  +
          ylab("R^2") +
          xlab("q") +
          scale_color_manual(values = color_codes)+
          ylim(zoom_r2_lower, 1) +
          theme_bw()  +
          theme(legend.position="bottom", text = element_text(size = 20)) +
          labs(color = "")
        
        if (!is.null(dev.list())) dev.off()
        
        
        ggsave(p_r2_partition_zoom,
               filename = paste0(out_folder_R2,
                                 paste0("R2_zoom_", partition_out, 
                                        "_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
               width = 8, height = 6 )
        
        if (!is.null(dev.list())) dev.off()
        
        ggsave(p_r2_partition_zoom,
               filename = paste0(out_folder_R2,
                                 paste0("R2_zoom_", partition_out, 
                                        "_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
               width = 8, height = 6 )
        
        
        
        
        if (do_rough_r2) {
          
          p_r2_both_rough <- ggplot(summ_all_r2_long %>% 
                                      filter(beta.num == beta_num, 
                                             nComp == n_comps_loop,
                                             partition == partition_out),
                                    aes(x = q, y = r2, color = method)) +
            geom_line(linewidth = 1, alpha = 0.8)  +
            ylab("R^2") +
            xlab("q") +
            scale_color_manual(values = color_codes)+
            ylim(0, 1) +
            theme_bw()  +
            theme(legend.position="bottom", text = element_text(size = 20)) +
            labs(color = "")
          
          if (!is.null(dev.list())) dev.off()
          
          
          ggsave(p_r2_both_rough,
                 filename = paste0(out_folder_R2,
                                   paste0("R2_rough_", partition_out, 
                                          "_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
                 width = 8, height = 6 )
          
          if (!is.null(dev.list())) dev.off()
          
          ggsave(p_r2_both_rough,
                 filename = paste0(out_folder_R2,
                                   paste0("R2_rough_", partition_out, 
                                          "_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
                 width = 8, height = 6 )
          
          p_r2_both_rough_zoom <- ggplot(summ_all_r2_long %>% 
                                           filter(beta.num == beta_num,
                                                  nComp == n_comps_loop,
                                                  partition == partition_out),
                                         aes(x = q, y = r2, color = method)) +
            geom_line(linewidth = 1, alpha = 0.8)  +
            ylab("R^2") +
            xlab("q") +
            scale_color_manual(values = color_codes)+
            ylim(zoom_r2_lower, 1) +
            theme_bw()  +
            theme(legend.position="bottom", text = element_text(size = 20)) +
            labs(color = "")
          
          if (!is.null(dev.list())) dev.off()
          
          
          ggsave(p_r2_both_rough_zoom,
                 filename = paste0(out_folder_R2,
                                   paste0("R2_rough_", partition_out, 
                                          "_zoom_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".png")  ),
                 width = 8, height = 6 )
          
          if (!is.null(dev.list())) dev.off()
          
          ggsave(p_r2_both_rough_zoom,
                 filename = paste0(out_folder_R2,
                                   paste0("R2_rough_", partition_out, 
                                          "_zoom_beta_", beta_num_to_text(beta_num), "_nComp", n_comps_loop,".pdf")  ),
                 width = 8, height = 6 )
          
        } # if do_rough_r2
        
        
        
      } # loop in partition
    } # loop R^2
    
  } # loop ncomponents
  
  
  # Betas -------------------------------------------------------------------
  
  df_true_betas <- all_betas %>% 
    ungroup() %>% 
    filter(method == "pFFPLS") %>% 
    mutate(z = z_true, method = "True Beta") %>% 
    dplyr::select(-z_true)
  
  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_OB",
                                              "pFFPLS",
                                              "pFFR_I",
                                              "pFFR_RS",
                                              "True Beta") )) %>%
    dplyr::select(-z_true) %>% 
    rbind(df_true_betas) %>% 
    group_by(method, beta.num, nComp, p, q) %>%
    summarise(mean_z = mean(z))
  
  
  plot_mean_betas <- function(summ_all_betas,
                              beta_num = 3,
                              n.Comp = 4,
                              bins = NA,
                              bin.width = 0.5,
                              line.size = 1,
                              num_col_wrap = 2) {
    
    betas_limits <-  summ_all_betas %>%
      filter(beta.num == beta_num) %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    plot_data <- summ_all_betas %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = p, y = q, z = mean_z)) +
      geom_raster(aes(fill = mean_z)) +
      facet_wrap( ~ method, ncol = num_col_wrap) +
      scale_fill_viridis() +
      theme_bw()+
      coord_fixed()+
      labs(fill = "") + xlab("p") + ylab("q") +
      theme(text = element_text(size = 20))
    # print(betas_plot)
    
    betas_plot_cont <- ggplot(data = plot_data,
                              mapping = aes(x = p, y = q, z = mean_z)) +
      geom_contour(aes(color=after_stat(level)), linewidth = line.size,
                   bins = bins, binwidth = bin.width) +
      # facet_grid( ~ method) +
      facet_wrap( ~ method, ncol = num_col_wrap) +
      scale_color_viridis() +
      theme_bw() +
      coord_fixed()+xlab("p") + ylab("q") +
      theme(text = element_text(size = 20))
    # print(betas_plot_cont)
    
    return(list(
      p_both_fill = betas_plot,
      p_both_cont = betas_plot_cont
    ))
    
  }
  
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    
    for (n.Beta in unique(summ_all_betas$beta.num)) {
      
      out_folder_mean_betas <- paste0(out_folder, "mean_beta_", beta_num_to_text(n.Beta), "/")
      
      if (!dir.exists(out_folder_mean_betas)) {
        dir.create(out_folder_mean_betas)
      }
      
      
      p_mean_fill <- plot_mean_betas(summ_all_betas ,
                                     beta_num = n.Beta,
                                     n.Comp = n.Comp,
                                     bins = NULL,
                                     bin.width = 0.1,
                                     line.size = 0.8,
                                     num_col_wrap = 3)[["p_both_fill"]]
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_mean_fill,
             path = out_folder_mean_betas,
             filename = paste0("2d_mean_beta_", n.Beta,
                               "_nComp_", n.Comp,
                               ".png"),
             width = 16, height = 8 )
      
      
      
      if (!is.null(dev.list())) dev.off()
      
      ggsave(p_mean_fill,
             path = out_folder_mean_betas,
             filename = paste0("2d_mean_beta_", n.Beta,
                               "_nComp_", n.Comp,
                               ".pdf"),
             width = 16, height = 8 )
      
      
      
    } # loop nComp
  } # loop beta.num
  
  
  
  
  ## 3D betas ---------------------------
  
  
  
  df_true_betas <- all_betas %>% 
    ungroup() %>% 
    filter(method == "pFFPLS", nComp == 1, rep_num == 1) %>% 
    mutate(method = "True Beta", z = z_true) %>% 
    dplyr::select(-z_true)
  
  
  
  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_OB",
                                              "pFFPLS",
                                              "pFFR_I",
                                              "pFFR_RS",
                                              "True Beta") )) %>%
    dplyr::select(-z_true) %>% 
    group_by(method, beta.num, nComp, p, q) %>%
    summarise(mean_z = mean(z))
  
  
  
  
  ## 3D beta as 2D -----
  
  plot_3D_betas_as_2D <- function(summ_all_betas,
                                  df_true_betas,
                                  beta_num = 1,
                                  n.Comp = 4,
                                  path = "3D_beta_plots/",
                                  theta = 30,   # Angle for viewing (rotation)
                                  phi = 30) {   # Angle for viewing (tilt)
    
    # Ensure the output directory exists
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    
    # Define True Beta:
    beta_true <- df_true_betas %>% 
      filter(beta.num == beta_num) %>% 
      .[["z"]]
    
    x_true <-  df_true_betas %>% 
      filter(method == "True Beta") %>% 
      filter(beta.num == beta_num) %>% 
      .[["p"]] %>% 
      unique()
    
    y_true <- df_true_betas %>% 
      filter(method == "True Beta") %>% 
      filter(beta.num == beta_num) %>% 
      .[["q"]] %>% 
      unique()
    
    z_true <- beta_true %>% matrix( nrow = length(x_true), ncol = length(y_true) )
    
    # Get plot limits out of estimations:
    betas_limits <- summ_all_betas %>%
      filter(beta.num == beta_num) %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    # Restrict to the actual beta we're studying:
    plot_data <- summ_all_betas %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    # Create list of all estimated betas, for all methods:
    estimated_betas <- list()
    estimated_betas_matrix <- list()
    
    for (unique_method in unique(plot_data[["method"]])) {
      
      plot_data_unique_method <- plot_data %>% 
        filter(method == unique_method)
      
      # Get the x and y values out of each estimation
      x <- plot_data_unique_method %>% 
        .[["p"]] %>% 
        unique()
      
      y <- plot_data_unique_method %>% 
        .[["q"]] %>% 
        unique()
      
      # Vectorized:
      estimated_betas[[unique_method]] <- plot_data_unique_method %>% 
        .[["mean_z"]]
      
      # Matrix form:
      estimated_betas_matrix[[unique_method]] <- matrix( 
        estimated_betas[[unique_method]],  
        nrow = length(x), 
        ncol = length(y) )
      
    }
    
    # Get the range of the estimated betas, including truth:
    zs_scale <- range(betas_limits, beta_true)
    
    
    # Iterate over each unique method:
    for (unique_method in names(estimated_betas_matrix)) {
      
      x <- plot_data %>% 
        filter(method == unique_method) %>% 
        .[["p"]] %>% 
        unique()
      
      y <- plot_data %>% 
        filter(method == unique_method) %>% 
        .[["q"]] %>% 
        unique()
      
      # File paths for saving
      pdf_file <- paste0(path, unique_method, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".pdf")
      eps_file <- paste0(path, unique_method, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".eps")
      png_file <- paste0(path, unique_method, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".png")
      
      # Save as PDF
      pdf(pdf_file, width = 7, height = 5)
      persp(x = x,
            y = y,
            z = estimated_betas_matrix[[unique_method]],
            col = "white",
            xlab = "p",
            ylab = "q",
            zlab = "z",
            zlim = zs_scale,
            theta = theta,
            phi = phi,
            expand = 0.5,
            shade = 0.5,
            ticktype = "detailed")
      dev.off()
      
      # Save as EPS
      postscript(eps_file, width = 7, height = 5, horizontal = FALSE, paper = "special")
      persp(x = x, 
            y = y, 
            z = estimated_betas_matrix[[unique_method]], 
            col = "white",
            xlab = "p", 
            ylab = "q", 
            zlab = "z", 
            zlim = zs_scale, 
            theta = theta, 
            phi = phi,
            expand = 0.5, 
            shade = 0.5, 
            ticktype = "detailed")
      dev.off()
      
      # Save as PNG
      png(png_file, width = 800, height = 600, res = 100)
      persp(x = x, 
            y = y, 
            z = estimated_betas_matrix[[unique_method]], 
            col = "white",
            xlab = "p", 
            ylab = "q", 
            zlab = "z", 
            zlim = zs_scale, 
            theta = theta, 
            phi = phi,
            expand = 0.5, 
            shade = 0.5, 
            ticktype = "detailed")
      dev.off()
      
      print(paste("Saved 3D plot for", unique_method, "as PDF, EPS, and PNG."))
    } # end iterate over each method
    
    # Now the true beta
    
    pdf_file <- paste0(path, "True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".pdf")
    eps_file <- paste0(path, "True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".eps")
    png_file <- paste0(path, "True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".png")
    
    # Save as PDF
    pdf(pdf_file, width = 7, height = 5)
    persp(x = x_true,
          y = y_true,
          z = z_true,
          col = "white",
          xlab = "p",
          ylab = "q",
          zlab = "z",
          zlim = zs_scale,
          theta = theta,
          phi = phi,
          expand = 0.5,
          shade = 0.5,
          ticktype = "detailed")
    dev.off()
    
    # Save as EPS
    postscript(eps_file, width = 7, height = 5, horizontal = FALSE, paper = "special")
    persp(x = x_true,
          y = y_true,
          z = z_true,
          col = "white",
          xlab = "p", 
          ylab = "q", 
          zlab = "z", 
          zlim = zs_scale, 
          theta = theta, 
          phi = phi,
          expand = 0.5, 
          shade = 0.5, 
          ticktype = "detailed")
    dev.off()
    
    # Save as PNG
    png(png_file, width = 800, height = 600, res = 100)
    persp(x = x_true,
          y = y_true,
          z = z_true,
          col = "white",
          xlab = "p", 
          ylab = "q", 
          zlab = "z", 
          zlim = zs_scale, 
          theta = theta, 
          phi = phi,
          expand = 0.5, 
          shade = 0.5, 
          ticktype = "detailed")
    dev.off()
    
    
    
    
  } # ends function 3D betas as 2D
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    for (n.Beta in unique(summ_all_betas$beta.num)) {
      
      out_folder_mean_betas2 <- paste0(out_folder, "3DBeta2D_", beta_num_to_text(n.Beta), "/")
      
      if (!dir.exists(out_folder_mean_betas2)) {
        dir.create(out_folder_mean_betas2)
      }
      
      plot_3D_betas_as_2D(summ_all_betas,
                          df_true_betas,
                          beta_num = n.Beta,
                          n.Comp = n.Comp,
                          path = out_folder_mean_betas2,
                          theta = 40,  # You can adjust the angle here
                          phi = 25)
    }
  }
  
  
}# end function

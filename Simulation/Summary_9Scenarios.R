###----------------------------------------------------###
###        Summary for Nine Simulation Scenarios       ###
###----------------------------------------------------###

rm(list=ls())

library(ggplot2)
library(reshape2)
library(gridExtra)

Meth <- c("HER-p", "HER-np", "BER", "UNO", "ITS")

make_cell <- function(mean, sd, digits=3){
  paste0(formatC(mean, format="f", digits=digits),
         " (", formatC(sd, format="f", digits=digits), ")")
}

safe_method_names <- function(x){
  x[x=="HER"] <- "HER-p"
  return(x)
}

### Plots
method_order <- c("HER-p", "HER-np", "BER", "UNO", "ITS")

roc_cols <- c("HER-p"  = "black",
              "HER-np" = "black",
              "BER"    = "red",
              "UNO"    = "blue",
              "ITS"    = "blue")

roc_ltys <- c("HER-p"  = "solid",
              "HER-np" = "dashed",
              "BER"    = "dashed",
              "UNO"    = "solid",
              "ITS"    = "dotted")

Perf_tab <- NULL
Phi_tab <- NULL
Runtime_tab <- NULL

for(s in 1:9){

  load(paste0("Result-s", s, ".RData"))

  AUC_by_rep <- apply(AUC, c(1,3), mean, na.rm=TRUE)
  F1_by_rep  <- apply(F1,  c(1,3), mean, na.rm=TRUE)
  DE_by_rep  <- apply(DE,  c(1,3), mean, na.rm=TRUE)

  colnames(AUC_by_rep) <- safe_method_names(colnames(AUC_by_rep))
  colnames(F1_by_rep)  <- safe_method_names(colnames(F1_by_rep))
  colnames(DE_by_rep)  <- safe_method_names(colnames(DE_by_rep))

  tmp_perf <- data.frame(
    Scenario = s,
    Method = colnames(AUC_by_rep),
    AUC_mean = as.numeric(apply(AUC_by_rep, 2, mean, na.rm=TRUE)),
    AUC_sd   = as.numeric(apply(AUC_by_rep, 2, sd,   na.rm=TRUE)),
    F1_mean  = as.numeric(apply(F1_by_rep,  2, mean, na.rm=TRUE)),
    F1_sd    = as.numeric(apply(F1_by_rep,  2, sd,   na.rm=TRUE)),
    Density_Error_mean = as.numeric(apply(DE_by_rep, 2, mean, na.rm=TRUE)),
    Density_Error_sd   = as.numeric(apply(DE_by_rep, 2, sd,   na.rm=TRUE))
  )
  Perf_tab <- rbind(Perf_tab, tmp_perf)

  tmp_runtime <- data.frame(
    Scenario = s,
    Method = colnames(Runtime),
    Runtime_mean_sec = as.numeric(apply(Runtime, 2, mean, na.rm=TRUE)),
    Runtime_sd_sec   = as.numeric(apply(Runtime, 2, sd,   na.rm=TRUE))
  )
  tmp_runtime$Method <- safe_method_names(tmp_runtime$Method)
  Runtime_tab <- rbind(Runtime_tab, tmp_runtime)

  tmp_phi <- data.frame(
    Scenario = s,
    Method = c("HER-p", "HER-np", "BER"),
    Abs_Bias = as.numeric(apply(abs(AB), 3, mean, na.rm=TRUE)),
    Coverage = as.numeric(apply(CP, 3, mean, na.rm=TRUE)),
    Avg_Length = as.numeric(apply(AL, 3, mean, na.rm=TRUE))
  )
  Phi_tab <- rbind(Phi_tab, tmp_phi)
}

write.csv(Perf_tab, file="Performance_Table_AUC_F1_Density.csv", row.names=FALSE)
write.csv(Runtime_tab, file="Runtime_Table.csv", row.names=FALSE)
write.csv(Phi_tab, file="Phi_Table_AB_CP_AL.csv", row.names=FALSE)


###----------------------------------------------------###
### Separate AUC boxplot for each scenario             ###
###----------------------------------------------------###

for(s in 1:9){

  load(paste0("Result-s", s, ".RData"))

  AUC_by_rep <- apply(AUC, c(1,3), mean, na.rm=TRUE)
  colnames(AUC_by_rep) <- safe_method_names(colnames(AUC_by_rep))
  AUC_by_rep <- AUC_by_rep[, method_order, drop=FALSE]

  df <- as.data.frame(AUC_by_rep)
  df_long <- melt(df, id.vars=NULL, variable.name="Method", value.name="AUC")
  df_long$Method <- factor(df_long$Method, levels=method_order)

  p_auc <- ggplot(df_long, aes(x=Method, y=AUC, fill=Method)) +
    geom_boxplot(coef=10) +
    theme_minimal() +
    labs(title=paste0("Scenario ", s), x="", y="AUC") +
    scale_fill_brewer(palette="Set1") +
    guides(fill="none") +
    theme(
      plot.title=element_text(hjust=0.5, size=15, face="bold"),
      axis.text.x=element_text(size=12, angle=30, hjust=1),
      axis.title.y=element_text(size=13)
    )

  ggsave(paste0("AUC_boxplot_scenario_", s, ".pdf"),
         p_auc, width=6, height=4.5)
}

###----------------------------------------------------###
### F1 boxplots for 9 scenarios                        ###
###----------------------------------------------------###

F1_result <- list()
pp_f1 <- list()

for(s in 1:9){
  load(paste0("Result-s", s, ".RData"))
  F1_result[[s]] <- apply(F1, c(1,3), mean, na.rm=TRUE)
  colnames(F1_result[[s]]) <- safe_method_names(colnames(F1_result[[s]]))
  F1_result[[s]] <- F1_result[[s]][, method_order, drop=FALSE]

  df <- as.data.frame(F1_result[[s]])
  df_long <- melt(df, id.vars=NULL, variable.name="Method", value.name="F1")
  df_long$Method <- factor(df_long$Method, levels=method_order)

  pp_f1[[s]] <- ggplot(df_long, aes(x=Method, y=F1, fill=Method)) +
    geom_boxplot(coef=10) +
    theme_minimal() +
    labs(title=paste0("Scenario ", s), x="", y="F1 score") +
    scale_fill_brewer(palette="Set1") +
    guides(fill="none") +
    theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
          axis.text.x=element_text(size=9, angle=30, hjust=1))
}
combined_f1_plot <- grid.arrange(grobs=pp_f1, nrow=3)
ggsave("F1_boxplots_9_scenarios.pdf", combined_f1_plot, width=12, height=9)

###----------------------------------------------------###
### Density error boxplots for 9 scenarios             ###
###----------------------------------------------------###

DE_result <- list()
pp_de <- list()

for(s in 1:9){
  load(paste0("Result-s", s, ".RData"))
  DE_result[[s]] <- apply(DE, c(1,3), mean, na.rm=TRUE)
  colnames(DE_result[[s]]) <- safe_method_names(colnames(DE_result[[s]]))
  DE_result[[s]] <- DE_result[[s]][, method_order, drop=FALSE]

  df <- as.data.frame(DE_result[[s]])
  df_long <- melt(df, id.vars=NULL, variable.name="Method", value.name="Density_Error")
  df_long$Method <- factor(df_long$Method, levels=method_order)

  pp_de[[s]] <- ggplot(df_long, aes(x=Method, y=Density_Error, fill=Method)) +
    geom_boxplot(coef=10) +
    theme_minimal() +
    labs(title=paste0("Scenario ", s), x="", y="Density error") +
    scale_fill_brewer(palette="Set1") +
    guides(fill="none") +
    theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
          axis.text.x=element_text(size=9, angle=30, hjust=1))
}
combined_de_plot <- grid.arrange(grobs=pp_de, nrow=3)
ggsave("Density_Error_boxplots_9_scenarios.pdf", combined_de_plot, width=12, height=9)


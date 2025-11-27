rm(list=ls())
library(ggplot2)
library(reshape2)
library(gridExtra)

## Figure (boxplot)
AUC_result <- list()

for(s in 1:4){
  load( paste0("Result-s", s, ".RData") )
  AUC_result[[s]] <- apply(AUC, c(1,3), mean)
  dimnames(AUC_result[[s]])[[2]][1] <- "HER-p"
}


pp <- list()
for(s in 1:4){
  df <- as.data.frame(AUC_result[[s]])
  df_long <- melt(df, id.vars=NULL, variable.name="Category", value.name="Value")
  
  pp[[s]] <- ggplot(df_long, aes(x=Category, y=Value, fill=Category)) +
    geom_boxplot(coef=10) +  
    theme_minimal() +
    labs(title=paste0("Scenario ", s), x="", y="AUC") +
    scale_fill_brewer(palette="Set1") + 
    guides(fill="none") + 
    theme(
      plot.title=element_text(hjust=0.5, size=15, face="bold"),
      axis.text.x=element_text(size=12)
    )
}


combined_plot <- grid.arrange(pp[[1]], pp[[2]], pp[[3]], pp[[4]], nrow=2) 

# save 
ggsave("AUC-plot.pdf", combined_plot, width=12, height=7)



## Table 
tab <- NULL

load("Result-s1.RData")
tab <- rbind(tab, c(apply(abs(AB), 3, mean), apply(CP, 3, mean), apply(AL, 3, mean)))

load("Result-s2.RData")
tab <- rbind(tab, c(apply(abs(AB), 3, mean), apply(CP, 3, mean), apply(AL, 3, mean)))

load("Result-s3.RData")
tab <- rbind(tab, c(apply(abs(AB), 3, mean), apply(CP, 3, mean), apply(AL, 3, mean)))

load("Result-s4.RData")
tab <- rbind(tab, c(apply(abs(AB), 3, mean), apply(CP, 3, mean), apply(AL, 3, mean)))


tab <- tab[,c(1,4,7,2,5,8,3,6,9)]
write.csv(100*tab, file="tab.csv")






library(tidyverse)
library(latex2exp)

#Generate Importance plot of threshold
#####
load('Results/Control_Strategies/urchin_threshold_GSA.RData')
threshold_GSA <-
  data.frame(Parameter = rownames(RF$importance),
             Importance = RF$importance[, 1])
threshold_GSA$Parameter <-
  factor(threshold_GSA$Parameter, levels = threshold_GSA$Parameter[order(threshold_GSA$Importance,decreasing = TRUE)])
threshold_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16)
  )+scale_x_discrete(labels=unname(TeX(c("$\\gamma_A$","$\\sigma_A$","$\\delta_U$","$a_0$","$\\delta_A$","$a_A$","$\\gamma_S$","$\\beta$","$L$","$R$","$\\gamma_G$","$\\epsilon$","$a_U$"))))

#Generate Importance plot of spread
#####
load('Results/Control_Strategies/spread_rate_GSA.RData')
spread_GSA <-
  data.frame(Parameter = rownames(RF$importance),
             Importance = RF$importance[, 1])
spread_GSA$Parameter <-
  factor(spread_GSA$Parameter, levels = spread_GSA$Parameter[order(spread_GSA$Importance,decreasing = TRUE)])
spread_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16)
  )+scale_x_discrete(labels=unname(TeX(c("$\\sigma_A$","$\\delta_A$","$a_0$","$u_0$","$\\gamma_A$","$a_A$","$L$","$\\mu_u$","$\\delta_U$","$R$","$\\gamma_S$","$\\mu_A$","$a_U$","$\\epsilon$","$\\beta$","$\\eta$","$\\mu_S$","$A_c$","$\\gamma_G$"))))
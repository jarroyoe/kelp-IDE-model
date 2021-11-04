library(tidyverse)
library(ggridges)
library(here)
library(Rmisc)

load(paste(sep='',here::here(),'/Results/ABC/final_distributions.RData'))

#Plot densities for all parameters
allEstimates %>% mutate(Region=ifelse(Region==1,'Little River','Timber Cove')) %>%
  ggplot(aes(x = NormalizedValue, y = Parameter, fill = Parameter)) +
  geom_density_ridges() + facet_wrap( ~ Region) + theme_bw() + theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size = 18),
    legend.position = "none"
  )
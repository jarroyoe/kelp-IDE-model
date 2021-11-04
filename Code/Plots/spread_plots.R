library(tidyverse)
library(Rmisc)
library(ggpubr)
library(plotly)

#Generate plot of threshold with single strategies
#####
load('Results/Control_Strategies/critical_U_control_strategies_1_1.Rdata')
singleStrategies <- c(2, 3, 4, 5, 9, 13, 17, 33, 49)
basethreshold <- criticalUs[1]
threshold_with_replantation <-
  data.frame(Sigma=rep('+10%',4),Intensity = c(0,weak, mid, strong),
                                               Threshold = c(basethreshold,criticalUs[singleStrategies[7:9]]))
load('Results/Control_Strategies/critical_U_control_strategies.Rdata')
singleStrategies <- c(2, 3, 4, 5, 9, 13, 17, 33, 49)
basethreshold <- criticalUs[1]
threshold_with_replantation <-
  rbind(threshold_with_replantation,data.frame(Sigma=rep('Baseline',4),Intensity = c(0,weak, mid, strong),
             Threshold = c(basethreshold,criticalUs[singleStrategies[7:9]])))
load('Results/Control_Strategies/critical_U_control_strategies_0_9.Rdata')
singleStrategies <- c(2, 3, 4, 5, 9, 13, 17, 33, 49)
basethreshold <- criticalUs[1]
threshold_with_replantation <-
  rbind(threshold_with_replantation,data.frame(Sigma=rep('-10%',4),Intensity = c(0,weak, mid, strong),
             Threshold = c(basethreshold,criticalUs[singleStrategies[7:9]])))
threshold_with_replantation %>% ggplot(aes(Intensity, Threshold/60,colour=Sigma)) +
  geom_line() + geom_point() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),legend.position = c(0.9,0.2),
    text= element_text(size=14)
  )+xlab('Kelp Outplanting Intensity')+ylab("Threshold urchin density (urchins/m²)")

#Generate single strategies versus spread rate
#####
load('Results/Control_Strategies/spread_rate_calculator_eta_1.Rdata')
singleStrategies <- c(1,2, 3, 4,1, 5, 9, 13,1, 17, 33, 49)
#singleStrategies <- c(1,5,9,13)
baseSpread <- spreadRates[1]
singleStrategies_spread <- data.frame(Spatial_Extent=rep('Short',12),Strategy=rep(c('Urchin harvesting','Seed outplanting','Kelp outplanting'),each=4),
                                      Intensity=c(0,rep(c(weak,mid,strong),n=3)),Spread_Rate=spreadRates[singleStrategies])

#singleStrategies_spread <- data.frame(Spatial_Extent=rep('Short',4),Strategy=rep(c('Seed outplanting'),each=4),
#                                      Intensity=c(0,rep(c(weak,mid,strong),n=1)),Spread_Rate=spreadRates[singleStrategies])
load('Results/Control_Strategies/spread_rate_calculator_eta_10.Rdata')
#singleStrategies_spread <- rbind(singleStrategies_spread,data.frame(Spatial_Extent=rep('Wide',4),Strategy=rep(c('Seed outplanting'),each=4),
#                                      Intensity=c(0,rep(c(weak,mid,strong),n=1)),Spread_Rate=spreadRates[singleStrategies]))
singleStrategies_spread <- rbind(singleStrategies_spread,data.frame(Spatial_Extent=rep('Wide',12),Strategy=rep(c('Urchin harvesting','Seed outplanting','Kelp outplanting'),each=4),
                                      Intensity=c(0,rep(c(weak,mid,strong),n=3)),Spread_Rate=spreadRates[singleStrategies]))

p1 <- singleStrategies_spread %>% filter(Spatial_Extent=='Short') %>% ggplot(aes(Intensity,Spread_Rate,colour=Strategy,shape=Strategy))+
  geom_point()+geom_line() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank(),
    text=element_text(size=12)
  )+ylab("Spread Rate")+xlab("Ongoing restoration intensity")+ggtitle("a) Control near the oasis")+ylim(NA,1.1*max(singleStrategies_spread$Spread_Rate))

  p2 <- singleStrategies_spread %>% filter(Spatial_Extent=='Wide') %>% ggplot(aes(Intensity,Spread_Rate,colour=Strategy,shape=Strategy))+
    geom_point()+geom_line()+theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.y = element_blank(),
      axis.title.y=element_blank(),
      axis.ticks.y = element_blank(),
      plot.background=element_blank(),
      text=element_text(size=12)
    )+ggtitle("b) Control farther of the oasis")+xlab("Ongoing restoration intensity")+ylim(NA,1.1*max(singleStrategies_spread$Spread_Rate))
  ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")

#####  
#Generate 2d plot of strategies versus spread
load('Results/Control_Strategies/spread_rate_calculator_eta_1.RData')
combinedUrchin_Kelp <- c(18,19,20,34,35,36,50,51,52)
combinedUrchin_Seed <- c(6,7,8,10,11,12,14,15,16)
combinedSeed_Kelp <- c(21,25,29,37,41,45,53,57,61)

urchin_kelp <- data.frame(Urchin_Harvesting=unlist(controlStrategies[combinedUrchin_Kelp])[which(1:27 %% 3==1)],
                          Kelp_Outplanting=unlist(controlStrategies[combinedUrchin_Kelp])[which(1:27 %% 3==0)],
                          Spread_Rate=spreadRates[combinedUrchin_Kelp])
urchin_kelp %>% ggplot(aes(x=Urchin_Harvesting,y=Kelp_Outplanting))+geom_tile(aes(fill=Spread_Rate))+
  scale_fill_distiller(palette = "YlGnBu")+theme_bw()+theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+xlab('Urchin Harvesting Intensity')+ylab('Kelp Outplanting Intensity')+labs(fill='Spread Rate')+ggtitle('a)')

urchin_seed <- data.frame(Urchin_Harvesting=unlist(controlStrategies[combinedUrchin_Seed])[which(1:27 %% 3==1)],
                          Seed_Outplanting=unlist(controlStrategies[combinedUrchin_Seed])[which(1:27 %% 3==2)],
                          Spread_Rate=spreadRates[combinedUrchin_Seed])
urchin_seed %>% ggplot(aes(x=Urchin_Harvesting,y=Seed_Outplanting))+geom_tile(aes(fill=Spread_Rate))+
  scale_fill_distiller(palette = "YlGnBu")+theme_bw()+theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+xlab('Urchin Harvesting Intensity')+ylab('Seed Outplanting Intensity')+labs(fill='Spread Rate')+ggtitle('b)')

seed_kelp <- data.frame(Seed_Outplanting=unlist(controlStrategies[combinedSeed_Kelp])[which(1:27 %% 3==2)],
                          Kelp_Outplanting=unlist(controlStrategies[combinedSeed_Kelp])[which(1:27 %% 3==0)],
                          Spread_Rate=spreadRates[combinedSeed_Kelp])
seed_kelp %>% ggplot(aes(x=Seed_Outplanting,y=Kelp_Outplanting))+geom_tile(aes(fill=Spread_Rate))+
  scale_fill_distiller(palette = "YlGnBu")+theme_bw()+theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+xlab('Seed Outplanting Intensity')+ylab('Kelp Outplanting Intensity')+labs(fill='Spread Rate')+ggtitle('c)')

#Generate 3d plot of strategies versus spread  
#####
  unlistedStrategies<-unlist(controlStrategies)
  x <- unlistedStrategies[which(1:192%%3 ==1)]
  y <- unlistedStrategies[which(1:192%%3 ==2)]
  z <- unlistedStrategies[which(1:192%%3 ==0)]
  plot_ly(x=x,y=y,z=z,type='isosurface',value=spreadRates) %>% 
    layout(title='Spread rate of kelps at different control strategies intensities',
           scene=list(xaxis=list(title='Urchin harvesting'),
                      yaxis=list(title='Seed plantation'),
                      zaxis=list(title='Kelp replantation')),
           legend=list(title=list(text='Spread rate')))
#####
#Generate spread vs initial conditions plots
load("Results/Control_Strategies/spread_rate_kelp_densities.Rdata")
kelpDensitiesSpread <-
  data.frame(
    A0 = kelpDensities,
    spread = spreadRates,
    Control = rep("No ongoing control",length(kelpDensities))
  )
load("Results/Control_Strategies/spread_rate_kelp_densities_strong.Rdata")
kelpDensitiesSpread <-
  rbind(kelpDensitiesSpread,data.frame(
    A0 = kelpDensities,
    spread = spreadRates,
    Control = rep("Strong ongoing control",length(kelpDensities))
  ))
load("Results/Control_Strategies/spread_rate_kelp_densities_weak.Rdata")
kelpDensitiesSpread <-
  rbind(kelpDensitiesSpread,data.frame(
    A0 = kelpDensities,
    spread = spreadRates,
    Control = rep("Weak ongoing control",length(kelpDensities))
  ))

p1 <- ggplot(kelpDensitiesSpread,aes(A0,spread,colour=Control,shape=Control))+geom_line()+geom_point()+theme_bw()+ theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  plot.background = element_blank(),text=element_text(size=12)
)+ylab("Spread Rate")+xlab("Initial kelp density")+ylab("Spread Rate")+ggtitle("a)")+ylim(0,1.1*max(kelpDensitiesSpread$spread))

load("Results/Control_Strategies/spread_rate_urchin_densities.Rdata")
urchinDensitiesSpread <-
  data.frame(
    U0 = uthresholds,
    spread = spreadRates,
    Control = rep("No ongoing control",length(uthresholds))
  )
load("Results/Control_Strategies/spread_rate_urchin_densities_strong.Rdata")
urchinDensitiesSpread <-
  rbind(urchinDensitiesSpread,data.frame(
    U0 = uthresholds,
    spread = spreadRates,
    Control = rep("Strong ongoing control",length(uthresholds))
  ))
load("Results/Control_Strategies/spread_rate_urchin_densities_weak.Rdata")
urchinDensitiesSpread <-
  rbind(urchinDensitiesSpread,data.frame(
    U0 = uthresholds,
    spread = spreadRates,
    Control = rep("Weak ongoing control",length(uthresholds))
  ))
p2 <- ggplot(urchinDensitiesSpread,aes(U0,spread,colour=Control,shape=Control))+geom_line()+geom_point()+theme_bw()+theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text.y = element_blank(),
  axis.title.y=element_blank(),
  axis.ticks.y = element_blank(),
  plot.background=element_blank(),
  text=element_text(size=12)
)+ylab("Spread Rate")+ggtitle("b)")+xlab("Initial urchin density")+ylab("Spread Rate")+ylim(0,1.1*max(kelpDensitiesSpread$spread))

ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")
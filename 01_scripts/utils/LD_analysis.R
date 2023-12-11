# Script produced by Florent Sylvestre and adapted from https://stats.stackexchange.com/questions/437477/calculate-accelerated-bootstrap-interval-in-r
# Script estimating the Correction bootstrap with acceleration 

library(data.table)
library(tidyverse)
library(viridis)
library(patchwork)


bootstrap <- function(data, N, func){
  bootmetric <- rep(NA, N)
  print(length(data))
  if(length(data) > 150000){data <- sample(data, 150000, replace =F)}
  for(i in 1:N){
    bootmetric[i] <- func(sample(data, length(data),replace = TRUE))
  }
  return(bootmetric)
}


BCB <- function(data, N.boot, func, alpha = 0.05){
  n <- length(data)
  u <- c(alpha/2, 1 - alpha/2)
  theta_hat <- func(data)
  theta_boot <- bootstrap(data, N.boot, func)
  z0 <- qnorm(mean(theta_boot <= theta_hat))
  zu <- qnorm(u)
  
  return(quantile(theta_boot, pnorm(z0 + (z0+zu))))
}


## File import
SV_SNP_LD <- fread("01_ld_files/SV_SNPs_concat_100_allchr.ld", h=T) %>% as.data.frame
SNP_SNP_LD <- fread("01_ld_files/SNPs_SNPs_concat_100_allchr.ld", h= T) %>% as.data.frame
n.boot <- 1000

## Main script
# Ajusting needed column:

SV_SNP_LD %<>% select(-c(SNP_A, SNP_B,POS,CAND_SVLEN)) %>%
               rename("TYPE" = CAND_SVTYPE) %>%
               mutate(across(c(BP_A, BP_B, R2), as.numeric)) %>%
               mutate(DIST =  pmin(abs(BP_A - BP_B),abs(BP_A - END))) %>%
               select(-END)

SNP_SNP_LD %<>% select(-c(SNP_A, SNP_B)) %>%
                mutate(across(c(BP_A, BP_B, R2), as.numeric)) %>%
                mutate(TYPE = "SNP", DIST = abs(BP_A - BP_B))


data <- rbind(SV_SNP_LD,SNP_SNP_LD) %>% as.data.frame
rm(SV_SNP_LD)
rm(SNP_SNP_LD)


R2_large <- list()
for(type in c("DEL", "INS","SNP")){
  mean = tapply(data$R2[data$TYPE == type],
                data$DIST[data$TYPE == type]%/%100 ,
                median)
  print("Start IC")
  IC =  tapply(data$R2[data$TYPE == type],
               data$DIST[data$TYPE == type]%/%100,
               BCB,
               n.boot,
               median)
  print("IC_done")
  IC <- do.call(rbind, IC)


  R2_large[[type]] = data.frame(TYPE = type,
                          DIST = as.numeric(names(mean))*100,
                          MEAN = mean,
                          LOWBOUND = IC[,1],
                          UPBOUND = IC[,2])
}
R2_large <- do.call(rbind, R2_large)
print("Zoom")
R2_zoom <- list()
for(type in c("DEL", "INS","SNP")){
  print(type)
  sub <- data %>% filter(DIST <=1000)
  mean = tapply(sub$R2[sub$TYPE == type],
                sub$DIST[sub$TYPE == type]%/%5 ,
                median)

  IC =  tapply(sub$R2[sub$TYPE == type],
               sub$DIST[sub$TYPE == type]%/%5,
               BCB,
               n.boot,
               median)

  IC <- do.call(rbind, IC)


  R2_zoom[[type]] = data.frame(TYPE = type,
                          DIST = as.numeric(names(mean))*5,
                          MEAN = mean,
                          LOWBOUND = IC[,1],
                          UPBOUND = IC[,2])
}
R2_zoom <- do.call(rbind, R2_zoom)

print("plots")
###Plots

## Assign a color to each svtype in a named vector
cols_svtypes <- c(DEL = "#440154FF", DUP = "#3B528BFF", INS = "#21908CFF", INV = "#5DC863FF", SNP = "#FDE725FF")
large <- R2_large %>%
  ggplot(aes(x = DIST)) +
  geom_ribbon(aes(ymin = UPBOUND, ymax = LOWBOUND, fill = TYPE), col = NA, alpha = 0.3) +
  geom_line(aes( y = MEAN, col = TYPE)) +
  theme_bw() +
  scale_color_manual(values = cols_svtypes) +
  scale_fill_manual(values = cols_svtypes)+
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_y_continuous(limits = c(0,0.1))+
  xlab("Distance (bp)") +
  ylab(expression(paste("Median ", r^2)))

zoom <- R2_zoom %>%
  ggplot(aes(x = DIST)) +
  geom_ribbon(aes(ymin = UPBOUND, ymax = LOWBOUND, fill = TYPE),col = NA, alpha = 0.3) +
  geom_line(aes( y = MEAN, col = TYPE)) +
  theme_bw() +
  scale_color_manual(values = cols_svtypes) +
  scale_fill_manual(values = cols_svtypes) +
  scale_y_continuous(limits = c(0,0.1)) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))+
  xlab("Distance (bp)") +
  ylab(expression(paste("Median ", r^2)))

p <- large + inset_element(zoom, 0.6,0.65,.99,.99, align_to = "plot")+
  plot_layout(guides = 'collect')

##output:
write.table(R2_large,
            "02_results/Large_r2_summary_table",
            row.names = F,
            quote = F)
write.table(R2_zoom,
            "02_results/Zoom_r2_summary_table",
            row.names = F,
            quote = F)

ggsave(p,
       file = "02_results/LD_decay.pdf",
       device = "pdf",
       dpi = 700)


ggsave(p,
       file = "02_results/LD_decay.png",
       device = "png",
       dpi = 700)







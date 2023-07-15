# Data from Table 3.4.1. of Andrew F. Siegel, 
# in Practical Business Statistics (Seventh Edition), 2016
rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(coda)
library(MASS)
library(tidyverse)
library(readxl)
library(gridExtra)
library(latex2exp)
library(parallel)
library(HDInterval)
color_scheme_set("brightblue")

## data import
df1 <- read.csv("Deposits_of_Banks_and_Savings_Institutions.csv")
colnames(df1)[2] <- "Deposits($ billions)"

KDE_list <- density(df1$`Deposits($ billions)`,
                    bw = "SJ",
                    kernel = "gaussian")

df_KDE <- data.frame(x = KDE_list$x, y = KDE_list$y)
p1 <- df_KDE %>%
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    scale_x_continuous(breaks = seq(0,900,100)) +
    xlim(0,max(KDE_list$x)) +
    ylab("Density") +
    xlab("Bank deposits ($billions)") +
    theme_bw() +
    geom_vline(xintercept=mean(df1$`Deposits($ billions)`),
               color="blue",
               linetype="solid") + 
    geom_vline(xintercept=median(df1$`Deposits($ billions)`),
               color="red",
               linetype="dashed") +
    geom_vline(xintercept=KDE_list$x[which.max(KDE_list$y)],
               color="#FFA500",
               linetype="dotted")

p1

## some calculations
mean(df1$`Deposits($ billions)`)
mean(df1$`Deposits($ billions)` < mean(df1$`Deposits($ billions)`))
median(df1$`Deposits($ billions)`)

## data for stan
y <- df1$`Deposits($ billions)`
N <- length(y)
dat <- list(N = N,
            y = y)
stan_file <- "../DTP_t_iid.stan"
stan_mod <- cmdstan_model(stan_file)

## pdf of DTP-student-t
dDTP <- function(x,theta,sigma1,sigma2,delta1,delta2){
  epsilon <- sigma1*dt(0,delta2)/(sigma1*dt(0,delta2)+sigma2*dt(0,delta1))
  p1 <- 2*epsilon/sigma1*dt((x-theta)/sigma1,delta1)*(x<theta)
  p2 <- 2*(1-epsilon)/sigma2*dt((x-theta)/sigma2,delta2)*(x>=theta)
  return(p1+p2)
}

## MCMC
stan_fit <- stan_mod$sample(
  data = dat,
  seed = 100,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  iter_warmup = 10000,
  save_warmup = FALSE,
  iter_sampling = 20000,
)

## parameter estimation 
par.est <- stan_fit$summary()
print(par.est)

## traceplot
df_post <- stan_fit$draws(variables = c("theta","sigma1","sigma2",
                                        "delta1","delta2"), 
                          format = "df")
pdf("traceplot.pdf",height = 4,width = 8)
mcmc_trace(df_post)
dev.off()

## density estimation plot

xaxis <- seq(0,max(df1$`Deposits($ billions)`),length.out=1000)

df_mix <- data.frame(x=xaxis,
                     y=dDTP(xaxis,
                            par.est$mean[6],
                            par.est$mean[2],
                            par.est$mean[3],
                            par.est$mean[4],
                            par.est$mean[5]
                            ))
df_M <- data.frame(x = c(par.est$mean[6],
                         median(df1$`Deposits($ billions)`),
                         mean(df1$`Deposits($ billions)`)),
                   type = c("mode","median","mean"))
df_M$type <- factor(df_M$type,levels = df_M$type)
colnames(df_M)[2] <- "measures of central tendency"

p_out <- df_mix %>%
  ggplot(aes(x=x,y=y)) +
  geom_histogram(data = data.frame(y = df1$`Deposits($ billions)`),
                 aes(x = y, y = after_stat(density)),
                 fill = rgb(0,0,0,0.1),
                 color = rgb(0,0,0,0.3),
                 bins = 85) +
  geom_line(color = "black") +
  theme_bw()+
  ylab("density") +
  xlab("Bank deposits ($billions)") +
  geom_vline(data = df_M,
             aes(xintercept = x,
                 color = `measures of central tendency`,
                 linetype = `measures of central tendency`),
             size = 0.8) +
  scale_linetype_manual(values = c(2,4,1)) +
  scale_color_manual(values = c("red","orange","blue")) +
  scale_x_continuous(breaks = seq(0,800,100),
                     expand = c(0.005,0.005),
                     sec.axis = sec_axis(~.,
                                         breaks = c(round(par.est$mean[6],1),
                                                    round(median(df1$`Deposits($ billions)`),1),
                                                    round(mean(df1$`Deposits($ billions)`),1)
                                                    ))) +
  scale_y_continuous(breaks = seq(0,0.03,0.005)) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(size = 5),
        legend.position = "bottom")

p_out

ggsave("bank_data.pdf",p_out,width = 7,height = 3)


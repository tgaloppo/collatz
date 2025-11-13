# plots the autocorrelation data produced by
# "correlation.cpp"
#
# $ ./correlation > acorr.csv
#
library(ggplot2)
library(dplyr)

ac <- read.csv("acorr.csv")

plt <- ac %>% 
  mutate(nbits=factor(nbits)) %>% 
  group_by(nbits, lag) %>% 
  summarize(mu=mean(autocorrelation), low=mean(autocorrelation)-2*sd(autocorrelation), high=mean(autocorrelation)+2*sd(autocorrelation)) %>% 
  ggplot() + 
    geom_line(aes(x=lag, y=mu, color="Mean")) +
    geom_line(aes(x=lag, y=low, color="95% CI")) +
    geom_line(aes(x=lag, y=high, color="95% CI")) +
    scale_color_manual(name = "Key",
                       values = c("Mean" = "red", "95% CI" = "gray")) +
    scale_y_continuous(limits=c(-2,2)) + 
    facet_wrap(~nbits, scale="free_x")

print(plt)

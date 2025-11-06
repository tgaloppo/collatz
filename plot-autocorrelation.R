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
  summarize(mu=mean(autocorrelation)) %>% 
  ggplot(aes(x=lag, y=mu)) + 
    geom_line() + 
    facet_wrap(~nbits, scale="free_x")

print(plt)

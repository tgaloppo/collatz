# plots the fractional positions of iterates
# as computed with "fractional-position.cpp"
# and overlays Benford distribution
#
# $ ./fractional-position > fracpos.csv
#
library(ggplot2)
library(dplyr)

dbenford <- function(a) 1 / (a * log(2))

fp <- read.csv("fractional-position.csv")

plt <- fp %>% mutate(nbits=factor(nbits)) %>% 
  ggplot(aes(x=alpha)) + 
    geom_histogram(aes(y = after_stat(density))) + 
    stat_function(fun=dbenford, color="red", linewidth=0.5) + 
    facet_wrap(~ nbits)

ggsave("benford-dist.eps", plot = plt)

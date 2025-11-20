# create semi-log plot of LSB Shift
# (Delta L) frequencies observed
# during single 1M bit trajectory.
#
# $ ./delta-L-freq > delta-L-freq.csv
#
library(ggplot2)
library(dplyr)

data <- read.csv("delta-L-freq.csv")

plt <- data %>% 
        ggplot(aes(x=dL, y=Count)) + 
        geom_line() + 
        geom_point() + 
        scale_y_log10() +
        xlab("LSB Shift (Delta L)") +
        ylab("log10( Count )")

ggsave("delta-L-freq.eps", plot = plt)

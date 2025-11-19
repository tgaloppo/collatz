# plots the Ulam estimate of the spectral
# gap of the transfer operator of the map
# U(a) for different input sizes.
#
# $ ./ulam-gap > ulam-gap-out.csv
#
library(ggplot2)
library(dplyr)

data <- read.csv("ulam-gap-out.csv", header=F)
names(data) <- c("InputSizeLog10", "Lambda0", "Lambda1", "Lambda2", "Gap", "GapLog10")

plt <- data %>% 
        ggplot(aes(x=InputSizeLog10, y=GapLog10)) + 
          geom_point() +
          geom_line() +
          scale_x_continuous(limits=c(1, 7)) +
          scale_y_continuous(limits=c(-7, -1)) +
          xlab("Input Size (log10)") +
          ylab("Spectral Gap (log10)")

ggsave("spectral-gap.eps", plot = plt)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

df = read.table(args[1], header=T)
df$chrom = factor(df$chrom, levels=c(1:22,"X","Y"))
library(ggplot2)
ggplot(df, aes(y=counts, x=chrom, fill=sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~sv, nrow=3)


library(dplyr)
df %>% group_by(sample, sv) %>% summarize(total=sum(counts))

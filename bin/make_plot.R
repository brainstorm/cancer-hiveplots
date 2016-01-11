df = read.table("~/repos/pipelines/glasgow/sv/res/lumpy.tsv", header=T)
df$chrom = factor(df$chrom, levels=c(1:22,"X","Y"))
library(ggplot2)
ggplot(df, aes(y=counts, x=chrom, fill=sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~sv, nrow=3)


library(dplyr)
df %>% group_by(sample, sv) %>% summarize(total=sum(counts))

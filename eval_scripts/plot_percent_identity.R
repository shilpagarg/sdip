library(tidyverse)

options(scipen=5)

contigs <- read_tsv("/project/heller_d-data/shilpa/CHM13/regions/contigs/pooled.t5.b5.d2.polished.fa.fai", col_names = FALSE)
us.hg38 <- read_tsv("/project/heller_d-data/shilpa/CHM13/regions/eval/t5.b5.d2/tables/polished.to.hg38.tbl")
us.bacs <- read_tsv("/project/heller_d-data/shilpa/CHM13/regions/eval/t5.b5.d2/tables/polished.to.bacs.tbl")

sda.hg38 <- read_tsv("/project/heller_d-data/shilpa/CHM13/regions/eval/sda/tables/polished.to.hg38.tbl")
sda.bacs <- read_tsv("/project/heller_d-data/shilpa/CHM13/regions/eval/sda/tables/polished.to.bacs.tbl")

resolved_hg38 <- hg38[hg38$perID_by_all >= 99.8 & hg38$query_start == 0 & (hg38$query_length - hg38$query_end == 0), "query_name"]$query_name
resolved_bacs <- bacs[bacs$perID_by_all >= 99.8 & (bacs$query_start == 0 | bacs$reference_start == 0) & (bacs$query_length - bacs$query_end == 0 | bacs$reference_length - bacs$reference_end == 0), "query_name"]$query_name
length(resolved_hg38)
sum(hg38[hg38$perID_by_all >= 99.8 & hg38$query_start == 0 & (hg38$query_length - hg38$query_end == 0), "query_length"])
length(setdiff(contigs$X1, resolved_hg38))
length(setdiff(resolved_bacs, resolved_hg38))
length(union(resolved_bacs, resolved_hg38))

aligned_hg38 <- hg38[hg38$perID_by_all >= 97.5 & hg38$query_start == 0 & (hg38$query_length - hg38$query_end == 0), "query_name"]$query_name
aligned_bacs <- bacs[bacs$perID_by_all >= 97.5 & (bacs$query_start == 0 | bacs$reference_start == 0) & (bacs$query_length - bacs$query_end == 0 | bacs$reference_length - bacs$reference_end == 0), "query_name"]$query_name
length(union(aligned_hg38, aligned_bacs))

data <- data.frame(ids = seq(97.5, 100, 0.05))
cumlength.us.hg38 <- c()
cumlength.sda.hg38 <- c()
cumlength.us.bacs <- c()
cumlength.sda.bacs <- c()
for (id in data$ids) {
  cumlength.us.hg38 <- c(cumlength.us.hg38, sum(us.hg38[us.hg38$perID_by_all <= id, "query_length"]))
  cumlength.sda.hg38 <- c(cumlength.sda.hg38, sum(sda.hg38[sda.hg38$perID_by_all <= id, "query_length"]))
  cumlength.us.bacs <- c(cumlength.us.bacs, sum(us.bacs[us.bacs$perID_by_all <= id, "query_length"]))
  cumlength.sda.bacs <- c(cumlength.sda.bacs, sum(sda.bacs[sda.bacs$perID_by_all <= id, "query_length"]))
}
data$cumlength.us.hg38 <- cumlength.us.hg38
data$cumlength.sda.hg38 <- cumlength.sda.hg38
data$cumlength.us.bacs <- cumlength.us.bacs
data$cumlength.sda.bacs <- cumlength.sda.bacs
data$matched <- data$ids > 99

data %>%
  pivot_longer(c(cumlength.us.hg38, cumlength.sda.hg38, cumlength.us.bacs, cumlength.sda.bacs), names_to = "key", values_to = "value") %>%
  ggplot(aes(x=ids, y=(value / 1000000), color=key)) +
    geom_line(size=1.5) +
    labs(x = "Best percent identity match", y = "Cumulative contig size in Mb", color = "Contigs / Reference")


ggplot(data, aes(x=ids)) +
  geom_area(aes(y=cumlength.us.bacs), fill="grey") +
  geom_area(aes(y=cumlength.sda.bacs), fill="black")

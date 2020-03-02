library(tidyverse)
u <- read_tsv("/project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.ref.bed", col_names=F)
r <- read_tsv("/project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.covered.ref.bed", col_names=F)
s <- read_tsv("/project/heller_d-data/shilpa/NA12878/resolved_ucsc/unmerged/segdups.uncovered.covered.sdip.bed", col_names=F)

u$set <- "Unresolved by base assembly"
r$set <- "Resolved by base assembly"
s$set <- "Resolved by SDip"
all <- rbind(u, s) %>% mutate(length=X3-X2, identity=X6*100)
all$set <- factor(all$set, levels=c("Unresolved by base assembly", "Resolved by SDip"))

png("/project/pacbiosv/code/sdip/eval_scripts/pngs/fig2a.png", width=1000, height=1000, res=200)

all %>% ggplot(aes(x=length/1000, y=identity, color=set)) +
geom_density_2d() +
scale_x_log10() +
labs(x="Segmental duplication length (kbps)", y="Percent sequence identity", colour="Subset") +
facet_grid(. ~ set) +
scale_color_discrete(guide=FALSE)


dev.off()

library(VennDiagram)
df <- read.delim("/project/heller_d-data/shilpa/NA12878/regions/svim/merged.sdip.sda.bacs.txt", sep="")
venn.plot <- venn.diagram(list(SDip=which(df[,1]==1), SDA=which(df[,2]==1), BACs=which(df[,3]==1)),
# Circles
lwd = 2,
lty = 'blank',
fill = c("red", "blue" ,"green"),
alpha = c(0.5, 0.5, 0.5),
# Numbers
cex = 2,
fontface = "bold",
fontfamily = "sans",
# Set names
cat.cex = 2,
cat.fontface = "bold",
cat.default.pos = "outer",
cat.pos = c(-27, 27, 135),
cat.dist = c(0.055, 0.055, 0.085),
cat.fontfamily = "sans",
rotation = 1,
filename = NULL);

png("/project/pacbiosv/code/sdip/eval_scripts/pngs/fig2b.png", width=1000, height=1000, pointsize=24)
grid.draw(venn.plot)
dev.off()

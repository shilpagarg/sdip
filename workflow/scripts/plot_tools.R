library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "parameters", "metric", "value"))
res$caller = factor(res$caller, levels=c('sniffles', 'sdip'), labels=c('Sniffles', 'SDip'))
res$metric = factor(res$metric, levels=c('precision', 'recall'), labels=c('precision', 'recall'))


res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    ggplot(aes(recall, precision, color=caller)) +
      geom_point(size=1.0) +
      geom_path() +
      scale_shape_manual(values=c(16,17)) +
      scale_color_manual(values=c("goldenrod2", "firebrick2")) +
      labs(y = "Precision", x = "Recall", color = "Tool") +
      lims(x=c(0,100), y=c(0,100)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))


ggsave(args[2], width=5, height=4)

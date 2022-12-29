rm(list=ls())

library(here)
library(choros)
library(ggplot2)

data_dir <- file.path(here(), "data", "meydan_2020")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "disome", "disome_codon_corr.Rda"))

# generate plot -----------------------------------------------------------

disome_codon_corr <- data.frame(position=c(names(disome_codon_corr$count),
                                           names(disome_codon_corr$correct_250)),
                                delta_corr=c(disome_codon_corr$count,
                                             disome_codon_corr$correct_250),
                                type=c(rep("Raw counts", length(disome_codon_corr$count)),
                                       rep("Corrected counts"), length(disome_codon_corr$correct_250)))
disome_codon_corr$position <- sub("^n", "-", disome_codon_corr$position)
disome_codon_corr$position <- sub("^p", "", disome_codon_corr$position)
disome_codon_corr$position <- factor(disome_codon_corr$position,
                                     levels=unique(disome_codon_corr$position))
disome_codon_corr$fill <- "other"
disome_codon_corr$fill[grepl("^A", as.character(disome_codon_corr$position))] <- "A"
disome_codon_corr$fill[grepl("^P", as.character(disome_codon_corr$position))] <- "P"
disome_codon_corr$fill[grepl("^E", as.character(disome_codon_corr$position))] <- "E"
disome_codon_corr$fill[disome_codon_corr$position %in% c("-16", "-15", "3", "4")] <- "bias"
disome_codon_corr$type <- factor(disome_codon_corr$type, levels=unique(disome_codon_corr$type))

figure_6D <- ggplot(disome_codon_corr, aes(x=position, y=delta_corr, fill=fill)) +
  geom_col() + theme_classic(base_size=8) + facet_wrap(vars(type), nrow=2) +
  scale_fill_manual(values=setNames(c(RColorBrewer::brewer.pal(4, "Set1"), "grey"),
                                    c("bias", "E", "P", "A", "other"))) +
  theme(legend.position="none") +
  xlab("codon position") + ylab(expression(Delta*" correlation"))

ggsave(filename=file.path(figures_dir, "figure_7D.pdf"),
       plot=figure_6D, device="pdf", width=6.25, height=2.5)
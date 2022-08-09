rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data", "simulated_data")
figures_dir <- file.path(here(), "figures")

expts <- c("no_bias", "f3_bias", "f5_bias", "f3_f5_bias")

# load data ---------------------------------------------------------------

for(expt in expts) {
  load(file.path(data_dir, expt, paste0(expt, "_codon_corr.Rda")))
}

# make plot ---------------------------------------------------------------

raw_codon_corr <- data.frame(position=unlist(lapply(expts,
                                                    function(x) {
                                                      names(get(paste0(x, "_codon_corr"))[[1]])
                                                    })),
                             codon_corr=unlist(lapply(expts,
                                                      function(x) {
                                                        get(paste0(x, "_codon_corr"))[[1]]
                                                      })),
                             expt=unlist(lapply(expts,
                                                function(x) {
                                                  rep(x, length(get(paste0(x, "_codon_corr"))[[1]]))
                                                })))
raw_codon_corr$type <- "Raw counts"
raw_codon_corr$expt <- factor(raw_codon_corr$expt, levels=expts)
levels(raw_codon_corr$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")
raw_codon_corr$position <- sub("n", "-", raw_codon_corr$position)
raw_codon_corr$position <- sub("p", "", raw_codon_corr$position)
raw_codon_corr$position <- factor(raw_codon_corr$position, 
                                  levels=unique(raw_codon_corr$position))
raw_codon_corr$label <- as.character(raw_codon_corr$position)
raw_codon_corr$label[raw_codon_corr$label %in% c("-5", "-4", "3", "4")] <- "bias"
raw_codon_corr$label[!(raw_codon_corr$label %in% c("A", "P", "E", "bias"))] <- "other"
fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")
names(fill_colors) <- c("bias", "E", "P", "A", "other")

figure_2B <- ggplot(raw_codon_corr, aes(x=position, y=codon_corr, fill=label)) + 
  geom_col() + scale_fill_manual(values=fill_colors) + theme_classic(base_size=8) + 
  theme(legend.position="none", panel.spacing=unit(0.25, "in"),
        axis.text.x=element_text(size=5)) + 
  facet_grid_sc(rows=vars(expt), cols=vars(type)) + 
  xlab("codon position") + ylab(expression(Delta*" correlation"))

ggsave(filename=file.path(figures_dir, "figure_2B.pdf"),
       plot=figure_2B, device="pdf", width=1.5, height=4, units="in")

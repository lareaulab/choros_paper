rm(list=ls())

library(here)
library(choros)
library(simRiboSeq)
library(facetscales)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

expts <- c("fixedLinker_fixedPrimer", "randomLinker_randomPrimer")
names(expts) <- c("fixed linker\nfixed primer", "random linker\nrandom primer")

for(expt in expts) {
  load(file.path(data_dir, "lecanda_2016", expt,
                 paste0(expt, "_codon_corr.Rda")))
}

# make plot ---------------------------------------------------------------

scales_y <- lapply(expts, 
                   function(x) {
                     scale_y_continuous(limits=c(0, max(unlist(get(paste0(x, "_codon_corr"))))))
                   })

plot_max <- max(sapply(expts, function(x) { max(unlist(get(paste0(x, "_codon_corr")))) }))

corrected_codon_corr <- data.frame(position=unlist(lapply(expts,
                                                          function(x) {
                                                            names(get(paste0(x, "_codon_corr"))[[2]])
                                                          })),
                                   codon_corr=unlist(lapply(expts,
                                                            function(x) {
                                                              get(paste0(x, "_codon_corr"))[[2]]
                                                            })),
                                   expt=unlist(lapply(expts,
                                                      function(x) {
                                                        rep(x, length(get(paste0(x, "_codon_corr"))[[2]]))
                                                      })))
corrected_codon_corr$type <- "Corrected counts"
corrected_codon_corr$expt <- factor(corrected_codon_corr$exp, levels=expts)
levels(corrected_codon_corr$expt) <- names(expts)
corrected_codon_corr$position <- sub("n", "-", corrected_codon_corr$position)
corrected_codon_corr$position <- sub("p", "", corrected_codon_corr$position)
corrected_codon_corr$position <- factor(corrected_codon_corr$position, 
                                        levels=unique(corrected_codon_corr$position))
corrected_codon_corr$label <- as.character(corrected_codon_corr$position)
corrected_codon_corr$label[corrected_codon_corr$label %in% c("-5", "-4", "3", "4")] <- "bias"
corrected_codon_corr$label[!(corrected_codon_corr$label %in% c("A", "P", "E", "bias"))] <- "other"
fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")
names(fill_colors) <- c("bias", "E", "P", "A", "other")

figure_4C <- ggplot(corrected_codon_corr, aes(x=position, y=codon_corr, fill=label)) + 
  geom_col() + scale_fill_manual(values=fill_colors) + 
  theme_classic(base_size=6) + theme(legend.position="none") + 
  # facet_grid_sc(rows=vars(expt), cols=vars(type),
  #               scales=list(y=scales_y)) + 
  facet_grid(expt~"Corrected counts") + 
  xlab("codon position") + ylab(expression(Delta*" correlation")) + 
  coord_cartesian(ylim=c(0, plot_max))

ggsave(filename=file.path(figures_dir, "figure_4C.pdf"),
       plot=figure_4C, device="pdf", width=1.5, height=2.5, units="in")

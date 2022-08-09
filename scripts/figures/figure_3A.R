rm(list=ls())

library(here)
library(choros)
library(facetscales)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

tunney_codon_corr_fname <- file.path(data_dir, "tunney_2018", "tunney",
                                     "tunney_codon_corr.Rda")
load(tunney_codon_corr_fname)
schuller_codon_corr_fname <- file.path(data_dir, "schuller_2017", "schuller",
                                       "schuller_codon_corr.Rda")
load(schuller_codon_corr_fname)
weinberg_codon_corr_fname <- file.path(data_dir, "weinberg_2016", "weinberg",
                                       "weinberg_codon_corr.Rda")
load(weinberg_codon_corr_fname)

# make plot ---------------------------------------------------------------

tunney_max <- max(unlist(tunney_codon_corr))
schuller_max <- max(unlist(schuller_codon_corr))
weinberg_max <- max(unlist(weinberg_codon_corr))

scales_y <- list(
  "Tunney" = scale_y_continuous(limits=c(0, tunney_max)),
  "Schuller" = scale_y_continuous(limits=c(0, schuller_max)),
  "Weinberg" = scale_y_continuous(limits=c(0, weinberg_max))
)

expts <- c("tunney", "schuller", "weinberg")
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
raw_codon_corr$expt <- factor(raw_codon_corr$exp, levels=expts)
levels(raw_codon_corr$expt) <- c("Tunney", "Schuller", "Weinberg")
raw_codon_corr$position <- sub("n", "-", raw_codon_corr$position)
raw_codon_corr$position <- sub("p", "", raw_codon_corr$position)
raw_codon_corr$position <- factor(raw_codon_corr$position, 
                                  levels=unique(raw_codon_corr$position))
raw_codon_corr$label <- as.character(raw_codon_corr$position)
raw_codon_corr$label[raw_codon_corr$label %in% c("-5", "-4", "3", "4")] <- "bias"
raw_codon_corr$label[!(raw_codon_corr$label %in% c("A", "P", "E", "bias"))] <- "other"
fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")
names(fill_colors) <- c("bias", "E", "P", "A", "other")

figure_3A <- ggplot(raw_codon_corr, aes(x=position, y=codon_corr, fill=label)) + 
  geom_col() + scale_fill_manual(values=fill_colors) + theme_classic(base_size=8) + 
  theme(legend.position="none", panel.spacing=unit(0.25, "in")) + 
  facet_grid_sc(rows=vars(expt), cols=vars(type),
                scales=list(y=scales_y)) + 
  xlab("codon position") + ylab(expression(Delta*" correlation"))

ggsave(filename=file.path(figures_dir, "figure_3A.pdf"),
       plot=figure_3A, device="pdf", width=1.75, height=3, units="in")

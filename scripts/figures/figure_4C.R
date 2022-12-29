rm(list=ls())

library(here)
library(choros)
library(simRiboSeq)
library(facetscales)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

expts <- c("fixedLinker_fixedPrimer", "randomLinker_randomPrimer")
names(expts) <- c("fixed linker\nstandard primer", "random linker\nrandom primer")

for(expt in expts) {
  load(file.path(data_dir, "lecanda_2016", expt,
                 paste0(expt, "_coef.Rda")))
  assign(paste0(expt, "_coef"),
         subset(get(paste0(expt, "_coef")),
                group %in% c("f5", "f3", "d5:f5", "d3:f3")))
}

# aggregate data ----------------------------------------------------------

all_terms <- fixedLinker_fixedPrimer_coef[, c("group", "term")]
all_terms$fixed_estimate <- fixedLinker_fixedPrimer_coef$estimate
all_terms$fixed_stderr <- fixedLinker_fixedPrimer_coef$std_error
all_terms$random_estimate <- randomLinker_randomPrimer_coef$estimate
all_terms$random_stderr <- randomLinker_randomPrimer_coef$std_error
all_terms$group <- factor(all_terms$group, levels=c("f5", "d5:f5", "f3", "d3:f3"))

figure_4D <- ggplot(all_terms, aes(x=fixed_estimate, y=random_estimate, col=group)) + 
  geom_point(size=0.1) + geom_abline(slope=1, intercept=0, color="grey25", size=0.25) + 
  facet_wrap(~group) + theme_classic(base_size=8) + 
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
  xlab("fixed adapter / standard RT primer") + ylab("random adapter\nrandom RT primer")

ggsave(filename=file.path(figures_dir, "figure_4C.pdf"),
       plot=figure_4D, device="pdf", width=3, height=2, units="in")

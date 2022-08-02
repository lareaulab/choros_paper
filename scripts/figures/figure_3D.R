rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

tunney_coef_fname <- file.path(data_dir, "tunney_2018", "tunney",
                               "tunney_coef.Rda")
load(tunney_coef_fname)
schuller_coef_fname <- file.path(data_dir, "schuller_2017", "schuller",
                                 "schuller_coef.Rda")
load(schuller_coef_fname)
weinberg_coef_fname <- file.path(data_dir, "weinberg_2016", "weinberg",
                                 "weinberg_coef.Rda")
load(weinberg_coef_fname)

yeast_codon_properties <- read.table(file.path(here(), "reference_data", 
                                               "yeast_codon_properties.txt"), 
                                     header=T)

# aggregate data ----------------------------------------------------------

tunney_A <- within(subset(tunney_coef, group=="A"), {
  expt <- "Tunney"
  tAI <- yeast_codon_properties$tAI[match(term, yeast_codon_properties$codon)]
})
schuller_A <- within(subset(schuller_coef, group=="A"), {
  expt <- "Schuller"
  tAI <- yeast_codon_properties$tAI[match(term, yeast_codon_properties$codon)]
})
weinberg_A <- within(subset(weinberg_coef, group=="A"), {
  expt <- "Weinberg"
  tAI <- yeast_codon_properties$tAI[match(term, yeast_codon_properties$codon)]
})
all_A <- rbind(tunney_A, schuller_A, weinberg_A)
all_A$expt <- factor(all_A$expt, levels=c("Tunney", "Schuller", "Weinberg"))

plot_text <- data.frame(label=sapply(levels(all_A$expt),
                                     function(x) {
                                       tmp_cor <- round(with(subset(all_A, expt==x),
                                                             cor(tAI, exp(estimate))), digits=3)
                                       # as.character(expression(rho*as.character(paste("=", tmp_cor))))
                                       paste("cor =", tmp_cor)
                                     }),
                        tAI=1, 
                        estimate=log(5.8),
                        expt=factor(levels(all_A$expt), levels=levels(all_A$expt)))

figure_3D <- ggplot(all_A, aes(x=tAI, y=exp(estimate))) + 
  geom_point(size=0.5, color="grey25") + 
  geom_errorbar(aes(ymin=exp(estimate + qnorm(0.025)*std_error),
                    ymax=exp(estimate + qnorm(0.975)*std_error))) + 
  geom_smooth(method="lm", formula=y~x) + 
  theme_classic(base_size=6) + facet_grid(~expt) + 
  ylab(expression("exp("*beta^A*")")) + 
  geom_text(data=plot_text, aes(label=label), hjust=1, size=1.5) 

ggsave(filename=file.path(figures_dir, "figure_3D.pdf"),
       plot=figure_3D, device="pdf", width=3, height=1.5, units="in")

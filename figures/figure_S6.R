rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data", "lecanda_2016")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

expts <- c("fixedLinker_fixedPrimer", "randomLinker_randomPrimer")

for(expt in expts) {
  load(file.path(data_dir, expt, paste0(expt, "_coef.Rda")))
}

# generate plot -----------------------------------------------------------

codons <- subset(fixedLinker_fixedPrimer_coef, group=="A")$term

codon_coef <- lapply(c("A", "P", "E"),
                     function(x) {
                       fixed_coef <- subset(fixedLinker_fixedPrimer_coef, group==x)
                       fixed_coef <- fixed_coef[match(codons, fixed_coef$term),]
                       random_coef <- subset(randomLinker_randomPrimer_coef, group==x)
                       random_coef <- random_coef[match(codons, random_coef$term),]
                       data.frame(codon=codons,
                                  site=x,
                                  fixed_est=fixed_coef$estimate,
                                  fixed_stderr=fixed_coef$std_error,
                                  random_est=random_coef$estimate,
                                  random_stderr=random_coef$std_error)
                     })
codon_coef <- do.call(rbind, codon_coef)
codon_coef$site <- factor(codon_coef$site, levels=c("E", "P", "A"))

figure_S6 <- ggplot(codon_coef, aes(x=fixed_est, y=random_est)) + 
  geom_point(size=0.5, color="grey25") + 
  geom_abline(slope=1, intercept=0, color="blue") + 
  theme_classic(base_size=6) + facet_grid(~site) + 
  xlab(expression("fixed adapter, fixed RT primer "*beta)) + 
  ylab(expression("random adapter, random RT primer "*beta))

ggsave(filename=file.path(figures_dir, "figure_S6.pdf"),
       plot=figure_S6, device="pdf", width=6.5, height=2, units="in")

rm(list=ls())

library(here)
library(choros)
library(facetscales)

data_dir <- file.path(here(), "data", "wu_2019")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

expts <- c("WT", "3AT")

for(expt in expts) {
  load(file.path(data_dir, paste0("CHXTIG_", expt),
                 paste0("CHXTIG_", expt, "_coef.Rda")))
}

# generate plot -----------------------------------------------------------

coef_WT <- subset(CHXTIG_WT_coef, group=="A")
coef_3AT <- subset(CHXTIG_3AT_coef, group=="A")
codons <- coef_WT$term
coef_A <- data.frame(codons=codons,
                     est_WT=coef_WT$estimate[match(codons, coef_WT$term)],
                     stderr_WT=coef_WT$std_error[match(codons, coef_WT$term)],
                     est_3AT=coef_3AT$estimate[match(codons, coef_3AT$term)],
                     stderr_3AT=coef_3AT$std_error[match(codons, coef_3AT$term)])

figure_5A <- ggplot(coef_A,
                    aes(x=est_WT, y=est_3AT,
                        xmin=est_WT + qnorm(0.025)*stderr_WT,
                        xmax=est_WT + qnorm(0.975)*stderr_WT,
                        ymin=est_3AT + qnorm(0.025)*stderr_3AT,
                        ymax=est_3AT + qnorm(0.975)*stderr_3AT)) +
  theme_classic(base_size=8) + geom_abline(slope=1, intercept=0, color="blue") +
  geom_point(size=0.5, color="grey25") +
  geom_errorbar(color="grey25", alpha=0.5) + geom_errorbarh(color="grey25", alpha=0.5) +
  geom_text(data=subset(coef_A, abs(est_WT)>2 | abs(est_3AT)>2), color="grey25",
            aes(label=codons), position=position_nudge(x=-0.35, y=0.2), size=2) +
  xlab(expression("WT "*beta^A)) + ylab(expression("3AT "*beta^A))

ggsave(filename=file.path(figures_dir, "figure_5A.pdf"),
       plot=figure_5A, device="pdf", width=3.2, height=2, units="in")

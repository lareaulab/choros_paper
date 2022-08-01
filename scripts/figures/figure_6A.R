rm(list=ls())

library(here)
library(choros)
library(facetscales)

data_dir <- "~/footprint-bias/expts/wu_2019"
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

expts <- c("WT", "3AT")

for(expt in expts) {
  load(file.path(data_dir, paste0("CHX_TIG_", expt),
                 paste0("CHX_TIG_", expt, "_fit_150.Rda")))
  assign(paste0("CHXTIG_", expt, "_coef"),
         parse_coefs(get(paste0("CHX_TIG_", expt, "_fit_150"))))
}

# generate plot -----------------------------------------------------------

coef_A <- data.frame(codon=subset(CHXTIG_WT_coef, group=="A")$term,
                     coef_WT=subset(CHXTIG_WT_coef, group=="A")$estimate,
                     stderr_WT=subset(CHXTIG_WT_coef, group=="A")$std_error,
                     coef_3AT=subset(CHXTIG_3AT_coef, group=="A")$estimate,
                     stderr_3AT=subset(CHXTIG_3AT_coef, group=="A")$std_error)

figure_6A <- ggplot(coef_A,
                    aes(x=coef_WT, y=coef_3AT,
                        xmin=coef_WT + qnorm(0.025)*stderr_WT,
                        xmax=coef_WT + qnorm(0.975)*stderr_WT,
                        ymin=coef_3AT + qnorm(0.025)*stderr_3AT,
                        ymax=coef_3AT + qnorm(0.975)*stderr_3AT)) +
  geom_point(size=0.5, color="grey25") +
  geom_errorbar(color="grey25", alpha=0.5) +
  geom_errorbarh(color="grey25", alpha=0.5) +
  theme_classic(base_size=6) + geom_abline(slope=1, intercept=0, color="blue") +
  geom_text(data=subset(coef_A, abs(coef_WT)>2 | abs(coef_3AT)>2), color="grey25",
            aes(label=codon), position=position_nudge(x=0.2, y=0.2), size=2) +
  xlab(expression("WT "*beta^A)) + ylab(expression("3AT "*beta^A))

ggsave(filename=file.path(figures_dir, "figure_6A.pdf"),
       plot=figure_6A, device="pdf", width=3.25, height=2, units="in")

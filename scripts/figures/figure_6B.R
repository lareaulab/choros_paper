rm(list=ls())

library(here)
library(choros)
library(facetscales)

data_dir <- "~/footprint-bias/expts/wu_2019"
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "CHXTIG_WT_vs_3AT", "WT_3AT_fit.Rda"))

# generate plot -----------------------------------------------------------

model_coef <- parse_coefs(WT_3AT_fit)

## 3AT was reference level for expt --> need to take negative value

model_coef <- subset(model_coef, group=="A:expt")

figure_6B <- ggplot(model_coef,
                    aes(x=-estimate, y=-log10(p+1e-25),
                        xmin=-(estimate + qnorm(0.025)*std_error),
                        xmax=-(estimate + qnorm(0.975)*std_error))) +
  geom_point(size=0.5, color="grey25") + # geom_errorbarh(color="grey25") +
  coord_cartesian(xlim=c(-2,4)) + theme_classic(base_size=6) +
  xlab(expression(beta^{"A:expt"})) + ylab("-log10(p)") +
  geom_hline(yintercept=0, color="grey25", size=0.5) +
  geom_vline(xintercept=0, color="grey25", size=0.5) +
  geom_text(data=subset(model_coef, -log10(p) > 20), size=2, color="grey25",
            aes(label=term), position=position_nudge(x=-0.2))

ggsave(filename=file.path(figures_dir, "figure_6B.pdf"),
       plot=figure_6B, device="pdf", width=3.25, height=2, units="in")

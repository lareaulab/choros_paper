rm(list=ls())

library(here)
library(choros)
library(facetscales)

data_dir <- file.path(here(), "data", "wu_2019")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "CHXTIG_WT_v_3AT", "CHXTIG_WT_v_3AT_coef.Rda"))

# generate plot -----------------------------------------------------------


CHXTIG_WT_v_3AT_coef <- subset(CHXTIG_WT_v_3AT_coef, group=="A:expt")

figure_5B <- ggplot(CHXTIG_WT_v_3AT_coef,
                    aes(x=estimate, y=-log10(p+1e-25),
                        xmin=-(estimate + qnorm(0.025)*std_error),
                        xmax=-(estimate + qnorm(0.975)*std_error))) +
  geom_point(size=0.5, color="grey25") + # geom_errorbarh(color="grey25") +
  coord_cartesian(xlim=c(-1.5,5)) + theme_classic(base_size=8) +
  xlab(expression(beta^{"A:expt"})) + ylab("-log10(p)") +
  geom_hline(yintercept=0, color="grey25", size=0.5) +
  geom_vline(xintercept=0, color="grey25", size=0.5) +
  geom_text(data=subset(CHXTIG_WT_v_3AT_coef, -log10(p) > 12), size=2, color="grey25",
            aes(label=term), position=position_nudge(x=-0.27)) + 
  geom_text(data=subset(CHXTIG_WT_v_3AT_coef, term %in% c("CGA", "CGG")), size=2, color="grey25",
            aes(label=term), position=position_nudge(x=-0.27))

ggsave(filename=file.path(figures_dir, "figure_5B.pdf"),
       plot=figure_5B, device="pdf", width=3.2, height=2, units="in")

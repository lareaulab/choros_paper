rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data", "simulated_data")
figures_dir <- file.path(here(), "figures")

expts <- c("no_bias", "f3_bias", "f5_bias", "f3_f5_bias")

# load data ---------------------------------------------------------------

for(expt in expts) {
  load(file.path(data_dir, expt, paste0(expt, "_coef.Rda")))
}

data(f5_bias, f3_bias, codon_TE, package="simRiboSeq")
f5_bias <- f5_bias / f5_bias["AAA"]
f3_bias <- f3_bias / f3_bias["AAA"]
codon_TE <- codon_TE / codon_TE["AAA"]

# 5' terms  ---------------------------------------------------------------

f5_coefs <- lapply(expts,
                   function(expt) {
                     tmp_coef <- within(subset(get(paste0(expt, "_coef")),
                                               group=="f5"),
                                        expt <- expt)
                     if(expt %in% c("no_bias", "f3_bias")) {
                       tmp_coef$sim_prob <- 1
                     } else {
                       tmp_coef$sim_prob <- f5_bias[as.character(tmp_coef$term)]
                     }
                     return(tmp_coef)
                   })
f5_coefs <- do.call(rbind, f5_coefs)

# 3' terms ----------------------------------------------------------------

f3_coefs <- lapply(expts,
                   function(expt) {
                     tmp_coef <- within(subset(get(paste0(expt, "_coef")),
                                               group=="f3"),
                                        expt <- expt)
                     if(expt %in% c("no_bias", "f5_bias")) {
                       tmp_coef$sim_prob <- 1
                     } else {
                       tmp_coef$sim_prob <- f3_bias[as.character(tmp_coef$term)]
                     }
                     return(tmp_coef)
                   })
f3_coefs <- do.call(rbind, f3_coefs)

# A-site terms ------------------------------------------------------------

A_coefs <- lapply(expts,
                  function(expt) {
                    within(subset(get(paste0(expt, "_coef")),
                                  group=="A"), {
                                    expt <- expt
                                    sim_prob <- codon_TE[as.character(term)]
                                  })
                  })
A_coefs <- do.call(rbind, A_coefs)

# generate plot -----------------------------------------------------------

all_coefs <- rbind(f5_coefs, f3_coefs, A_coefs)
all_coefs$expt <- factor(all_coefs$expt, levels=expts)
levels(all_coefs$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")
all_coefs$group <- factor(all_coefs$group, levels=c("f5", "f3", "A"))
levels(all_coefs$group) <- c("5' recovery", "3' recovery", "A-site codon")

figure_2A <- ggplot(all_coefs, aes(x=sim_prob, y=exp(estimate))) +
  geom_abline(slope=1, intercept=0, col="blue", alpha=0.5) +
  geom_point(col="grey25", size=0.5) + theme_classic(base_size=8) +
  facet_grid(expt ~ group) +
  xlab("simulation probability") + ylab(expression("exp("*beta*")")) +
  theme(panel.spacing=unit(0.25, "in"), axis.text.x=element_text(size=6))

ggsave(filename=file.path(figures_dir, "figure_2A.pdf"),
       plot=figure_2A, device="pdf", width=3, height=4, units="in")

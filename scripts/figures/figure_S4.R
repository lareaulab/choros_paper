rm(list=ls())

library(here)
library(choros)
library(simRiboSeq)

data_dir <- file.path(here(), "data", "simulated_data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

expts <- c("no_bias", "f3_bias", "f5_bias", "f3_f5_bias")
for(expt in expts) {
  load(file.path(data_dir, expt, paste0(expt, "_coef.Rda")))
}

data(d5_bias, d3_bias, package="simRiboSeq")

# aggregate data ----------------------------------------------------------

d5_bias <- d5_bias/d5_bias["15"]
d3_bias <- d3_bias/d3_bias["10"]

all_digest_coefs <- lapply(expts,
                           function(x) {
                             coef_obj <- paste0(x, "_coef")
                             f5_coef <- within(subset(get(coef_obj), group=="d5"), {
                               expt <- x
                               sim_prob <- d5_bias[match(term, names(d5_bias))]
                             })
                             f3_coef <- within(subset(get(coef_obj), group=="d3"), {
                               expt <- x
                               sim_prob <- d3_bias[match(term, names(d3_bias))]
                             })
                             return(rbind(f5_coef, f3_coef))
                           })
all_digest_coefs <- do.call(rbind, all_digest_coefs)
all_digest_coefs$expt <- factor(all_digest_coefs$expt, levels=expts)
levels(all_digest_coefs$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")

# make plot ---------------------------------------------------------------

figure_S4 <- ggplot(all_digest_coefs, aes(x=sim_prob, y=exp(estimate))) + 
  geom_point() + theme_classic(base_size=6) + 
  facet_grid(expt ~ group) + 
  geom_abline(slope=1, intercept=0, color="blue") + 
  xlab("relative simulation probability") + 
  ylab(expression("exp("*beta*")"))

ggsave(filename=file.path(figures_dir, "figure_S4.pdf"),
       plot=figure_S4, device="pdf", width=6, height=4)

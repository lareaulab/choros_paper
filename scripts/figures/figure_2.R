rm(list=ls())

library(choros)
library(simRiboSeq)
library(facetscales)

setwd("~/footprint-bias/expts/simulated_data/")

expts <- c("noBias", "n3Bias", "p5Bias", "bothBias")

for(expt in expts) {
  load(file.path(expt, "model_fit_100genes.Rda"))
  assign(paste0(expt, "_coef"),
         parse_coefs(model_fit_100))
}

data(f5_bias, f3_bias, codon_TE, package="simRiboSeq")

# f5 plot -----------------------------------------------------------------

f5_bias <- data.frame(term=names(f5_bias), prob=f5_bias)
f5_bias$f5 <- substr(f5_bias$term, 1, 2)
f5_bias <- aggregate(prob ~ f5, f5_bias, mean)
f5_bias$prob_norm <- f5_bias$prob / f5_bias$prob[which(f5_bias$f5=="AA")]

f5_coefs <- lapply(expts,
                   function(x) {
                     within(subset(get(paste0(x, "_coef")),
                                   group=="f5"), {
                                     expt <- x
                                   })
                   })
f5_coefs <- do.call(rbind, f5_coefs)
f5_coefs$sim_prob <- f5_bias$prob_norm[match(f5_coefs$term, f5_bias$f5)]
f5_coefs$sim_prob[f5_coefs$expt %in% c("noBias", "n3Bias")] <- 1

f5_coefs$expt <- factor(f5_coefs$expt, levels=expts)
levels(f5_coefs$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")

f5_plot <- ggplot(f5_coefs,
                  aes(x=sim_prob, y=exp(estimate),
                      ymin=exp(estimate + qnorm(0.025)*std_error),
                      ymax=exp(estimate + qnorm(0.975)*std_error))) +
  geom_point(size=0.5, color="grey25") + geom_abline(slope=1, intercept=0, col="blue") +
  theme_classic(base_size=6) + facet_grid(expt~.) +
  xlab("relative 5' recovery probability") + ylab(expression("exp("*beta^{f^5}*")"))

ggsave(filename="~/choros_paper/figures/figure_2A.pdf",
       plot=f5_plot, device="pdf", width=1.15, height=4, units="in")

# f3 plot -----------------------------------------------------------------

f3_bias <- data.frame(term=names(f3_bias), prob=f3_bias)
f3_bias$prob_norm <- f3_bias$prob / f3_bias$prob[which(f3_bias$term == "AAA")]

f3_coefs <- lapply(expts,
                   function(x) {
                     within(subset(get(paste0(x, "_coef")),
                                   group=="f3"), {
                                     expt <- x
                                   })
                   })
f3_coefs <- do.call(rbind, f3_coefs)
f3_coefs$sim_prob <- f3_bias$prob_norm[match(f3_coefs$term, f3_bias$term)]
f3_coefs$sim_prob[f3_coefs$expt %in% c("noBias", "p5Bias")] <- 1

f3_coefs$expt <- factor(f3_coefs$expt, levels=expts)
levels(f3_coefs$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")

f3_plot <- ggplot(f3_coefs,
                  aes(x=sim_prob, y=exp(estimate),
                      ymin=exp(estimate + qnorm(0.025)*std_error),
                      ymax=exp(estimate + qnorm(0.975)*std_error))) +
  geom_point(size=0.5, color="grey25") + geom_abline(slope=1, intercept=0, col="blue") +
  theme_classic(base_size=6) + facet_grid(expt~.) +
  xlab("relative 3' recovery probability") + ylab(expression("exp("*beta^{f^3}*")"))

ggsave(filename="~/choros_paper/figures/figure_2B.pdf",
       plot=f3_plot, device="pdf", width=1.15, height=4, units="in")

# A plot ------------------------------------------------------------------

codon_TE <- codon_TE/codon_TE["AAA"]

A_coefs <- lapply(expts,
                  function(x) {
                    within(subset(get(paste0(x, "_coef")),
                                  group=="A"), {
                                    expt <- x
                                    codon_TE <- codon_TE[match(term, names(codon_TE))]
                                  })
                  })
A_coefs <- do.call(rbind, A_coefs)

A_coefs$expt <- factor(A_coefs$expt, levels=expts)
levels(A_coefs$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")

A_plot <- ggplot(A_coefs,
                  aes(x=codon_TE, y=exp(estimate),
                      ymin=exp(estimate + qnorm(0.025)*std_error),
                      ymax=exp(estimate + qnorm(0.975)*std_error))) +
  geom_point(size=0.5, color="grey25") + geom_abline(slope=1, intercept=0, col="blue") +
  theme_classic(base_size=6) + facet_grid(expt~.) +
  xlab("A-site codon weight") + ylab(expression("exp("*beta^A*")"))

ggsave(filename="~/choros_paper/figures/figure_2C.pdf",
       plot=A_plot, device="pdf", width=1.15, height=4, units="in")

# position importance plot: raw counts ------------------------------------

raw_codon_corr <- lapply(expts,
                         function(x) {
                           load(file.path(x, "bias_codon.Rda"))
                           tmp <- bias_codon[[1]]
                           data.frame(position=names(tmp),
                                      codon_corr=tmp,
                                      expt=x)
                         })
raw_codon_corr <- do.call(rbind, raw_codon_corr)
raw_codon_corr$type <- "Raw counts"
raw_codon_corr$expt <- factor(raw_codon_corr$exp, levels=expts)
levels(raw_codon_corr$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")
raw_codon_corr$position <- sub("n", "-", raw_codon_corr$position)
raw_codon_corr$position <- sub("p", "", raw_codon_corr$position)
raw_codon_corr$position <- factor(raw_codon_corr$position,
                                  levels=unique(raw_codon_corr$position))
raw_codon_corr$label <- as.character(raw_codon_corr$position)
raw_codon_corr$label[raw_codon_corr$label %in% c("-5", "-4", "3", "4")] <- "bias"
raw_codon_corr$label[!(raw_codon_corr$label %in% c("A", "P", "E", "bias"))] <- "other"
fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")
names(fill_colors) <- c("bias", "E", "P", "A", "other")

scales_y <- lapply(expts,
                   function(x) {
                     load(file.path(x, "bias_codon.Rda"))
                     tmp <- c(bias_codon[[1]], bias_codon[[2]])
                     return(scale_y_continuous(limits=c(0, max(tmp))))
                   })
names(scales_y) <- levels(raw_codon_corr$expt)

figure_2D <- ggplot(raw_codon_corr, aes(x=position, y=codon_corr, fill=label)) +
  geom_col() + scale_fill_manual(values=fill_colors) +
  theme_classic(base_size=6) + theme(legend.position="none") +
  facet_grid_sc(rows=vars(expt), cols=vars(type),
                scales=list(y=scales_y)) +
  xlab("codon position") + ylab(expression(Delta*" correlation"))

ggsave(filename="~/choros_paper/figures/figure_2D.pdf",
       plot=figure_2D, device="pdf", width=1.5, height=4, units="in")

# position importance plots: corrected counts -----------------------------

corrected_codon_corr <- lapply(expts,
                         function(x) {
                           load(file.path(x, "bias_codon.Rda"))
                           tmp <- bias_codon[[2]]
                           data.frame(position=names(tmp),
                                      codon_corr=tmp,
                                      expt=x)
                         })
corrected_codon_corr <- do.call(rbind, corrected_codon_corr)
corrected_codon_corr$type <- "corrected counts"
corrected_codon_corr$expt <- factor(corrected_codon_corr$exp, levels=expts)
levels(corrected_codon_corr$expt) <- c("no bias", "3' bias", "5' bias", "3' and 5' bias")
corrected_codon_corr$position <- sub("n", "-", corrected_codon_corr$position)
corrected_codon_corr$position <- sub("p", "", corrected_codon_corr$position)
corrected_codon_corr$position <- factor(corrected_codon_corr$position,
                                  levels=unique(corrected_codon_corr$position))
corrected_codon_corr$label <- as.character(corrected_codon_corr$position)
corrected_codon_corr$label[corrected_codon_corr$label %in% c("-5", "-4", "3", "4")] <- "bias"
corrected_codon_corr$label[!(corrected_codon_corr$label %in% c("A", "P", "E", "bias"))] <- "other"
fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")
names(fill_colors) <- c("bias", "E", "P", "A", "other")

scales_y <- lapply(expts,
                   function(x) {
                     load(file.path(x, "bias_codon.Rda"))
                     tmp <- c(bias_codon[[1]], bias_codon[[2]])
                     return(scale_y_continuous(limits=c(0, max(tmp))))
                   })
names(scales_y) <- levels(corrected_codon_corr$expt)

figure_2E <- ggplot(corrected_codon_corr, aes(x=position, y=codon_corr, fill=label)) +
  geom_col() + scale_fill_manual(values=fill_colors) +
  theme_classic(base_size=6) + theme(legend.position="none") +
  facet_grid_sc(rows=vars(expt), cols=vars(type),
                scales=list(y=scales_y)) +
  xlab("codon position") + ylab(expression(Delta*" correlation"))

ggsave(filename="~/choros_paper/figures/figure_2E.pdf",
       plot=figure_2E, device="pdf", width=1.5, height=4, units="in")

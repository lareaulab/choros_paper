rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data", "chou_2017", "processed_data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "codon_cts.Rda"))

# compute mean and variance -----------------------------------------------

mean_variance <- codon_cts[, c("transcript", "cod_idx")]
mean_variance$mean <- rowMeans(codon_cts[, grepl("^SRR", colnames(codon_cts))])
mean_variance$variance <- sapply(seq(nrow(mean_variance)),
                                 function(x) {
                                   var(unlist(codon_cts[x, grepl("^SRR", colnames(codon_cts))]))
                                 })
save(mean_variance, file=file.path(data_dir, "mean_variance.Rda"))

# compute mean-variance fits ----------------------------------------------

mean_variance$mean_sq <- mean_variance$mean^2

poisson_fit <- lm(variance ~ 0 + offset(mean), mean_variance)
quasipoisson_fit <- lm(variance ~ 0 + mean, mean_variance)
negbin_fit <- lm(variance ~ 0 + offset(mean) + mean_sq, mean_variance)

fit_prediction <- data.frame(mean=seq(1, 1000, length.out=500))
fit_prediction$mean_sq <- fit_prediction$mean^2
fit_prediction$poisson <- predict(poisson_fit, newdata=fit_prediction)
fit_prediction$quasipoisson <- predict(quasipoisson_fit, newdata=fit_prediction)
fit_prediction$negbin <- predict(negbin_fit, newdata=fit_prediction)

# compute goodness of fit -------------------------------------------------

AIC(poisson_fit, quasipoisson_fit)
BIC(poisson_fit, quasipoisson_fit)

AIC(poisson_fit, negbin_fit)
BIC(poisson_fit, negbin_fit)

# generate plot -----------------------------------------------------------

figure_1C <- ggplot(subset(mean_variance, mean < 1000), aes(x=mean, y=variance)) +
  theme_classic(base_size=6) + geom_point(size=0.5, alpha=0.1) +
  # scale_fill_gradient(low="grey", high="blue") + labs(fill="") +
  geom_line(data=fit_prediction, aes(y=poisson), col="purple") +
  geom_line(data=fit_prediction, aes(y=quasipoisson), col="green") +
  geom_line(data=fit_prediction, aes(y=negbin), col="red")

ggsave(filename=file.path(figures_dir, "figure_1C.pdf"),
       plot=figure_1C, device="pdf", width=2, height=2, units="in")
ggsave(filename=file.path(figures_dir, "figure_1C.png"),
       plot=figure_1C, device="png", width=2, height=2, units="in", dpi="print")

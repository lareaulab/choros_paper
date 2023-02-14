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

# load qPCR data ----------------------------------------------------------

se <- function(x) { sd(x)/sqrt(length(x)) }

oligos <- c("ATA", "TCC", "CCA", "CGT", "GAC", "GGG")
names(oligos) <- c(1:3, 7:9)

# read in qPCR data
qpcr_fname <- file.path(here(), "reference_data", "circligase_qpcr.csv")
raw_data <- read.csv(qpcr_fname, stringsAsFactors=F)
raw_data <- subset(raw_data, as.character(raw_data$Oligo) %in% names(oligos))
raw_data$end_seq <- oligos[as.character(raw_data$Oligo)]
raw_data$primer_type <- ifelse(raw_data$PRIMERS == "NM827_NM828", "circ", "cont")

# calculate efficiency per primer
eff <- aggregate(eff ~ primer_type, data=raw_data, FUN=mean)

# calcaulate amplification
raw_data$primer_eff <- eff$eff[match(raw_data$primer_type, eff$primer)]
raw_data$amp <- raw_data$primer_eff^raw_data$Cy0

# average across replicates
qpcr_data <- aggregate(amp ~ end_seq + primer_type + CL, data=raw_data, FUN=mean)
qpcr_data$amp_se <- aggregate(amp ~ end_seq + primer_type + CL, data=raw_data, FUN=se)$amp
qpcr_data <- reshape(qpcr_data, v.names=c("amp", "amp_se"), idvar=c("end_seq", "CL"),
                     timevar="primer_type", direction="wide")

# calculate ratios
qpcr_data$ratio <- with(qpcr_data, amp.cont/amp.circ)
qpcr_data$ratio <- qpcr_data$ratio / max(qpcr_data$ratio)

# error propagation
qpcr_data$se_percent <- with(qpcr_data,
                             sqrt( (amp_se.circ / amp.circ)^2 + (amp_se.cont / amp.cont)^2 ))
qpcr_data$se_abs <- with(qpcr_data, ratio * se_percent)

# aggregate data ----------------------------------------------------------

tunney_f5 <- within(subset(tunney_coef, group=="f5" & term %in% qpcr_data$end_seq), {
  expt <- "Tunney"
  rel_eff <- qpcr_data$ratio[match(term, qpcr_data$end_seq)]
})
schuller_f5 <- within(subset(schuller_coef, group=="f5" & term %in% qpcr_data$end_seq), {
  expt <- "Schuller"
  rel_eff <- qpcr_data$ratio[match(term, qpcr_data$end_seq)]
})
weinberg_f5 <- within(subset(weinberg_coef, group=="f5" & term %in% qpcr_data$end_seq), {
  expt <- "Weinberg"
  rel_eff <- qpcr_data$ratio[match(term, qpcr_data$end_seq)]
})
all_f5 <- rbind(tunney_f5, schuller_f5, weinberg_f5)
all_f5$expt <- factor(all_f5$expt, levels=c("Tunney", "Schuller", "Weinberg"))

plot_text <- data.frame(label=sapply(levels(all_f5$expt),
                                     function(x) {
                                       tmp_cor <- round(with(subset(all_f5, expt==x),
                                                             cor(rel_eff, exp(estimate))), digits=3)
                                       # as.character(expression(rho*as.character(paste("=", tmp_cor))))
                                       paste("cor =", tmp_cor)
                                     }),
                        rel_eff=1, 
                        estimate=log(5),
                        expt=factor(levels(all_f5$expt), levels=levels(all_f5$expt)))

figure_3C <- ggplot(all_f5, aes(x=rel_eff, y=exp(estimate))) + 
  geom_smooth(method="lm", formula=y~x) + 
  geom_point(size=0.5, color="grey25") + 
  geom_errorbar(aes(ymin=exp(estimate + qnorm(0.025)*std_error),
                    ymax=exp(estimate + qnorm(0.975)*std_error)),
                col="grey25", size=0.5) + 
  theme_classic(base_size=8) + facet_grid(~expt) + 
  xlab("CircLigase in vitro ligation efficiency") + ylab(expression("exp("*beta^f^5*")")) + 
  geom_text(data=plot_text, aes(label=label), hjust=1, size=1.5) + 
  geom_text(data=data.frame(rel_eff=1, estimate=log(3.5), 
                            label="does not use\nCircLigase", 
                            expt=factor("Weinberg", levels=levels(all_f5$expt))),
            color="red", aes(label=label), size=1.5, hjust=1) + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename=file.path(figures_dir, "figure_3D.pdf"),
       plot=figure_3C, device="pdf", width=3, height=1.5, units="in")

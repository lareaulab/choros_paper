rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")
ref_dir <- file.path(here(), "reference_data")

# load data ---------------------------------------------------------------

yor1_id <- "YGR281W"

scer_lengths_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- load_lengths(scer_lengths_fname)
scer_lengths$num_codons <- with(scer_lengths, cds_length/3)

weinberg_bam_fname <- "/mnt/lareaulab/amok/choros_paper/data/weinberg_2016/weinberg/weinberg_bam.Rda"
load(weinberg_bam_fname)
codon_cts <- aggregate(cbind(count, correct_250) ~ transcript + cod_idx,
                       data=subset(weinberg_bam, transcript==yor1_id),
                       FUN=sum, na.rm=T, na.action=na.omit)
codon_cts$diff <- with(codon_cts, correct_250 - count)
codon_cts$new_diff <- 1.5*codon_cts$diff
codon_cts$new_correct_250 <- with(codon_cts, count + new_diff)
codon_cts$new_correct_250[codon_cts$cod_idx==510] <- codon_cts$new_correct_250[codon_cts$cod_idx==510] - 3

yor1 <- data.frame(transcript=yor1_id,
                   cod_idx=seq(scer_lengths$num_codons[scer_lengths$transcript==yor1_id]))
yor1 <- rbind(within(yor1,
                     {
                       count <- codon_cts$count[match(cod_idx, codon_cts$cod_idx)]
                       type <- "raw"
                     }),
              within(yor1,
                     {
                       count <- codon_cts$new_correct_250[match(cod_idx, codon_cts$cod_idx)]
                       type="corrected"
                     }))
yor1$count[is.na(yor1$count)] <- 0
yor1$type <- factor(yor1$type, levels=c("raw", "corrected"))

# generate figures --------------------------------------------------------

figure_1A_top <- ggplot(subset(yor1, cod_idx %in% seq(from=505, to=515) & type=="raw"),
                        aes(x=cod_idx, y=count, fill=type)) +
  geom_col(fill="grey45") + theme_classic(base_size=6) +
  xlab("codon position") + theme(legend.position="none") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab("position along gene") + ylab("true\nfootprint\ncount") +
  coord_cartesian(ylim=c(0, 23))

figure_1A_bottom <- ggplot(subset(yor1, cod_idx %in% seq(from=505, to=515) & type=="raw"),
                           aes(x=cod_idx, y=count)) +
  geom_bar(stat="identity", fill="grey45") + theme_classic(base_size=6) +
  xlab("codon position") + theme(legend.position="none") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab("position along gene") + ylab("observed\nfootprint\ncount") +
  coord_cartesian(ylim=c(0, 23)) +
  geom_col(data=subset(yor1, cod_idx %in% seq(from=505, to=515) & type=="corrected"),
           fill="transparent", linetype="dashed", col="black")

figure_1A <- figure_1A_top / figure_1A_bottom

ggsave(filename=file.path(figures_dir, "figure_1A.pdf"),
       plot=figure_1A, device="pdf", width=3, height=1.5, units="in")

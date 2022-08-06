rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

yor1_id <- "YGR281W"

scer_lengths_fname <- "/mnt/lareaulab/amok/choros_paper/reference_data/scer.transcripts.20cds20.lengths.txt"
scer_lengths <- load_lengths(scer_lengths_fname)
scer_lengths$num_codons <- with(scer_lengths, cds_length/3)

weinberg_bam_fname <- "/mnt/lareaulab/amok/choros_paper/data/weinberg_2016/weinberg/weinberg_bam.Rda"
load(weinberg_bam_fname)
codon_cts <- aggregate(cbind(count, correct_250) ~ transcript + cod_idx,
                       data=subset(weinberg_bam, transcript==yor1_id),
                       FUN=sum, na.rm=T, na.action=na.omit)

yor1 <- data.frame(transcript=yor1_id,
                   cod_idx=seq(scer_lengths$num_codons[scer_lengths$transcript==yor1_id]))
yor1 <- rbind(within(yor1,
                     {
                       count <- codon_cts$count[match(cod_idx, codon_cts$cod_idx)]
                       type <- "raw"
                     }),
              within(yor1,
                     {
                       count <- codon_cts$correct_250[match(cod_idx, codon_cts$cod_idx)]
                       type="corrected"
                     }))
yor1$count[is.na(yor1$count)] <- 0
yor1$type <- factor(yor1$type, levels=c("raw", "corrected"))

# generate figures --------------------------------------------------------

figure_1A_top <- ggplot(subset(yor1, cod_idx %in% seq(from=505, to=515) & type=="raw"),
                        aes(x=cod_idx, y=count, fill=type)) +
  geom_col(position="identity", alpha=0.25) + theme_classic(base_size=6) +
  xlab("codon position") + theme(legend.position="none") +
  scale_fill_manual(values=c("raw"="grey45", "corrected"="blue")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab("") + ylab("footprint count") + coord_cartesian(ylim=c(0, 22.5))

figure_1A_bottom <- ggplot(subset(yor1, cod_idx %in% seq(from=505, to=515)),
                           aes(x=cod_idx, y=count, fill=type)) +
  geom_col(position="identity", alpha=0.25) + theme_classic(base_size=6) +
  xlab("codon position") + theme(legend.position="none") +
  scale_fill_manual(values=c("raw"="grey45", "corrected"="blue")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab("") + ylab("footprint count") + coord_cartesian(ylim=c(0, 22.5))

figure_1A <- figure_1A_top / figure_1A_bottom

ggsave(filename=file.path(figures_dir, "figure_1A.pdf"),
       plot=figure_1A, device="pdf", width=3, height=1.25, units="in")

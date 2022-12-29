rm(list=ls())

library(here)
library(choros)
library(ggplot2)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

expts <- c("tunney", "schuller", "lecanda")

# load datasets -----------------------------------------------------------

tunney_bam_fname <- file.path(data_dir, "tunney_2018", "tunney",
                              "tunney_bam.Rda")
load(tunney_bam_fname)

schuller_bam_fname <- file.path(data_dir, "schuller_2017", "schuller",
                                "schuller_bam.Rda")
load(schuller_bam_fname)

lecanda_bam_fname <- file.path(data_dir, "lecanda_2016", "randomLinker_randomPrimer",
                               "randomLinker_randomPrimer_bam.Rda")
load(lecanda_bam_fname)

# generate ribogrid data at start codon -----------------------------------

generate_ribogrid <- function(bam_dat, grid_start = -20, grid_end = 40) {
  # bam_dat: data.frame; output from load_bam()
  ## `correct_250` column: corrected counts
  # grid_start: integer; 5'-most position of ribogrid plot
  # grid_end: integer; 3'-most position of ribogrid plot
  # 1. annotate RPF 5' nucleotide position
  bam_dat$nt_idx <- with(bam_dat, cod_idx * 3 - d5 - 3)
  # 2. annotate RPF length
  bam_dat$rpf_length <- with(bam_dat, d5 + d3 + 3)
  # 3. subset data to ribogrid coordinates
  bam_dat <- subset(bam_dat, nt_idx >= grid_start & nt_idx <= grid_end)
  # 4. aggregate RPF counts by nt_idx and rpf_length
  bam_dat <- aggregate(cbind(count, correct_250) ~ nt_idx + rpf_length,
                       data = bam_dat, FUN = sum)
  # 5. subset to most abundant rpf_length at nt_idx == 0
  ## corresponds to RPF 5' coordinate at A of AUG start codon
  which_length <- subset(bam_dat, nt_idx == 0)
  which_length <- which_length$rpf_length[which.max(which_length$count)]
  bam_dat <- subset(bam_dat, rpf_length == which_length)
  # 6. flesh out counts for unobserved nt_idx
  unobserved_pos <- seq.int(grid_start, grid_end)
  unobserved_pos <- unobserved_pos[!(unobserved_pos %in% bam_dat$nt_idx)]
  bam_dat <- rbind(data.frame(nt_idx = unobserved_pos,
                              rpf_length = which_length,
                              count = 0,
                              correct_250 = 0),
                   bam_dat)
  # 7. convert to long format
  bam_dat <- reshape::melt(bam_dat, measure.vars = c("count", "correct_250"),
                           variable_name = "count_type")
  levels(bam_dat$count_type) <- c("uncorrected", "corrected")
  return(bam_dat)
}

tunney_ribogrid <- generate_ribogrid(tunney_bam)
schuller_ribogrid <- generate_ribogrid(schuller_bam)
lecanda_ribogrid <- generate_ribogrid(randomLinker_randomPrimer_bam)

all_ribogrid <- rbind(data.frame(schuller_ribogrid, expt = "Schuller"),
                      data.frame(tunney_ribogrid, expt = "Tunney"),
                      data.frame(lecanda_ribogrid, expt = "Lecanda"))

# generate plot -----------------------------------------------------------

figure_5 <- ggplot(all_ribogrid,
                   aes(x = nt_idx, y = value / 1000, col = count_type)) +
  geom_line() + theme_classic() + facet_grid(expt ~ ., scales = "free_y") +
  scale_color_manual(values = c("uncorrected" = "red",
                                "corrected" = "grey40")) +
  xlab("position of footprint 5' end (nt)") +
  ylab("footprint count, thousands") +
  labs(col = "") + theme(legend.position = "bottom")

ggsave(filename = file.path(figures_dir, "figure_5.pdf"),
       plot = figure_5, device = "pdf", width = 3.5, height = 4.5)
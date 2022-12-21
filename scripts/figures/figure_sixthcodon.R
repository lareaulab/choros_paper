rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

pdfname <- file.path(figures_dir, "figure_sixthcodon.pdf")

# load data ---------------------------------------------------------------
tunney_bam_fname <- file.path(data_dir, "tunney_2018", "tunney",
                                     "tunney_bam.Rda")
load(tunney_bam_fname)
schuller_bam_fname <- file.path(data_dir, "schuller_2017", "schuller",
                                       "schuller_bam.Rda")
load(schuller_bam_fname)
lecanda_bam_fname <- file.path(data_dir, "lecanda_2016", "randomLinker_randomPrimer",
                                       "randomLinker_randomPrimer_bam.Rda")
load(lecanda_bam_fname)


from <- -20
to <- 40


# sum counts by nucleotide position and footprint length -------------------------------------------
meta_by <- function( bam ) {
  
  # convert the footprint counts to nucleotide indexes
  bam$nt_idx <- with( bam, cod_idx * 3 - d5 - 3)
  bam$length <- with( bam, d5 + d3 + 3)
  
  meta <- by( bam, list( bam$nt_idx, bam$length), function(x){ 

      c( nt_idx = unique(x$nt_idx),
      length = unique(x$length),
      cod_idx = unique(x$cod_idx),
      count = sum( x$count, na.rm=T ),
      corrected = sum( x$correct_250, na.rm=T ) )
  })
  
  meta <- data.frame( do.call( rbind, meta ))

  meta[ order( meta$cod_idx, meta$nt_idx, meta$length ), ]
}

# make a ribogrid metagene for the start region  -------------------------------------------
make_grid <- function( meta_df, from, to, column ) {
  
  combos <- expand.grid( nt_idx = from:to, length = min(meta_df$length) : max(meta_df$length) )
  
  meta <- merge( meta_df, combos, by = names(combos), all.x = F, all.y = T, sort = T)
  
  meta[ is.na(meta) ] = 0
  
  matrix( meta[ , column ], 
          nrow = length( min(meta$length) : max(meta$length) ),
          ncol = length( from:to ),
          dimnames = list( min(meta$length) : max(meta$length), from:to ))
}

plot_meta_start <- function( uncor_grid, cor_grid, from, to, name, row_idx ) {

  # plot just the footprint length that is most abundant with its 5' end at nt position 0 (A of the AUG)

  fp_size_at_0 <- row.names(uncor_grid)[ row_idx ]
  
  plot( from:to, uncor_grid[ row_idx, ]/1000,
        type = "l", col = "gray40",
        xlab = NA, ylab = NA,
#        xlab = "position of footprint 5' end (nt)",
#        ylab = "footprint count, thousands",
        axes = F, bty = "n" )
  axis( 1, tick = T, lwd = 0, lwd.ticks = 1 )
  axis( 2, tick = T, lwd = 0, lwd.ticks = 1 )
  lines( from:to, cor_grid[ row_idx, ]/1000, col = "red" )
  text( par("usr")[2], par("usr")[4], pos = 2, labels = paste0( name, " (", fp_size_at_0, " nt)" ))
}

generate_plot <- function( sample, name ) {
  
  sampledata <- get( paste0( sample, "_bam" ) )
  
  sample_meta <- meta_by( sampledata )
  uncor_grid <-  make_grid( sample_meta, from, to, "count" )
  cor_grid <-  make_grid( sample_meta, from, to, "corrected" )
  
  # figure out the footprint size that is most abundant with its 5' end at the A of the AUG, 
  # so we can see if the "bump" at that position is an artifact corrected well by choros
  at_0 = sample_meta[ sample_meta$nt_idx == 0, ]
  fp_size_at_0 = at_0$length[ which.max(at_0$count) ]
  row_idx <- which( row.names(uncor_grid) == fp_size_at_0 )
  
  plot_meta_start( uncor_grid, cor_grid, from, to, name, row_idx )
}

# plot the three datasets and save -------------------------------------------
pdf( pdfname, height = 3.5, width = 4.5)
par( xpd = NA )
par( mfrow = c(3,1) )
par( mar = c(2.6, 4.1, 1.1, 2.1) ) # default: c(5.1, 4.1, 4.1, 2.1)
par( oma = c(3, 1, 0, 0) )
generate_plot( "schuller", "Schuller" )
generate_plot( "tunney", "Tunney" )
title( ylab = "footprint count, thousands" )
generate_plot( "randomLinker_randomPrimer", "Lecanda")
title( xlab = "position of footprint 5' end (nt)" )

legend( "bottomright", legend = c("uncorrected","corrected"), 
        col = c("black", "red"), pch = 15, bty = "n", inset=c(0,-1.5), xpd = NA)

dev.off()
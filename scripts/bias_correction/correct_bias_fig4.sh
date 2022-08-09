#!/bin/bash

R CMD BATCH correct_bias_lecanda_process_data.R

R CMD BATCH correct_bias_lecanda_fixedLinker_fixedPrimer.R

R CMD BATCH correct_bias_lecanda_randomLinker_randomPrimer.R

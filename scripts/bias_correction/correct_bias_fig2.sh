#!/bin/bash

R CMD BATCH correct_bias_simulated_no_bias.R

R CMD BATCH correct_bias_simulated_f3_bias.R

R CMD BATCH correct_bias_simulated_f5_bias.R

R CMD BATCH correct_bias_simulated_f3_f5_bias.R

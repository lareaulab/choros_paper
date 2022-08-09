#!/bin/bash

R CMD BATCH correct_bias_wu_process_data.R

R CMD BATCH correct_bias_wu_WT.R

R CMD BATCH correct_bias_wu_3AT.R

R CMD BATCH correct_bias_wu_comparison.R

#!/bin/bash

# hicue extract test_out/ps_loops test_data/positions/SRP_dd_genes.bed test_data/matrices/Control.mcool::resolutions/1000 --circulars NC_014500.1 --windows 100000 --detrending ps --diag_mask 2000 --flip -t 8 --save_tmp --no-loci --display_strand --contact_range 5000 1000000 100000
hicue extract test_out/bed2d_test test_data/positions/plants_loops.bed2d test_data/matrices/Control.mcool::resolutions/1000 --circulars NC_014500.1 --windows 100000 --detrending ps --flip -t 8 --save_tmp --display_strand
# hicue extract test_out/human test_data/positions/some_human_genes.bed test_data/matrices/human_PL69.mcool::resolutions/1000 --windows 100000 --detrending ps --flip -t 8 --save_tmp --no-loci --display_strand
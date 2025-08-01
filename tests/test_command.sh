#!/bin/bash

hicue extract_opti test_out/bacteria test_data/SRP_dd_genes.bed test_data/Control.mcool::resolutions/1000 --loci --circulars NC_014500.1 --windows 100000 --detrending patch --flip
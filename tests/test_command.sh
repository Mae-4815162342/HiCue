#!/bin/bash

# hicue extract test_out/ps_loops test_data/positions/SRP_dd_genes.bed test_data/matrices/Control.mcool::resolutions/1000 --circulars NC_014500.1 --windows 100000 --detrending patch --flip -t 8 --save_tmp --no-loci --batch --display_strand
# hicue extract test_out/bed2d_test test_data/positions/plants_loops.bed2d test_data/matrices/Control.mcool::resolutions/1000 --circulars NC_014500.1 --windows 100000 --detrending patch --flip -t 8 --save_tmp --display_strand --nb_pos 3 --random_jitter 100
# hicue extract test_out/human test_data/positions/some_human_genes.bed test_data/others/human_PL69.mcool::resolutions/1000 --windows 100000 --detrending patch --flip -t 8 --save_tmp --no-loci --batch --display_strand 
# hicue extract test_out/SecB_pileups test_data/positions/SRP_dd_genes.bed test_data/others/banks.txt --windows 50000,100000 --circulars NC_014500.1 --detrending patch --flip -t 8 --no-loci --no-batch --display_strand --nb_pos 5
# hicue extract test_out/SecB_loops test_data/positions/SRP_dd_genes.bed test_data/others/banks.txt --windows 50000,100000 --circulars NC_014500.1 --detrending ps --flip -t 8 --no-loci --no-batch --display_strand --loops
# hicue extract test_out/SecB_mcool test_data/positions/SRP_dd_genes.bed test_data/others/banks_mcool.txt --binnings 5000,10000 --windows 50000,100000 --circulars NC_014500.1 --detrending patch --flip -t 8 --no-loci --no-batch --display_strand 

# hicue extract /data/Maelys/D_dadantii_analysis/directionnality/HEG_5perc /data/Maelys/D_dadantii_analysis/directionnality/HEG_5perc.gff /data/Maelys/D_dadantii_analysis/directionnality/dd_library.txt --circulars NC_014500.1 --save_tmp --windows 100000 --detrending patch --method median --flip --nb_pos 2 --display_strand
# hicue extract /data/Maelys/D_dadantii_analysis/directionnality/HEG_5perc_loop /data/Maelys/D_dadantii_analysis/directionnality/HEG_5perc.gff /data/Maelys/D_dadantii_analysis/directionnality/dd_library.txt --circulars NC_014500.1 --windows 100000 --detrending patch --method median --loops --separate_by direction --save_tmp

# hicue extract --method mean --min_dist 0 --windows 20000 \
#     /home/sardine/Bureau/example_agglo_cohesin/new_version_hicue \
#     /home/sardine/Bureau/example_agglo_cohesin/pairs_peaks_cohesins3.txt.bg2.10kb.50kb.2.bed2d \
#     /home/sardine/Bureau/example_agglo_cohesin/valid_idx_pcrfree.pairs.2000.cool \
#     --save_tmp \
#     --detrending patch \
#     --cmap_color seismic \
#     --cmap_limits -0.2 0.2 \
#     --nb_pos 10


# chromosight quantify \
#     /home/sardine/Bureau/example_agglo_cohesin/pairs_peaks_cohesins3.txt.bg2.10kb.50kb.2.bed2d \
#     /home/sardine/Bureau/example_agglo_cohesin/valid_idx_pcrfree.pairs.2000.cool \
#     /home/sardine/Bureau/example_agglo_cohesin/chromosight_test

# hicue tracks test_out/tracks_test test_data/tracks/WT.bw test_data/matrices/Control.mcool::resolutions/1000 --circulars NC_014500.1 --windows 50000 --detrending patch --nb_pos 1 --flip -t min 1500 --batch --loops --save_tmp --display_strand
# hicue tracks test_out/tracks_test/positions test_data/tracks/WT.bw test_data/matrices/Control.mcool::resolutions/1000 --positions test_data/positions/plants_loops.bed2d --circulars NC_014500.1 --windows 50000 --detrending patch --flip -p high 100 --save_tmp --no-loci --batch --display_strand
# hicue tracks test_out/tracks_test/positions test_data/tracks/Endive4_2025.bw test_data/matrices/Endive4_2025.mcool::resolutions/5000 --positions test_data/positions/plants_loops.bed2d --circulars NC_014500.1 --windows 50000 --detrending patch --flip -p high 100 --save_tmp --no-loci --batch --display_strand --cmap_limits -0.5 0.5
# hicue tracks test_out/tracks_test/positions test_data/tracks/Endive4_2025.bw test_data/matrices/Endive4_2025.mcool::resolutions/5000 --positions test_data/positions/plants_loops.bed --min_dist 100000 --loops --circulars NC_014500.1 --windows 50000 --detrending patch --flip -p high 100 --save_tmp --no-loci --batch --display_strand

# hicue regions test_out/regions_test_none test_data/regions/3D7_only_chrom7.bed test_data/regions/MicroC-Bartfai-Control.mcool::resolutions/5000 --loops --separate_by chroms --min_region_size 5000 --padding 1.0 --batch --no-display_log --format png,pdf --save_tmp --detrending none -e 20,51 # --cmap_limits -0.2 0.2
hicue regions test_out/regions_test_none test_data/regions/3D7_only_chrom7.bed test_data/regions/MicroC-Bartfai-Control.mcool::resolutions/5000 --diag_mask 5000 --separate_by chroms --min_region_size 5000 --padding 1.0 --batch --no-display_log --format png,pdf --save_tmp --detrending none -e 20,51 # --cmap_limits -0.2 0.2
# hicue regions test_out/regions_test_ps test_data/regions/3D7_only_chrom7.bed test_data/regions/MicroC-Bartfai-Control.mcool::resolutions/5000 --loops --separate_by chroms --min_region_size 5000 --padding 1.0 --batch --no-loci --no-display_log --format png,pdf --save_tmp --detrending ps -e 20,51 #--cmap_limits 0 2
# hicue regions test_out/regions_test_ps test_data/regions/3D7_only_chrom7.bed test_data/regions/MicroC-Bartfai-Control.mcool::resolutions/5000 --separate_by chroms --min_region_size 5000 --padding 1.0 --batch --no-loci --no-display_log --format png,pdf --save_tmp --detrending ps -e 20,51 #--cmap_limits 0 2
# hicue regions test_out/regions_test_patch test_data/regions/3D7_only_chrom7.bed test_data/regions/MicroC-Bartfai-Control.mcool::resolutions/5000 --loops --separate_by chroms --min_region_size 5000 --padding 1.0 --nb_pos 10 --batch --no-display_log --format png,pdf --save_tmp --detrending patch -e 20,51 #--cmap_limits 0 2
# hicue regions test_out/regions_test_patch test_data/regions/3D7_only_chrom7.bed test_data/regions/MicroC-Bartfai-Control.mcool::resolutions/5000 --separate_by chroms --min_region_size 5000 --padding 1.0 --nb_pos 10 --batch --no-display_log --format png,pdf --save_tmp --detrending patch -e 20,51 #--cmap_limits 0 2
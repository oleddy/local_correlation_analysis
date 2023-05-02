# Local correlation analysis
This repository contains microscopy data analysis scripts used in the manuscript "Immunopeptidomics reveals determinants of Mycobacterium tuberculosis antigen presentation on MHC class I." 

Preprint: https://www.biorxiv.org/content/10.1101/2022.08.30.505882v1

Published paper:

Owen Leddy, Forest M White, Bryan D Bryson (2023) Immunopeptidomics reveals determinants of Mycobacterium tuberculosis antigen presentation on MHC class I eLife 12:e84070. https://doi.org/10.7554/eLife.84070

The script correlation_analysis_multiprocess.py segments an intracellular compartment of interest in one channel (e.g., bacterium-containing phagosomes indicated by GFP expressed by the bacteria) and measures the local correlation (within a sliding window) between two channels within that compartment, averaging over all pixels within the compartment. It reads in the target directory of images for each experimental condition from a spreadsheet (see example_conditions_table.csv). This computation is done in parallel across the specified directories. The output is a CSV file for each experimental condition listing the size (in pixels) and average local correlation value for each object, which will be generated in the image directory corresponding to that condition. 

The script threshold_percents.py convert these CSV files of raw (size, average correlation) pairs for each given condition into a proportion of objects exceeding a specified threshold for local correlation (i.e., percent of objects co-localized). A threshold can also be set for the minimum size of a given object to help emove noise. 

The script pos_area_per_cell.py computes the average area per cell that exceeds a given threshold in one image channel.

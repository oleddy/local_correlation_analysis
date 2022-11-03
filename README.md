# Local correlation analysis
This repository contains microscopy data analysis scripts used in the manuscript "Immunopeptidomics reveals determinants of Mycobacterium tuberculosis antigen presentation on MHC class I" (https://www.biorxiv.org/content/10.1101/2022.08.30.505882v1).

The script correlation_analysis_multiprocess.py segments an intracellular compartment of interest in one channel (e.g., bacterium-containing phagosomes indicated by GFP expressed by the bacteria) and measures the local correlation (within a sliding window) between two channels within that compartment, averaging over all pixels within the compartment. It reads in the target directory of images for each experimental condition from a spreadsheet (see example_conditions_table.csv). This computation is done in parallel across the specified directories.

The script pos_area_per_cell.py computes the average area per cell that exceeds a given threshold in one image channel.

# README
Last updated 22 December 2021.

# CONTACT
Contact Benjamin Blonder (benjamin.blonder@berkeley.edu) with any questions.

# INSTRUCTIONS
All the analyses in the study can be replicated by running each of the companion R files in order:

+ analyze neon - 0 prep data.R
+ analyze neon - 1 prep topography.R
+ analyze neon - 2 analysis.R
+ analyze neon - 3 rf.R
+ analyze neon - 4 spectral separation confusion.R
+ analyze neon - 5 si damage by year.R
+ analyze neon - 6 aerial vs ground damage.R

The scripts require access to numerous datafiles that are available either in repositories we have made available publicly, or which are available on other public websites. The files are not included with this code due to the multi-GB size such an archive would require.

# DATASET DOWNLOADS
Required datasets that must be downloaded and renamed prior to analysis are:

+ Download from ESS-DIVE: https://data.ess-dive.lbl.gov/view/doi:10.15485/1779716
	+ 'layers/min_phase_cover.tif' (aspen cover)
	+ 'layers/min_phase_cytotype_medfilt-seived.tif' (cytotype predictions)

+ Download from ESS-DIVE: https://data.ess-dive.lbl.gov/view/doi%3A10.15485%2F1618131
	+ 'layers/min_phase_shade_tch_tiled.tif' (canopy height)
	+ 'layers/dtm_mosaic_min_phase_me.tif' (elevation)
	+ 'neon_wavelengths.txt' (spectral bands)

+ Download from ESS-DIVE: https://data.ess-dive.lbl.gov/view/doi:10.15485/1618133
	+ 'layers/N.tif' (nitrogen)
	+ 'layers/d13C.tif' (d13C)
	+ 'layers/N_unc.tif' (nitrogen uncertainty)
	+ 'layers/d13C_unc.tif' (d13C uncertainty)
	+ 'layers/min_phase_wtrl_tiled.tif' (cwc)

+ Download from ESS-DIVE: https://data.ess-dive.lbl.gov/view/ess-dive-b81a5f7de6ad201-20210413T173256168
	+ 'layers/ht_change_2015_2019_1m_masked.tif' (LiDAR height change)

+ Download from USFS: https://www.fs.usda.gov/detail/r2/forest-grasslandhealth/?cid=fsbdev3_041629 - download ZIP files for years 2000-2018, rename as appropriately labeled in script:
	+ 'usfs r2/*.shp' (USFS damage)

+ Download from Will Petry: https://github.com/wpetry/USTreeAtlas
	+ 'Little range map/poputrem.shp' (Range map of aspen)

+ Download from Dryad: https://doi.org/10.6078/D1898W
	+ 'aspen data site-level processed 30 Mar 2020.csv'

+ Available in this archive, but originally generated from https://data.ess-dive.lbl.gov/view/doi:10.15485/1779716:
	+ 'confusion matrix/*' (confusion matrix statistics)
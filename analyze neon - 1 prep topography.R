library(raster)
library(rgeos)
library(rgdal)

r_elev = raster('layers/dtm_mosaic_min_phase_me.tif')
r_aspect = terrain(r_elev,opt='aspect',unit='radians')
r_cos_aspect = cos(r_aspect)
r_slope = terrain(r_elev,opt='slope',unit='degrees')

writeRaster(r_cos_aspect, 'layers/r_cos_aspect.tif',options="COMPRESS=LZW",overwrite=TRUE)
writeRaster(r_slope, 'layers/r_slope.tif',options="COMPRESS=LZW",overwrite=TRUE)


r_ws = aggregate(r_elev,20) > -Inf
s_ws = rasterToPolygons(r_ws, dissolve=TRUE)
s_ws_simple = gSimplify(s_ws,tol=0.001)
shapefile(s_ws_simple,file='watershed.shp',overwrite=TRUE)





# clip traits
r_n = raster('layers/N.tif')
r_n_unc = raster('layers/N_unc.tif')
r_n[r_n_unc > 0.15 | r_n < 0 | r_n > 10] = NA
r_n = crop(r_n, extent(r_elev))
r_n = mask(r_n, r_elev)
writeRaster(r_n,file='r_trait_N_processed.tif',options="COMPRESS=LZW",overwrite=TRUE)

r_d13C = raster('layers/d13C.tif')
r_d13C_unc = raster('layers/d13C_unc.tif')
r_d13C[r_d13C_unc>0.8 | r_d13C < -30 | r_d13C > -22] = NA
r_d13C = crop(r_d13C, extent(r_elev))
r_d13C = mask(r_d13C, r_elev)
writeRaster(r_d13C,file='r_trait_d13C_processed.tif',options="COMPRESS=LZW",overwrite=TRUE)


#r_cwc = raster('layers/min_phase_wtrl_tiled.tif')
#r_lwc = raster('layers/LWC.tif')
#r_lwc_unc = raster('layers/LWC_unc.tif')
#r_lwc[r_lwc_unc>2 | r_lwc < 0 | r_lwc > 200] = NA
#r_lwc = crop(r_lwc, extent(r_elev))
#r_lwc = mask(r_lwc, r_elev)
#writeRaster(r_lwc,file='r_trait_lwc_processed.tif',options="COMPRESS=LZW",overwrite=TRUE)


# clip the damage mask
r_damage_masked = raster('r_damage_by_first_year.tif')
r_damage_masked = mask(r_damage_masked, r_elev)
writeRaster(r_damage_masked, 'r_damage_masked.tif',options="COMPRESS=LZW",overwrite=TRUE)


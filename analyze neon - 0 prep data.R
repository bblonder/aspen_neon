library(terra)
library(sf)
library(fasterize)

terraOptions(verbose=TRUE,todisk=TRUE,memfrac=0.8)

r_aspen_classification = rast('layers/min_phase_cover.tif')
r_cytotype = rast('layers/min_phase_cytotype_medfilt-seived.tif')
r_shade = rast('layers/min_phase_shade_tch_tiled.tif')

r_mask = (r_aspen_classification==0 & r_shade==1)
r_mask_na = clamp(r_mask, lower=0.9, values=FALSE)

r_cytotype_masked = mask(r_cytotype, r_mask_na)
r_cytotype_masked_clamped = clamp(r_cytotype_masked, lower=-0.0001, values=FALSE)
r_cytotype_masked_thresholded = (r_cytotype_masked_clamped > 0.5)

# get 5 m guess
r_cytotype_masked_thresholded_5m = disaggregate(aggregate(r_cytotype_masked_thresholded, 5, fun='median', na.rm=T), 5, method='near')
# replace NA values based on the guess
r_cytotype_masked_thresholded_covered = cover(r_cytotype_masked_thresholded, r_cytotype_masked_thresholded_5m)

# write raw output
writeRaster(r_cytotype_masked_thresholded, file='r_cytotype_masked.tif',overwrite=TRUE)
# write 5 m buffer output
writeRaster(r_cytotype_masked_thresholded_covered, file='r_cytotype_masked_buffer_5m.tif',overwrite=TRUE)


# get pixel coordinates of non-NA cytotypes
indices_cytotype =  which(!is.na(values(r_cytotype_masked_thresholded)))
write.csv(indices_cytotype, file='indices_cytotype.csv',row.names=F)
# get pixel coordinates of non-NA cytotypes 5 m buffer
indices_cytotype_buffer_5m =  which(!is.na(values(r_cytotype_masked_thresholded_covered)))
write.csv(indices_cytotype_buffer_5m, file='indices_cytotype_buffer_5m.csv',row.names=F)



# do ian mortality
r_ian_canopyheight = rast('layers/ht_change_2015_2019_1m_masked.tif')
r_ian_canopyheight_projected = project(x=r_ian_canopyheight, y=r_aspen_classification)
#r_ian_canopyheight_projected_median = focal(r_ian_canopyheight_projected, w=matrix(1/25, nc=5, nr=5), fun=median, na.rm=T)
# aggregate
r_ian_canopyheight_projected_clamped = clamp(r_ian_canopyheight_projected, lower=-100, values=FALSE)

writeRaster(r_ian_canopyheight_projected_clamped, file='r_ian_canopyheight_projected_clamped.tif',overwrite=TRUE)




# DO DAMAGE


rasterize_damage <- function(input)
{
  damage_projected = st_transform(input, crs(r_aspen_classification))
  damage_projected = st_cast(damage_projected, "MULTIPOLYGON")
  r_damage = fasterize(damage_projected, r_aspen_classification)
  #r_damage[is.na(r_damage[])] <- 0 
  #r_damage_masked = mask(r_damage, r_mask_na)  
  
  return(r_damage)
}


# load in all rasters

d_00 = st_read('usfs r2/2000 - stelprdb5183421/r200_dmg.shp')
d_00a = d_00[d_00$HOST1==746,]

d_01 = st_read('usfs r2/2001 - stelprdb5183422/r201_dmg.shp')
d_01a = d_01[d_01$HOST1==746,]

d_02 = st_read('usfs r2/2002 - stelprdb5183428/r202_dmg.shp')
d_02a = d_02[d_02$HOST1==746,]

d_03 = st_read('usfs r2/2003 - stelprdb5183433/r203_dmg.shp')
d_03a = d_03[d_03$HOST1==746,]

d_04 = st_read('usfs r2/2004 - stelprdb5183437/r204_dmg.shp')
d_04a = d_04[d_04$HOST1==746,]

d_05 = st_read('usfs r2/2005 - stelprdb5183443/r205_dmg.shp')
d_05a = d_05[d_05$HOST1==746,]

d_06 = st_read('usfs r2/2006 - stelprdb5183445/r206_dmg.shp')
d_06a = d_06[d_06$HOST1==746,]

d_07 = st_read('usfs r2/2007 - stelprdb5183448/r207_dmg.shp')
d_07a = d_07[d_07$HOST1==746,]

d_08 = st_read('usfs r2/2008 - stelprdb5183454/r208_dmg.shp')
d_08a = d_08[d_08$HOST1==746,]

d_09 = st_read('usfs r2/2009 - stelprdb5183456/r209_dmg.shp')
d_09a = d_09[d_09$HOST1==746,] # codes change this year

d_10 = st_read('usfs r2/2010 - stelprdb5247576/r210_dmg.shp')
d_10a = d_10[d_10$HOST1==746,]

d_11 = st_read('usfs r2/2011 - stelprdb5349026/r211_dmg.shp')
d_11a = d_11[d_11$HOST1==746,]

d_12 = st_read('usfs r2/2012 - stelprdb5409983/r212_dmg.shp')
d_12a = d_12[d_12$HOST1==746,]

d_13 = st_read('usfs r2/2013 - stelprdb5447183/r213_dmg.shp')
d_13a = d_13[d_13$HOST1==746,]

d_14 = st_read('usfs r2/2014 - stelprd3829398/r214_dmg.shp')
d_14a = d_14[d_14$HOST1==746,]

d_15 = st_read('usfs r2/2015 - fseprd490660/r215_dmg.shp')
d_15a = d_15[d_15$HOST1==746,]

d_16 = st_read('usfs r2/2016 - fseprd533949/r216_dmg.shp')
d_16a = d_16[d_16$HOST1==746,]

d_17 = st_read('usfs r2/2017 - fseprd574289/r217data2.gdb',layer='r217_dmg')
d_17a = d_17[d_17$HOST_CODE==746,]

d_18 = st_read('usfs r2/2018 - fseprd604735/R218_DMG.gdb',layer='r218_dmg')
d_18a = d_18[d_18$HOST_CODE==746,]; d_18a = d_18a[!is.na(d_18a$Area),]

# note we are not using other damage type codes (24032 vs 80001 vs 12900...)


damage_all = list(d_00a, d_01a, d_02a, d_03a, d_04a, d_05a, d_06a, d_07a, d_08a, d_09a, d_10a, d_11a, d_12a, d_13a, d_14a, d_15a, d_16a, d_17a, d_18a)
names(damage_all) = paste("Year",2000:2018,sep=".")
# save some memory
rm(list=ls(pattern="^d_[0-9]*"))

for (i in 1:length(damage_all))
{
  print(i)
  
  r_out = rasterize_damage(damage_all[[i]])
  writeRaster(r_out,file=sprintf('r_damage_%s.tif',names(damage_all)[i]),overwrite=TRUE, datatype="INT2U", options=c("COMPRESS=LZW"))
  
  png(width=800,height=1200,file=sprintf('r_damage_%s.png',names(damage_all)[i]))
  plot(r_out,axes=F)
  dev.off()
  
  rm(r_out)
}




# now combine all rasters
stack_damage = rast(dir(path='.',pattern='r*Year\\.[0-9][0-9][0-9][0-9]\\.tif$'))
damage_summed = app(stack_damage, fun=sum, na.rm=T)

writeRaster(damage_summed, 'r_damage_summed.tif',options="COMPRESS=LZW", overwrite=TRUE)
png(width=800,height=1200,file=sprintf('r_damage_summed.png'))
plot(damage_summed,axes=F)
dev.off()

# now modify to include 
stack_damage_copy = stack_damage
for (i in 1:nlyr(stack_damage_copy))
{
  year_this = as.numeric(gsub("r_damage_Year.","",names(stack_damage_copy)[i],fixed=TRUE))
  print(year_this)
  # convert 1 values to the year
  stack_damage_copy[[i]] = year_this * stack_damage_copy[[i]]
}

damage_by_first_year = app(stack_damage_copy, fun=min, na.rm=TRUE)

writeRaster(damage_by_first_year, 'r_damage_by_first_year.tif',options="COMPRESS=LZW", overwrite=TRUE)
png(width=800,height=1200,file=sprintf('r_damage_by_first_year.png'))
plot(damage_by_first_year,axes=F)
dev.off()










# OLD

#rollup <- st_read('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2018/usfs data/national rollup/Rollup_all.shp')
#ads_2018 = st_read('mortality data usfs/fseprd604735/R218_DMG.gdb',layer='r218_dmg')
#damage_aspen = ads_2018[(!is.na(ads_2018$HOST_CODE) & ads_2018$HOST_CODE==746) | (!is.na(ads_2018$HOST_CODE2) & ads_2018$HOST_CODE2==746),] 
#writeRaster(r_rollup_masked, file='r_rollup_masked.tif')

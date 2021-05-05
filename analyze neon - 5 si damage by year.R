library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(rasterVis)
library(ggpubr)
library(MASS)
library(DHARMa)
library(visreg)
library(ggExtra)
library(ggspatial)
library(ggsn)
library(RStoolbox)
library(wesanderson)
library(maps)
library(ggrepel)
library(rgdal)
library(mgcv)
library(terra)
library(RColorBrewer)
library(viridis)
library(pals)
library(mgcViz)
library(tidyr)

# load damage


r_lidar = rast('layers/ht_change_2015_2019_1m_masked.tif')

r_topo = rast('layers/r_cos_aspect.tif')

#r_lidar_proj = project(r_lidar, r_damage_usfs_masked)
#r_damage_usfs_masked_clipped_to_lidar = mask(r_damage_usfs_masked, r_lidar_proj)

result = NULL
for (year in 2000:2018)
{
  print(year)
  
  r_damage_usfs_masked = rast('r_damage_masked.tif')
  r_aspen = rast('layers/min_phase_cover.tif')==0
  
  r_this_year = (r_damage_usfs_masked==year)
  
  r_damage_not_in_neon = r_this_year & r_aspen
  r_damage_in_neon = r_this_year & !r_aspen
  
  count_damage_not_in_neon = freq(r_damage_not_in_neon, value=1)[1,"count"]
  count_damage_in_neon = freq(r_damage_in_neon, value=1)[1,"count"]
  result_this = data.frame(year=year, count_damage_not_in_neon = count_damage_not_in_neon, count_damage_in_neon=count_damage_in_neon)
  result = rbind(result, result_this)
  
  # memory management for terra (huge temp files)
  rm(r_this_year)
  rm(r_damage_not_in_neon)
  rm(r_damage_in_neon)
  tmpFiles(remove=TRUE)
}


result_short = result %>% gather(Variable, Value, 2:3)



#freq_damage_usfs = freq(r_damage_usfs_masked)
#freq_damage_usfs_clipped_to_lidar = freq(r_damage_usfs_masked_clipped_to_lidar)

freq_all = rbind(data.frame(freq_damage_usfs, Region='All'),data.frame(freq_damage_usfs_clipped_to_lidar, Region='LiDAR only'))

g_damage_by_year_usfs = ggplot(result_short,aes(x=year,y=Value,fill=Variable)) +
  scale_fill_brewer(palette='Set1', name='2018 aspen\ncover map',labels=c('Inside','Outside')) +
  geom_bar(stat='identity') +
  theme_bw() + 
  xlab("Year") + 
  ylab(expression(paste("Area impacted (km"^2,")"))) +
  scale_x_continuous(breaks=seq(2000,2018,by=1),limits=c(2000,2018)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

ggsave(g_damage_by_year_usfs, file='g_damage_by_year_usfs.png',width=7,height=4)

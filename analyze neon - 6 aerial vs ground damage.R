library(dplyr)
library(terra)
library(ggplot2)
library(ggpubr)
library(reshape2)

data_ground = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')
rasters_usfs_2015_2018 = lapply(c('r_damage_Year.2018.tif','r_damage_Year.2017.tif','r_damage_Year.2016.tif','r_damage_Year.2015.tif'), function(fn) {
  r = rast(fn)
  r[is.na(r)] = 0
  return(r)
  })
rasters_usfs = sum(do.call("c",rasters_usfs_2015_2018))

rasters_lidar = (rast('r_ian_canopyheight_projected_clamped.tif') < -3)

data_damage = data.frame(data_ground %>% dplyr::select(Site_Code, Count.adult.dead), 
                         damage_usfs=terra::extract(rasters_usfs, data_ground %>% dplyr::select(X.UTM, Y.UTM))[,2],
                         damage_lidar=terra::extract(rasters_lidar, data_ground %>% dplyr::select(X.UTM, Y.UTM))[,2]) %>%
  rename(site_code=Site_Code)

data_damage_melted = melt(data_damage, 
                          id.vars=c('site_code','Count.adult.dead')) %>%
                      na.omit

g_final_counts = ggplot(data_damage_melted, aes(x=Count.adult.dead/11,fill=factor(value),color=factor(value),group=factor(value))) + 
  geom_density() +
  facet_wrap(~variable,labeller=as_labeller(c(damage_usfs='USFS damage (2015-2018)',damage_lidar='Lidar damage (2015-2019)'))) + 
  xlab('Fraction dead trees (ground-based census)') +
  ylab('Probability density') +
  theme_bw() +
  scale_color_manual(values=c('darkblue','red'),name='Remotely-sensed\ndamage status') +
  scale_fill_manual(values=c(rgb(0,0,0.8,0.5),rgb(1,0,0,0.5)),name='Remotely-sensed\ndamage status') +
  xlim(0,1) +
  theme(legend.position = 'bottom')

ggsave(g_final_counts, file='g_final_counts.png',width=7,height=4)
ggsave(g_final_counts, file='g_final_counts.pdf',width=7,height=4)

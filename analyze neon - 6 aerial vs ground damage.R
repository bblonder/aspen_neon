library(dplyr)
library(terra)
library(ggplot2)
library(ggpubr)

data_ground = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')
rasters_usfs_2015_2018 = lapply(c('r_damage_Year.2018.tif','r_damage_Year.2017.tif','r_damage_Year.2016.tif','r_damage_Year.2015.tif'), function(fn) {
  r = rast(fn)
  r[is.na(r)] = 0
  return(r)
  })
rasters_usfs = sum(do.call("c",rasters_usfs_2016_2018))

rasters_lidar = (rast('r_ian_canopyheight_projected_clamped.tif') < -3)


data_damage = data.frame(data_ground %>% dplyr::select(Site_Code, Count.adult.dead, Count.adult.damaged, Unhealthy_Site), 
                         damage_usfs=terra::extract(rasters_usfs, data_ground %>% dplyr::select(X.UTM, Y.UTM))[,2],
                         damage_lidar=terra::extract(rasters_lidar, data_ground %>% dplyr::select(X.UTM, Y.UTM))[,2]) %>%
  rename(site_code=Site_Code)


g_lidar = ggplot(data_damage %>% filter(!is.nan(damage_lidar)),aes(x=Count.adult.dead,col=factor(damage_lidar))) +
  geom_density() +
  theme_bw() + 
  xlab("Number of dead trees") +
  ggtitle("Lidar, 2015-2018") +
  scale_color_manual(values=c("blue","red"),name='Canopy mortality')

g_usfs = ggplot(data_damage,aes(x=Count.adult.dead,col=factor(damage_usfs))) +
  geom_density() +
  theme_bw() + 
  xlab("Number of dead trees") +
  ggtitle("USFS 2015-2018") +
  scale_color_manual(values=c("blue","red"),name='Canopy mortality')

ggsave(ggarrange(g_usfs, g_lidar,common.legend = TRUE, align='hv',legend = 'right'),
       file='g_validation of remotely sensed mortality.png',width=8,height=5)

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


# load topography rasters
r_elev = rast('layers/dtm_mosaic_min_phase_me.tif')
r_cos_aspect = rast('layers/r_cos_aspect.tif')
r_slope = rast('layers/r_slope.tif')

# load trait maps (pre masked and clipped)
r_trait_n = rast('r_trait_N_processed.tif')
r_trait_d13C = rast('r_trait_d13C_processed.tif')
r_trait_lwc = rast('r_trait_lwc_processed.tif')
r_trait_cwc = rast('layers/min_phase_wtrl_tiled.tif')

# load height
r_height = rast('layers/tch_mosaic_min_phase_me.tif')

# load cytotype
r_cytotype_masked = rast('r_cytotype_masked.tif') # 0 = triploid, 1=diploid
r_cytotype_masked_buffer_5m = rast('r_cytotype_masked_buffer_5m.tif') # 0 = triploid, 1=diploid

# make aspen mask
r_aspen_cover = !is.na(r_cytotype_masked)

# load damage
r_damage_usfs_masked = rast('r_damage_masked.tif')
r_damage_lidar_masked = rast('r_ian_canopyheight_projected_clamped.tif')

r_damage_usfs_masked_binary = !is.na(r_damage_usfs_masked)
r_damage_lidar_masked_binary = r_damage_lidar_masked < -3 # set 3m cutoff

# load in watershed boundaries
s_ws = shapefile('watershed.shp')
s_ws_df = fortify(s_ws)


# count area
#v_elev = getValues(r_elev)
#length(which(!is.na(v_elev)))

# stack up predictors
s_all = c(r_aspen_cover,
          r_cytotype_masked,
          r_cytotype_masked_buffer_5m,
          #r_rollup_masked, 
          r_damage_usfs_masked_binary,
          r_damage_lidar_masked_binary,
          r_elev, 
          r_cos_aspect, 
          r_slope, 
          r_trait_n, 
          #r_lma_cropped, 
          r_trait_d13C,
          #r_trait_lwc,
          r_trait_cwc,
          r_height)
names(s_all) = c("Aspen.Cover",
                 "Cytotype",
                 "Cytotype.buffer.5m",
                 "Damage.USFS",
                 "Damage.Lidar",
                 "Elevation",
                 "Cos.aspect",
                 "Slope",
                 "Trait.N",
                 #"Trait.LMA",
                 "Trait.d13C",
                 #"Trait.LWC",
                 "Trait.CWC",
                 "Height.Canopy")

s_all_10m = aggregate(s_all, fact=10,fun='mean',na.rm=TRUE)

indices_cytotype_10m = which(s_all_10m[["Aspen.Cover"]][]>0)
# pick a sample of cells that do contain aspen
indices_cytotype = read.csv('indices_cytotype.csv')
#indices_cytotype_sample = indices_cytotype$x
#indices_cytotype_buffer_5m = read.csv('indices_cytotype_buffer_5m.csv')
#indices_cytotype_buffer_5m_sample = indices_cytotype_buffer_5m$x


generate_df <- function(s_all_this, indices)
{
  df = s_all_this[indices]
  # add x-y coordinates
  df = cbind(df, xyFromCell(s_all_this, cell=indices))
  
  #df = df %>% rename(Cytotype.nobuffer = Cytotype)
  
  df = df %>%
    rename(Aspen.Cover.fraction=Aspen.Cover,
           Cytotype.fractionDiploid=Cytotype, 
           Damage.USFS.fraction=Damage.USFS,
           Damage.Lidar.fraction=Damage.Lidar) %>%
    dplyr::select(-Cytotype.buffer.5m)
  
  # %>%
    #mutate(Cytotype.nobuffer=factor(ifelse(Cytotype.nobuffer > 0.5,"diploid","triploid"))) %>%
    #mutate(Cytotype.buffer.5m=factor(ifelse(Cytotype.buffer.5m > 0.5,"diploid","triploid"))) %>%
    #mutate(Damage.Year = Damage.USFS) %>%
    #mutate(Damage.USFS = as.numeric(!is.na(Damage.USFS)))
  
  return(df)
}

# extract values
df_all_10m = generate_df(s_all_10m, indices_cytotype_10m)
#df_all_10m = df_all_10m %>% dplyr::select(-Cytotype.buffer.5m) %>% rename(Cytotype = Cytotype.nobuffer)

df_all = generate_df(s_all, indices_cytotype$x)
#df_all = df_all %>% dplyr::select(-Cytotype.buffer.5m) %>% rename(Cytotype = Cytotype.nobuffer)
#df_all_buffer_5m = generate_df(indices_cytotype_buffer_5m_sample)
#df_all_buffer_5m = df_all_buffer_5m %>% dplyr::select(-Cytotype.nobuffer) %>% rename(Cytotype = Cytotype.buffer.5m)



# count values
#df_all$Cytotype %>% table / sum(!is.na(df_all$Cytotype))
#df_all_10m$Cytotype %>% table / sum(!is.na(df_all_10m$Cytotype))



#save(df_all,file='df_all.Rdata')
save(df_all_10m,file='df_all_10m.Rdata')






# FIGURE 1
palette = colorRampPalette(brewer.pal(9,"Spectral"))(100)
palette2 = colorRampPalette(magma(100))(100)
palette_trait = colorRampPalette(parula(100))(100)

g_map_elev = ggR(raster(r_elev),geom_raster=TRUE) + 
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_gradientn(colors=palette,na.value='white') +
  
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm")) +
  labs(fill='Elevation (m)') 





s_ws_ll = spTransform(s_ws,  "+init=epsg:4121 +proj=longlat +ellps=GRS80")
s_ws_ll_df = fortify(s_ws_ll)

USA <- map_data("state")

aspen = shapefile('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2020/aspen site selection for roots/Little range map/poputrem.shp')
crs(aspen) = "+init=epsg:4121 +proj=longlat +ellps=GRS80"
aspen_df = broom::tidy(aspen)

ws_df = broom::tidy(ws_boundaries)

palette_map = wes_palette("Chevalier1")

data_cities = us.cities %>% 
  filter(country.etc %in% "CO") %>%
  mutate(name=gsub(" CO","",name)) %>%
  filter(name %in% c("Denver", "Fort Collins", "Pueblo"))

NAmap <- ggplot() + geom_polygon(data = USA, 
                                 aes(x=long, y = lat, group = group), 
                                 fill = 'white', 
                                 color=palette_map[4]) +
  geom_polygon(data = aspen_df %>% filter(hole==FALSE), aes(x=long, y=lat, group=group),alpha=0.9,color=NA,fill='lightgreen') +
  geom_polygon(data = aspen_df %>% filter(hole==TRUE), aes(x=long, y=lat, group=group),color=NA,fill='white') +
  geom_label_repel(data=data_cities,aes(x=long,y=lat,label=name),
                   size = 2,
                   nudge_x      = 0.15,
                   direction    = "y",
                   hjust        = 0,
                   label.padding=0.05,
                   alpha=0.5) +
  geom_polygon(data = s_ws_ll_df, aes(x=long, y=lat, group=group),fill='purple',col='purple',size=0.01) +
  geom_point(data=data_cities,aes(x=long,y=lat),size=0.1) +
  theme_bw() +
  coord_equal(ylim=c(34.5,42),xlim = c(-109,-102)) +
  xlab("Longitude (°)") +
  ylab("Latitude (°)")







data_ground = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')


inset_x = 327571
inset_y = 4309454
delta_2 = 200
df_inset1 = data.frame(x=inset_x,y=inset_y,xmin=inset_x-delta_2, xmax=inset_x+delta_2,ymin=inset_y-delta_2,ymax=inset_y+delta_2)



g_map_cytotype = gplot(raster(r_cytotype_masked),maxpixels=1e6) + 
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA) +
  geom_tile(aes(fill=ifelse(value>0.5,"diploid","triploid"))) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_manual(values=c("blue","red"),name='Cytotype',na.translate=FALSE) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm")) +
  geom_rect(data=df_inset1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE)

r_cytotype_inset1 = crop(r_cytotype_masked,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))

g_map_cytotype_inset1 = gplot(raster(r_cytotype_inset1),maxpixels=1e6) + 
  geom_tile(aes(fill=ifelse(value>0.5,"diploid","triploid"))) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_manual(values=c("blue","red"),name='Cytotype',na.translate=FALSE) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm"))

g_map_cytotype_all = ggarrange(NAmap, cowplot::ggdraw(), g_map_cytotype, g_map_cytotype_inset1,
                               nrow=2,ncol=2,labels=c('(a)','','(b)','(c)'),align='hv',common.legend = TRUE,legend='bottom')

#ggsave(g_map_cytotype_all, file='g_map_cytotype.pdf',width=12,height=12)
ggsave(g_map_cytotype_all, file='g_map_cytotype.png',width=6,height=7)
rm(g_map_cytotype_all)


# see https://datacarpentry.org/r-raster-vector-geospatial/02-raster-plot/
# https://rdrr.io/github/statnmap/cartomisc/src/R/gplot_data.R
gplot_data <- function(y, maxpixels = 50000)  {
  #x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  x <- spatSample(y, maxpixels, as.raster = TRUE)
  #coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  coords <- xyFromCell(x, seq_len(ncell(x)))
  ## Extract values
  #dat <- utils::stack(as.data.frame(raster::getValues(x)))
  dat <- utils::stack(as.data.frame(values(x)))
  
  names(dat) <- c('value', 'variable')
  # If only one variable
  #if (dat$variable[1] == "raster::getValues(x)") {
  #if (dat$variable[1] == "values(x)") {
  #  dat$variable <- names(x)
  #}
  
  dat <- dplyr::as_tibble(data.frame(coords, dat))
  
  #if (!is.null(levels(x))) {
  #  dat <- dplyr::left_join(dat, levels(x)[[1]],
  #                          by = c("value" = "ID"))
  #}
  dat
}


r_cytotype_masked_cropped = crop(r_cytotype_masked,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))
df_aspen_for_map = gplot_data(r_cytotype_masked,maxpixels=2e7)
df_aspen_for_map_cropped = gplot_data(r_cytotype_masked_cropped,maxpixels=1e7)

df_damage_usfs_for_map = gplot_data(r_damage_usfs_masked,maxpixels=1e7)

r_damage_usfs_masked_inset = crop(r_damage_usfs_masked,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))
df_damage_usfs_for_map_inset = gplot_data(r_damage_usfs_masked_inset,maxpixels=1e7)


r_damage_lidar_clipped = mask(r_damage_lidar_masked < -3, r_cytotype_masked)#_buffer_5m) # use the buffer for visuals
df_damage_lidar_for_map = gplot_data(r_damage_lidar_clipped,maxpixels=1e7) # pick 3m cutoff

r_damage_lidar_inset = crop(r_damage_lidar_clipped,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))
df_damage_lidar_inset = gplot_data(r_damage_lidar_inset,maxpixels=1e7)


g_damage = ggplot(df_aspen_for_map,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray') + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_damage_usfs_for_map,aes(x=x,y=y,fill=value),inherit.aes = FALSE, alpha=0.75) +
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA,inherit.aes=FALSE) +
  scale_fill_gradientn(colors=palette,na.value=NA,name='Damage year',breaks=c(2006,2008,2010,2012,2014,2016,2018)) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  geom_rect(data=df_inset1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE)


g_damage_inset = ggplot(df_aspen_for_map_cropped,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray') + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_damage_usfs_for_map_inset,aes(x=x,y=y,fill=value),inherit.aes = FALSE, alpha=0.75) +
  scale_fill_gradientn(colors=palette,na.value=NA,name='Damage year',breaks=c(2006,2008,2010,2012,2014,2016,2018)) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm"))

g_lidar = ggplot(df_aspen_for_map,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray') + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_damage_lidar_for_map,aes(x=x,y=y,fill=value),inherit.aes = FALSE) + # choose 3m cutoff
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA,inherit.aes=FALSE) +
  scale_fill_gradientn(colors=c('lightblue','red'),na.value=NA,name='Damage',breaks=c(0,1)) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  geom_rect(data=df_inset1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE)

g_lidar_inset = ggplot(df_aspen_for_map_cropped,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray') + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_damage_lidar_inset,aes(x=x,y=y,fill=value),inherit.aes = FALSE, alpha=0.75) +
  scale_fill_gradientn(colors=c('lightblue','red'),na.value=NA,name='Damage',breaks=c(0,1)) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm"))

# note line below needs the plotting raster in place (this gets run below...)
df_height_for_map = gplot_data(r_height_for_plotting,maxpixels=1e7)
df_height_for_map$value[is.nan(df_height_for_map$value)] <- NA

g_map_height = ggplot() + 
  geom_tile(data=df_height_for_map,aes(x=x,y=y,fill=value)) +
  scale_fill_viridis(name='Canopy height (m)',na.value=NA) +
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA,inherit.aes=FALSE) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm"))

g_damage_final = ggarrange(g_map_elev, g_map_height, g_damage, g_lidar, g_damage_inset, g_lidar_inset, nrow=3,ncol=2,labels=c('(a)','(b)','(c)','(d)','(e)','(f)'),align='hv')
ggsave(g_damage_final, file='g_damage_final.pdf',width=8,height=10)
ggsave(g_damage_final, file='g_damage_final.png',width=8,height=10)
rm(g_damage_final)





# show ground based plots relative to remotely sensed data



inset_A_x = 328091
inset_A_y = 4309700
inset_B_x = 322848
inset_B_y = 4311324  
delta_inset = 700
df_A = data.frame(xmin=inset_A_x-delta_inset,xmax=inset_A_x+delta_inset,ymin=inset_A_y-delta_inset,ymax=inset_A_y+delta_inset,x=inset_A_x,y=inset_A_y,label='b')
df_B = data.frame(xmin=inset_B_x-delta_inset,xmax=inset_B_x+delta_inset,ymin=inset_B_y-delta_inset,ymax=inset_B_y+delta_inset,x=inset_B_x,y=inset_B_y,label='c')

g_map_cytotype_ground = ggplot(df_aspen_for_map,aes(x=x,y=y,fill=value)) + 
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA) +
  geom_tile(aes(fill=ifelse(value>0.5,"diploid","triploid"))) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_manual(values=c("blue","red"),name='Cytotype (predicted)',na.translate=FALSE) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  geom_point(data=data_ground, aes(x=X.UTM,y=Y.UTM,color=tolower(Ploidy_level)),inherit.aes=FALSE,alpha=0.6,size=1) +
  scale_color_manual(values=c("purple","orange"),name='Cytotype (ground truth)',na.translate=FALSE) +  
  geom_rect(data=df_A,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE) +
  geom_text(data=df_A,aes(x=x,y=y,label=label),color='black',inherit.aes = FALSE) +
  geom_rect(data=df_B,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE) +
  geom_text(data=df_B,aes(x=x,y=y,label=label),color='black',inherit.aes = FALSE)

r_cytotype_groundtruth_cropped_A = crop(r_cytotype_masked,extent(c(inset_A_x-delta_inset, inset_A_x+delta_inset,inset_A_y-delta_inset,inset_A_y+delta_inset)))
df_cytotype_for_map_A = gplot_data(r_cytotype_groundtruth_cropped_A,maxpixels=5e6)
data_ground_A = data_ground %>%
  filter(X.UTM >= inset_A_x-delta_inset & X.UTM <= inset_A_x+delta_inset) %>%
  filter(Y.UTM >= inset_A_y-delta_inset & Y.UTM <= inset_A_y+delta_inset)

g_map_cytotype_ground_A = ggplot(df_cytotype_for_map_A,aes(x=x,y=y,fill=value)) + 
  geom_tile(aes(fill=ifelse(value>0.5,"diploid","triploid"))) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_manual(values=c("blue","red"),name='Cytotype (predicted)',na.translate=FALSE) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  geom_point(data=data_ground_A, aes(x=X.UTM,y=Y.UTM,color=tolower(Ploidy_level)),inherit.aes=FALSE,alpha=0.6,size=1) +
  scale_color_manual(values=c("purple","orange"),name='Cytotype (ground truth)',na.translate=FALSE) +
  theme(legend.position = 'none')


r_cytotype_groundtruth_cropped_B = crop(r_cytotype_masked,extent(c(inset_B_x-delta_inset, inset_B_x+delta_inset,inset_B_y-delta_inset,inset_B_y+delta_inset)))
df_cytotype_for_map_B = gplot_data(r_cytotype_groundtruth_cropped_B,maxpixels=5e6)
data_ground_B = data_ground %>%
  filter(X.UTM >= inset_B_x-delta_inset & X.UTM <= inset_B_x+delta_inset) %>%
  filter(Y.UTM >= inset_B_y-delta_inset & Y.UTM <= inset_B_y+delta_inset)

g_map_cytotype_ground_B = ggplot(df_cytotype_for_map_B,aes(x=x,y=y,fill=value)) + 
  geom_tile(aes(fill=ifelse(value>0.5,"diploid","triploid"))) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_manual(values=c("blue","red"),name='Cytotype (predicted)',na.translate=FALSE) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  geom_point(data=data_ground_B, aes(x=X.UTM,y=Y.UTM,color=tolower(Ploidy_level)),inherit.aes=FALSE,alpha=0.6,size=1) +
  scale_color_manual(values=c("purple","orange"),name='Cytotype (ground truth)',na.translate=FALSE) +
  theme(legend.position = 'none')


g_map_cytotye_ground_all = ggarrange(g_map_cytotype_ground, ggarrange(g_map_cytotype_ground_A, g_map_cytotype_ground_B,nrow=2,ncol=1,labels=c('b','c')), 
                                     nrow=1,ncol=2,widths = c(2,1),labels='a',common.legend = TRUE,legend='bottom')

ggsave(g_map_cytotye_ground_all, file='g_fig_groundtruth.pdf', width=8,height=8)
ggsave(g_map_cytotye_ground_all, file='g_fig_groundtruth.png', width=8,height=8,dpi=600)






warning('redo at 10m?')
# look at cytotype distributions across mortality classes
xt.usfs = as.data.frame(xtabs(~Cytotype.fractionDiploid + Damage.USFS.fraction, data=df_all)) %>% 
  #filter(Damage==TRUE) %>%
  group_by(Damage.USFS.fraction) %>%
  mutate(NormFreq = Freq/sum(Freq))
xt.usfs.totals = xt.usfs %>% group_by(Damage.USFS.fraction) %>% summarize(NumPixels = sum(Freq))
# use 3m cutoff
xt.lidar = as.data.frame(xtabs(~Cytotype.fractionDiploid + Damage.Lidar.fraction,data=df_all)) %>% 
  #filter(Damage==TRUE) %>%
  group_by(Damage.Lidar.fraction) %>%
  mutate(NormFreq = Freq/sum(Freq))
xt.lidar.totals = xt.lidar %>% group_by(Damage.Lidar.fraction) %>% summarize(NumPixels = sum(Freq))

g_damage_usfs_by_cytotype = ggplot(xt.usfs,aes(x=(Damage.USFS.fraction==1),y=NormFreq,fill=Cytotype.fractionDiploid)) + 
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(values=c("blue","red"),name='Cytotype',breaks=0:1,labels=c('Diploid','Triploid')) +
  ylab("Fraction of area with damage") +
  xlab("Aerially-surveyed damage (2000-2018)") +
  geom_text(data=xt.usfs.totals, aes(label=NumPixels,x=(Damage.USFS.fraction==1),y=0.05),inherit.aes=FALSE) +
  theme_bw()# +
#theme(legend.position = 'none')
g_damage_lidar_by_cytotype = ggplot(xt.lidar,aes(x=Damage.Lidar.fraction,y=NormFreq,fill=Cytotype.fractionDiploid)) + 
  geom_bar(stat='identity',position='stack') +
  scale_fill_manual(values=c("blue","red"),name='Cytotype',breaks=0:1,labels=c('Diploid','Triploid')) +
  ylab("Fraction of area with damage") +
  xlab("Lidar-detected damage (2015-2019)") +
  geom_text(data=xt.lidar.totals, aes(label=NumPixels,x=Damage.Lidar.fraction,y=0.05),inherit.aes=FALSE) +
  theme_bw()

g_inset_xt = ggarrange(g_damage_usfs_by_cytotype, g_damage_lidar_by_cytotype,
                        nrow=1,ncol=2,common.legend = TRUE, legend='bottom',labels=c('(a)','(b)'))
ggsave(g_inset_xt,file='g_inset_xt.pdf',width=8,height=5)
ggsave(g_inset_xt,file='g_inset_xt.png',width=8,height=5)

# do chi square tests
chisq.test(xtabs(~Cytotype.fractionDiploid + Damage.USFS.fraction, data=df_all))
chisq.test(xtabs(~Cytotype.fractionDiploid + Damage.Lidar.fraction, data=df_all))
# get a smaller subset for effective model fitting
#set.seed(1)
#df_all_ss = df_all %>%
#  filter(Elevation >= quantile(Elevation,0.01,na.rm=T) & 
#           Elevation <= quantile(Elevation,0.99,na.rm=T) & 
#           Slope >= quantile(Slope,0.01,na.rm=T) & 
#           Slope <= quantile(Slope,0.99,na.rm=T) & 
#           Cos.aspect >= quantile(Cos.aspect,0.01,na.rm=T) & 
#           Cos.aspect <= quantile(Cos.aspect,0.99,na.rm=T) & 
#           x >= quantile(x,0.01,na.rm=T) & 
#           x <= quantile(x,0.99,na.rm=T) &
#           y >= quantile(y,0.01,na.rm=T) & 
#           y <= quantile(y,0.99,na.rm=T))


# clear some memory
rm(g_damage)
rm(g_lidar)
rm(g_map_cytotype_all)
rm(g_map_cytotype_ground)

#df_all_buffer_5m_ss = df_all_buffer_5m %>% 
#  dplyr::select(-Damage.Year) %>% 
#  na.omit

#m_damage = bam(factor(Damage.USFS)~Cytotype + s(Elevation,by=Cytotype) + s(Slope,by=Cytotype) + s(Cos.aspect,by=Cytotype) + s(x,y),
#               family='binomial',
#               data=df_all_ss,
#               control = gam.control(trace = TRUE),
#               discrete=TRUE)

#https://converged.yt/talks/creemcrackers-splines/talk.pdf




# make trait maps
r_height_for_plotting = mask(r_height, r_cytotype_masked)
r_n_for_plotting = mask(r_trait_n, r_cytotype_masked)
r_d13C_for_plotting = mask(r_trait_d13C, r_cytotype_masked)
r_cwc_for_plotting = mask(r_trait_cwc, r_cytotype_masked)

q_h = quantile(spatSample(r_height_for_plotting,10000), c(0.01,0.99), na.rm=T)
q_n = quantile(spatSample(r_n_for_plotting,10000), c(0.01,0.99), na.rm=T)
q_d13C = quantile(spatSample(r_d13C_for_plotting,10000, na.rm=TRUE), c(0.01,0.99), na.rm=T)
q_cwc = quantile(spatSample(r_cwc_for_plotting,10000, na.rm=TRUE), c(0.01,0.99), na.rm=T)

df_n_for_map = gplot_data(r_n_for_plotting,maxpixels=1e7)
df_d13c_for_map = gplot_data(r_d13C_for_plotting,maxpixels=1e7)
df_cwc_for_map = gplot_data(r_cwc_for_plotting,maxpixels=1e7)



g_map_trait_N = ggplot(df_n_for_map,aes(x=x,y=y,fill=value)) + 
  geom_tile() + 
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA,inherit.aes=FALSE) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_gradientn(colors=palette_trait,na.value=NA,limits=q_n) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm")) +
  labs(fill='Nitrogen (%)') +
  theme(legend.position='bottom') +
  geom_rect(data=df_inset1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE)

g_map_trait_d13C = ggplot(df_d13c_for_map,aes(x=x,y=y,fill=value)) + 
  geom_tile() + 
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA,inherit.aes=FALSE) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_gradientn(colors=palette_trait,na.value=NA,limits=q_d13C) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm")) +
  labs(fill=expression(paste(delta^{13},'C (per mil)')))  +
  theme(legend.position='bottom') +
  geom_rect(data=df_inset1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE)

g_map_trait_cwc = ggplot(df_cwc_for_map,aes(x=x,y=y,fill=value)) + 
  geom_tile() + 
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='gray',fill=NA,inherit.aes=FALSE) +
  xlab("") + ylab("") +
  coord_equal() +
  scale_fill_gradientn(colors=palette_trait,na.value=NA,limits=q_cwc) +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm")) +
  labs(fill=expression(paste('Canopy water content (mL m'^{-2},')')))  +
  theme(legend.position='bottom') +
  geom_rect(data=df_inset1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=NA,color='black',inherit.aes = FALSE)




r_cwc_inset = crop(r_cwc_for_plotting,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))
df_cwc_inset = gplot_data(r_cwc_inset,maxpixels=1e7)

g_cwc_inset = ggplot(df_aspen_for_map_cropped,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray',show.legend = FALSE) + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_cwc_inset,aes(x=x,y=y,fill=value),inherit.aes = FALSE, alpha=0.75) +
  scale_fill_gradientn(colors=palette_trait,na.value=NA,limits=q_cwc) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  labs(fill=expression(paste('Canopy water content (mL m'^{-2},')')))  +
  theme(legend.position='bottom')

r_d13c_inset = crop(r_d13C_for_plotting,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))
df_d13c_inset = gplot_data(r_d13c_inset,maxpixels=1e7)

g_d13c_inset = ggplot(df_aspen_for_map_cropped,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray',show.legend = FALSE) + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_d13c_inset,aes(x=x,y=y,fill=value),inherit.aes = FALSE, alpha=0.75) +
  scale_fill_gradientn(colors=palette_trait,na.value=NA,limits=q_d13C) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  labs(fill=expression(paste(delta^{13},'C (per mil)')))  +
  theme(legend.position='bottom')


r_n_inset = crop(r_n_for_plotting,extent(c(inset_x-delta_2, inset_x+delta_2,inset_y-delta_2,inset_y+delta_2)))
df_n_inset = gplot_data(r_n_inset,maxpixels=1e7)

g_n_inset = ggplot(df_aspen_for_map_cropped,aes(x=x,y=y,alpha=factor(1-as.numeric(is.na(value))))) + 
  geom_tile(fill='gray',show.legend = FALSE) + 
  scale_alpha_discrete(na.value=1,range=c(0,1),name='Aspen cover') +
  geom_tile(data=df_n_inset,aes(x=x,y=y,fill=value),inherit.aes = FALSE, alpha=0.75) +
  scale_fill_gradientn(colors=palette_trait,na.value=NA,limits=q_n) +
  theme_bw() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  labs(fill='Nitrogen (%)') +
  theme(legend.position='bottom')




g_map_traits = ggarrange(g_map_trait_cwc, g_map_trait_d13C, g_map_trait_N, 
                         g_cwc_inset, g_d13c_inset, g_n_inset,
                         nrow=2,ncol=3,align='hv',labels=c('(a)','(b)','(c)','(d)','(e)','(f)'))
ggsave(g_map_traits,file='g_map_traits.pdf',width=10,height=7)
ggsave(g_map_traits,file='g_map_traits.png',width=10,height=7)








# do cytotype cuts
df_all_cut = df_all_10m %>% 
  mutate(Elevation.cut=cut(Elevation,breaks=seq(2600,3700,by=200),dig.lab=10)) %>%
  mutate(Slope.cut=cut(Slope,breaks=seq(0,50,by=10),dig.lab=10)) %>%
  mutate(Cos.aspect.cut=cut(Cos.aspect,breaks=seq(-1,1,by=0.5),dig.lab=10)) %>%
  mutate(Cytotype = ifelse(Cytotype.fractionDiploid>0.5,'diploid','triploid'))

summaries_cut_elevation = df_all_cut %>% 
  ungroup %>%
  group_by(Cytotype, Elevation.cut) %>% 
  summarize(area=n())

summaries_cut_slope = df_all_cut %>% 
  ungroup %>%
  group_by(Cytotype, Slope.cut) %>% 
  summarize(area=n())

summaries_cut_cos_aspect = df_all_cut %>% 
  ungroup %>%
  group_by(Cytotype, Cos.aspect.cut) %>% 
  summarize(area=n())

g_cut_elevation = ggplot(summaries_cut_elevation, aes(x=Cytotype,y=area,fill=Elevation.cut)) + 
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_brewer(palette='GnBu',name='Elevation (m)',na.translate=FALSE) +
  ylab(expression(paste("Area (m"^2,")")))

g_cut_slope = ggplot(summaries_cut_slope, aes(x=Cytotype,y=area,fill=Slope.cut)) + 
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_brewer(palette='GnBu',name='Slope (°)',na.translate=FALSE) +
  ylab(expression(paste("Area (m"^2,")")))

g_cut_cos_aspect = ggplot(summaries_cut_cos_aspect, aes(x=Cytotype,y=area,fill=Cos.aspect.cut)) + 
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_brewer(palette='GnBu',name='Cosine aspect\n(dimensionless)',na.translate=FALSE) +
  ylab(expression(paste("Area (m"^2,")")))

g_cut = ggarrange(g_cut_elevation, g_cut_slope, g_cut_cos_aspect,nrow=2,ncol=2,align='hv',labels=c('(a)','(b)','(c)'))
ggsave(g_cut, file='g_cut.pdf',width=12,height=12)
ggsave(g_cut, file='g_cut.png',width=12,height=12)








# write rasters for GEE
writeRaster(r_height_for_plotting, file='gee_canopy_height.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_n_for_plotting, file='gee_leaf_n.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_d13C_for_plotting, file='gee_leaf_d13C.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_cwc_for_plotting, file='gee_leaf_cwc.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_lwc_for_plotting, file='gee_leaf_lwc.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_slope, file='gee_topo_slope.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_cos_aspect, file='gee_topo_cos_aspect.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_elev, file='gee_topo_elevation.tif', options="COMPRESS=LZW", overwrite=TRUE)
writeRaster(r_cytotype_masked, file='gee_cytotype.tif', options="COMPRESS=LZW", overwrite=TRUE, datatype='INT1S')
writeRaster(r_damage_usfs_masked, file='gee_damage_usfs.tif', options="COMPRESS=LZW", overwrite=TRUE, datatype='INT1S')
writeRaster(r_damage_lidar_masked, file='gee_damage_lidar.tif', options="COMPRESS=LZW", overwrite=TRUE, datatype='INT1S')
            
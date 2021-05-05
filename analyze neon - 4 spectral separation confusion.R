library(ggplot2)
library(tidyr)
library(dplyr)
library(raster)
library(dplyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(ggpubr)

refl = read.csv('combined_clone_spectral_data_buffer_5_roi.csv')
refl = refl %>% dplyr::select(starts_with("refl_band"), site_code, ploidy_level, clone, x, y, wtrl, wtrv, shade, tch, cover)
# trim bad bands
bad_bands = paste("refl_band",1+c(0:8, 192:205, 284:327, 417:425),sep="_")
refl[,bad_bands] = NA

wl = read.csv('neon_wavelengths.txt',header=F)$V1
names(refl)[grep("refl_band*",names(refl))] = paste("W",wl,sep="")


# add other metadata
# sex
#data_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex.csv') %>%
#  rename(site_code=Site_Code)

# aerial damage
data_ground = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')
map_damage_usfs = raster('r_damage_by_first_year.tif')
map_damage_lidar = (raster('r_ian_canopyheight_projected_clamped.tif') < -3)
data_damage = data.frame(data_ground %>% dplyr::select(Site_Code), 
                         damage_usfs=raster::extract(map_damage_usfs, data_ground %>% dplyr::select(X.UTM, Y.UTM)),
                         damage_lidar=as.logical(raster::extract(map_damage_lidar, data_ground %>% dplyr::select(X.UTM, Y.UTM)))) %>%
                rename(site_code=Site_Code) %>%
                mutate(damage_usfs=!is.na(damage_usfs))

refl_all = refl %>%
#  left_join(data_sex, by='site_code') %>%
  left_join(data_damage, by='site_code')

# filter to aspen
refl_aspen_unshaded = refl_all %>%
  filter(shade==1 & cover==0)

data_combined_flat = refl_aspen_unshaded %>% 
#  dplyr::select(starts_with("W",ignore.case=FALSE),ploidy_level,geneticSexID, damage_usfs, damage_lidar, site_code) %>%
  dplyr::select(starts_with("W",ignore.case=FALSE),ploidy_level, damage_usfs, damage_lidar, site_code) %>%
#  gather(key, value,-ploidy_level, -geneticSexID, -damage_usfs, -damage_lidar, -site_code) %>%
  gather(key, value,-ploidy_level, -damage_usfs, -damage_lidar, -site_code) %>%
  mutate(wavelength=as.numeric(as.character(gsub("W","",key)))) %>%
  mutate(reflectance=as.numeric(value)/10000) %>%
  rename(Cytotype=ploidy_level) %>%
  dplyr::select(-key,-value)



# make summaries

data_summary_cytotype = data_combined_flat %>% 
  filter(!is.na(Cytotype)) %>% 
  group_by(wavelength, Cytotype) %>% 
  summarize(m=mean(reflectance),s=sd(reflectance)) %>%
  arrange(Cytotype,wavelength) %>%
  filter(!(Cytotype %in% c("nan",'')))

#data_summary_geneticSexID = data_combined_flat %>% 
#  filter(!is.na(geneticSexID)) %>% 
#  group_by(wavelength, geneticSexID) %>% 
#  summarize(m=mean(reflectance),s=sd(reflectance)) %>%
#  arrange(geneticSexID,wavelength) %>%
#  filter(!(geneticSexID %in% c("nan",'')))

data_summary_damage_usfs = data_combined_flat %>% 
  filter(!is.na(damage_usfs)) %>% 
  group_by(wavelength, damage_usfs) %>% 
  summarize(m=mean(reflectance),s=sd(reflectance)) %>%
  arrange(damage_usfs,wavelength) %>%
  filter(!(damage_usfs %in% c("nan",'')))

data_summary_damage_lidar = data_combined_flat %>% 
  filter(!is.na(damage_lidar)) %>% 
  group_by(wavelength, damage_lidar) %>% 
  summarize(m=mean(reflectance),s=sd(reflectance)) %>%
  arrange(damage_lidar,wavelength)

g_spectra_cytotype = ggplot(data_summary_cytotype, aes(x=wavelength,ymin=m-s,ymax=m+s,y=m)) + 
  geom_ribbon(alpha=0.1,aes(fill=Cytotype)) +
  geom_line(aes(col=Cytotype)) +
  scale_color_manual(values=c("blue","red")) +
  scale_fill_manual(values=c("blue","red")) +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance")

#g_spectra_geneticSexID = ggplot(data_summary_geneticSexID, aes(x=wavelength,ymin=m-s,ymax=m+s,y=m)) + 
#  geom_ribbon(alpha=0.1,aes(fill=geneticSexID)) +
#  geom_line(aes(col=geneticSexID)) +
#  scale_color_manual(values=c("blue","red")) +
#  scale_fill_manual(values=c("blue","red")) +
#  theme_bw() +
#  xlab("Wavelength (nm)") +
#  ylab("Reflectance")

g_spectra_damage_usfs = ggplot(data_summary_damage_usfs, aes(x=wavelength,ymin=m-s,ymax=m+s,y=m)) + 
  geom_ribbon(alpha=0.1,aes(fill=damage_usfs)) +
  geom_line(aes(col=damage_usfs)) +
  scale_color_manual(values=c("blue","red"),name='USFS damage') +
  scale_fill_manual(values=c("blue","red"),name='USFS damage') +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance")

g_spectra_damage_lidar = ggplot(data_summary_damage_lidar, aes(x=wavelength,ymin=m-s,ymax=m+s,y=m)) + 
  geom_ribbon(alpha=0.1,aes(fill=damage_lidar)) +
  geom_line(aes(col=damage_lidar)) +
  scale_color_manual(values=c("blue","red"),name='Lidar damage') +
  scale_fill_manual(values=c("blue","red"),name='Lidar damage') +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance")


# write summary
write.csv(data_summary_cytotype,file='data_summary_cytotype.csv',row.names=F)
#write.csv(data_combined_flat %>% dplyr::select(-geneticSexID),file='data_combined_flat.csv',row.names=F)
write.csv(data_combined_flat,file='data_combined_flat.csv',row.names=F)











# confusion matrices


cm_cytotype_train = read.csv('confusion matrix/confusion_cytotype_train.csv')
row.names(cm_cytotype_train) = cm_cytotype_train[,1]
#cm_cytotype_train = cm_cytotype_train[,-1]
cm_cytotype_train_melted = melt(cm_cytotype_train,id.var="X") %>%
  rename(predicted=X,observed=variable) %>%
  mutate(predicted=gsub("predicted ","",predicted)) %>%
  mutate(observed=gsub("observed\\.","",observed))

g_cm_cytotype_train = ggplot(cm_cytotype_train_melted, aes(x=predicted,y=observed,fill=value,label=value)) +
  geom_tile() +
  geom_text(color='white') +
  scale_fill_gradient(low='midnightblue',high='limegreen',name='Count') +
  theme_minimal() + 
  coord_equal()

cm_cytotype_test = read.csv('confusion matrix/confusion_cytotype_test.csv')
row.names(cm_cytotype_test) = cm_cytotype_test[,1]
#cm_cytotype_test = cm_cytotype_test[,-1]
cm_cytotype_test_melted = melt(cm_cytotype_test,id.var="X") %>%
  rename(predicted=X,observed=variable) %>%
  mutate(predicted=gsub("predicted ","",predicted)) %>%
  mutate(observed=gsub("observed\\.","",observed))

g_cm_cytotype_test = ggplot(cm_cytotype_test_melted, aes(x=predicted,y=observed,fill=value,label=value)) +
  geom_tile() +
  geom_text(color='white') +
  scale_fill_gradient(low='midnightblue',high='limegreen',name='Count') +
  theme_minimal() + 
  coord_equal()



tp_test = cm_cytotype_test["predicted triploid","observed.triploid"]
tn_test = cm_cytotype_test["predicted diploid","observed.diploid"]
fp_test = cm_cytotype_test["predicted triploid","observed.diploid"]
fn_test = cm_cytotype_test["predicted diploid","observed.triploid"]

sens_test = tp_test/(tp_test+fn_test)
spec_test = tn_test/(tn_test+fp_test)






cm_cover_train = read.csv('confusion matrix/extract_1120_v2_nl_6_dr_0.4_nn_200_class_breakout_train.csv',header=FALSE)
dn = c("aspen","built environment","conifer","dry meadow","mesic meadow","misc bare","snow","water","woody riparian")
names(cm_cover_train) = dn
cm_cover_train = cbind(observed=dn,cm_cover_train)
cm_cover_train_melted = melt(cm_cover_train, id.var='observed') %>%
  rename(predicted=variable)

g_cm_cover_train = ggplot(cm_cover_train_melted, aes(x=predicted,y=observed,fill=value,label=value)) +
  geom_tile() +
  geom_text(color='white') +
  scale_fill_gradient(low='midnightblue',high='limegreen',name='Count') +
  theme_minimal() + 
  coord_equal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cm_cover_test = read.csv('confusion matrix/extract_1120_v2_nl_6_dr_0.4_nn_200_class_breakout_test.csv',header=FALSE)
dn = c("aspen","built environment","conifer","dry meadow","mesic meadow","misc bare","snow","water","woody riparian")
names(cm_cover_test) = dn
cm_cover_test = cbind(observed=dn,cm_cover_test)
cm_cover_test_melted = melt(cm_cover_test, id.var='observed') %>%
  rename(predicted=variable)

g_cm_cover_test = ggplot(cm_cover_test_melted, aes(x=predicted,y=observed,fill=value,label=value)) +
  geom_tile() +
  geom_text(color='white') +
  scale_fill_gradient(low='midnightblue',high='limegreen',name='Count') +
  theme_minimal() + 
  coord_equal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




tp_test = cm_cover_test_melted %>% filter(predicted=="aspen" & observed=="aspen") %>% dplyr::select(value) %>% summarize(sum(value)) %>% as.numeric
tn_test = cm_cover_test_melted %>% filter(predicted!="aspen" & observed!="aspen") %>% dplyr::select(value) %>% summarize(sum(value)) %>% as.numeric
fp_test = cm_cover_test_melted %>% filter(predicted=="aspen" & observed!="aspen") %>% dplyr::select(value) %>% summarize(sum(value)) %>% as.numeric
fn_test = cm_cover_test_melted %>% filter(predicted!="aspen" & observed=="aspen") %>% dplyr::select(value) %>% summarize(sum(value)) %>% as.numeric

sens_test = tp_test/(tp_test+fn_test)
spec_test = tn_test/(tn_test+fp_test)







g_full_cm_spectra = ggarrange(g_spectra_cytotype, 
          ggarrange(g_cm_cover_train, g_cm_cover_test, g_cm_cytotype_train, g_cm_cytotype_test,nrow=2,ncol=2,labels=c('(b)','(c)','(d)','(e)')),
          nrow=2,ncol=1,heights=c(1,2.5), labels=c('(a)',''))
ggsave(g_full_cm_spectra, file='g_full_cm_spectra.pdf',width=10,height=12)
ggsave(g_full_cm_spectra, file='g_full_cm_spectra.png',width=10,height=12)


g_cm_cytotype = ggarrange(g_cm_cytotype_train, g_cm_cytotype_test, align='hv',labels=c("(a)","(b)"))
g_cm_cover = ggarrange(g_cm_cover_train, g_cm_cover_test, align='hv',labels=c("(a)","(b)"))
ggsave(g_cm_cover, file='g_cm_cover.pdf',width=12,height=5)

ggsave(g_spectra_cytotype, file='g_spectra_cytotype.png',width=7,height=4)
#ggsave(g_spectra_geneticSexID, file='g_spectra_sex.png',width=7,height=4)

g_spectra_damage_all = ggarrange(g_spectra_damage_usfs, g_spectra_damage_lidar,nrow=2,ncol=1,labels=c("(a)","(b)"))
ggsave(g_spectra_damage_all, file='g_spectra_damage_all.png',width=7,height=6)


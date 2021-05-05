library(dplyr)
library(ranger)
library(randomForest)
library(pdp)
library(ggplot2)
library(caret)
library(viridis)
library(ggpubr)

load('df_all_10m.Rdata')

xvars_list = c('cyto','topo','trait','height')

xvars_all = unlist(sapply(1:length(xvars_list), function(i) { apply(combn(xvars_list,m=i),2,paste,collapse="+") }))

#q = 0.05

df_all_gridded = df_all_10m %>% 
  # do the spatial gridding
  mutate(x_grid = cut(x, breaks=seq(min(x),max(x),by=100))) %>%
  mutate(y_grid = cut(y, breaks=seq(min(y),max(y),by=100))) %>%
  mutate(grid_id = as.numeric(factor(paste(x_grid, y_grid)))) %>%
  # filter to above a certain threshold
  filter(Aspen.Cover.fraction > 0.5) %>%
  # make a clipped USFS Damage
  mutate(Damage.USFS.fraction.clipped = ifelse(is.na(Damage.Lidar.fraction), NA, Damage.USFS.fraction))


  #filter((Elevation > quantile(Elevation,q/2, na.rm=T)) & (Elevation < quantile(Elevation,1-q/2, na.rm=T)) &
  #  (Slope > quantile(Slope,q/2, na.rm=T)) & (Slope < quantile(Slope,1-q/2, na.rm=T)) &
  #  (Cos.aspect > quantile(Cos.aspect,q/2, na.rm=T)) & (Cos.aspect < quantile(Cos.aspect,1-q/2, na.rm=T)) &
  #  (Elevation > quantile(Elevation,q/2, na.rm=T)) & (Elevation < quantile(Elevation,1-q/2, na.rm=T)))


split_train_test <- function(df, fraction=0.8, xvar, yvar)
{
  n_cells = max(df$grid_id)
  train_ids = sample(1:n_cells, fraction*n_cells, replace=FALSE)
  test_ids = setdiff(1:n_cells, train_ids)
  
  train=df %>% filter(grid_id %in% train_ids) %>% dplyr::select(any_of(c(xvar, yvar))) %>% na.omit
  test=df %>% filter(grid_id %in% test_ids) %>% dplyr::select(any_of(c(xvar, yvar))) %>% na.omit
  
  train_x = as.data.frame(train[,xvar])
  names(train_x) = xvar
  train_y = as.numeric(train[,yvar,drop=TRUE])
  test_x = as.data.frame(test[,xvar])
  names(test_x) = xvar
  test_y = as.numeric(test[,yvar,drop=TRUE])
  
  return(list(train_x=train_x, train_y=train_y, test_x=test_x, test_y=test_y))
}

train_test_model <- function(tt, xvar, categorical)
{
  print('training model')
  if (categorical==TRUE)
  {
    m_rf = ranger(x=tt$train_x, y=tt$train_y,
                  max.depth = 10,
                  num.trees = 1000,
                  classification=TRUE,
                  verbose=TRUE,
                  probability = TRUE,
                  write.forest=TRUE,
                  importance = "impurity")
  }
  else
  {
    m_rf = ranger(x=tt$train_x, y=tt$train_y,
                  max.depth = 10,
                  num.trees = 1000,
                  classification=FALSE,
                  verbose=TRUE,
                  probability = FALSE,
                  write.forest=TRUE,
                  importance = "impurity")    
  }
  
  print('testing model')
  labels_predicted = predict(m_rf, tt$test_x)$predictions
  
  print('summarizing model')
  if (categorical==TRUE)
  {
    # apply threshold to probabilistic model
    labels_predicted = as.numeric(factor(labels_predicted[,1] < 0.5))
    labels_true = as.numeric(factor(tt$test_y))
    
    cm = confusionMatrix(factor(labels_predicted), factor(labels_true))$byClass
    result_this = data.frame(t(cm), xvars=paste(xvar,collapse="*"), np_train = nrow(tt$train_x),nv=ncol(tt$train_x))
  }
  else
  {
    labels_true = tt$test_y
    
    stats = postResample(labels_true, labels_predicted)
    
    result_this = data.frame(rmse=stats["RMSE"], r2=stats["Rsquared"], xvars=paste(xvar,collapse="*"), np_train = nrow(tt$train_x),nv=ncol(tt$train_x))
  }
  
  return(list(result=result_this,model=m_rf))
}


try_models <- function(yvar, categorical=TRUE, iter=1:5, np=NULL, xvartype)
{
  if (is.null(np))
  {
    if (categorical==TRUE)
    {
      np = df_all_gridded %>%
        group_by_(yvar) %>% 
        summarize(count=n()) %>%
        pull %>%
        min
      
    }
    else
    {
      np = nrow(df_all_gridded)
    }
    print(np)
  }
  
  params = expand.grid(xvartype=xvartype,
                       iter=iter, 
                       yvar=yvar, 
                       np=np, 
                       stringsAsFactors = FALSE)
  
  result = NULL
  models = vector(mode="list",length=nrow(params))
  datasets = vector(mode="list",length=nrow(params))
  for (i in 1:nrow(params))
  {
    start_time <- Sys.time()
    
    print('sampling data')
    if (categorical==TRUE)
    {
      # balance response variable
      df_all_ss = df_all_gridded %>%
        group_by_(params$yvar[i]) %>%
        sample_n(floor(params$np[i]/2)) # only works for binary variables - won't work for more classes
    }
    else
    {
      df_all_ss = df_all_gridded %>%
        sample_n(params$np[i])
    }
    
    print('splitting data')
    
    xvar = NULL
    if (length(grep("cyto",params$xvartype[i]))==1)
    {
      xvar=c(xvar, 'Cytotype.fractionDiploid')
    }
    if (length(grep("topo",params$xvartype[i]))==1)
    {
      xvar=c(xvar, 'Elevation','Cos.aspect','Slope')
    }
    if (length(grep("trait",params$xvartype[i]))==1)
    {
      xvar=c(xvar, 'Trait.N','Trait.d13C','Trait.CWC')
    }
    if (length(grep("height",params$xvartype[i]))==1)
    {
      xvar=c(xvar, 'Height.Canopy')
    }
    print(xvar)
    
    tt = split_train_test(df_all_ss, xvar=xvar, yvar=params$yvar[i])
  
    output = train_test_model(tt, xvar, categorical=categorical)
    
    # copy over outputs
    result_this = output$result
    models[[i]] = output$model
    datasets[[i]] = tt
    
    # report last stats
    end_time <- Sys.time()
    result_this$runtime = end_time - start_time
    result_this$np = params$np[i]
    result_this$xvartype = params$xvartype[i]
    result_this$yvar = params$yvar[i]
    
    result = rbind(result, result_this)
    
    print(i/nrow(params))
  }
  
  return(list(result=result,models=models,datasets=datasets))
}






# run models
#np_run = c(1e3, 1e4, 1e5)
models_results_cytotype = try_models(yvar=c("Cytotype.fractionDiploid"),xvartype=xvars_all[ !grepl("cyto",xvars_all) & !grepl("trait",xvars_all) ],np=NULL,iter = 1:5, categorical=FALSE)
save(models_results_cytotype, file='models_results_cytotype.Rdata')

models_results_traits = try_models(yvar=c("Trait.d13C","Trait.N","Trait.CWC"),xvartype =xvars_all[!grepl("trait",xvars_all) ],np=NULL,iter = 1:5, categorical=FALSE)
save(models_results_traits, file='models_results_traits.Rdata')

# need separate lines because resampling to balance response has diff # of observations for each damage variable
models_results_damage_usfs = try_models(yvar=c("Damage.USFS.fraction"),xvartype =xvars_all, np=NULL, iter = 1:5, categorical=FALSE)
save(models_results_damage_usfs, file='models_results_damage_usfs.Rdata')

models_results_damage_lidar = try_models(yvar=c("Damage.Lidar.fraction"),xvartype =xvars_all, np=NULL, iter = 1:5, categorical=FALSE)
save(models_results_damage_lidar, file='models_results_damage_lidar.Rdata')

models_results_damage_usfs_clipped = try_models(yvar=c("Damage.USFS.fraction.clipped"),xvartype =xvars_all, np=NULL, iter = 1:5, categorical=FALSE)
save(models_results_damage_usfs_clipped, file='models_results_damage_usfs_clipped.Rdata')




quantile_criteria <- function(x) {
  x<quantile(x,0.99,na.rm=TRUE) & x > quantile(x,0.01,na.rm=TRUE)
  }



df_train_ss = df_all_gridded %>% 
  dplyr::select(-Damage.USFS.fraction,-Damage.Lidar.fraction) %>% 
  na.omit
df_train_ss = df_train_ss[quantile_criteria(df_train_ss$Trait.CWC) & 
                            quantile_criteria(df_train_ss$Trait.d13C) &
                            quantile_criteria(df_train_ss$Trait.N) &
                            quantile_criteria(df_train_ss$Elevation) &
                            quantile_criteria(df_train_ss$Cos.aspect) &
                            quantile_criteria(df_train_ss$Slope) &
                            quantile_criteria(df_train_ss$Trait.N) &
                            quantile_criteria(df_train_ss$Height.Canopy)
                            ,]

## GET RESULTS FOR DAMAGE
do_pdp <- function(model_list, ids, pred.vars, grid.resolution, categorical)
{
  result = lapply(ids, function(id) {
    print(ids)
    
    result_this = partial(model_list$models[[id]], pred.var=pred.vars, 
            grid.resolution=grid.resolution,
            train=df_train_ss,#model_list$datasets[[id]]$train_x, # use the same training data throughout
            prob=categorical,
            progress='text')
    result_this$replicate=id
    return(result_this)
    
    })  
  
  result_all = do.call("rbind",result)
  return(result_all)
}


summarize_pdp <- function(pdp)
{
  # need to do approximate groupings as the partial picks different predictor values based on samples...
  pdp_list = pdp %>% 
    as.tbl %>% 
    group_by(replicate) %>% 
    filter(Cytotype.fractionDiploid %in% c(0,1)) %>%
    #mutate(Cytotype.fractionDiploid=as.numeric(Cytotype.fractionDiploid>0.5)) %>% 
    group_split
  
  # stack up all the predictions by replicate into an array
  pdp_list = lapply(pdp_list, as.matrix)
  pdp_array = array(unlist(pdp_list), dim = c(nrow(pdp_list[[1]]), ncol(pdp_list[[1]]), length(pdp_list)))
  
  # calculate elementwise means
  means = apply(pdp_array, c(1,2), mean) %>% 
    as.data.frame
  names(means) = names(pdp)
  
  sds = apply(pdp_array, c(1,2), sd) %>% 
    as.data.frame
  names(sds) = names(pdp)
  
  sds = sds %>% 
    dplyr::select(-replicate)# %>%
  # put cytotype back
  #mutate(Cytotype.fractionDiploid=factor(Cytotype.fractionDiploid,labels=levels(pdp$Cytotype.fractionDiploid)))
  
  # lose the replicate column
  final = means %>% 
    dplyr::select(-replicate) %>%
    mutate(yhat.mean = yhat) %>%
    mutate(yhat.sd = sds$yhat) %>%
    mutate(Cytotype.fractionDiploid = factor(Cytotype.fractionDiploid,levels=c(0,1),labels=c("triploid","diploid"),ordered=TRUE))# %>%
    # put cytotype back
    #mutate(Cytotype.fractionDiploid=factor(Cytotype.fractionDiploid,labels=levels(pdp$Cytotype.fractionDiploid)))

  return(final)
}






## DAMAGE MODELS
best_ids_damage_usfs = models_results_damage_usfs$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_damage_usfs = do_pdp(models_results_damage_usfs, 
                          ids=best_ids_damage_usfs,
                          pred.vars=c('Cytotype.fractionDiploid','Elevation','Trait.CWC'),
                          categorical=FALSE,
                          grid.resolution = 5)
save(pdp_damage_usfs,file='pdp_damage_usfs.Rdata')

# DAMAGE lidar
best_ids_damage_lidar = models_results_damage_lidar$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_damage_lidar = do_pdp(models_results_damage_lidar, 
                                 ids=best_ids_damage_lidar,
                                 pred.vars=c('Cytotype.fractionDiploid','Elevation','Trait.CWC'),
                                 categorical=FALSE,
                                 grid.resolution = 5)
save(pdp_damage_lidar,file='pdp_damage_lidar.Rdata')

# DAMAGE usfs.clipped
best_ids_damage_usfs_clipped = models_results_damage_usfs_clipped$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_damage_usfs_clipped = do_pdp(models_results_damage_usfs_clipped, 
                                 ids=best_ids_damage_usfs_clipped,
                                 pred.vars=c('Cytotype.fractionDiploid','Elevation','Trait.CWC'),
                                 categorical=FALSE,
                                 grid.resolution = 5)
save(pdp_damage_usfs_clipped,file='pdp_damage_usfs_clipped.Rdata')










### CYTOTYPE
best_ids_cytotype = models_results_cytotype$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_cytotype = do_pdp(models_results_cytotype, 
                      ids=best_ids_cytotype,
                      pred.vars=c('Elevation','Cos.aspect'),
                      categorical=FALSE,
                      grid.resolution = 5)
save(pdp_cytotype,file='pdp_cytotype.Rdata')

# TRAIT CWC
best_ids_trait_cwc = models_results_traits$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(yvar=="Trait.CWC" & np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_trait_cwc = do_pdp(models_results_traits, 
                       ids=best_ids_trait_cwc,
                       pred.vars=c('Cytotype.fractionDiploid','Cos.aspect','Height.Canopy'),
                       categorical=FALSE,
                       grid.resolution = 5)
save(pdp_trait_cwc,file='pdp_trait_cwc.Rdata')

# TRAIT d13c
best_ids_trait_d13c = models_results_traits$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(yvar=="Trait.d13C" & np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_trait_d13c = do_pdp(models_results_traits, 
                        ids=best_ids_trait_d13c,
                        pred.vars=c('Cytotype.fractionDiploid','Cos.aspect','Height.Canopy'),
                        categorical=FALSE,
                        grid.resolution = 5)
save(pdp_trait_d13c,file='pdp_trait_d13c.Rdata')

# TRAIT n
best_ids_trait_n = models_results_traits$result %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(yvar=="Trait.N" & np == max(np) & xvartype==names(which.max(sapply(xvartype,nchar))) ) %>%
  pull(rowid)

pdp_trait_n = do_pdp(models_results_traits, 
                     ids=best_ids_trait_n,
                     pred.vars=c('Cytotype.fractionDiploid','Cos.aspect','Height.Canopy'),
                     categorical=FALSE,
                     grid.resolution = 5)
save(pdp_trait_n,file='pdp_trait_n.Rdata')


## or alternatively load precomputed results
#load('pdp_damage_usfs.Rdata')
#load('pdp_damage_lidar.Rdata')
#load('pdp_damage_usfs_clipped.Rdata')
#load('pdp_cytotype.Rdata')
#load('pdp_trait_cwc.Rdata')
#load('pdp_trait_d13c.Rdata')
#load('pdp_trait_n.Rdata')



######## PLOT PDPS ########

# DAMAGE USFS
pdp_summary_damage_usfs = summarize_pdp(pdp_damage_usfs)

g_pdp_damage_usfs.mean = ggplot(pdp_summary_damage_usfs, aes(x=Elevation,y=Trait.CWC,fill=yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle('Damage fraction - USFS') +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

g_pdp_damage_usfs.cv = ggplot(pdp_summary_damage_usfs, aes(x=Elevation,y=Trait.CWC,fill=yhat.sd/yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

g_pdp_damage_usfs = ggarrange(g_pdp_damage_usfs.mean, g_pdp_damage_usfs.cv,
                              nrow=1,ncol=2,labels=c('(a)','(b)'), align='hv')
#ggsave(g_pdp_damage_usfs,file='g_pdp_damage_usfs.png',width=7,height=7)




# DAMAGE LIDAR
pdp_summary_damage_lidar = summarize_pdp(pdp_damage_lidar)

g_pdp_damage_lidar.mean = ggplot(pdp_summary_damage_lidar, aes(x=Elevation,y=Trait.CWC,fill=yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle('Damage fraction - LiDAR') +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

g_pdp_damage_lidar.cv = ggplot(pdp_summary_damage_lidar, aes(x=Elevation,y=Trait.CWC,fill=yhat.sd/yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

g_pdp_damage_lidar = ggarrange(g_pdp_damage_lidar.mean, g_pdp_damage_lidar.cv,
                               nrow=1,ncol=2,labels=c('(c)','(d)'), align='hv')
#ggsave(g_pdp_damage_lidar,file='g_pdp_damage_lidar.png',width=7,height=7)
ggsave(ggarrange(g_pdp_damage_usfs, g_pdp_damage_lidar, 
                 nrow=2,ncol=1), file='g_pdp_damage.png',width=8,height=5)

# DAMAGE USFS clipped
pdp_summary_damage_usfs_clipped = summarize_pdp(pdp_damage_usfs_clipped)

g_pdp_damage_usfs_clipped.mean = ggplot(pdp_summary_damage_usfs_clipped, aes(x=Elevation,y=Trait.CWC,fill=yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle('Damage fraction - USFS (clipped)') +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

g_pdp_damage_usfs_clipped.cv = ggplot(pdp_summary_damage_usfs_clipped, aes(x=Elevation,y=Trait.CWC,fill=yhat.sd/yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

g_pdp_damage_usfs_clipped = ggarrange(g_pdp_damage_usfs_clipped.mean, g_pdp_damage_usfs_clipped.cv,
          nrow=2,ncol=1,labels=c('(a)','(b)'), align='hv')
ggsave(g_pdp_damage_usfs_clipped,file='g_pdp_damage_usfs_clipped.png',width=7,height=7)



# CYTOTYPE
pdp_summary_cytotype = summarize_pdp(pdp_cytotype %>% as_tibble %>% mutate(Cytotype.fractionDiploid=0)) # make a dummy column for cytotype just for plotting

g_pdp_cytotype.mean = ggplot(pdp_summary_cytotype, aes(x=Elevation,y=Cos.aspect,fill=yhat)) + 
  geom_tile() + 
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle('Prob(Diploid cytotype)') +
  xlab('Elevation (m)') + ylab(expression(paste("Cosine aspect")))

g_pdp_cytotype.cv = ggplot(pdp_summary_cytotype, aes(x=Elevation,y=Cos.aspect,fill=yhat.sd/yhat)) + 
  geom_tile() + 
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Elevation (m)') + ylab(expression(paste("Cosine aspect")))

g_pdp_cytotype = ggarrange(g_pdp_cytotype.mean, g_pdp_cytotype.cv,
          nrow=2,ncol=1,labels=c('(a)','(b)'), align='hv')
ggsave(g_pdp_cytotype,file='g_pdp_cytotype.png',width=4.5,height=7)





# TRAIT CWC
pdp_summary_trait_cwc = summarize_pdp(pdp_trait_cwc)

g_pdp_trait_cwc.mean = ggplot(pdp_summary_trait_cwc, aes(x=Cos.aspect,y=Height.Canopy,fill=yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle(expression(paste('Canopy water content (CWC, mL m'^{-2},')'))) +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy height (m)")))

g_pdp_trait_cwc.cv = ggplot(pdp_summary_trait_cwc, aes(x=Cos.aspect,y=Height.Canopy,fill=yhat.sd/yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy height (m)")))

g_pdp_trait_cwc = ggarrange(g_pdp_trait_cwc.mean, g_pdp_trait_cwc.cv,
          nrow=1,ncol=2,labels=c('(a)','(b)'),align='hv')


# TRAIT d13c
pdp_summary_trait_d13c = summarize_pdp(pdp_trait_d13c)

g_pdp_trait_d13c.mean = ggplot(pdp_summary_trait_d13c, aes(x=Cos.aspect,y=Height.Canopy,fill=yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle(expression(paste('Carbon isotope shift (',delta^13,'C, per mil)'))) +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy height (m)")))

g_pdp_trait_d13c.cv = ggplot(pdp_summary_trait_d13c, aes(x=Cos.aspect,y=Height.Canopy,fill=abs(yhat.sd/yhat))) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy height (m)")))

g_pdp_trait_d13c = ggarrange(g_pdp_trait_d13c.mean, g_pdp_trait_d13c.cv,
          nrow=1,ncol=2,labels=c('(c)','(d)'), align='hv')



# TRAIT N
pdp_summary_trait_n = summarize_pdp(pdp_trait_n)

g_pdp_trait_n.mean = ggplot(pdp_summary_trait_n, aes(x=Cos.aspect,y=Height.Canopy,fill=yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_viridis(name=expression(mu),option='magma') +
  theme_bw() +
  ggtitle('Nitrogen content (N, %)') +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy height (m)")))

g_pdp_trait_n.cv = ggplot(pdp_summary_trait_n, aes(x=Cos.aspect,y=Height.Canopy,fill=yhat.sd/yhat)) + 
  geom_tile() + 
  facet_wrap(~Cytotype.fractionDiploid) +
  scale_fill_gradient(name=expression(paste(sigma,"/",mu)),low='lightgreen',high='magenta', limits=c(0,0.2)) +
  theme_bw() +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy height (m)")))

g_pdp_trait_n = ggarrange(g_pdp_trait_n.mean, g_pdp_trait_n.cv,
          nrow=1,ncol=2,labels=c('(e)','(f)'), align='hv')

# TRAITS ALL
ggsave(ggarrange(g_pdp_trait_cwc, g_pdp_trait_d13c, g_pdp_trait_n,nrow=3,ncol=1),
       file='g_pdp_trait.png',width=9,height=7)





















## PLOT PERFORMANCE OF MODELS


assign_levels_xvars <- function(xv)
{
  xv = factor(xv, levels=xvars_all, ordered=TRUE)
  return(xv)
}




g_cytotype_performance = ggplot(models_results_cytotype$result %>% mutate(xvartype=assign_levels_xvars(xvartype)),
                                aes(x=xvartype,y=r2)) + 
  facet_wrap(~yvar) +
  geom_violin(draw_quantiles=0.5) +
  theme_bw() +
  xlab('Predictor variables') +
  ylim(0,0.5) + ylab(expression(paste("R"^2))) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))

ggsave(g_cytotype_performance, file='g_cytotype_performance.png',width=7,height=4)


g_trait_performance = ggplot(models_results_traits$result %>% mutate(xvartype=assign_levels_xvars(xvartype)),
                             aes(x=xvartype,y=r2)) + 
  facet_wrap(~yvar) +
  geom_violin(draw_quantiles=0.5) +
  theme_bw() +
  xlab('Predictor variables') + ylab(expression(paste("R"^2))) + 
  ylim(0,0.5) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))

ggsave(g_trait_performance, file='g_trait_performance.png',width=7,height=4)


g_damage_usfs_performance = ggplot(models_results_damage_usfs$result %>% mutate(xvartype=assign_levels_xvars(xvartype)),
                                                    aes(x=xvartype,y=r2)) + 
  facet_wrap(~yvar) +
  geom_violin(draw_quantiles=0.5) +
  theme_bw() +
  xlab('Predictor variables') +
  ylim(0,0.5) + ylab(expression(paste("R"^2))) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))

#ggsave(g_damage_usfs_performance, file='g_damage_usfs_performance.png',width=7,height=4)

g_damage_lidar_performance = ggplot(models_results_damage_lidar$result %>% mutate(xvartype=assign_levels_xvars(xvartype)),
                                   aes(x=xvartype,y=r2)) + 
  facet_wrap(~yvar) +
  geom_violin(draw_quantiles=0.5) +
  theme_bw() +
  xlab('Predictor variables') +
  ylim(0,0.5) + ylab(expression(paste("R"^2))) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))

#ggsave(g_damage_lidar_performance, file='g_damage_lidar_performance.png',width=7,height=4)


g_damage_usfs_clipped_performance = ggplot(models_results_damage_usfs_clipped$result %>% mutate(xvartype=assign_levels_xvars(xvartype)),
                                    aes(x=xvartype,y=r2)) + 
  facet_wrap(~yvar) +
  geom_violin(draw_quantiles=0.5) +
  theme_bw() +
  xlab('Predictor variables') +
  ylim(0,0.5) + ylab(expression(paste("R"^2))) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))


ggsave(ggarrange(g_damage_usfs_performance, g_damage_lidar_performance, g_damage_usfs_clipped_performance, 
                 nrow=3,ncol=1, 
                 labels=c('(a)','(b)','(c)'),
                 common.legend = TRUE,
                 align='hv'),
  file='g_damage_performance.png',width=7,height=8)
#ggsave(g_damage_usfs_clipped_performance, file='g_damage_usfs_clipped_performance.png',width=7,height=4)




# summarized PDP stats
# cytotype
models_results_cytotype$result %>% 
  slice(best_ids_cytotype) %>% 
  summarize(r2.mean=mean(r2),r2.sd=sd(r2))

# trait cwc
models_results_traits$result %>% 
  slice(best_ids_trait_cwc) %>% 
  summarize(vars=xvartype[1],r2.mean=mean(r2),r2.sd=sd(r2))

# trait d13c
models_results_traits$result %>% 
  slice(best_ids_trait_d13c) %>% 
  summarize(vars=xvartype[1],r2.mean=mean(r2),r2.sd=sd(r2))

# trait n
models_results_traits$result %>% 
  slice(best_ids_trait_n) %>% 
  summarize(vars=xvartype[1],r2.mean=mean(r2),r2.sd=sd(r2))



# damage usfs
models_results_damage_usfs$result %>% 
  slice(best_ids_damage_usfs) %>% 
  summarize(vars=xvartype[1],r2.mean=mean(r2),r2.sd=sd(r2))

# damage lidar
models_results_damage_lidar$result %>% 
  slice(best_ids_damage_lidar) %>% 
  summarize(vars=xvartype[1],r2.mean=mean(r2),r2.sd=sd(r2))





# also assess variable prevalence for DAMAGE
cuts_elev = seq(min(df_train_ss$Elevation),max(df_train_ss$Elevation),length.out = 10)
cuts_cwc = seq(min(df_train_ss$Trait.CWC),max(df_train_ss$Trait.CWC),length.out = 10)
counts_for_damage = df_all_10m %>% 
  select(Cytotype.fractionDiploid, Elevation, Trait.CWC) %>%
  mutate(Cytotype = factor(Cytotype.fractionDiploid>0.5,levels=c(FALSE,TRUE),labels=c("triploid","diploid"))) %>%
  mutate(Elevation.cut = cut(Elevation,breaks=cuts_elev)) %>%
  mutate(Trait.CWC.cut = cut(Trait.CWC,breaks=cuts_cwc)) %>%
  group_by(Cytotype, Elevation.cut, Trait.CWC.cut) %>%
  summarize(count=n()) %>%
  mutate(Elevation.cut = cuts_elev[Elevation.cut]) %>%
  mutate(Trait.CWC.cut = cuts_cwc[Trait.CWC.cut])

g_counts_damage = ggplot(counts_for_damage, aes(x=Elevation.cut,y=Trait.CWC.cut,fill=sqrt(count / sum(count)))) + 
  geom_tile() + 
  facet_wrap(~Cytotype) +
  theme_bw() +
  scale_fill_viridis(option='C',name='√Prevalence') +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

ggsave(g_counts_damage, file='g_counts_damage.png', width=8, height=4)



# also assess variable prevalence for TRAITS
cuts_cos_aspect = seq(min(df_train_ss$Cos.aspect),max(df_train_ss$Cos.aspect),length.out = 10)
cuts_height_canopy = seq(min(df_train_ss$Height.Canopy),max(df_train_ss$Height.Canopy),length.out = 10)
counts_for_traits = df_all_10m %>% 
  select(Cytotype.fractionDiploid, Cos.aspect, Height.Canopy) %>%
  mutate(Cytotype = factor(Cytotype.fractionDiploid>0.5,levels=c(FALSE,TRUE),labels=c("triploid","diploid"))) %>%
  mutate(Cos.aspect.cut = cut(Cos.aspect,breaks=cuts_cos_aspect)) %>%
  mutate(Height.canopy.cut = cut(Height.Canopy,breaks=cuts_height_canopy)) %>%
  group_by(Cytotype, Cos.aspect.cut, Height.canopy.cut) %>%
  summarize(count=n()) %>%
  mutate(Cos.aspect.cut = cuts_cos_aspect[Cos.aspect.cut]) %>%
  mutate(Height.canopy.cut = cuts_height_canopy[Height.canopy.cut])

g_counts_traits = ggplot(counts_for_traits, aes(x=Cos.aspect.cut,y=Height.canopy.cut,fill=sqrt(count / sum(count)))) + 
  geom_tile() + 
  facet_wrap(~Cytotype) +
  theme_bw() +
  scale_fill_viridis(option='C',name='√Prevalence') +
  xlab('Cosine aspect') + ylab(expression(paste("Canopy heigth (m)")))

ggsave(g_counts_traits, file='g_counts_traits.png', width=8, height=4)



# also assess variable prevalence for CYTOTYPE
cuts_elev = seq(min(df_train_ss$Elevation),max(df_train_ss$Elevation),length.out = 10)
cuts_cwc = seq(min(df_train_ss$Trait.CWC),max(df_train_ss$Trait.CWC),length.out = 10)
counts_for_cytotype = df_all_10m %>% 
  select(Elevation, Trait.CWC) %>%
  mutate(Elevation.cut = cut(Elevation,breaks=cuts_elev)) %>%
  mutate(Trait.CWC.cut = cut(Trait.CWC,breaks=cuts_cwc)) %>%
  group_by(Elevation.cut, Trait.CWC.cut) %>%
  summarize(count=n()) %>%
  mutate(Elevation.cut = cuts_elev[Elevation.cut]) %>%
  mutate(Trait.CWC.cut = cuts_cwc[Trait.CWC.cut])

g_counts_cytotype = ggplot(counts_for_cytotype, aes(x=Elevation.cut,y=Trait.CWC.cut,fill=sqrt(count / sum(count)))) + 
  geom_tile() + 
  theme_bw() +
  scale_fill_viridis(option='C',name='√Prevalence') +
  xlab('Elevation (m)') + ylab(expression(paste("Canopy water content (mL m"^{-2},")")))

ggsave(g_counts_cytotype, file='g_counts_cytotype.png', width=5,height=4)

                                

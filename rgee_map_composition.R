## 1. Pacotes a serem utilizados ----

x <- c('rgdal','tmap','dplyr','raster','rgeos','ggplot2','gstat',
       'openair','sf','sp', 'hydroGOF','Metrics','geobr', 'nasapower',
       'googledrive','geojsonio','mapview','googleCloudStorageR', 'RStoolbox')

# Função para carregar todos os pacotes, ou instalar e carregar caso não não tenha algum pacote
lapply(x, function(i){
  if(i %in% rownames(installed.packages())){
    library(i,character.only = T)
  }else{
    install.packages(i)
    library(i,character.only = T)
  }
}
)

remotes::install_github("r-spatial/rgee", force = TRUE)
library(rgee)


### 2. Corpo ----


#ee_install()# ------------- INSTALAR TERMINAL GEE

#ee_check()

ee_Initialize(
  email = "gabriel.ferraz@agrotools.com.br",
  drive = TRUE,
  gcs = F
) 

crs = '+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs'

embargo_nugeo <- readOGR("Arquivos vetoriais/nugeo.shp")

embargo_ibama <- readOGR("Arquivos vetoriais/ibama_500.shp")

poligonos <- c(embargo_ibama,embargo_nugeo)

for (poly in poligonos) {
  poly <- spTransform(poly,crs(crs))
}


centroide <- as.data.frame(coordinates(gCentroid(embargo_nugeo, byid = T))) # definindo centróide

x = mean(centroide$x)
y = mean(centroide$y)

point <- ee$Geometry$Point( x, y) 

extensao_f <- fortify(as(extent(embargo_nugeo), 'SpatialPolygons'))
lpi_extensao <- list()
lpi_extensao <- lapply(1:nrow(extensao_f), function(i){
  a <- as.numeric(extensao_f$long[i])
  b <- as.numeric(extensao_f$lat[i])
  c <- append(a,b)
  return(c)
}
)

extensao <- ee$Geometry$Polygon(lpi_extensao)

datas <-  seq(from = as.Date('2013-04-12'), to = Sys.Date(), by = 'quarter')


dir.create(file.path(getwd(), '/ndvi_map'), showWarnings = FALSE)

dir.create(file.path(getwd(), '/rgb_map'), showWarnings = FALSE)

dir.create(file.path(getwd(), '/rgb_img'), showWarnings = FALSE)

range_datas <- seq_along(datas)[-length(seq_along(datas))]

################################################################################

for (dia in range_datas) {
  
  start = ee$Date(as.character(datas[dia]))
  
  finish = ee$Date(as.character(datas[dia+1]))

################################################################################    
  
  if (datas[dia-1] < '2012-05-05'){
    imagecol <- ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$  
      filterBounds(point)$ 
      filterDate(start, finish)$ 
      sort("CLOUD_COVER", F) # 
    
    image <- imagecol$first()
    
    rgb <- ee$Image(image)$
      select(c("B3", "B4", "B5"))
    
    getNDVI <- function(image) {
        image$normalizedDifference(c('B4', 'B3'))
    }
    
    ndvi1 <- getNDVI(image)
    
    clipped_rgb <- rgb$clip(extensao)
    
    clipped_ndvi <- ndvi1$clip(extensao)
    
    data_da_imagem <- ee_get_date_img(image, time_end = FALSE)
    data_e_hora <- data_da_imagem$time_start
    data_da_imagem_str <- gsub("-","",substr(as.character(data_da_imagem$time_start),1,10))
    
    
    rgb_raster <- ee_as_raster(clipped_rgb,region = ee$Geometry$Polygon(lpi_extensao),
                               via = 'drive', dsn = paste0("rgb_img/",data_da_imagem_str,".tif"))
    
    ndvi_raster <- ee_as_raster(clipped_ndvi,region = ee$Geometry$Polygon(lpi_extensao),
                                via = 'drive')
    new_rgb_raster <- rgb_raster
    
    for (i in  1:3) {
      new_rgb_raster[[i]] <- stretch(rgb_raster[[i]], maxv = 255, minv = 0, minq = 0.1, maxq = 0.9)  
    } 
    
    mapa_rgb <- tm_shape(new_rgb_raster)+tm_rgb(r=1, g=2, b=3, max.value = max(maxValue(new_rgb_raster)))+
      tm_shape(embargo_nugeo)+tm_borders(lwd = 2, col = 'red')+
      tm_shape(embargo_ibama)+tm_borders(lwd = 2, col = 'blue')+
      tm_polygons(alpha=0)+tm_scale_bar(position=c(0.08, 0.90),width = 0.1)+ # Barra de escala
      tm_compass(type="arrow", position=c(0.10,0.07), show.labels = 1,size=2.5,fontsize = 0.6)+
      tm_grid(labels.margin.y = -0.3,alpha=0.1,n.x=5,n.y=5)+ # Grid
      tm_layout(main.title = paste0("RGB ", data_e_hora), main.title.position = 'center')
    
    mapa_ndvi <- tm_shape(ndvi_raster)+tm_raster(title=paste0('NDVI ', data_e_hora), style = 'fixed', breaks = seq(-1, 1, by = 0.2) )+
      tm_layout(legend.outside = T)+
      tm_shape(embargo_nugeo)+tm_borders(lwd = 2, col = 'red')+
      tm_shape(embargo_ibama)+tm_borders(lwd = 2, col = 'blue')+
      tm_scale_bar(position=c(0.08, 0.90),width = 0.1)+ # Barra de escala
      tm_compass(type="arrow", position=c(0.10,0.07), show.labels = 1,size=2.5,fontsize = 0.6)#+ # Norte
    
    
    tmap_save(tm = mapa_rgb, file =  paste0('rgb_map/mapa',data_da_imagem_str,'.png'))
    tmap_save(tm = mapa_ndvi, file =  paste0('ndvi_map/mapa',data_da_imagem_str,'.png'))
    
################################################################################   
    
  } else if (datas[dia] > '2013-04-11') {
    imagecol <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")$ 
      filterBounds(point)$
      filterDate(start, finish)$ 
      sort("CLOUD_COVER", F) 
    
    image <- imagecol$first()
    
    rgb <- ee$Image(image)$
      select(c("B4", "B5", "B6"))
    
    getNDVI <- function(i) {
        image$normalizedDifference(c('B5', 'B4'))
      } 
    
    
    ndvi1 <- getNDVI(image)
    
    clipped_rgb <- rgb$clip(extensao)
    
    clipped_ndvi <- ndvi1$clip(extensao)
    
    data_da_imagem <- ee_get_date_img(image, time_end = FALSE)
    data_e_hora <- data_da_imagem$time_start
    data_da_imagem_str <- gsub("-","",substr(as.character(data_da_imagem$time_start),1,10))
    
    
    rgb_raster <- ee_as_raster(clipped_rgb,region = ee$Geometry$Polygon(lpi_extensao),
                               via = 'drive', dsn = paste0("rgb_img/",data_da_imagem_str,".tif"))
    ndvi_raster <- ee_as_raster(clipped_ndvi,region = ee$Geometry$Polygon(lpi_extensao),
                via = 'drive')
    new_rgb_raster <- rgb_raster
    
    for (i in  1:3) {
      new_rgb_raster[[i]] <- stretch(rgb_raster[[i]], maxv = 255, minv = 0, minq = 0.1, maxq = 0.9)  
    } 
    
    mapa_rgb <- tm_shape(new_rgb_raster)+tm_rgb(r=1, g=2, b=3, max.value = max(maxValue(new_rgb_raster)))+
      tm_shape(embargo_nugeo)+tm_borders(lwd = 2, col = 'red')+
      tm_shape(embargo_ibama)+tm_borders(lwd = 2, col = 'blue')+
      tm_polygons(alpha=0)+tm_scale_bar(position=c(0.08, 0.90),width = 0.1)+ # Barra de escala
      tm_compass(type="arrow", position=c(0.10,0.07), show.labels = 1,size=2.5,fontsize = 0.6)+
      tm_grid(labels.margin.y = -0.3,alpha=0.1,n.x=5,n.y=5)+ # Grid
      tm_layout(main.title = paste0("RGB ", data_e_hora), main.title.position = 'center')
    
    mapa_ndvi <- tm_shape(ndvi_raster)+tm_raster(title=paste0('NDVI ', data_e_hora), style = 'fixed', breaks = seq(-1, 1, by = 0.2) )+
      tm_layout(legend.outside = T)+
      tm_shape(embargo_nugeo)+tm_borders(lwd = 2, col = 'red')+
      tm_shape(embargo_ibama)+tm_borders(lwd = 2, col = 'blue')+
      tm_scale_bar(position=c(0.08, 0.90),width = 0.1)+ # Barra de escala
      tm_compass(type="arrow", position=c(0.10,0.07), show.labels = 1,size=2.5,fontsize = 0.6)#+ # Norte
  
    
    tmap_save(tm = mapa_rgb, file =  paste0('rgb_map/mapa',data_da_imagem_str,'.png'))
    tmap_save(tm = mapa_ndvi, file =  paste0('ndvi_map/mapa',data_da_imagem_str,'.png'))
 
 
###############################################################################    
        
  } else {
    next
  }
  
}
  






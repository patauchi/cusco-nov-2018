########################### CONSTRUIR MODELOS KU ############################
###                       A.T. Peterson and M. E. Cobos
##                                 Mayo 2018
# 1. Definir directorio de trabajo
# 2. Calibrar modelos
# 3. Evaluacion de parametros del modelo
# 4. Evaluacion del modelo final
##

# Datos de ocurrencia
#library(rgbif)
library(spThin)
#library(red)
library(dismo)
library(dplyr)
library(readr)
#library(rebird)
library(spocc)
library(kuenm)

####################################################################
#####################  PARTE 1. Definir el directorio de trabajo ####
####################################################################
scientificName = "Patagioenas subvinacea"
sp.name <- 'Patagioenas_subvinacea'


# carpeta principal de trabajo
main.getwd <- c("/Users/p.joseratauchi/Documents/Project-R/KU_models")
setwd(main.getwd)
# Directorio de trabajo para la especie 1
# Crear carpeta principal de la especie

dir.create(sp.name)
#unlink(paste0(getwd(),"/", sp.name), recursive = TRUE) # Elimina la carpeta

setwd(paste0(getwd(),'/',sp.name))

getwd() # Nuevo directorio (debe conincidir con la especie de estudio)
dir.create('G_variables')
dir.create('M_variables')

##### Para correr distintos set de variables
dir.create('M_variables/set1')
#dir.create('M_variables/set2')
#dir.create('M_variables/set3')

# Variables para la transferencia
dir.create('G_variables/set1')
#dir.create('G_variables/set2')
#dir.create('G_variables/set3')

# Agregar distintos GCMs al modelo de transferencia
# Futuro
dir.create('G_variables/set1/fut_he45')
dir.create('G_variables/set1/fut_he26')
dir.create('G_variables/set1/fut_he60')
dir.create('G_variables/set1/fut_he85')
dir.create('G_variables/set1/fut_mc45')
dir.create('G_variables/set1/fut_cc45')


####################################################################################
################################     PARTE 2. Datos de ocurrencia   ################
####################################################################################



#######  datos de ocurrencia
(data_occs <- occ(query = scientificName, 
                  from = c('gbif','bison','inat','ebird','ecoengine','vertnet','idigbio','obis','ala')))

data_occs <- occ2df(data_occs)

data_occs <- data_occs %>%
  na.omit() %>%
  select(long = longitude, lat = latitude)


####
data_occss <- data.frame(long=as.numeric(as.character(data_occs$long)),
                         lat=as.numeric(as.character(data_occs$lat)))
####
data_occ <- data_occss

#
data_occ <- data.frame(SPEC=sp.name, Longitude = data_occ$long, Latitude=data_occ$lat)

###### Rango altitudinal de los puntos
altitud <- raster('/Volumes/SEAGATE/clima/altitud/altitud.tif',
                  full.names=TRUE)
rango.d <- raster::extract(altitud, data_occ[,2:3], method='bilinear')
hist(rango.d)


data_occ <- data.frame(data_occ, rango.d)
rm(altitud)
data_occ <- data_occ %>%
  filter(rango.d > 2000 & rango.d < 1800)

hist(data_occ$rango.d)


# Ajuste de ocurrencias
data_occ <- spThin::thin.algorithm(data_occ[ ,2:3], thin.par = 25, reps = 1)
data_occ <- data_occ[[1]]


# Raster plot
env <- raster('/Volumes/SEAGATE/clima/pre/5km/bio_1.tif')

## crear extent
lco <- (SpatialPoints(data_occ))

coords = matrix(c(extent(lco)[1] - 0.5, extent(lco)[3] - 0.5,
                  extent(lco)[2] + 0.5, extent(lco)[3] - 0.5,
                  extent(lco)[2] + 0.5, extent(lco)[4] + 0.5,
                  extent(lco)[1] - 0.5, extent(lco)[4] + 0.5), 
                ncol = 2, byrow = TRUE)


pM = Polygon(coords)

pM = SpatialPolygons(list(Polygons(list(pM), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#pM = SpatialPolygonsDataFrame(list(Polygons(list(pM), ID = "a")))
#plot(pM)

env <- crop(env,pM)
env <- mask(env,pM)

### Create background points
bg.pts <- randomPoints(env, 500 )

plot(env)
#points(bg.pts)
points(data_occ, pch=21, bg=2)


### Checkerboard1 partitioning method
library(ENMeval)
chk1.pts <- get.randomkfold(data_occ, bg.pts, 4)
plot(env)
points(data_occ, pch=23, bg=chk1.pts$occ.grp)

length(chk1.pts$occ.grp[chk1.pts$occ.grp==1])
length(chk1.pts$occ.grp[chk1.pts$occ.grp==2])
length(chk1.pts$occ.grp[chk1.pts$occ.grp==3])
length(chk1.pts$occ.grp[chk1.pts$occ.grp==4])

data_occ <- data.frame(data_occ, fac=chk1.pts$occ.grp)


### 
# data_occ <- data_occ[-c(6,48),]
##

train <- data_occ %>%
  mutate(Specie = sp.name) %>%
  filter(fac != 1) %>%
  #filter(!fac %in% c(2,3)) %>%
  select(Specie, Longitude, Latitude) %>%
  data.frame()

test <- data_occ %>%
  mutate(Specie = sp.name) %>%
  filter(fac == 1) %>%
  #filter(fac %in% c(2,3)) %>%
  select(Specie, Longitude, Latitude) 

joint <- data_occ %>%
  mutate(Specie = sp.name) %>%
  select(Specie, Longitude, Latitude) 


write.csv(train, 'sp_train.csv', row.names = F)
write.csv(test, 'sp_test.csv', row.names = F)
write.csv(joint, 'sp_joint.csv', row.names = F)

##################################################################
###########    PARTE 3. Crop Environmental   ####################
#################################################################

# Models

# Variables a 5 km

pre <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/pre/5km',
                                         full.names=TRUE, pattern='.tif')))


f1 <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/fut/5km/cc50/cc45',
                                        full.names=TRUE, pattern='.tif')))
f2 <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/fut/5km/he50/he26',
                                        full.names=TRUE, pattern='.tif')))
f3 <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/fut/5km/he50/he45',
                                        full.names=TRUE, pattern='.tif')))
f4 <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/fut/5km/he50/he60',
                                        full.names=TRUE, pattern='.tif')))
f5 <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/fut/5km/he50/he85',
                                        full.names=TRUE, pattern='.tif')))
f6 <- suppressWarnings(stack(list.files('/Volumes/SEAGATE/clima/fut/5km/mc50/mc45',
                                        full.names=TRUE, pattern='.tif')))


# remover
fut<-c(f1,f2,f3,f4,f5,f6)
rm(f1,f2,f3,f4,f5,f6)
library(readr)
sp.m <- read_csv("sp_joint.csv"); sp.m <- as.data.frame(sp.m)


library(sdStaf)
library(rgdal)
bgeo <- readOGR('/Volumes/SEAGATE/2017/RECURSOS/biogeo/Lowenberg_Neto_2014.shp')

## Mediante zona buffer
mask.M <- stim.M(sp.m[,2:3], radio = 500, bgeo = bgeo)

points(sp.m[,2:3], pch=21, bg=2,cex=0.5)
crs(mask.M) <- c('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

r.dc <- reduce.env(env = pre, transfer = fut, occ_data = sp.m[,2:3], mask = mask.M)
rm(mask.M,bgeo)
#rm(f,f1,f2,f3,f4,f5,f10,f11,f12,f21,f22, transfer)


plot(r.dc@cropa$bio_1)
points(sp.m[,2:3], pch=16,col='blue', cex=0.5)

var.reduce <- c('bio_8','bio_19','bio_9','bio_18',
                'bio_1','bio_2','bio_3','bio_10','bio_11','bio_16','bio_12','bio_17','bio_6','bio_15'); cor.show(r.dc, rm=TRUE, var.rm = var.reduce)

pn <- dropLayer(r.dc@cropa, var.reduce)


for(i in 1:nlayers(pn)){
  writeRaster(pn[[i]],
              paste0('M_variables', "/set1/", names(pn[[i]]),".asc",sep=''),overwrite=TRUE)
}



f1n <- dropLayer(r.dc@project[[1]], var.reduce)
for(i in 1:nlayers(f1n)){
  writeRaster(f1n[[i]],
              paste0('G_variables', "/set1/",'fut_cc45/', names(f1n[[i]]),".asc",sep=''),overwrite=TRUE)
}


f2n <- dropLayer(r.dc@project[[2]], var.reduce)
for(i in 1:nlayers(f2n)){
  writeRaster(f2n[[i]],
              paste0('G_variables', "/set1/",'fut_he26/', names(f2n[[i]]),".asc",sep=''),overwrite=TRUE)
}

f3n <- dropLayer(r.dc@project[[3]], var.reduce)
for(i in 1:nlayers(f3n)){
  writeRaster(f3n[[i]],
              paste0('G_variables', "/set1/",'fut_he45/', names(f3n[[i]]),".asc",sep=''),overwrite=TRUE)
}


f4n <- dropLayer(r.dc@project[[4]], var.reduce)
for(i in 1:nlayers(f4n)){
  writeRaster(f4n[[i]],
              paste0('G_variables', "/set1/",'fut_he60/', names(f4n[[i]]),".asc",sep=''),overwrite=TRUE)
}


f5n <- dropLayer(r.dc@project[[5]], var.reduce)
for(i in 1:nlayers(f5n)){
  writeRaster(f5n[[i]],
              paste0('G_variables', "/set1/",'fut_he85/', names(f5n[[i]]),".asc",sep=''),overwrite=TRUE)
}

f6n <- dropLayer(r.dc@project[[6]], var.reduce)
for(i in 1:nlayers(f6n)){
  writeRaster(f5n[[i]],
              paste0('G_variables', "/set1/",'fut_mc45/', names(f6n[[i]]),".asc",sep=''),overwrite=TRUE)
}

rm(pn,f1n, f2n, f3n, f4n, f5n, f6n, pre)

########################################################################
##################     PARTE 4. Construcción de modelos   ###############
#########################################################################

# NOTA: ESCRIBIR LAS VARIABLES EN EL NUEVAS CARPETAS QUE SERAN UTILIZADAS
# ver (cortar variables ambientales)

library(kuenm)
library(ENMeval)
# 4.1. Calibracion de modelos ####
# Variables with information to be used as arguments. Change "YOUR/DIRECTORY" by your actual directory.
occ_joint <- "sp_joint.csv"
occ_tra <- "sp_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.5, 2, 0.5), seq(2.2, 6, 0.5))

f_clas <- c("l", "lq","h","lp" ,"lqp","lqh", "lqph", "lqpt", "lqpth",
            "lpth","lqth",'lpt','qph')

background <- 10000
maxent_path <- paste0(getwd())
wait <- FALSE
#invisible <- TRUE # Windows (solo para windows)
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, background = background,
          maxent.path = maxent_path, wait = wait, run = run)

# 4.2. Evaluacion de parametros del modelo ####
occ_test <- "sp_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 100
kept <- TRUE
selection <- "OR_AICc"
# Note, some of the variables used here as arguments were already created for previous function

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = FALSE)

# 4.3. Modelo final y transferencia ####

batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
G_var_dir <- "G_variables"
out_format <- "logistic"
project <- TRUE
ext_type <- "no_ext"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)

###########################################################
###########   PARTE 5. Modelos de estabilidad y dispersion   ##############
##########################################################

# Lectura del modelo final
#
file_model <- '/M_2.35_F_lpt_set1_NE/'
m.pre<- raster(paste0(batch_fin, file_model ,sp.name,'_median.asc'))
#cur<- raster(paste0(batch_fin, file_model ,sp.name,'_current_median.asc'))
e1<- raster(paste0(batch_fin, file_model,sp.name,'_fut_he26_median.asc'))
e2<- raster(paste0(batch_fin, file_model,sp.name,'_fut_he45_median.asc'))
e3<- raster(paste0(batch_fin, file_model,sp.name,'_fut_he60_median.asc'))
e4<- raster(paste0(batch_fin, file_model,sp.name,'_fut_he85_median.asc'))
e5<- raster(paste0(batch_fin, file_model,sp.name,'_fut_mc45_median.asc'))
e6<- raster(paste0(batch_fin, file_model,sp.name,'_fut_cc45_median.asc'))

## Analysis of species

## Variables de Estabilidad
stab.fut <- stack(e2, e5,e6)

## Variables - Mapas de dispersion
disp.vars <- stack(e1, e2, e3, e4)

joint <- read.csv('sp_joint.csv')
# Umbrales a la media
library(readr)
library(dplyr)
maxentResults <- read_csv(paste0(batch_fin,file_model,"maxentResults.csv"))

maxentResults[nrow(maxentResults), ] %>%
  select(Fixed_5 =`Fixed cumulative value 5 Logistic threshold`,
         Fixed_10 =`Fixed cumulative value 10 Logistic threshold`,
         MTP = `Minimum training presence Logistic threshold`,
         Perc_10 = `10 percentile training presence Logistic threshold`,
         MTSPS= `Maximum training sensitivity plus specificity Logistic threshold`,
         ETSPS = `Equal training sensitivity and specificity Logistic threshold`)


#### Stability  manual
Th_err = 0.1# Porcentaje de error
#
vn<-extract(m.pre, joint[,2:3]); vns<-length(vn) * Th_err; round(vns)
sort(vn,decreasing = T)[length(vn)-vns] - 0.01


#
library(sdStaf)
t <- 0.258


##### Current distribution 
ts <- c(0, t, 0, t, 1, 1)
mr <- matrix(ts, ncol = 3, byrow = TRUE)

current_dis <- reclassify(m.pre, mr)

#p1 <- reclassify(sp@proj[[2]], mr)
#p2 <- reclassify(s4, mr)
plot(current_dis)
points(joint[,2:3], pch=21, bg=2, cex=0.5)

trs <- maxentResults[nrow(maxentResults), ] %>%
  select(`10 percentile training presence Logistic threshold`) %>%
  as.numeric() * 1000
storage.mode(trs) = "integer"

sp.stab.fut <- stability(current = m.pre, future = stab.fut, thr.value = t)
plot(sp.stab.fut@map)


### Mapa continuo con umbral
t1 <- c(0, t, 0)
mr <- matrix(t1, ncol = 3, byrow = TRUE)

current_dis_con <- reclassify(m.pre, mr)
plot(current_dis_con)
#points(joint[,2:3], pch=21, bg=2, cex=0.7)

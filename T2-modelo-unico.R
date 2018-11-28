########################### CONSTRUIR MODELOS KU ############################
###                      P.J. Atauchi, A.T. Peterson & M. E. Cobos
##                                 Mayo 2018
#-------------------------------------------------------------------------
# 1. Definir directorio de trabajo
# 2. Calibrar modelos
# 3. Evaluacion de parametros del modelo
# 4. Evaluacion del modelo final
##

# Datos de ocurrencia
library(kuenm)
library(sdStaf)
library(ENMeval)
#### 1. Definir directorio de trabajo
setwd('/Users/p.joseratauchi/Desktop/Synallaxis_zimmeri') #Ruta de la carpeta Synallaxis_zimeri


### 2. Calibrar modelos
occ_joint <- "sp_joint.csv"
joint <- read.csv('sp_joint.csv')
occ_tra <- "sp_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
background <- 10000
maxent_path <- paste0(getwd())
wait <- TRUE
#invisible <- TRUE # Windows (solo para windows)
run <- TRUE

reg_mult <- c(seq(0.5, 6, 0.8))

f_clas <- "basic"

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, background = background,
          maxent.path = maxent_path, wait = wait, run = run)

### 3. Evalución de los modelos
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

# 4. Modelo final y transferencia ####
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

### 4. Analisis de datos
file_model <- '/M_0.5_F_lq_set1_NE/'
presente <- raster(paste0(batch_fin, file_model ,'Synallaxis_zimmeri','_median.asc'))
futuro45 <- raster(paste0(batch_fin, file_model,'Synallaxis_zimmeri','_fut_he45_median.asc'))
futuro85 <- raster(paste0(batch_fin, file_model,'Synallaxis_zimmeri','_fut_he85_median.asc'))

### Analisis de datos

presente

plot(presente, main='Distribución Actual')
points(joint[,2:3], pch=21, bg=2,cex=0.7)

plot(futuro45, main='Distribución Futuro optimista')
points(joint[,2:3], pch=21, bg=2,cex=0.7)

plot(futuro85, main='Distribución Futro pesimista')
points(joint[,2:3], pch=21, bg=2,cex=0.7)

### Analisis estadistico

valores_presente <- rasterToPoints(presente, spatial = T)@data$Synallaxis_zimmeri_median

valores_futuro45 <- rasterToPoints(futuro45, spatial = T)@data$Synallaxis_zimmeri_fut_he45_median

valores_futuro85 <- rasterToPoints(futuro85, spatial = T)@data$Synallaxis_zimmeri_fut_he85_median


#### Test de Normailidad
require(nortest)


## Test de significancia


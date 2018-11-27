x <- c(1:10)
class(x)

# Ejemplo: Clase lapiz
Lapiz <- setClass("Lapiz",
                     slots = c(tamaño="numeric",
                               color = "character",
                               punta="character"))

Tecnico_2B <- new('Lapiz', tamaño = 10,
                 color = 'verde oscuro',
                 punta = '0.5 mm')

Tecnico_3B <- new('Lapiz', tamaño = 8,
                  color = 'verde oscuro',
                  punta = '0.7 mm')

Lapicero <- new('Lapiz', tamaño = 13,
                  color = 'Azul oscuro',
                  punta = '0.4 mm')



## Clases usados en SIG
# Librerias
library(sdStaf)
library(dismo)
library(raster)

predictor <- stack(list.files(path=paste(system.file(package="dismo"),'/ex', sep=''),
                              pattern='grd', full.names=TRUE ))

# Mostrar clase del objeto 'predictor'
class(predictor)

# Graficar los objetos de 'predictor'
plot(predictor)

# Mostrar clase del objeto predictor$bio1
class(predictor$bio1)

# Graficar el objeto
class(predictor$bio1)







### Ejemplo 

library(raster)

# RasterLayer con parametros por defecto
x <- raster()

# Con otros parametros
x <- raster(ncol=36, nrow=18, xmn=-1000, xmx=1000, ymn=-100, ymx=900)

# que puede cambiar?
res(x)

# change resolution
res(x) <- 100

res(x)

ncol(x)

# Cambiar el numero de columnas (afecta la resolucion)
ncol(x) <- 18
ncol(x)

# Definir la projeccion
projection(x) <- "+proj=utm +zone=18 +datum=WGS84"







#### Class puntos
library(sdStaf)
library(raster)

# Cargar datos Phytotoma raimondii
data('phytotoma')

# Crear la clase puntos
phytotoma_puntos <- SpatialPoints(phytotoma[,2:3])

# Mostrar la clase del objeto 'phytotoma_puntos'
class(phytotoma_puntos)

# Graficar el objeto 'phytotoma_puntos'
plot(phytotoma_puntos)











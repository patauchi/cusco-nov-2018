########################################################################
#
#                         Usando R como GIS
#                         P. Joser Atauchi
#                  Museo de Historia Natural de Cusco
#                          28 de nov 2018
#----------------------------------------------------------------------
#
# Instalar Librerias
list.of.packages <- c("sdStaf", "dismo", "spThin", "spocc", "devtools",
                      "tidyverse", "ENMeval", "SDMTools", "raster", "rgdal", "rgeos", "readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Instalar librera modelos de nicho
devtools::install_github('marlonecobos/kuenm')


#
#

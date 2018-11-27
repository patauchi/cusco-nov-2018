# Curso de R aplicado a GIS

Usando sistemas de información geográfica(SIG) aplicado al modelamiento de nicho ecológico (ENM).

## Pasos:

1. Descargar el archivo "Instalar-paquetes-gis.R".

2. Abrir el archivo con RStudio.

3. Ejecutar todas las lineas de comando. Estas lineas permite installar todos los paquetes necesarios para elborar los modelos de nicho ecológico y usar R como GIS.

4. Descargar e instalar RTools.exe desde [https://cran.r-project.org/bin/windows/Rtools/]

5. Ejecutar la siguiente linea de código para probar si los paquetes fueron instalados correctamente

Importante ...!!!

``` r
require(sdStaf)
require(kuenm)
require(SDMTools)
require(spThin)
require(raster)
require(dismo)

```

6. Instalar Java JDK software [https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html]. Asegurarse instalar versión x86 o X64.

7. Instalar GDAL [https://trac.osgeo.org/osgeo4w/]

8. Instalar rJava package y comprobar

``` r
install.packages("rJava")
require(rJava)
```
9. Descargar MaxEnt desde [https://biodiversityinformatics.amnh.org/open_source/maxent/]

10. Estamos listos para trabajar!.


NOTA: En caso error contactar a *patauchi@gmail.com*


## Entorno de trabajo
El modelamiento de nicho ecológico trabaja mediante el siguiente esquema de trabajo.

![Alt text](https://github.com/patauchi/cusco-nov-2018/blob/master/tb.png "Frameworks")

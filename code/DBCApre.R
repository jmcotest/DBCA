rm(list=ls(all=TRUE))
library(writexl)
#Generando estructuras de carpetas
# Eliminar la carpeta salidas y todo su contenido
if (dir.exists(file.path(".", "salidas"))) {
  unlink(file.path(".", "salidas"), recursive = TRUE)
}

# Crear nuevamente las carpetas
for (carpeta in c("graficos", "excel", "pdf", "html")) {
  dir.create(
    file.path(".", "salidas", carpeta),
    recursive = TRUE,
    showWarnings = FALSE
  )
}
dir.create(file.path(".","script"),showWarnings = FALSE, recursive=TRUE)
dir.create(file.path(".","rmd"),showWarnings = FALSE, recursive=TRUE)
dir.create(file.path(".","qmd"),showWarnings = FALSE, recursive=TRUE)
dir.create(file.path(".","BD"),showWarnings = FALSE, recursive=TRUE)

# Datos de ejemplo
#' 
#' # Ejercicio
#' Hoppe registró observaciones con el objeto de comparar 7 fungicidas para semillas y un control sin tratamiento respecto a la emergencia de plántulas de maíz infectado con _Diplodia spp._ El experimento se llevó a cabo en un invernadero con seis bloques en un diseño en bloques completos al azar. Cada unidad experimental consistía en 25 semillas. Los datos y tratamientos codificados se dan en la tabla.
#' 
#' Los tratamientos corresponden a:
#' 
#' - A control sin aplicación.
#' - B y C fungicidas mercúricos.
#' - D y H funguicida no mercurico compañia 1
#' - E, F y G funguicida no mercurico compañia 2
#' - F y G son nuevas formulaciones de E.
#' 
#' (Tomado de: Bioestadística Steel, R. G. D., & Torrie, J. H. (1988). Bioestadística: principios y procedimientos (2.ª ed.). McGraw-Hill. p.198-199)
#' 
#' ## Datos Experimentales
#' 
## ----echo=FALSE---------------------------------------------------------------
datos <- data.frame(
  Bloque = factor(rep(1:6, times = 8)),
  Tratamiento = factor(rep(c("A", "B", "C", "D", "E", "F", "G", "H"), each = 6)),
  Emergencia = c(
    8, 8, 9, 7, 7, 5,       # A
    16, 19, 24, 22, 19, 19, # B
    14, 16, 14, 13, 14, 13, # C
    10, 11, 12, 8, 7, 3,    # D
    8, 7, 1, 1, 3, 2,       # E
    8, 8, 3, 3, 3, 7,       # F
    7, 6, 6, 6, 4, 4,       # G
    12, 19, 9, 11, 9, 5     # H
  )
)

#####################################
# 1. ANÁLISIS GRÁFICO DE LOS DATOS  #
#####################################
source(file.path(".","code","funciones.R"))

#Estadistica descriptiva
descriptiva<-aggregate(Emergencia~Tratamiento,data=datos,resumen)
# Operaciona previas a la gráfica descriptiva
descriptiva<-data.frame(Tratamiento=descriptiva$Tratamiento,descriptiva$Emergencia)

# Uso de la funcion graficos
graficos_personalizado(
  data_g = descriptiva,formula_g=Tratamiento~Media+Mediana+Inferior+Superior,
  datos=datos,formula_m=Emergencia~Tratamiento,
  ordenar="Media", typeplot="lineas-violines", prefixjpg="inicial",
  xlab="Fungicidas para semillas", 
  ylab="Plántulas emergidas por 25 semillas")
  #Para el de lineas
  #Fig 1. Gráfico descriptivo para el contenido de nitróneno en 6 tratamientos de cepas de Rhizoboum.
  #.      La barras indican la media del tratamiento la cual también está resaltada en negro. La Mediana
  #       se resalta en color naranja. Las observaciones originales están en color azul y las barras 
  #.      de error indican los límites del intervalo de confianza de 95%.

#####################################
# 2. MODELACIÓN ESTADÍSTICA         #
#####################################
modelo<-formula(Emergencia~Tratamiento+Bloque)

m <- lm(modelo, data = datos)
summary(m)
anova(m)

library(MCMCglmm)
m1 <- MCMCglmm(Emergencia ~ Tratamiento, random = ~ Bloque,
               data = datos,
               rcov = ~ idh(Tratamiento):units,
               nitt = 1e6,
               verbose = FALSE)

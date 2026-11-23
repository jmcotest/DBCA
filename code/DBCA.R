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


###############################
# 2.1 ANALISIS CLÁSICO        #
###############################

######################################################
# 2.1.1. Las varianzas homogeneas o hereogéneas      #

# NOTA IMPORTANTE:
#   ESTO SÓLO SE HACE SÍ POR TEORÍA TENEMOS DUDAS 
#   DE LA HOMOGENEIDAD/HETEROGENEIDAD DE LAS VARIANZAS.

#   ASÍ A CONTINUACIÓN SE PRESENTA EL CÓDIGO PARA
#   DETERMINAR HOMOGENEIDAD DE VARIANZAS 
#   POR RAZONES EMPÍRICAS.

# RECORDAR: 
#   LAS HIPÓTESIS NULAS JAMÁS SE PUEDEN PROBAR, 
#   SOLO SE PUEDEN RECHAZAR.
####################################################

#Es posible hacerlo pero no es tan fácil. Shukla (1972,1973)

######################################################
# 2.1.2. ANOVA varianza homogenea     #
  m <- lm(modelo, data = datos)
  summary(m)

  anova_m<-anova(m)
  anova_m

# Prueba de Shapiro-Wilk
sh<-shapiro.test(residuals(m))
sh

######################################################
# 2.1.3. Cálculo de los predichos obtenidos por mínimos cuadrados - LSMEAN     #
# Deben usarse en el caso de desbalance por repeticiones
media <- emmeans::emmeans(m, ~ Tratamiento)

#######################################################
# 2.1.4. Prueba de comparacion múltiple     #
#Opcion 1: Usando la función TukeyHSD
tukey <- TukeyHSD(aov(m), "Tratamiento")
tukey

library(multcomp)
library(multcompView)
letras <- multcomp::cld(media,adjust = "tukey", Letters = letters, 
  decreasing = TRUE, alpha = 0.05)
letras

#####################################
# 3. COMPLEMENTAR Y AJUSTAR GRÁFICO #
#####################################
####################################################################################
# Conversión del resultado a un data.frame
tabla_letras <- as.data.frame(letras)
# Eliminación de espacios adicionales en la columna de letras
tabla_letras$.group <- trimws(tabla_letras$.group)
# Visualización de la tabla final
tabla_letras

analisis_post_modelo <- merge(descriptiva,tabla_letras, by.x = "Tratamiento",
 by.y = "Tratamiento", all.x = TRUE)

analisis_post_modelo[["LS Mean"]] <- analisis_post_modelo[["emmean"]] 

graficos_personalizado(
  data_g = analisis_post_modelo,formula_g=Tratamiento~`LS Mean`+Mediana+lower.CL+upper.CL,
  datos=datos,formula_m=Emergencia~Tratamiento,
  ordenar="LS Mean", typeplot="lineas-violines",
  prefixjpg="clasico_homogeneo",
  xlab=expression("Cepas de " * italic("Rhizobium") * " sp."), 
  ylab="Nitrógeno (mg)", 
  add = function(g) {
    with(g,{
      text(x=bp, y=upper.CL, labels=.group, col = "#1B5E20", cex=0.9, offset=0.5,pos=3L)
    })
  }
)
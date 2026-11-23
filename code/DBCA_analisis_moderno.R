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
#' Hoppe registró observaciones con el objeto de comparar 7 fungicidas para semillas y un control sin 
#' tratamiento respecto a la emergencia de plántulas de maíz infectado con _Diplodia spp._ El experimento se llevó a cabo en un invernadero con seis bloques en un diseño en bloques completos al azar. Cada unidad experimental consistía en 25 semillas. Los datos y tratamientos codificados se dan en la tabla.
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
# 2.2 ANALISIS MODERNO        #
###############################

######################################################
# 2.1.1. Las varianzas homogeneas o hereogéneas      #
# RECORDAR:
#  Esto solo se hace si es necesario, es decir, si por teoría no se puede deducir 
#  que las varianzas sean homogéneas o heterogéneas
library(MCMCglmm)

##Prueba de homogenedida de varianzas aprox. Lo cierto es que prueba plausibliidad de modelos.
set.seed(458796)
m0<-MCMCglmm(Emergencia~Tratamiento, random=~Bloque, 
    data = datos, nitt=1e6L,verbose=FALSE)

k <- nlevels(datos$Tratamiento)
prior1 <- list(
   R = list(
     V = diag(k),nu = k + 2),
   G = list(
     G1 = list(
       V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)
   )
 )

m1<-MCMCglmm(Emergencia~Tratamiento, random=~Bloque, 
  rcov=~idh(Tratamiento):units,prior = prior1,
  data = datos, nitt=1e6L,verbose=FALSE)

m0$DIC
m1$DIC
m0$DIC-m1$DIC

#' 
#' Interpretación:
#' 
#' - Menor DIC: mejor ajuste considerando complejidad.
#' - $\Delta$DIC < 2: evidencia débil a favor de m1.
#' - $\Delta$DIC entre 3 y 7: evidencia moderada a favor de m1.
#' - $\Delta$DIC > 10: evidencia fuerte a favor de m1.

# Modelo escogido con varianzas heterogeneas
m<-MCMCglmm(Emergencia~Tratamiento, random=~Bloque, 
  rcov=~idh(Tratamiento):units,prior = prior1,
  data = datos, nitt=1e6L,verbose=FALSE)

## Cálculo de los predichos obtenidos Bayes - Bayesian Linear Prediction
# Predicciones de todas las iteraciones MCMC
tmp<- m[["X"]] %*% t(m[["Sol"]])
grupo <- datos[["Tratamiento"]]
indices <- split(seq_len(nrow(tmp)), grupo)

MCMC_trt <- lapply(indices, function(i) colMeans(tmp[i, , drop = FALSE]))

#Resumen de distribucion *a posteriori*
resumenMCMC<-lapply(MCMC_trt,function(x) c(mean(x),HPDinterval(as.mcmc(x),prob=0.90)))
resumenMCMC<-do.call(rbind,resumenMCMC)
colnames(resumenMCMC)<-c("Est Bayes", "lower.HPD", "upper.HPD")

# Analisis de la distribución *a posteriori* para cada tratamiento: 
# Probabilidad de ser el mejor y el peor tratamiento
# Esto es mejor que las pruebas de comparación multiples en cuanto a una toma 
# de decisiones más clara.
# El autor quien realizó la prouesta es "Julian Besag, David Higdon, 
# Bayesian Analysis of Agricultural Field Experiments, 
# Journal of the Royal Statistical Society Series B: Statistical Methodology, 
# Volume 61, Issue 4, November 1999, Pages 691–746, 
# https://doi.org/10.1111/1467-9868.00201"

MatrixMCMC<-as.matrix(do.call(cbind,MCMC_trt))
rank_MatrixMCMC<-apply(MatrixMCMC,1L,rank)
pbest<-rowSums(rank_MatrixMCMC==nlevels(grupo))/ncol(rank_MatrixMCMC)
pworst<-rowSums(rank_MatrixMCMC==1)/ncol(rank_MatrixMCMC)

pbest_c<-format(round(pbest,4L), nsmall=4L) 
pworst_c<-format(round(pworst,4L), nsmall=4L)

#Nuevo probabilidad de superar un valor de referencia
p_ref<-sapply(MCMC_trt,function(x) sum(x > 15L) / length(x))
p_ref_c<-sapply(p_ref,function(x) format(round(x,4L), nsmall=4L))

#Nuevo probabilidad de superar el testigo Compuesto
p_check<-sapply(MCMC_trt,function(x) sum(x > MCMC_trt[["Compuesto"]]) / length(x))
p_check_c<-sapply(p_check,function(x) format(round(x,4L), nsmall=4L))
p_check_c["A"]<-""
#Juntando complementos
complementos<-data.frame(Tratamiento=rownames(resumenMCMC), resumenMCMC, pbest=pbest_c, pworst=pworst_c, pref=p_ref_c, pcheck=p_check_c)

#####################################
# 3. COMPLEMENTAR Y AJUSTAR GRÁFICO #
#####################################
analisis_post_modelo <- merge(descriptiva,complementos, by.x = "Tratamiento",
 by.y = "Tratamiento", all.x = TRUE)

analisis_post_modelo[["E Bayes"]] <- analisis_post_modelo[["Est.Bayes"]] 

graficos_personalizado(
  data_g = analisis_post_modelo,formula_g=Tratamiento~`E Bayes`+Mediana+lower.HPD+upper.HPD,
  datos=datos,formula_m=Emergencia~Tratamiento,
  typeplot="lineas-violines",
  ordenar="E Bayes",
  prefixjpg="moderno_heterogeneo", 
  xlab="Fungicidas para semillas", 
  ylab="Plántulas emergidas por 25 semillas",
  add = function(g) {
    with(g, {
      text(x=bp, y=lower.HPD, labels=pworst, col = "#000000", cex=0.6,offset=0.2,pos=1L)
      text(x=bp, y=upper.HPD, labels=pbest, col = "#000000", cex=0.6,offset=0.2,pos=3L)
      #abline(h=15,lty=2,lwd=2,col = "#8B0000")
      #text(x=bp, y=15, labels=pref, col = "#8B0000", cex=0.6,adj=0,offset=0.2,pos=3)
      #text(x=bp, y=`E Bayes`, labels=pcheck, col = "#E65100", cex=0.6,adj=0,offset=0.4,pos=2)
    })
  }
)

# 
# Gráfica resumen más completo basado en el resumen de la distribución *a posteriori*
# Fig n. Número de plantas emergidas en diferentes tratamiento de control de Diplodia spp._
#        Las barras de error indican los límites del intervalo de mas alta densidad a posterior de 90%. En la
#        parte superior de límite se indica la probabilidad de que ese tratamiento sea el mejor mientras en el
#        límite inferior se indica la probabilidad de ser el peor. En rojo se indica la probabilidad de que ese
#        tratamiento supere emergencia de 15 plantas. Al lado izquiero de la media se indica la probabilidad de
#        superar el testigo A. La línea punteada indica el valor de referencia de 15 plantas.
library(tsibble)
library(tidyverse)
library(sf)
# library(leaflet)
# library(deldir)
library(spatstat)
# library(tigris)
library(spatstat.geom)
library(spatstat.linnet)
library(spNetwork)
library(sp)
library(sfnetworks)
# library(tmap)
library(gridExtra)
library(grid)
# library(jpeg)
# library(raster)

# Smulacion --------------------------------------------


set.seed(10)
# Parámetros
lambda_parents <- 20  # Intensidad del proceso Poisson de los padres
mu_offspring <- 10    # Media de la cantidad de hijos por padre
R <- 0.1              # Radio del círculo alrededor del padre
par(mfrow = c(1, 3), mar=c(5, 4, 1, 1), oma=c(5, 5, 5, 5)) 
parents <- rpoispp(lambda_parents, win=owin(c(0, 1), c(0, 2)))
plot(parents, main="Puntos Padres", pch=19, cols = "blue")
#Generar puntos hijos alrededor de cada padre
offspring <- NULL
plot(parents, main="Discos de radio R", pch=16, xlab="", ylab="", asp=1, cols = "blue")
for (i in 1:parents$n) {
  num_offspring <- rpois(1, mu_offspring)
  if (num_offspring > 0) {
    angles <- runif(num_offspring, 0, 2 * pi)
    radii <- R * sqrt(runif(num_offspring))
    xoffs <- parents$x[i] + radii * cos(angles)
    yoffs <- parents$y[i] + radii * sin(angles)
    offspring <- rbind(offspring, data.frame(x = xoffs, y = yoffs))
  }
# Dibujar los discos
  symbols(parents$x[i], parents$y[i], circles = R, add = TRUE, inches = FALSE, fg = "red")
}

#  Dibujar sólo los puntos hijos dentro de los círculos,sin los puntos padres
offspring_points <- ppp(offspring$x, offspring$y, window = owin(c(0, 1), c(0, 2)))
plot(offspring_points, main="Proceso cluster de Matérn", pch=16, xlab="", ylab="", asp=1)
for (i in 1:parents$n) {
  symbols(parents$x[i], parents$y[i], circles = R, add = TRUE, inches = FALSE, fg = "red")
}



# Conteo por cuadrantes y test CRS--------------------------------------------

win <- owin(c(0, 1), c(0, 2))
offspring_points
plot(quadratcount(offspring_points, nx = 8, ny = 8))
quadrat.test(offspring_points, nx = 8, ny = 8)

# Calculo de intensidades y funcion K no homogenea--------------------------------------------

par(mfrow = c(2, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))

# Graficar cada densidad en su respectivo lugar
plot(density.ppp(offspring_points, sigma = 0.01, edge = TRUE), main = "Ancho de Banda: 0.01")
plot(density.ppp(offspring_points, sigma = 0.07, edge = TRUE), main = "Ancho de Banda: 0.07")
plot(density.ppp(offspring_points, sigma = 0.1, edge = TRUE), main = "Ancho de Banda: 0.1")
plot(density.ppp(offspring_points, sigma = 0.2, edge = TRUE), main = "Ancho de Banda: 0.2")

K_functioninom <- Kinhom(offspring_points,r = seq(0,0.15,by=0.01),correction = "border")
summary(K_functioninom)
plot(K_functioninom, main = "")

# Carga y tratamiento de datos reales--------------------------------------------

BASE<-read.csv("obs_data_112Valencia (1).csv", sep=",", header=TRUE) 
puntos<-BASE|>filter(year<=2019) |> 
  st_as_sf(coords=c("crime_lon","crime_lat"),crs= 4258)
ventana <- st_read("shape files-20221031T083303Z-001/shape files/valencia_outline.shp",quiet=TRUE)
calles <- st_read("shape files-20221031T083303Z-001/shape files/valencia_road.shp",quiet=TRUE)
calles<-st_transform(calles, crs = st_crs(puntos))
barrio<-st_read("barris-barrios/barris-barrios.shp",quiet=TRUE)
a<-"112" #codigo de barrio el Cabanyal
barrio<-st_read("barris-barrios/barris-barrios.shp",quiet=TRUE)
barrio<-barrio |> filter(coddistbar==a)
barrio<-st_transform(barrio, crs = st_crs(puntos))
interseccion<-st_intersection(calles,barrio)
puntos <- st_transform(puntos,st_crs(barrio))
puntos_barrio<-st_intersection(puntos,barrio)
calles_barrio<- interseccion[st_geometry_type(interseccion) == "LINESTRING", ]

# Construccion del proceso puntual. Ventana: El Cabanyal --------------------------------------------

ventana<- st_read("barris-barrios/barris-barrios.shp",quiet=TRUE)
ventana<-ventana |> filter(coddistbar==a)
ventana1 <- st_transform(ventana, crs = st_crs("EPSG:32630"))
ventana2<-st_union(ventana1[1])
window<-st_coordinates(st_reverse(ventana2[1]))
a<-st_coordinates(st_reverse(ventana2[1]))
window<-owin(poly=a[,c("X","Y")])
target_crs <- st_crs("+proj=utm +zone=30 +north +datum=WGS84 +units=m +no_defs")
a_sf_reprojected <- st_transform(puntos, target_crs)
app <- ppp( x = st_coordinates(a_sf_reprojected)[,1],
            y = st_coordinates(a_sf_reprojected)[,2],
            marks = c(1:3221),
            window = window)

unique(app, rule="deldir")
options(repr.plot.width = 50, repr.plot.height = 30)
plot(quadratcount(app, nx = 13, ny = 8))
quadrat.test(app, nx = 13, ny = 8)


# Construccion de la grilla --------------------------------------------

grid <- ventana%>%
  st_make_grid(n = c(13,8)) %>%
  st_intersection(ventana) %>% 
  st_cast("MULTIPOLYGON") %>%
  st_sf()
plot(grid)
puntos_barrio <- st_transform(puntos_barrio, crs = st_crs(grid))

inter <- st_intersects(grid, puntos_barrio)
grid$count <- lengths(inter)
ggplot(grid) + geom_sf(aes(fill = count))


ggplot(grid) + geom_sf(aes(fill = count))+
  scale_fill_continuous(type = "viridis")

grid.selec <- grid
# Crear el gráfico
grilla <- ggplot(grid.selec) + 
  geom_sf(aes(fill = count)) +
  geom_sf_text(aes(label = count), size = 3, color = "white") +  # Ajusta el tamaño aquí
  scale_fill_continuous(type = "viridis") +
  theme(
    axis.text = element_text(size = 12, face = "bold"),  # Agrega negrita aquí
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) +
  theme(
    plot.margin = margin(2, 2, 2, 2, "mm"),
    panel.spacing = unit(0, "lines")
  )


# Calculo de intensidades, ventana: El Cabanyal --------------------------------------------
plot(density.ppp(app, sigma = 0.5, edge = TRUE), main = "Ancho de Banda: 0.5")
plot(density.ppp(app, sigma = 15, edge = TRUE), main = "Ancho de Banda: 15")
plot(density.ppp(app, sigma = 55, edge = TRUE), main = "Ancho de Banda: 55")
plot(density.ppp(app, sigma = 120, edge = TRUE), main = "Ancho de Banda: 120")


# Calculo de K-Ripley, ventana: El Cabanyal --------------------------------------------
K_functioninom <- Kinhom(app,r = seq(0,180,by=10))
summary(K_functioninom)



# Construccion del proceso puntual en red lineal.--------------------------------------------

red_vial_ln<- as_sfnetwork(calles_barrio)
plot(red_vial_ln, col = '#017308', pch = 18, lwd = 1, cex = 2,graticule = TRUE, axes = TRUE)
calles_barrio <- st_transform(calles_barrio, crs = "+proj=utm +zone=30 +datum=WGS84")
puntos_barrio <- st_transform(puntos_barrio, crs = "+proj=utm +zone=30 +datum=WGS84")
red_vial_ln <- st_transform(red_vial_ln, crs = "+proj=utm +zone=30 +datum=WGS84")
red_vial_ln <- as.linnet(red_vial_ln)
lpp_VALE <- lpp(app,red_vial_ln)# correccion mediante proyeccion de los puntos.

# Calculo de la intensidad y de K sobre la red lineal.--------------------------------------------

lambda <- density.lpp(lpp_VALE,sigma = 77.26966 ,distance = "euclidean")
par(mfrow = c(1, 1), mar = c(0, 0, 2, 0), oma = c(0, 0, 0, 0))
plot(lambda,main="Ancho de Banda: 77.27")
lpp_VALE<-unmark(lpp_VALE)
m1<-lppm(lpp_VALE ~ 1)
r <- seq(0, 300, by=50)
k_finom<-linearKinhom(lpp_VALE,lambda = m1,r)
plot(k_finom,main="Función K no homogénea")
plot(lpp_VALE, pch=19,main=" ",cex = 1.2,cols = "red")
plot(red_vial_ln)


# Construccion de variables complementarias espaciales.--------------------------------------------

par(mfrow = c(4, 3), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), mai = c(0, 0, 0.2, 0))
pnight <- puntos_barrio$atm_dist
Datm <- density(app, weights=pnight, sigma=77.26966)
plot(Datm,main="Atm (Datm)")

pnight <- puntos_barrio$bank_dist
Dbank <- density(app, weights=pnight, sigma=77.26966)
plot(Dbank,main="Banco (Dbank)")

pnight <- puntos_barrio$bar_dist
Dbar <- density(app, weights=pnight, sigma=77.26966)
plot(Dbar,main="Bar (Dbar)")

pnight <- puntos_barrio$cafe_dist
Dcafe <- density(app, weights=pnight, sigma=77.26966)
plot(Dcafe,main="Cafetería (Dcafe)")

pnight <- puntos_barrio$industrial_dist
Dindustrial <- density(app, weights=pnight, sigma=77.26966)
plot(Dindustrial,main="Zona Industrial (Dindustrial)")

pnight <- puntos_barrio$market_dist
Dmarket <- density(app, weights=pnight, sigma=77.26966)
plot(Dmarket,main="Mercado(Dmarket)")

pnight <- puntos_barrio$nightclub_dist
Dnight <- density(app, weights=pnight, sigma=77.26966)
plot(Dnight,main="Nightclub (Dnight)")

pnight <- puntos_barrio$police_dist
Dpolice <- density(app, weights=pnight, sigma=77.26966)
plot(Dpolice,main="Dep Policia (Dpolice)")


pnight <- puntos_barrio$pub_dist
Dpub <- density(app, weights=pnight, sigma=77.26966)
plot(Dpub,main="Pub (Dpub)")

pnight <- puntos_barrio$restaurant_dist
Drestaurant <- density(app, weights=pnight, sigma=77.26966)
plot(Drestaurant,main="Restaurante (Drestaurant)")


pnight <- puntos_barrio$taxi_dist
Dtaxi <- density(app, weights=pnight, sigma=77.26966)
plot(Dtaxi,main="Estacionamiento Taxis (Dtaxi)")

# Modelo inicial seleccionado.--------------------------------------------

m3<-lppm(lpp_VALE ~ Datm+Dbank+Dbar+Dcafe+Dindustrial+Dmarket+Dnight+Dpolice)

# Proceso de cluster de matern sobre una red lineal--------------------------------------------

# Centers = (x,y)-coordenadas de los puntos padres
# R = El  parametro R del proceso de cluster de  Matern 
# alpha = El parametro alfa del proceso de cluster de Matern 
# LL = La red lineal sobre la que se simula cada patron de puntos.
rMatClustlpp <- function(Centers, R, alpha, LL) {
  X <- array(0,0)
  Y <- array(0,0)
  for(p in 1:length(Centers$data$x)) {
    BBCOutD <- disc(radius=R, centre=c(Centers$data$x[p],
                                       Centers$data$y[p]),npoly = 32)
    BBCD <- intersect.owin(LL$window, BBCOutD)
    if(volume(LL[BBCOutD])>0) {
      Xp <- rpoislpp(alpha/volume(LL[BBCOutD]), L=LL[BBCD])
      X <- append(X, as.numeric(Xp$data$x))
      Y <- append(Y, as.numeric(Xp$data$y))
    }
  }
  lpp(cbind(X,Y), LL)
}


# Algoritmo de hot spots--------------------------------------------
set.seed(2024)
##  Prediccion de la intensidad del modelo inicial ajustado----
EIP <- predict(m3)
##  Ajuste de dominio para posibles valores resultantes de R y alpha----
r <- seq(0, 300, by=50)
nsim <- 5
valpha <- seq(5, 100, by=20)
vR <- seq(50, 800, by=100)
k_finom<-linearKinhom(lpp_VALE,lambda = m3,r)

## Proceso de determinacion de valores R y alpha adecuados(puede tardar al rededor de 2 horas)----

Contrast <- array(0, c(length(valpha), length(vR)))
for(i in 1:length(valpha)) {
  for(j in 1:length(vR)) {
    # calculo de K media mediante simulaciones.
    KMC <- array(0,length(r))
    for(s in 1:nsim) {
      # Centros de un proceso de Poisson
      Centers <- rpoislpp(EIP/valpha[i], L=lpp_VALE[['domain']])
      XX <- rMatClustlpp(Centers, vR[j], valpha[i], lpp_VALE[['domain']])
      KMC <- KMC + linearKinhom(XX, lambda = m3, r=r)$est
    }
    # Calcular la diferencia entre la estimación y la media de K
    Contrast[i,j] <- sqrt(sum((k_finom$est-KMC/nsim)^2))
  }
}
plot(as.im(Contrast), main="Valores de contraste",cex.main = 0.8)
# Encontrar el valor mínimo del contraste
id <- which(Contrast == min(Contrast), arr.ind = TRUE)
alpha <- valpha[id[,1]]
R <- vR[id[,2]]
## Valores elegidos----
R
alpha
## Mostrando resultados obtenidos----
# Centros de un proceso de Poisson
Centers <- rpoislpp(EIP/alpha, L=lpp_VALE[['domain']])
plot(Centers, col = 'red', pch = 16, lwd = 1, cex = 1,graticule = TRUE, axes = TRUE)
XX <- rMatClustlpp(Centers, R, alpha, lpp_VALE[['domain']])
plot(XX, pch=19,main="Realización del proceso",cex = 0.5,cols = "red")
## Calculo de la intensidad y de K ----
lambda_modelo <- density.lpp(XX,sigma = 70 ,distance = "euclidean")
plot(lambda_modelo)
m_modelo<-lppm(XX ~ 1)
r <- seq(0, 300, by=50)
k_finom_modelo<-linearKinhom(XX,lambda = m_modelo,r)
plot(k_finom_modelo)


install.packages("flexclust")
install.packages("geometry")
install.packages("ggbump")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
install.packages("../../../narrative/diagram/narrativecharts_0.1.0/narrativecharts_0.1.0.tar.gz", type="source")

library(flexclust)
library(cluster)
library(tidyverse)
library(stringi)
library(FactoMineR)
library(factoextra)
library(geometry)
library(stringi)
library(plotly)
library(narrativecharts)
library(strucchange)
library(ggbump)
library(jpeg)
library(htmlwidgets)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#### 0. Functions####
angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

llenararith<-function(x){
  for(i in 1:length(x)){
    k=1
    if(!is.na(x[i])|(i==length(x))|(i==1)){
      x[i]<-x[i]
    }
    else{
      while((is.na(x[i+k]))&(i+k<length(x))){
        k=k+1
      }
      x[i]<-x[i-1]+((x[i+k]-x[i-1])/(k+1))
    }
  }
  return(x)
}
llenararith2<-function(x){
  for(i in 1:length(x)){
    k=1
    if(!is.na(x[i])|(i==length(x))|(i==1)){
      x[i]<-x[i]
    }
    else if(is.na(x[i-1])){
      while((is.na(x[i+k]))&(i+k<length(x)+1)){
        k=k+1
      }
      x[i]<-x[i+k]
    }
    else{
      while((is.na(x[i+k]))&(i+k<length(x))){
        k=k+1
      }
      if(!is.na(x[i+k])){
        x[i]<-x[i-1]+((x[i+k]-x[i-1])/(k+1))
      }
      else{
        x[i]<-x[i-1]
      }
    }
  }
  for(i in 1:length(x)){
    if(i==1&is.na(x[i])){
      x[i]<-x[i+1]
    } else if(i==length(x)&is.na(x[i])){
      x[i]<-x[i-1]
    } else{
      x[i]<-x[i]
    }
  }
  return(x)
}
avgsil<-function(df, k){
  km.res<-kmeans(df, centers=k, iter.max = 20, nstart = 25)
  ss<-silhouette(km.res$cluster, dist(df))
  return(mean(ss[,3]))
}
####Datos####
facultades<-read_csv("facultades.csv")

####Analisis ####
####Descriptivo####
jpeg("descriptivo1.jpg", width=800, height=800)
ggplot(facultades, aes(x=year, y=n, fill=general))+
  geom_col(position="stack", show.legend = F)+
  theme_minimal()+labs(title="Número de estudiantes por año", x="Año", y="Número de estudiantes", fill="Entidad")
dev.off()
jpeg("descriptivo2.jpg", width=800, height=800)
ggplot(facultades, aes(x=year, y=docen_totalf, fill=general))+
  geom_col(position="stack", show.legend = F)+
  theme_minimal()+labs(title="Número de académicos por año", x="Año", y="Número de académicos", fill="Entidad")
dev.off()
#####Tablas generales#####
total2<-facultades %>% group_by(year) %>% nest()

set.seed(1917)
tableU<-data.frame()
tableV<-data.frame()
tableL<-data.frame()
for (i in 1:nrow(total2)){
  df<-total2$data[[i]] %>% as.matrix()
  df[is.na(df)]<-0
  fredom<-nrow(df)-1
  rownames(df)<-df[,1]
  ye<-total2$year[i]
  mat<-df
  mat<-mat[, -1] 
  m<-as.data.frame(mat) %>% mutate_if(is.character, as.numeric)
  ca<-CA(m, graph=F, ncp=4)
  rowcoords<-ca$row$coord %>% as.data.frame()
  scale<-c()
  for(j in 1:4){ #if you want to scale the coordinates
    scale[j]<-sqrt(ca$eig[j,1])
  }
  for(j in 1:4){
    rowcoords[,j]<-rowcoords[,j]/scale[j]
  }
  rowcoords$point<-rownames(rowcoords)
  rowcoords$w<-ca$call$marge.row
  colnames(rowcoords)<-c("a", "b", "c", "d", "point", "w")
  colcoords<-ca$col$coord %>% as.data.frame()
  for(j in 1:4){ #if you want to scale the coordinates
    colcoords[,j]<-colcoords[,j]/scale[j]
  }
  colcoords$point<-rownames(colcoords)
  colcoords$w<-ca$call$marge.col
  colnames(colcoords)<-c("a", "b", "c", "d", "point", "w")
  rowcontrib<-ca$row$contrib %>% as.data.frame()
  rowcontrib$point<-rownames(rowcontrib)
  colnames(rowcontrib)<-c("conta", "contb", "contc", "contd", "point")
  colcontrib<-ca$col$contrib %>% as.data.frame()
  colcontrib$point<-rownames(colcontrib)
  colnames(colcontrib)<-c("conta", "contb", "contc", "contd", "point")
  rows<-left_join(rowcoords, rowcontrib)
  cols<-left_join(colcoords, colcontrib)
  pval<-pchisq(sum(ca$eig[,1][1:4])*ca$call$N, fredom, lower.tail = F)
  dimensions<-data.frame(dim=1:4, eigenvalue=ca$eig[,1][1:4], pervariance=ca$eig[,2][1:4], pval=pval)
  rows$year<-ye
  cols$year<-ye
  dimensions$year<-ye
  tableU<-rbind(tableU, rows)
  tableV<-rbind(tableV, cols)
  tableL<-rbind(tableL, dimensions)
}
write_csv(tableL, "tableL.csv")
write_csv(tableU, "tableU.csv")
write_csv(tableV, "tableV.csv")
#####Distancias#####
tableU<-read_csv("tableU.csv")
tableV<-read_csv("tableV.csv")
tableU<- tableU %>% group_by(year) %>% nest()
for(i in 1:nrow(tableU)){
  for(j in 1:nrow(tableU$data[[i]])){
    dt1<-tableU$data[[i]][-j,] %>% subset(select=c(a, b, c, d, w))
    dt2<-tableU$data[[i]][j,] %>% subset(select=c(a, b, c, d, w))
    dt1$dist<-sqrt((dt2$a)^2+(dt2$b)^2+(dt2$c)^2+(dt2$d)^2)
    dt1$distw<-dt1$dist*dt1$w
    dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
    dt1$wdist<-dt1$udist*dt1$w
    tableU$data[[i]]$dist[j]<-mean(dt1$dist)
    tableU$data[[i]]$distw[j]<-mean(dt1$distw)
    tableU$data[[i]]$udist[j]<-mean(dt1$udist)
    tableU$data[[i]]$wdist[j]<-mean(dt1$wdist)
  }
}
tableU<-tableU %>% unnest(cols=data)
dumd<-tableU %>% group_by(year) %>% 
  summarise(dist=mean(dist),
            distw=mean(distw),
            udist=mean(udist), 
            wdist=mean(wdist), .groups="keep")

tableV<- tableV %>% group_by(year) %>% nest()
for(i in 1:nrow(tableV)){
  for(j in 1:nrow(tableV$data[[i]])){
    dt1<-tableV$data[[i]][-j,] %>% subset(select=c(a, b, c, d, w))
    dt2<-tableV$data[[i]][j,] %>% subset(select=c(a, b, c, d, w))
    dt1$dist<-sqrt((dt2$a)^2+(dt2$b)^2+(dt2$c)^2+(dt2$d)^2)
    dt1$distw<-dt1$dist*dt1$w
    dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
    dt1$wdist<-dt1$udist*dt1$w
    tableV$data[[i]]$dist[j]<-mean(dt1$dist)
    tableV$data[[i]]$distw[j]<-mean(dt1$distw)
    tableV$data[[i]]$udist[j]<-mean(dt1$udist)
    tableV$data[[i]]$wdist[j]<-mean(dt1$wdist)
  }
}
tableV<-tableV %>% unnest(cols=data)
vdumd<-tableV %>% group_by(year) %>% 
  summarise(dist=mean(dist),
            distw=mean(distw),
            udist=mean(udist), 
            wdist=mean(wdist), .groups="keep")
write_csv(tableU, "tableU.csv")
write_csv(tableV, "tableV.csv")
write_csv(dumd, "dumd.csv")
write_csv(vdumd, "vdumd.csv")

#####Clusters U#####
tableU<-read_csv("tableU.csv")
tableU<- tableU %>% group_by(year) %>% nest()
ucentroids<-data.frame()
for(i in 1:nrow(tableU)){
  dum<-tableU$data[[i]]
  posiciones<-dum %>% subset(select=c(a,b,c,d))
  kmax<-(nrow(dum)/2) %>%floor()
  set.seed(1917)
  optclust<-data.frame(
    clust=fviz_nbclust(posiciones, kmeans, method="silhouette", k.max=kmax)$data$clusters %>% as.numeric, 
    sil=fviz_nbclust(posiciones,  kmeans, method="silhouette", k.max=kmax, nboot=15000)$data$y)
  oks<-fviz_nbclust(posiciones,  kmeans, method="silhouette", k.max=kmax, nboot=15000)$layers[[3]]$data$xintercept %>% as.numeric()
  if(oks==1){
    optclust$sild<-ifelse(((optclust$sil>lag(optclust$sil))&(optclust$sil>lead(optclust$sil)))|
                            (optclust$clust==max(optclust$clust)), 1, 0)
    sils<-subset(optclust, sild==1, select=c(clust, sil))
    oks<-sils$clust[1]  
  }
  dclusters<-kcca(posiciones, k=oks, weights = sqrt(dum$w), control = list(iter=15000))
  centroids<-parameters(dclusters) %>% as.data.frame()
  centroids$clusters<-rownames(centroids)
  dum$clusters<-dclusters@cluster
  tableU$data[[i]]<-dum
  centroids$year<-tableU$year[i]
  ucentroids<-rbind(ucentroids, centroids)
  }
tableU<-tableU %>% unnest(cols=data)
write_csv(tableU, "tableU.csv")
write_csv(ucentroids, "ucentroids.csv")

tableV<-read_csv("tableV.csv")
tableV<- tableV %>% group_by(year) %>% nest()
vcentroids<-data.frame()
vcentroidg<-data.frame()
for(i in 1:nrow(tableV)){
  dum<-tableV$data[[i]]
  posiciones<-dum %>% subset(select=c(a,b,c,d))
  kmax<-(nrow(dum)/2) %>%floor()
  set.seed(1917)
  optclust<-data.frame(
    clust=fviz_nbclust(posiciones, kmeans, method="silhouette", k.max=kmax)$data$clusters %>% as.numeric, 
    sil=fviz_nbclust(posiciones,  kmeans, method="silhouette", k.max=kmax, nboot=15000)$data$y)
  oks<-fviz_nbclust(posiciones,  kmeans, method="silhouette", k.max=kmax, nboot=15000)$layers[[3]]$data$xintercept %>% as.numeric()
  if(oks==1){
    optclust$sild<-ifelse(((optclust$sil>lag(optclust$sil))&(optclust$sil>lead(optclust$sil)))|
                            (optclust$clust==max(optclust$clust)), 1, 0)
    sils<-subset(optclust, sild==1, select=c(clust, sil))
    oks<-sils$clust[1]  
  }
  dclusters<-kcca(posiciones, k=oks, weights = sqrt(dum$w), control = list(iter=15000))
  centroids<-parameters(dclusters) %>% as.data.frame()
  centroids$clusters<-rownames(centroids)
  dum$clusters<-dclusters@cluster
  centroids$year<-tableV$year[i]
  vcentroids<-rbind(vcentroids, centroids)
  tableV$data[[i]]<-dum
}
tableV<-tableV %>% unnest(cols=data)
write_csv(tableV, "tableV.csv")
write_csv(vcentroids, "vcentroids.csv")
#####Contributions#####
tableU<-read_csv("tableU.csv")
tableV<-read_csv("tableV.csv")
tableL<-read_csv("tableL.csv")
tableV<-tableV %>% group_by(year) %>% nest()
for(i in 1:nrow(tableV)){
  tableV$data[[i]]$sumcont<-(tableV$data[[i]]$conta+tableV$data[[i]]$contb+tableV$data[[i]]$contc+tableV$data[[i]]$contd)/4
  tableV$data[[i]]$contdiso<-100*tableV$data[[i]]$w*(tableV$data[[i]]$dist^2)/sum(tableV$data[[i]]$w*(tableV$data[[i]]$dist^2))
  tableV$data[[i]]$contdisp<-100*tableV$data[[i]]$w*(tableV$data[[i]]$udist^2)/sum(tableV$data[[i]]$w*(tableV$data[[i]]$udist^2))
}
tableU<-tableU %>% group_by(year) %>% nest()
for(i in 1:nrow(tableU)){
  tableU$data[[i]]$sumcont<-(tableU$data[[i]]$conta+tableU$data[[i]]$contb+tableU$data[[i]]$contc+tableU$data[[i]]$contd)/4
  tableU$data[[i]]$contdiso<-100*tableU$data[[i]]$w*(tableU$data[[i]]$dist^2)/sum(tableU$data[[i]]$w*(tableU$data[[i]]$dist^2))
  tableU$data[[i]]$contdisp<-100*tableU$data[[i]]$w*(tableU$data[[i]]$udist^2)/sum(tableU$data[[i]]$w*(tableU$data[[i]]$udist^2))
}

tableV<-tableV %>% unnest(cols=data)
tableU<-tableU %>% unnest(cols=data)
meancontV<-tableV %>% group_by(point) %>% summarize(meancont=mean(sumcont),
                                                   contdiso=mean(contdiso), 
                                                 contdisp=mean(contdisp),.groups="keep")
meancontU<-tableU %>% group_by(point) %>% summarize(meancont=mean(sumcont),
                                                    contdiso=mean(contdiso), 
                                                    contdisp=mean(contdisp),.groups="keep")
write_csv(meancontV, "meancontV.csv")
write_csv(meancontU, "meancontU.csv")
####Narrativa####
#####estadistico#####
tableU<-read_csv("tableU.csv")
tableU$id<-tableU$point
tableU$year<-tableU$year %>% as.character() %>% as.numeric()
tableU$tipo<-ifelse(str_detect(tableU$point, "enep"), "ENEP", 
                    ifelse(str_detect(tableU$point, "^e"), "Escuela", 
                           ifelse(str_detect(tableU$point, "medicina|ingenier|derech|arqu"), "Fac Hist", 
                                  "Facultad")))
narrative1<-narrative.group(tableU, tableU$clusters, lim=2000, 
                            sub=tableU$tipo, phase=tableU$year,
                            type.init="distribution", type.calc="cross", sec.calc="level",
                            parallel=T, extraspace=0.2)
write_csv(narrative1, "narrative1.csv")

####Centroides####
######generales######
vcentroids<-read_csv("vcentroids.csv")
ucentroids<-read_csv("ucentroids.csv")
vcentroids<-vcentroids %>% group_by(year)%>% nest()
ucentroids<-ucentroids %>% group_by(year) %>% nest()
for(i in 1:nrow(vcentroids)){
  dum<-vcentroids$data[[i]]
  dum$maxangle<-NA
  dum$maxclus<-NA
  dum$minclus<-NA
  dum$minangle<-NA
  for(j in 1:nrow(dum)){
    angles<-data.frame()
    for(k in 1:nrow(ucentroids$data[[i]])){
      clus<-ucentroids$data[[i]]$clusters[k] %>% as.numeric()
      coord<-ucentroids$data[[i]][k,] %>% subset(select=c(a, b, c, d)) %>% as.numeric()
      angulo<-angle(dum[j,] %>% subset(select=c(a,b,c,d)) %>% as.numeric(), coord)
      angulos<-data.frame(clust=clus, angle=angulo)
      angles<-rbind(angles, angulos)
    }
    dum$minangle[j]<-min(angles$angle)
    dum$minclus[j]<-angles$clust[which.min(angles$angle)]
    dum$maxangle[j]<-max(angles$angle)
    dum$maxclus[j]<-angles$clust[which.max(angles$angle)]
  }
  vcentroids$data[[i]]<-dum
}
for(i in 1:nrow(ucentroids)){
  dum<-ucentroids$data[[i]]
  dum$maxangle<-NA
  dum$maxclus<-NA
  dum$minclus<-NA
  dum$minangle<-NA
  for(j in 1:nrow(dum)){
    angles<-data.frame()
    for(k in 1:nrow(vcentroids$data[[i]])){
      clus<-vcentroids$data[[i]]$clusters[k] %>% as.numeric()
      coord<-vcentroids$data[[i]][k,] %>% subset(select=c(a, b, c, d)) %>% as.numeric()
      angulo<-angle(dum[j,] %>% subset(select=c(a,b,c,d)) %>% as.numeric(), coord)
      angulos<-data.frame(clust=clus, angle=angulo)
      angles<-rbind(angles, angulos)
    }
    dum$minangle[j]<-min(angles$angle)
    dum$minclus[j]<-angles$clust[which.min(angles$angle)]
    dum$maxangle[j]<-max(angles$angle)
    dum$maxclus[j]<-angles$clust[which.max(angles$angle)]
  }
  ucentroids$data[[i]]<-dum
}
vcentroids<-vcentroids %>% unnest(cols=data)
ucentroids<-ucentroids %>% unnest(cols=data)
write_csv(vcentroids, "vcentroids.csv")
write_csv(ucentroids, "ucentroids.csv")

#####angulos#####
tableV<-read_csv("tableV.csv")
tableU<-read_csv("tableU.csv")
ucentroids<-read_csv("ucentroids.csv")
vcentroids<-read_csv("vcentroids.csv")
ucentroids<-ucentroids %>% subset(select=c(year, clusters, maxclus, maxangle, minclus, minangle))
vcentroids<-vcentroids %>% subset(select=c(year, clusters, maxclus, maxangle, minclus, minangle))
tableU<-tableU %>% left_join(ucentroids, by=c("year", "clusters"))
tableV<-tableV %>% left_join(vcentroids, by=c("year", "clusters"))
tableU<-tableU %>% group_by(year) %>% nest()
tableV<-tableV %>% group_by(year) %>% nest()
tableU<-tableU %>% unnest()
tableV<-tableV %>% unnest()
write_csv(tableU, "tableUangle.csv")
write_csv(tableV, "tableVangle.csv")

####Estabilidad#####
tableU<-read_csv("tableU.csv")
tableL<-read_csv("tableL.csv")
tablad<-tableL %>% group_by(year) %>% summarize(pervariance=sum(pervariance), .groups="keep")
dumd<-tableU %>% group_by(year) %>% 
  summarise(dist=mean(dist), distw=mean(distw),
            udist=mean(udist), wdist=mean(wdist),.groups="keep")
tablad$dist<-dumd$dist
tablad$distw<-dumd$distw
tablad$udist<-dumd$udist
tablad$wdist<-dumd$wdist
mediast<-data.frame(pervar=mean(tablad$pervariance), dist=mean(tablad$dist), distw=mean(tablad$distw), 
                    udist=mean(tablad$udist), wdist=mean(tablad$wdist))
facultades<-read_csv("facultades.csv")
total2<-facultades
######Con columnas#####
set.seed(1917)
tabladf<-data.frame()
for(k in 3:ncol(total2)){
  totaltemp<-total2[-k]
  totaltemp2<-totaltemp %>% group_by(year) %>% nest()
  tableU<-data.frame()
  tableL<-data.frame()
  for (i in 1:nrow(totaltemp2)){
    df<-totaltemp2$data[[i]] %>% as.matrix()
    df[is.na(df)]<-0
    fredom<-nrow(df)-1
    rownames(df)<-df[,1]
    ye<-total2$year[i]
    mat<-df
    mat<-mat[, -1] 
    m<-as.data.frame(mat) %>% mutate_if(is.character, as.numeric)
    ca<-CA(m, graph=F, ncp=4)
    rowcoords<-ca$row$coord %>% as.data.frame()
    scale<-c()
    for(j in 1:4){ #if you want to scale the coordinates
      scale[j]<-sqrt(ca$eig[j,1])
    }
    for(j in 1:4){
      rowcoords[,j]<-rowcoords[,j]/scale[j]
    }
    rowcoords$point<-rownames(rowcoords)
    colnames(rowcoords)<-c("a", "b", "c", "d", "point")
    rowcoords$w<-sqrt(ca$call$marge.row)
    pval<-pchisq(sum(ca$eig[,1][1:4])*ca$call$N, fredom, lower.tail = F)
    dimensions<-data.frame(dim=1:4, eigenvalue=ca$eig[,1][1:4], pervariance=ca$eig[,2][1:4], pval=pval)
    rowcoords$year<-ye
    dimensions$year<-ye
    tableU<-rbind(tableU, rowcoords)
    tableL<-rbind(tableL, dimensions)
  }
  tableU<- tableU %>% group_by(year) %>% nest()
  for(i in 1:nrow(tableU)){
    for(j in 1:nrow(tableU$data[[i]])){
      dt1<-tableU$data[[i]][-j,] %>% subset(select=c(a, b, c, d, w))
      dt2<-tableU$data[[i]][j,] %>% subset(select=c(a, b, c, d, w))
      dt1$dist<-sqrt((dt2$a)^2+(dt2$b)^2+(dt2$c)^2+(dt2$d)^2)
      dt1$distw<-dt1$dist*dt1$w
      dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
      dt1$wdist<-dt1$udist*dt1$w
      tableU$data[[i]]$dist[j]<-mean(dt1$dist)
      tableU$data[[i]]$distw[j]<-mean(dt1$distw)
      tableU$data[[i]]$udist[j]<-mean(dt1$udist)
      tableU$data[[i]]$wdist[j]<-mean(dt1$wdist)
    }
  }
  tableU<-tableU %>% unnest(cols=data)
  dumd<-tableU %>% group_by(year) %>% 
    summarise(dist=mean(dist),
              distw=mean(distw),
              udist=mean(udist), 
              wdist=mean(wdist), .groups="keep") %>%
    summarise(dist=mean(dist),
              distw=mean(distw),
              udist=mean(udist), 
              wdist=mean(wdist))
  dumv<-tableL %>% group_by(year) %>% 
    summarise(pervariance=sum(pervariance), .groups="keep")
  tablad<-data.frame(colus=colnames(total2)[k], 
                     pervar=mean(dumv$pervariance), 
                     dist=mean(tableU$dist), 
                     distw=mean(tableU$distw),
                     udist=mean(tableU$udist), 
                     wdist=mean(tableU$wdist))
  tabladf<-rbind(tabladf, tablad)
}
tablaparam<-tabladf
tablaparam$infdist<-abs(tablaparam$dist-mediast$dist)
tablaparam$infdistw<-abs(tablaparam$distw-mediast$distw)
tablaparam$infudist<-abs(tablaparam$udist-mediast$udist)
tablaparam$infwdist<-abs(tablaparam$wdist-mediast$wdist)
tablaparam$infvar<-abs(tablaparam$pervar-mediast$pervar)
tablaparam$infvarp<-100*tablaparam$infvar/abs(mediast$pervar)
tablaparam$infdisp<-100*tablaparam$infdist/mediast$dist
tablaparam$infdiswp<-100*tablaparam$infdistw/mediast$distw
tablaparam$infudisp<-100*tablaparam$infudist/mediast$udist
tablaparam$infwdisp<-100*tablaparam$infwdist/mediast$wdist
tablainfl1<-tablaparam %>% 
  select(colus, infdist, infdistw, infudist, infwdist, infvar, infdisp, infdiswp,infudisp, infwdisp, infvarp)
write_csv(tablainfl1, "tablainfl1.csv")
######Con entidades######
set.seed(1917)
tabladf<-data.frame()
entidades<-unique(total2$general)
for(k in 1:length(entidades)){
  totaltemp<-total2 %>% subset(general!=entidades[k])
  totaltemp2<-totaltemp %>% group_by(year) %>% nest()
  tableU<-data.frame()
  tableL<-data.frame()
  for (i in 1:nrow(totaltemp2)){
    df<-totaltemp2$data[[i]] %>% as.matrix()
    df[is.na(df)]<-0
    fredom<-nrow(df)-1
    rownames(df)<-df[,1]
    ye<-total2$year[i]
    mat<-df
    mat<-mat[, -1] 
    m<-as.data.frame(mat) %>% mutate_if(is.character, as.numeric)
    ca<-CA(m, graph=F, ncp=4)
    rowcoords<-ca$row$coord %>% as.data.frame()
    scale<-c()
    for(j in 1:4){ #if you want to scale the coordinates
      scale[j]<-sqrt(ca$eig[j,1])
    }
    for(j in 1:4){
      rowcoords[,j]<-rowcoords[,j]/scale[j]
    }
    rowcoords$point<-rownames(rowcoords)
    colnames(rowcoords)<-c("a", "b", "c", "d", "point")
    rowcoords$w<-sqrt(ca$call$marge.row)
    pval<-pchisq(sum(ca$eig[,1][1:4])*ca$call$N, fredom, lower.tail = F)
    dimensions<-data.frame(dim=1:4, eigenvalue=ca$eig[,1][1:4], pervariance=ca$eig[,2][1:4], pval=pval)
    rowcoords$year<-ye
    dimensions$year<-ye
    tableU<-rbind(tableU, rowcoords)
    tableL<-rbind(tableL, dimensions)
  }
  tableU<- tableU %>% group_by(year) %>% nest()
  for(i in 1:nrow(tableU)){
    for(j in 1:nrow(tableU$data[[i]])){
      dt1<-tableU$data[[i]][-j,] %>% subset(select=c(a, b, c, d, w))
      dt2<-tableU$data[[i]][j,] %>% subset(select=c(a, b, c, d, w))
      dt1$dist<-sqrt((dt2$a)^2+(dt2$b)^2+(dt2$c)^2+(dt2$d)^2)
      dt1$distw<-dt1$dist*dt1$w
      dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
      dt1$wdist<-dt1$udist*dt1$w
      tableU$data[[i]]$dist[j]<-mean(dt1$dist)
      tableU$data[[i]]$distw[j]<-mean(dt1$distw)
      tableU$data[[i]]$udist[j]<-mean(dt1$udist)
      tableU$data[[i]]$wdist[j]<-mean(dt1$wdist)
    }
  }
  tableU<-tableU %>% unnest(cols=data)
  dumd<-tableU %>% group_by(year) %>% 
    summarise(dist=mean(dist),
              distw=mean(distw),
              udist=mean(udist), 
              wdist=mean(wdist), .groups="keep") %>%
    summarise(dist=mean(dist),
              distw=mean(distw),
              udist=mean(udist), 
              wdist=mean(wdist))
  dumv<-tableL %>% group_by(year) %>% 
    summarise(pervariance=sum(pervariance), .groups="keep")
  tablad<-data.frame(entidades=entidades[k], 
                     pervar=mean(dumv$pervariance), 
                     dist=mean(tableU$dist), 
                     distw=mean(tableU$distw),
                     udist=mean(tableU$udist), 
                     wdist=mean(tableU$wdist))
  tabladf<-rbind(tabladf, tablad)
}
tablaparam<-tabladf
tablaparam$infudist<-abs(tablaparam$udist-mediast$udist)
tablaparam$infwdist<-abs(tablaparam$wdist-mediast$wdist)
tablaparam$infdist<-abs(tablaparam$dist-mediast$dist)
tablaparam$infdistw<-abs(tablaparam$distw-mediast$distw)
tablaparam$infvar<-abs(tablaparam$pervar-mediast$pervar)
tablaparam$infvarp<-100*tablaparam$infvar/abs(mediast$pervar)
tablaparam$infdisp<-100*tablaparam$infdist/mediast$dist
tablaparam$infdiswp<-100*tablaparam$infdistw/mediast$distw
tablaparam$infudisp<-100*tablaparam$infudist/mediast$udist
tablaparam$infwdisp<-100*tablaparam$infwdist/mediast$wdist
tablainfl2<-tablaparam %>% 
  select(entidades, infdist, infudist, infvar, infdisp, infdiswp,infudisp, infwdisp, infvarp)
write_csv(tablainfl2, "tablainfl2.csv")

####Gráficas####
#####narrative#####
narrative1<-read_csv("narrative1.csv")
narrative1$tipo<-ifelse(str_detect(narrative1$point, "enep"), "ENEP", 
                    ifelse(str_detect(narrative1$point, "^e"), "Escuela", 
                           ifelse(str_detect(narrative1$point, "medicina|ingenier|derech|arqu"), "Fac Hist", 
                                  "Facultad")))
narrative1<-correct.group(narrative1, narrative1$group, 100, type.calc = "cross", direction="forward")
clusters<-narrplot.group(df=narrative1, group=narrative1$clusters, 
                         phase=narrative1$year, pos=narrative1$pos, color=narrative1$point)+
  ggtitle("Clusters using Silhouette")
htmlclus<-ggplotly(clusters)
saveWidget(htmlclus, "clusters.html")
clusters2<-narrplot.group(narrative1, narrative1$clusters, narrative1$year, 
                          narrative1$pos, narrative1$tipo)+
  labs(title="Clusters using Silhouette")
jpeg("clusters2.jpg", width=1200, height=600, units="px")
clusters2
dev.off()

#####numero de clusters U#####
tableU<-read_csv("tableU.csv")
escuelas<-tableU %>% count(year)
gruposil<-tableU %>% group_by(year, clusters) %>% count()
maxsil<-gruposil %>% group_by(year) %>% summarise(maxs=max(n), means=mean(n), groups=length(n), .groups="keep")
escuelas<-escuelas %>% left_join(maxsil, by="year")
escuelas$mxms<-100*escuelas$maxs/escuelas$means
escuelas$maxsp<-100*escuelas$maxs/escuelas$n
escuelas$meansp<-100*escuelas$means/escuelas$n
escuelas2<-data.frame(year=seq(1959, 1979, 1)) %>% 
  left_join(subset(escuelas, select=c(year, mxms)), by="year")
escuelas2$mxms<-llenararith2(escuelas2$mxms)

rm(ci, br)
porcentajes<-ts(escuelas2$mxms, start = 1959)
br<-breakpoints(porcentajes~1)
ci <- confint(br)
jpeg("Ubreaksilhouette.jpg", width=800, height=600, units="px")
plot(porcentajes, las=2,xaxt="n",
     ylab="Max(n)/mean(n) (%)", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()

tableV<-read_csv("tableV.csv")
variables<-tableV %>% count(year)
gruposil<-tableV %>% group_by(year, clusters) %>% count()
maxsil<-gruposil %>% group_by(year) %>% summarise(maxs=max(n), means=mean(n), groups=length(n), .groups="keep")
variables<-variables %>% left_join(maxsil, by="year")
variables$mxms<-100*variables$maxs/variables$means
variables$maxsp<-100*variables$maxs/variables$n
variables$meansp<-100*variables$means/variables$n
variables2<-data.frame(year=seq(1959, 1979, 1)) %>% 
  left_join(subset(variables, select=c(year, mxms)), by="year")
variables2$mxms<-llenararith2(variables2$mxms)

rm(ci, br)
porcentajes<-ts(variables2$mxms, start = 1959)
br<-breakpoints(porcentajes~1)
jpeg("Vbreaksilhouette.jpg", width=800, height=600, units="px")
plot(porcentajes, las=2,xaxt="n",main="Concentración relativa clusters variables",sub="método silhouette", 
     ylab="Max(n)/mean(n) (%)", 
     xlab="Year")
lines(br)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()

#####Distancias#####
tableU<-read_csv("tableU.csv")
dumd<-read_csv("dumd.csv")
tableV<-read_csv("tableV.csv")
vdumd<-read_csv("vdumd.csv")
######Distancia media del origen, no ponderada######
dumd2<-data.frame(year=seq(1959, 1979, 1)) %>% 
  left_join(subset(dumd, select=c(year, dist, distw, udist, wdist)), by="year")
dumd2$dist<-llenararith2(dumd2$dist)
dumd2$distw<-llenararith2(dumd2$distw)
dumd2$udist<-llenararith2(dumd2$udist)
dumd2$wdist<-llenararith2(dumd2$wdist)
rm(ci, br)
distancias<-ts(dumd2$dist, start = 1959)
br<-breakpoints(distancias~1)
ci <- confint(br)
jpeg("breakdist.jpg", width=800, height=600, units="px")
plot(distancias, type="point", xaxt="n",ylab="Distancia del origen, no ponderada", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()
######Distancia del origen, ponderada######
rm(ci, br)
distancias<-ts(dumd2$distw, start = 1959)
br<-breakpoints(distancias~1)
ci <- confint(br)
jpeg("breakdistw.jpg", width=800, height=600, units="px")
plot(distancias, type="point", xaxt="n",ylab="Distancia del origen ponderada", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()
######Distancia media no ponderada######
rm(ci, br)
cdistancias<-ts(dumd2$udist, start = 1959)
br<-breakpoints(distancias~1)
ci <- confint(br)
jpeg("breakudist.jpg", width=800, height=600, units="px")
plot(distancias, type="point", las=2,xaxt="n",ylab="Distancia media de pares, no ponderada", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()

######Distancia media ponderada######
rm(ci, br)
distancias<-ts(dumd2$wdist, start = 1959)
br<-breakpoints(distancias~1, h=0.2)
ci <- confint(br)
jpeg("breakwdist.jpg", width=800, height=600, units="px")
plot(distancias, type="point", las=2,xaxt="n",ylab="Distancia media de pares, ponderada", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()
######Graficas######
plotd<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=dist, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=dist, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=dist), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=dist), color="black", se=F)+
  theme_minimal()+
  labs(x="Year", y="Distancia del origen, no ponderada")
plotdw<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=distw, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=distw, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=distw), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=distw), color="black", se=F)+
  theme_minimal()+
  labs(x="Year", y="Distancia del origen, ponderada")
plotwd<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=wdist, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=wdist, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=wdist), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=wdist), color="black", se=F)+
  theme_minimal()+
  labs(x="Year", y="Distancia media de pares, ponderada")
plotud<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=udist, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=udist, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=udist), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=udist), color="black", se=F)+
  theme_minimal()+
  labs(x="Year", y="Distancia media de pares, no ponderada")
saveWidget(ggplotly(plotd), "distance2.html")
saveWidget(ggplotly(plotdw), "distance2w.html")
saveWidget(ggplotly(plotud), "udistance2.html")
saveWidget(ggplotly(plotwd), "wdistance2.html")
pdf("plotdist.pdf")
plotd+theme(legend.position="none")
dev.off()
pdf("plotdistw.pdf")
plotdw+theme(legend.position="none")
dev.off()
pdf("plotudist.pdf")
plotud+theme(legend.position="none")
dev.off()
pdf("plotwdist.pdf")
plotwd+theme(legend.position="none")
dev.off()

distancias<-grid.arrange(plotd+theme(legend.position="none", title.position="none"), 
             plotdw+theme(legend.position="none", title.position="none"), 
             plotud+theme(legend.position="none", title.position="none"), 
             plotwd+theme(legend.position="none", title.position="none"), ncol=2)
jpeg("griddistanciam.jpg", width=800, height=1000, units="px")
distancias
dev.off()

library(jpeg)
plot1<-readJPEG("breakdist.jpg")
plot2<-readJPEG("breakdistw.jpg")
plot3<-readJPEG("breakudist.jpg")
plot4<-readJPEG("breakwdist.jpg")
plotfinal<-grid.arrange(rasterGrob(plot1), rasterGrob(plot3), rasterGrob(plot2), rasterGrob(plot4), ncol=2)
jpeg("griddistancias.jpg", width=800, height=1000, units="px")
plotfinal
dev.off()

#####tablas#####
meancontV<-read_csv("meancontV.csv")
tabinflV<-read_csv("tablainfl1.csv")
(mean( tabinflV$infdisp)+ mean(tabinflV$infdiswp)+ mean(tabinflV$infudisp)+ mean(tabinflV$infwdisp))/4
tabinflV$mean<-rowMeans(tabinflV %>% select(c(infdisp, infudisp)))
tabinflV<-tabinflV[order(tabinflV$mean, decreasing=T),]
tabinflV<-tabinflV[1:30,]
tabinflV<-subset(tabinflV, select=c(colus, infdisp, infudisp))
tablaiV<-left_join(tabinflV, meancontV, by=c("colus"="point"))
write_csv(tablaiV, "tablaiV.csv")

meancontU<-read_csv("meancontU.csv")
tabinflU<-read_csv("tablainfl2.csv")
(mean( tabinflU$infdisp)+ mean(tabinflU$infdiswp)+ mean(tabinflU$infudisp)+ mean(tabinflU$infwdisp))/4
tabinflU$mean<-rowMeans(tabinflU %>% select(c(infdisp, infudisp)))
tabinflU<-tabinflU[order(tabinflU$mean, decreasing=T),]
tabinflU<-subset(tabinflU, select=c(entidades, infdisp, infudisp))
tablaiU<-left_join(tabinflU, meancontU, by=c("entidades"="point"))
write_csv(tablaiU, "tablaiU.csv")

tableL<-read_csv("tableL.csv")
timedim<-ggplot(tableL, aes(x=year, y=pervariance, fill=dim))+
  geom_col(position = "stack")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="Percentage of variance", x="Year", y="Percentage of variance")
timedim
jpeg("plotdimensionstime.jpg", width=800, height=600, units="px")
print(timedim)
dev.off()

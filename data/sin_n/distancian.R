install.packages("flexclust")
install.packages("geometry")
install.packages("ggbump")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
install.packages("../../../../narrative/diagram/narrativecharts_0.1.0/narrativecharts_0.1.0.tar.gz", type="source")

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
facultades<-facultades %>% subset(select=colnames(facultades)[which(!str_detect(colnames(facultades), "num_|cap_"))])
total2<-facultades %>% group_by(year) %>% nest()
####Analisis ####
#####Tablas generales#####
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
  rowcoords$point<-rownames(rowcoords)
  rowcoords$w<-ca$call$marge.row
  colnames(rowcoords)<-c("a", "b", "c", "d", "point", "w")
  colcoords<-ca$col$coord %>% as.data.frame()
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
ggplot(tableL, aes(x=year, y=pervariance, fill=dim))+geom_col(position="stack")
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
    dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
    dt1$wdist<-dt1$udist*dt1$w
    tableU$data[[i]]$dist[j]<-mean(dt1$dist)
    tableU$data[[i]]$udist[j]<-mean(dt1$udist)
    tableU$data[[i]]$wdist[j]<-mean(dt1$wdist)
  }
}
tableU<-tableU %>% unnest(cols=data)
dumd<-tableU %>% group_by(year) %>% 
  summarise(dist=mean(dist),
            udist=mean(udist), 
            wdist=mean(wdist), .groups="keep")

tableV<- tableV %>% group_by(year) %>% nest()
for(i in 1:nrow(tableV)){
  for(j in 1:nrow(tableV$data[[i]])){
    dt1<-tableV$data[[i]][-j,] %>% subset(select=c(a, b, c, d, w))
    dt2<-tableV$data[[i]][j,] %>% subset(select=c(a, b, c, d, w))
    dt1$dist<-sqrt((dt2$a)^2+(dt2$b)^2+(dt2$c)^2+(dt2$d)^2)
    dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
    dt1$wdist<-dt1$udist*dt1$w
    tableV$data[[i]]$dist[j]<-mean(dt1$dist)
    tableV$data[[i]]$udist[j]<-mean(dt1$udist)
    tableV$data[[i]]$wdist[j]<-mean(dt1$wdist)
  }
}
tableV<-tableV %>% unnest(cols=data)
vdumd<-tableV %>% group_by(year) %>% 
  summarise(dist=mean(dist),
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
ucentroidg<-data.frame()
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
    optclust$gapd<-ifelse(((optclust$gap>lag(optclust$gap))&(optclust$gap>lead(optclust$gap)))|
                            (optclust$clust==max(optclust$clust)), 1, 0)
    gaps<-subset(optclust, gapd==1, select=c(clust, gap))
    oks<-gaps$clust[1]  
  }
  dclusters<-kcca(posiciones, k=oks, weights = dum$w, control = list(iter=15000))
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
write_csv(ucentroidg, "ucentroidg.csv")

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
    optclust$gapd<-ifelse(((optclust$gap>lag(optclust$gap))&(optclust$gap>lead(optclust$gap)))|
                            (optclust$clust==max(optclust$clust)), 1, 0)
    gaps<-subset(optclust, gapd==1, select=c(clust, gap))
    oks<-gaps$clust[1]  
  }
  dclusters<-kcca(posiciones, k=oks, weights = dum$w, control = list(iter=15000))
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
write_csv(vcentroidg, "vcentroidg.csv")
#####Contributions#####
tableV<-read_csv("tableV.csv")
contribs<-data.frame(year=unique(tableV$year), conta=NA, contb=NA, contc=NA, contd=NA)
tableV<-tableV %>% group_by(year) %>% nest()
for(i in 1:nrow(tableV)){
  contribs$conta[i]<-tableL$data[[i]]$point[which(tableL$data[[i]]$conta==max(tableL$data[[i]]$conta))]
  contribs$contb[i]<-tableL$data[[i]]$point[which(tableL$data[[i]]$contb==max(tableL$data[[i]]$contb))]
  contribs$contc[i]<-tableL$data[[i]]$point[which(tableL$data[[i]]$contc==max(tableL$data[[i]]$contc))]
  contribs$contd[i]<-tableL$data[[i]]$point[which(tableL$data[[i]]$contd==max(tableL$data[[i]]$contd))]
}

####Narrativa####
#####random#####
tableU<-read_csv("tableU.csv")
tableU$id<-tableU$point
tableU$year<-tableU$year %>% as.character() %>% as.numeric()
narrative1<-narrative.group(tableU, tableU$clusters, lim=200000, phase=tableU$year,type.init="distribution", type.calc="level", parallel = T)
narrative<-correct.group(narrative1, narrative1$clusters, 1000000, type.calc="level",direction = "forward", parallel = T)
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
tableU<-tableU %>% unnest(cols = c(data))
tableV<-tableV %>% unnest(cols = c(data))
write_csv(tableU, "tableUangle.csv")
write_csv(tableV, "tableVangle.csv")

####Estabilidad#####
tableU<-read_csv("tableU.csv")
tableL<-read_csv("tableL.csv")
tablad<-tableL %>% group_by(year) %>% summarize(pervariance=sum(pervariance), .groups="keep")
dumd<-tableU %>% group_by(year) %>% 
  summarise(dist=mean(dist), udist=mean(udist), .groups="keep")
tablad$dist<-dumd$dist
tablad$udist<-dumd$udist
mediast<-data.frame(pervar=mean(tablad$pervariance), dist=mean(tablad$dist), udist=mean(tablad$udist))
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
    rowcoords$point<-rownames(rowcoords)
    colnames(rowcoords)<-c("a", "b", "c", "d", "point")
    rowcoords$w<-ca$call$marge.row
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
      dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
      dt1$wdist<-dt1$udist*dt1$w
      tableU$data[[i]]$dist[j]<-mean(dt1$dist)
      tableU$data[[i]]$udist[j]<-mean(dt1$udist)
      tableU$data[[i]]$wdist[j]<-mean(dt1$wdist)
    }
  }
  tableU<-tableU %>% unnest(cols=data)
  dumd<-tableU %>% group_by(year) %>% 
    summarise(dist=mean(dist),
              udist=mean(udist), 
              wdist=mean(wdist), .groups="keep") %>%
    summarise(dist=mean(dist),
              udist=mean(udist), 
              wdist=mean(wdist))
  dumv<-tableL %>% group_by(year) %>% 
    summarise(pervariance=sum(pervariance), .groups="keep")
  tablad<-data.frame(colus=colnames(total2)[k], 
                     pervar=mean(dumv$pervariance), 
                     dist=mean(tableU$dist), 
                     udist=mean(tableU$udist))
  tabladf<-rbind(tabladf, tablad)
}
tablaparam<-tabladf
tablaparam$infudist<-abs(tablaparam$udist-mediast$udist)
tablaparam$infdist<-abs(tablaparam$dist-mediast$dist)
tablaparam$infvar<-abs(tablaparam$pervar-mediast$pervar)
tablaparam$infvarp<-100*tablaparam$infvar/abs(mediast$pervar)
tablaparam$infdisp<-100*tablaparam$infdist/mediast$dist
tablaparam$infudisp<-100*tablaparam$infudist/mediast$udist
tablainfl1<-tablaparam %>% 
  select(colus, infdist, infudist, infvar, infdisp, infudisp, infvarp)
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
    rowcoords$point<-rownames(rowcoords)
    colnames(rowcoords)<-c("a", "b", "c", "d", "point")
    rowcoords$w<-ca$call$marge.row
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
      dt1$udist<-sqrt((dt1$a-dt2$a)^2+(dt1$b-dt2$b)^2+(dt1$c-dt2$c)^2+(dt1$d-dt2$d)^2)
      dt1$wdist<-dt1$udist*dt1$w
      tableU$data[[i]]$dist[j]<-mean(dt1$dist)
      tableU$data[[i]]$udist[j]<-mean(dt1$udist)
      tableU$data[[i]]$wdist[j]<-mean(dt1$wdist)
    }
  }
  tableU<-tableU %>% unnest(cols=data)
  dumd<-tableU %>% group_by(year) %>% 
    summarise(dist=mean(dist),
              udist=mean(udist), 
              wdist=mean(wdist), .groups="keep") %>%
    summarise(dist=mean(dist),
              udist=mean(udist), 
              wdist=mean(wdist))
  dumv<-tableL %>% group_by(year) %>% 
    summarise(pervariance=sum(pervariance), .groups="keep")
  tablad<-data.frame(entidades=entidades[k], 
                     pervar=mean(dumv$pervariance), 
                     dist=mean(tableU$dist), 
                     udist=mean(tableU$udist))
  tabladf<-rbind(tabladf, tablad)
}
tablaparam<-tabladf
tablaparam$infudist<-abs(tablaparam$udist-mediast$udist)
tablaparam$infdist<-abs(tablaparam$dist-mediast$dist)
tablaparam$infvar<-abs(tablaparam$pervar-mediast$pervar)
tablaparam$infvarp<-100*tablaparam$infvar/abs(mediast$pervar)
tablaparam$infdisp<-100*tablaparam$infdist/mediast$dist
tablaparam$infudisp<-100*tablaparam$infudist/mediast$udist
tablainfl2<-tablaparam %>% 
  select(entidades, infdist, infudist, infvar, infdisp, infudisp, infvarp)
write_csv(tablainfl2, "tablainfl2.csv")

####Gráficas####
#####narrative#####
narrative1<-read_csv("narrative1.csv")
narrative1<-correct.group(narrative1, narrative1$group, 1000, type.calc = "level")
narrative0<-read_csv("narrative2.csv")
clusters<-narrplot.group(narrative1, narrative1$clusters, narrative1$year, narrative1$pos, narrative1$point)+
  ggtitle("Clusters using Silhouette")
htmlclus<-ggplotly(clusters)
saveWidget(htmlclus, "clusters.html")
pdf("plotclusters.pdf")
clusters+theme(legend.position="none")
dev.off()


narrative2<-read_csv("narrative2.csv")
clusters<-narrplot.group(narrative1, narrative1$clusters, narrative1$year, narrative1$pos, narrative1$point)+
  ggtitle("Clusters using Silhouette")
clusters+theme(legend.position="none")

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
  left_join(subset(escuelas, select=c(year, mxms)), by="year") %>%
  fill(mxms)

porcentajes<-ts(escuelas2$mxms, start = 1959)
br<-breakpoints(porcentajes~1)
ci <- confint(br)
jpeg("Ubreaksilhouette.jpg", width=800, height=600, units="px")
plot(porcentajes, type="line", las=2,xaxt="n",main="Concentración relativa entidades" ,sub="método silhouette", 
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
  left_join(subset(variables, select=c(year, mxms)), by="year") %>%
  fill(mxms)

rm(ci, br, porcentajes)
porcentajes<-ts(variables2$mxms, start = 1959)
br<-breakpoints(porcentajes~1)
#ci <- confint(br)
jpeg("Vbreaksilhouette.jpg", width=800, height=600, units="px")
plot(porcentajes, type="line", las=2,xaxt="n",main="Concentración relativa clusters variables",sub="método silhouette", 
     ylab="Max(n)/mean(n) (%)", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()

#####Distancias#####
tableU<-read_csv("tableU.csv")
dumd<-read_csv("dumd.csv")
tableV<-read_csv("tableV.csv")
vdumd<-read_csv("vdumd.csv")
######Distancia media del origen######
dumd2<-data.frame(year=seq(1959, 1979, 1)) %>% 
  left_join(subset(dumd, select=c(year, dist, udist, wdist)), by="year") %>%
  fill(dist, udist, wdist)
rm(ci, br)
distancias<-ts(dumd2$dist, start = 1959)
br<-breakpoints(distancias~1)
ci <- confint(br)
jpeg("breakdist.jpg", width=800, height=600, units="px")
plot(distancias, type="line", xaxt="n",main="Distancia media del origen", ylab="Distancia", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()
######Distancia media no ponderada######
cdistancias<-ts(dumd2$udist, start = 1959)
br<-breakpoints(distancias~1)
ci <- confint(br)
jpeg("breakudist.jpg", width=800, height=600, units="px")
plot(distancias, type="line", las=2,xaxt="n",main="Distancia media de pares", ylab="Distancia no ponderada", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()

######Distancia media ponderada######
rm(ci, br, porcentajes)
distancias<-ts(dumd2$wdist, start = 1959)
br<-breakpoints(distancias~1, h=0.2)
ci <- confint(br)
jpeg("breakwdist.jpg", width=800, height=600, units="px")
plot(distancias, type="line", las=2,xaxt="n",main="Distancia media de pares", ylab="Distancia ponderada", 
     xlab="Year")
lines(br)
lines(ci)
axis(side = 1, at = seq(1959, 1979, by=1), las=2)
dev.off()

plotd<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=dist, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=dist, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=dist), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=dist), color="black", se=F)+
  theme_minimal()+
  labs(title="Distancia media", x="Year", y="Distancia del origen")
plotwd<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=wdist, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=wdist, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=wdist), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=wdist), color="black", se=F)+
  theme_minimal()+
  labs(title="Distancia media", x="Year", y="Distancia media de pares, ponderada")
plotud<-ggplot()+
  geom_point(data=tableU, aes(x=year, y=udist, group=point, color=point))+ 
  geom_smooth(data=tableU, aes(x=year, y=udist, group=point, color=point), se=F)+ 
  geom_point(data=dumd, aes(x=year, y=udist), color="black")+
  geom_smooth(data=dumd, aes(x=year, y=udist), color="black", se=F)+
  theme_minimal()+
  labs(title="Distancia media", x="Year", y="Distancia media de pares, no ponderada")

saveWidget(ggplotly(plotd), "distance2.html")
saveWidget(ggplotly(plotud), "udistance2.html")
saveWidget(ggplotly(plotwd), "wdistance2.html")
pdf("plotdist.pdf")
plotd+theme(legend.position="none")
dev.off()
pdf("plotudist.pdf")
plotud+theme(legend.position="none")
dev.off()
pdf("plotwdist.pdf")
plotwd+theme(legend.position="none")
dev.off()
#####graficas sensibilidad#####
tabinfl<-read_csv("tablainfl1.csv")
varinf<-tabinfl %>% subset(select=c(colus, infdisp, infudisp)) %>%
  pivot_longer(cols=c(infdisp, infudisp), names_to="variable", values_to="influence") %>%
  group_by(colus) %>%
  summarize(mean=mean(influence), max=max(influence), min=min(influence),  .groups="keep")
varinf<-varinf[order(varinf$mean, decreasing = T),]
varinf<-varinf[1:30,]
varinf<-varinf[order(varinf$mean),]
infvariable<-varinf %>%
  ggplot(aes(x = fct_inorder(colus),y = mean)) +
  geom_errorbar(aes(ymin = min, ymax = max))+
  coord_flip() +
  labs(x="Variable", y="Percentage of influence", title="Influence of variables in distance")
tabinfl2<-read_csv("tablainfl2.csv")
entinf<-tabinfl2 %>% subset(select=c(entidades, infdisp, infudisp)) %>%
  pivot_longer(cols=c(infdisp, infudisp), names_to="variable", values_to="influence") %>%
  group_by(entidades) %>%
  summarize(mean=mean(influence), max=max(influence), min=min(influence),  .groups="keep")
entinf<-entinf[order(entinf$mean),]
entinf<-entinf %>% 
  ggplot(aes(x = fct_inorder(entidades),y = mean)) +
  geom_errorbar(aes(ymin = min, ymax = max))+
  coord_flip() +
  labs(x="Entity", y="Percentage of influence", title="Influence of entities in distance")
jpeg("influencevariable.jpg", width=800, height=600, units="px")
print(infvariable)
dev.off()
jpeg("influenceentity.jpg", width=800, height=600, units="px")
print(entinf)
dev.off()
tableL<-read_csv("tableL.csv")
jpeg("plotdimensionstime.jpg", width=800, height=600, units="px")
ggplot(tableL, aes(x=year, y=pervariance, fill=dim))+
  geom_col(position="stack")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="Percentage of variance", x="Dimension", y="Eigenvalue")
dev.off()

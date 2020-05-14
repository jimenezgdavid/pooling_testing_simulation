library(tidyverse)

# Diferentes parámetros de entrada
      # Se crea set inicial
      # Se crea un loop para ir desde p0 = infIni/nSus hasta po = 100
            # Con con cada uno de estos estados se calcula para 2:64 m:
                  # N, nInf, Ntest
                  # Se calcula esto para diferentes niveles de orden ninguno, Fam, Nei, Buil, p0
            # Se agregan los análisis empleando heurísticas para dichas mismas agrupaciones.
      # Se calcula la tasa de crecimiento para diferentes parámetros de entrada
# Se presentan valores encontrados en las simulaciones de forma ordenada.

createSus <- function(N = 10000, pI = 0.005){
      
      nIni<- max(round(N*pI),1)
      nSus <- N
      inf <- sample(1:nSus, nIni)
      nFam <- ceiling(nSus/5)
      nBuil <- ceiling(nFam/2)
      nNei <- ceiling(nFam/10)
      
      nei <- tibble(nei = 1:nNei,
                    lng = sample(1:nNei,nNei), 
                    lat = sample(1:nNei,nNei))
      
      buil <- tibble(buil = 1:nBuil,
                     nei = sample(1:nNei, size = nBuil , replace = T))%>%
            left_join(nei,by = "nei")
      
      fam <- tibble(fam = 1:nFam,
                    buil = sample(1:nBuil, size = nFam , replace = T))%>%
            left_join(buil, by ="buil")
      
      sus <- tibble(sus = 1:nSus,
                    fam = sample(1:nFam, size = nSus, replace = T),
                    inf = sus %in% inf)%>%
            left_join(fam, by ="fam")%>%
            group_by(fam)%>%
            mutate(infFam = sum(inf))%>%
            ungroup()%>%
            group_by(buil)%>%
            mutate(infBuil =  sum(inf))%>%
            ungroup()%>%
            group_by(nei)%>%
            mutate(infNei =  sum(inf))%>%
            ungroup()%>%
            mutate(infP = InfProb (infFam,infBuil,infNei,nIni/nSus))
            
      
      return(sus)
}

createInfPFun<- function(modelCoefFact = 0.5,
                         modelCoef=list(fam = log(2),
                                        buil = log(1.5), 
                                        neib = log(1.25),
                                        p0 = log(2))){
      
      InfProb <- function(fam, buil, nei, p0){
            mapply(
                  FUN = function(famX,builX,neiX, p0X){
                        
                        l<- modelCoefFact*(famX * modelCoef$fam + 
                                   builX * modelCoef$buil + 
                                   neiX * modelCoef$neib + 
                                   p0X * modelCoef$p0 )
                        
                        
                        p<- 1/(1+exp(-l))
                        
                  },
                  fam,buil,nei,p0
            )
            
      }
      
      return(InfProb)
      
}

nextIter<-function(sus, threshold, p0){
      sus<-sus%>%
            mutate(infP = InfProb (infFam,infBuil,infNei,p0),
                   inf = inf | infP > (threshold+runif(1,-0.05,+0.05)),
            ) %>%
            group_by(fam)%>%
            mutate(infFam = sum(inf))%>%
            ungroup()%>%
            group_by(buil)%>%
            mutate(infBuil =  sum(inf))%>%
            ungroup()%>%
            group_by(nei)%>%
            mutate(infNei =  sum(inf))%>%
            ungroup()
      
      return(sus)
      
}


MOptD<-function(p0){
      max(1,floor(-1/log(1-p0)))
}

MOptT<-function(p0){
      nMax<-10000
      tol<-0.4
      i<-1
      m0<-0.5
      m1<-sqrt( -1 / ( ( 1 - p0 ) ^ m0 * log( 1 - p0 ) ) ) 
      while(i < nMax & abs(m0-m1) < 0.5){
            m0<-m1
            m1<-sqrt( -1 / ( ( 1 - p0 ) ^ m0 * log( 1 - p0 ) ) )
            i<- i+1
      }
      max(1,floor(m1))
}

TGroup<- function(data, 
                  groupSize = 64, 
                  arrangeVar = "sus", 
                  p0Sup = NA, 
                  optFun = "T"){
      
      MOptFun<-ifelse(optFun == "T",MOptT,MOptD)
      
      groupSize<-ifelse(is.na(p0Sup), groupSize, MOptFun(p0Sup))
      
      nGroups<-ceiling(dim(data)[1]/groupSize)
      
      group<-rep_len(1:nGroups,length.out = length(data$inf)) %>% sort()
      
      nDesc<-
            data %>%
            arrange_(arrangeVar)%>%
            mutate(group = group)%>%
            group_by(group) %>%
            mutate(negTest = sum(inf) == 0) %>%
            ungroup() %>%
            summarise(nDesc = sum(negTest)) %>%
            unlist()
      
      tibble(arrangeVar = arrangeVar,
             groupSize = groupSize,
             totalSus = length(data$inf),
             nDesc = nDesc,
             nGroups = nGroups,
             opt = ifelse(!is.na(p0Sup),optFun,NA),
             totalTest = ifelse(groupSize == 1, dim(data)[1], (dim(data)[1] - nDesc) + nGroups ))
}

# set.seed(0)
# sus%>%
#     ggplot()+
#     geom_jitter(aes(lng,lat,col = inf, alpha= inf), width = 1, height = 1,)+
#     scale_color_manual(values= c("grey50","red"))+
#     scale_alpha_manual(values =c(0.3,0.80))+
#     theme_minimal()




nSim<- 100
pI <- sample( c(0.001, 0.005,0.01), nSim, replace = T) 
threshold <- runif(nSim,0.7,0.75)
n<- round(runif(nSim,500,3000))

# pI<-c(0.001,0.005,)
# threshold<-c(0.8,0.8)
# n<- c(1000,2000,3000)


simData<-
      mapply(
            FUN = function(pIX,thresholdX,nX){
                  #Se crea la función de probabilidad de contagio
                  InfProb<-createInfPFun()
                  
                  #Se crea la población inicial
                  sus<-createSus(N = nX, pI = pIX)
                  
                  #Se determina el p0 de la población en la iteración i
                  p0 <- sum(sus$inf)/dim(sus)[1]
                  i<- 1
                  
                  #Se crea la función para iterar en la poblacion
                  maxIter<-1000
                  
                  p0_n <- c(p0)
                  i_n  <- c(i)
                  
                  tGroupTest <- NULL
                  
                  while( p0 < 0.99 & i <= maxIter ){
                        
                        tGroupTest<-
                              bind_rows(
                                    tGroupTest,

                                    lapply(1:64,
                                           FUN = function(gSize){
                                                 bind_rows(
                                                       TGroup(sus, groupSize = gSize),
                                                       TGroup(sus, groupSize = gSize,arrangeVar = "fam"),
                                                       TGroup(sus, groupSize = gSize,arrangeVar = "buil"),
                                                       TGroup(sus, groupSize = gSize, arrangeVar ="nei"),
                                                       TGroup(sus, groupSize = gSize, arrangeVar = "infP")
                                                 )
                                           })%>%
                                          bind_rows()%>%
                                          group_by(arrangeVar)%>%
                                          mutate(opt = rank(totalTest,ties.method = "first"))%>%
                                          ungroup()%>%
                                          filter(opt == 1)%>%
                                          mutate(opt ="num",
                                                 i=i,
                                                 p0 = p0)
                                    # ,
                                    # 
                                    # bind_rows(
                                    #       TGroup(sus, p0 = p0),
                                    #       TGroup(sus, p0 = p0,arrangeVar = "fam"),
                                    #       TGroup(sus, p0 = p0,arrangeVar = "buil"),
                                    #       TGroup(sus, p0 = p0, arrangeVar ="nei"),
                                    #       TGroup(sus, p0 = p0, arrangeVar = "infP"),
                                    #       TGroup(sus, p0 = p0),
                                    #       TGroup(sus, p0 = p0,arrangeVar = "fam", optFun = "D"),
                                    #       TGroup(sus, p0 = p0,arrangeVar = "buil",optFun = "D"),
                                    #       TGroup(sus, p0 = p0, arrangeVar ="nei",optFun = "D"),
                                    #       TGroup(sus, p0 = p0, arrangeVar = "infP",optFun = "D")
                                    # )%>%
                                    #       mutate(i = i,
                                    #              p0 = p0)
                              )
                        
                        
                        sus <- nextIter(sus,threshold = thresholdX, p0 = p0)
                        p0 <- sum(sus$inf) / dim(sus)[1]
                        i<- i+1
                        
                        p0_n<- c(p0_n,p0)
                        i_n <- c(i_n,i)
                        
                        if(i == maxIter){
                              print("Stopped withoud converging to p0 = 1, reached iterations limit")
                        }
                  }
                  
                  # plot(i_n,p0_n)
                  # tibble(i = i_n, 
                  #        p0 = p0_n,
                  #        pI = pIX,
                  #        threshold = thresholdX,
                  #        n = nX )
                  # tGroupTest%>% ggplot()+geom_point(aes(p0,totalTest,col = arrangeVar))+facet_grid(opt~.)
                  tGroupTest%>%mutate(pI = pIX,threshold = thresholdX, n = nX )
                  
            } , pI, threshold, n, SIMPLIFY = F)

simData<-simData %>% bind_rows()

simData%>%
      write.csv("testing_sim/data/simData2.csv", row.names = F)


simData<-read.csv("testing_sim/data/simData1.csv")

plotBackG<- 
      simData %>%
      mutate(group = paste(threshold,n,pI,arrangeVar, sep = "_"))%>%
      select(-arrangeVar)
      
 simData %>% 
      mutate(totalTestR=totalTest/totalSus)%>%
      # .[[10]]%>% 
      ggplot()+
      geom_line(aes(p0,
                    totalTestR,
                    group = paste(threshold,arrangeVar),
                    col =arrangeVar),
                alpha = 0.1)+
       geom_line(data = plotBackG ,aes(p0,totalTestR, group = group))
      theme_minimal()

simData %>% 
       bind_rows() %>%
       mutate(totalTestR=totalTest/totalSus) %>%
       filter(totalTestR < 0.95) %>%
       # .[[10]]%>% 
       ggplot()+
       geom_smooth(aes(p0,totalTestR, col =arrangeVar ), formula = "y~x")+
       geom_point(aes(p0,totalTestR,col = arrangeVar))



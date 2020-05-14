# K variales subyacentes tags
# Misma familia o grupo social
# 
library(tidyverse)

# Dat set generation ####
nIni<- 3
nSus <- 1000
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
    ungroup()

susDist<-
    dist(sus%>%select(lng,lat),upper = T,diag = T)%>%
    as.matrix()
susInvDist<-1/(1+susDist)
maxInvDist<-mean(susInvDist %*% rep(1,nSus))


#Second generation
fact<-0.5
threshold<-0.85
coef<-list(fam = 1*fact, buil = 0.5*fact, neib = 0.25*fact, invDist = 2*fact)
b <-exp(1)
noiseSd<- 1

InfProb <- function(fam, buil, nei, invDist){
    mapply(
        FUN = function(famX,builX,neiX, invDistX){
            l<- famX * coef$fam + 
                builX * coef$buil + 
                neiX * coef$neib + 
                invDistX * coef$invDist +
                rnorm(sd = noiseSd,n = 1)*fact
            
            p<- 1/(1+b^(-l))
        },
        fam,buil,nei,invDist
    )
    
}

sus<-sus%>%    mutate(invDist = as.vector(susInvDist %*% inf / maxInvDist ))



sus<-sus%>%
    mutate(infP = InfProb (infFam,infBuil,infNei,invDist),
           inf = inf | infP > threshold,
           invDist =as.vector(susInvDist %*% inf / maxInvDist )
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

sum(sus$inf)/length(sus$inf)*100

p0<-sum(sus$inf)/length(sus$inf)

# dir<- "C:/Users/David Jimenez/OneDrive - Instituto Costarricense de Electricidad/03-OTROS/2020_COVID19/testing_sim"
# 
# write.csv(sus,file.path(dir,paste0(nIni,"_",nSus,"_","sim.csv")), row.names = F)

sus<-read.csv("testing_sim/1_1000_sim.csv")

# set.seed(0)
# sus%>%
#     ggplot()+
#     geom_jitter(aes(lng,lat,col = inf, alpha= inf), width = 1, height = 1,)+
#     scale_color_manual(values= c("grey50","red"))+
#     scale_alpha_manual(values =c(0.3,0.80))+
#     theme_minimal()
    

ggsave("out/inf_pop.png",width = 7,height = 7, dpi = "retina")


#Implementación del algoritmo de test

#Caso Base testear a todo el mundo:

sus%>%
    summarise(nTest = n(), inf = sum(inf))

#Caso base usando n tamaños y escogiendo de forma aleatoria:

TGroup<- function(data, 
                  groupSize = 64, 
                  arrangeVar = "sus", 
                  p0 = NA, 
                  MOptFun = MOptT){
    
    groupSize<-ifelse(is.na(p0), groupSize, MOptFun(p0))
    
    nGroups<-ceiling(dim(data)[1]/groupSize)
    
    nDesc<-
        data %>%
        arrange_(arrangeVar)%>%
        mutate(group = rep_len(1:nGroups,length.out = length(data$inf)) %>% sort())%>%
        group_by(group) %>%
        mutate(negTest = sum(inf) == 0) %>%
        ungroup() %>%
        summarise(nDesc = sum(negTest)) %>%
        unlist()
    
    tibble(arrangeVar = arrangeVar,
           groupSize = groupSize,
           nDesc = nDesc,
           nGroups = nGroups,
           totalTest = dim(data)[1] - nDesc + nGroups,
           MOptFun=as.character(substitute(package)))
}



TGroupSimp<- function(data, groupSize=length(inf)/64){
    # groupSize<-ifelse(is.na(p0), groupSize, floor(1/sqrt(p0/100)) )
    nGroups<-ceiling(dim(data)[1]/groupSize)
    nDesc<-data%>%
        mutate(group = rep_len(1:nGroups,length.out = length(data$inf)) %>% sort())%>%
        group_by(group)%>%
        mutate(negTest = sum(inf) == 0)%>%
        ungroup()%>%
        summarise(nDesc = sum(negTest)) %>%
        unlist()
    
    tibble(groupSize = groupSize,
           nDesc = nDesc,
           nGroups = nGroups,
           totalTest = dim(data)[1] - nDesc + nGroups)
}


groupSize<-2:64

tGroupTest<-
    lapply(groupSize,
       FUN = function(gSize){
           bind_rows(
               TGroup(sus, groupSize = gSize)%>% mutate(testType = "group"),
               TGroup(sus, groupSize = gSize, arrangeVar = c("fam"))%>% mutate(testType = "fam"),
               TGroup(sus, groupSize = gSize, arrangeVar = c("buil"))%>% mutate(testType = "buil"),
               TGroup(sus, groupSize = gSize, arrangeVar = c("nei"))%>% mutate(testType = "nei"),
               TGroup(sus, groupSize = gSize, arrangeVar = c("infP"))%>% mutate(testType = "prob")
           )
           })%>%
    bind_rows()

tGroupTest %>%
    ggplot()+
    geom_point(aes(nGroups, totalTest, col = testType))


sus%>%
    arrange(-infP)%>%
    mutate(ind = 1:length(sus))%>%
    ggplot()+
    geom_point(aes(ind,infP))



#Pruebas de algoritmos:
#

MOptD<-function(p0){
    max(1,floor(-1/log(1-p0)))
}

MOptT<-function(p0){
    nMax<-10000
    tol<-0.4
    i<-0
    m0<-0.5
    m1<-sqrt( -1 / ( ( 1 - p0 ) ^ m0 * log( 1 - p0 ) ) ) 
    while(i < nMax & abs(m0-m1) < 0.5){
        m0<-m1
        m1<-sqrt( -1 / ( ( 1 - p0 ) ^ m0 * log( 1 - p0 ) ) )
        i<- i+1
    }
    max(1,floor(m1))
}

tibble(p0 = seq(0.001,1,0.01),
       m1 = sapply(p0,FUN = MOptD),
       m2 = sapply(p0,FUN = MOptT))%>%
    pivot_longer(names_to = "var",
                 values_to = "val",
                 cols = c(m1,m2))%>%
    ggplot()+
    geom_point(aes(p0,val,col=var))+
    scale_y_log10()


# n<-1000
# p0<-p0
# pop<-tibble(inf=1:n %in% sample(1:n,n*p0))

TVarGS1<- 
    function(data,dir = 1, MOptFun = MOptT){
        p0<-sum(data$inf)/length(data$inf)
        
        data<-data %>% 
            arrange(-dir*infP)%>%
            mutate(pEst =infP/mean(infP)*p0,
                   m0 = sapply(pEst,FUN  = MOptFun))
        
        i<-1
        j<-1
        nMax<-1000
        group<-c()
        m<-data$m0
        
        while(j <= length(m) | nMax < i ){
            k<-min((j+m[j]-1),length(m))
            group<-c(group,(m[j:k]*0+1)*i)
            j<-k+1
            i<-i+1
        }
        
        nDesc<-data%>%
            mutate(group=group)%>%
            group_by(group)%>%
            mutate(negTest = sum(inf) == 0)%>%
            ungroup()%>%
            summarise(nDesc = sum(negTest)) %>%
            unlist()
        
        nOne<-data%>%
            mutate(group = group)%>%
            group_by(group)%>%
            summarise(one = n() == 1)%>% 
            ungroup()%>%
            summarise(nOne = sum(one))%>%
            unlist()
        
        tibble(nDesc = nDesc,
               nGroups = i-1,
               totalTest = dim(data)[1] - nDesc + nGroups - nOne)
               #totalTest = dim(data)[1] - nDesc + nGroups)
}


groupSize<-2:64

tGroupTest<-
    lapply(groupSize,
           FUN = function(gSize){
               bind_rows(
                   TGroup(sus, groupSize = gSize)%>% mutate(testType = "group"),
                   TGroup(sus, groupSize = gSize, arrangeVar = c("fam"))%>% 
                       mutate(testType = "fam"),
                   TGroup(sus, groupSize = gSize, arrangeVar = c("buil"))%>%
                       mutate(testType = "buil"),
                   TGroup(sus, groupSize = gSize, arrangeVar = c("nei"))%>% 
                       mutate(testType = "nei"),
                   TGroup(sus, groupSize = gSize, arrangeVar = c("infP"))%>% 
                       mutate(testType = "prob"),
                   TGroup(sus,p0 =  p0)%>%
                       mutate(testType = "optimT", one =T),
                   TGroup(sus,p0 =  p0,MOptFun = MOptD)%>%
                       mutate(testType = "optimD",one =T),
                   TGroup(sus,p0 =  p0, arrangeVar = c("infP"))%>%
                       mutate(testType = "optimT_p", one =T),
                   TGroup(sus,p0 =  p0,MOptFun = MOptD, arrangeVar = c("infP"))%>%
                       mutate(testType = "optimD_p",one =T),
                   TVarGS1(sus)%>%
                       mutate(testType = "TL2R",one =T),
                   TVarGS1(sus,dir = -1)%>%
                       mutate(testType = "TR2L",one =T),
                   TVarGS1(sus,MOptFun = MOptD)%>%
                       mutate(testType = "TL2Rs_D",one =T),
                   TVarGS1(sus,MOptFun = MOptD,dir = -1)%>%
                       mutate(testType = "TR2Ls_D",one =T)
               )
           })%>%
    bind_rows()%>%
    mutate(one = !is.na(one))

htmlwidgets::saveWidget(
    plotly::ggplotly(
        tGroupTest %>%
            ggplot()+
            geom_point(aes(nGroups, totalTest, col =testType, 
                           size =factor(one*1), shape = factor(one*1)))+
            scale_size_manual(values = c(1,5))+
            scale_shape_manual(values = c(19,10))+
            scale_color_viridis_d()
    ),file = "test.htm"
)

ggsave("out/graph_pop.png",width = 7,height = 7, dpi = "retina")

tGroupTest%>%
    filter(totalTest == min(totalTest))




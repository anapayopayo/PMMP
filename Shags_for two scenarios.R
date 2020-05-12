
###  Seasonally mobile metapopulation lambda and elasticities analitical solution  #### 
###  Author: Ana Payo Payo                                                         ####
###  Date:      07/08/2019                                                         ####

#### House keeping  #### 
#### Clean the workspace 
rm(list=ls())
ls()
#### Load libraries  ####
library(varhandle)
library(ggplot2)
library(popbio)
library(magick)
library(gridExtra)
library(grid)
library(lattice)
library("RColorBrewer")
library(viridis)
####          get/set working directory        #### 
getwd()
#### Model projection matrix #### 

#Mlife   matrix ####
A=expression(0,0,  0,0,-f1_R*s011_R*s021_R*(m01-1),     -f1_M*s011_M*s021_M*(m01-1),        0,0,0,0,0,0,
             0,0,  0,0,f1_R*m01*s011_R*s022_R*w0*w0,f1_M*m01*s011_M*s022_M*w0*w0,       0,0,0,0,0,0,
             
             s111_R*s121_R,0,0,0,0,0,                                                       0,0,0,0,0,0,
             0,s111_M*s122_M*w1*w1,0,0,0,0,                                                 0,0,0,0,0,0,
             
             0,0,s211_R*s221_R,0,s211_R*s221_R,0,                                           0,0,0,0,0,0,
             0,0,0,s211_M*s222_M*w2*w2,0,s211_M*s222_M*w2*w2,                               0,0,0,0,0,0,
             
             0,0,0,0,0,0,                                                                   0,0,0,0,-f2_R*s012_R*s022_R*(m02-1),-f2_M*s012_M*s022_M*(m02-1),
             0,0,0,0,0,0,                                                                   0,0,0,0,f2_R*m02*s012_R*s021_R*w0*w0,f2_M*m02*s012_M*s021_M*w0*w0,
             
             0,0,0,0,0,0,                                                                   s112_R*s122_R,0,0,0,0,0,
             0,0,0,0,0,0,                                                                   0,s112_M*s121_M*w1*w1,0,0,0,0,
             
             0,0,0,0,0,0,                                                                   0,0,s212_R*s222_R,0,s212_R*s222_R,0,
             0,0,0,0,0,0,                                                                   0,0,0,s212_M*s221_M*w2*w2,0,s212_M*s221_M*w2*w2)

#Dispersal matrix ####
D=expression(0,0,0,0,0,0,d,0,0,0,0,0,
             0,0,0,0,0,0,0,d,0,0,0,0,
             0,0,0,0,0,0,0,0,d,0,0,0,
             0,0,0,0,0,0,0,0,0,d,0,0,
             0,0,0,0,0,0,0,0,0,0,d,0,
             0,0,0,0,0,0,0,0,0,0,0,d,
             d,0,0,0,0,0,0,0,0,0,0,0,
             0,d,0,0,0,0,0,0,0,0,0,0,
             0,0,d,0,0,0,0,0,0,0,0,0,
             0,0,0,d,0,0,0,0,0,0,0,0,
             0,0,0,0,d,0,0,0,0,0,0,0,
             0,0,0,0,0,d,0,0,0,0,0,0)

### AA metapopulation projection matrix  AA = D M   ###### 

AA =expression(
  
  0,                           0,                      0,                           0, f1_R*s011_R*s021_R*(d - 1)*(m01 - 1), f1_M*s011_M*s021_M*(d - 1)*(m01 - 1),                      0,                           0,                      0,                           0,      -d*f2_R*s012_R*s022_R*(m02 - 1),      -d*f2_M*s012_M*s022_M*(m02 - 1),
  0,                           0,                      0,                           0, -f1_R*m01*s011_R*s022_R*w0^2*(d - 1), -f1_M*m01*s011_M*s022_M*w0^2*(d - 1),                      0,                           0,                      0,                           0,        d*f2_R*m02*s012_R*s021_R*w0^2,        d*f2_M*m02*s012_M*s021_M*w0^2,
  -s111_R*s121_R*(d - 1),                           0,                      0,                           0,                                    0,                                    0,        d*s112_R*s122_R,                           0,                      0,                           0,                                    0,                                    0,
  0, -s111_M*s122_M*w1^2*(d - 1),                      0,                           0,                                    0,                                    0,                      0,        d*s112_M*s121_M*w1^2,                      0,                           0,                                    0,                                    0,
  0,                           0, -s211_R*s221_R*(d - 1),                           0,               -s211_R*s221_R*(d - 1),                                    0,                      0,                           0,        d*s212_R*s222_R,                           0,                      d*s212_R*s222_R,                                    0,
  0,                           0,                      0, -s211_M*s222_M*w2^2*(d - 1),                                    0,          -s211_M*s222_M*w2^2*(d - 1),                      0,                           0,                      0,        d*s212_M*s221_M*w2^2,                                    0,                 d*s212_M*s221_M*w2^2,
  0,                           0,                      0,                           0,      -d*f1_R*s011_R*s021_R*(m01 - 1),      -d*f1_M*s011_M*s021_M*(m01 - 1),                      0,                           0,                      0,                           0, f2_R*s012_R*s022_R*(d - 1)*(m02 - 1), f2_M*s012_M*s022_M*(d - 1)*(m02 - 1),
  0,                           0,                      0,                           0,        d*f1_R*m01*s011_R*s022_R*w0^2,        d*f1_M*m01*s011_M*s022_M*w0^2,                      0,                           0,                      0,                           0, -f2_R*m02*s012_R*s021_R*w0^2*(d - 1), -f2_M*m02*s012_M*s021_M*w0^2*(d - 1),
  d*s111_R*s121_R,                           0,                      0,                           0,                                    0,                                    0, -s112_R*s122_R*(d - 1),                           0,                      0,                           0,                                    0,                                    0,
  0,        d*s111_M*s122_M*w1^2,                      0,                           0,                                    0,                                    0,                      0, -s112_M*s121_M*w1^2*(d - 1),                      0,                           0,                                    0,                                    0,
  0,                           0,        d*s211_R*s221_R,                           0,                      d*s211_R*s221_R,                                    0,                      0,                           0, -s212_R*s222_R*(d - 1),                           0,               -s212_R*s222_R*(d - 1),                                    0,
  0,                           0,                      0,        d*s211_M*s222_M*w2^2,                                    0,                 d*s211_M*s222_M*w2^2,                      0,                           0,                      0, -s212_M*s221_M*w2^2*(d - 1),                                    0,          -s212_M*s221_M*w2^2*(d - 1))

#### Model parametrization ####

##Fredericksen 2008 survival values  as a reference of age survival relationships.
#newborn survival 0.513
#juvenile survival 0.737
#adult survival 0.858 
## To derive these relationships, we used adult survival as a reference and forced newborn and juvenile survival to be proportions
#of adult survival. This can be achieved by first dividing all annual age class survival values  by annual adult survival  and, then,
#calculating its square root to transform it into seasonal relationships. 

    a0=round(sqrt(0.513/0.858),3)  #newborn survival 0.513
    a1=round(sqrt(0.737/0.858),3)  #juvenile survival 0.737
    a2=round(sqrt(0.858/0.858),3)  #adult survival 0.858 

    ## relationship is a2:a1:a0 ->  1:0.927:0.773

####Define parameters for each scenario #####
    nscn= 2
    m<-c(0.441,0.306)
    l<-matrix(0,ncol=nscn,nrow=1)
    sss<-matrix(0,ncol=nscn,nrow=12)
    e<-array(0, dim=c(34,12,nscn))
    lab<-c("For observed values  (m)", "Twice the observed values  (2m)")

dim(e)
for (scn in 1:2){
  #scenario 2 is 2*m 
  vr<- list(
    f1_R=2.05/2,          # Mean offspring IoM Resident
    f1_M=1.75/2,          # Mean offspring IoM Migrant
    
    f2_R=1.21/2,          # Mean offspring BoB Resident
    f2_M=0.97/2,          # Mean offspring BoB Migrant
    
    s011_R=0.984*a0,      # Breeding  survival newborn  IoM Resident
    s011_M=0.984*a0,      # Breeding  survival newborn  IoM Migrant 
    s012_R=0.993*a0,      # Breeding  survival newborn  BoB Resident 
    s012_M=0.993*a0,      # Breeding  survival newborn  BoB Migrant 
    
    s021_R=0.917 *a0,     # Non breeding  survival newborn   IoM Resident
    s021_M=0.918 *a0,     # Non breeding  survival newborn   IoM Migrant 
    s022_R=0.922 *a0,     # Non breeding  survival newborn   BoB Resident
    s022_M=0.912 *a0,     # Non breeding  survival newborn   BoB Migrant
    
    s111_R=0.984*a1,      # Breeding  survival adult    IoM Resident
    s111_M=0.984*a1,      # Breeding  survival adult    IoM Migrant 
    s112_R=0.993*a1,      # Breeding  survival adult    BoB Resident
    s112_M=0.993*a1,      # Breeding  survival adult    BoB Migrant
    
    s121_R=0.917*a1,      # Non breeding  survival juvenile and adult   IoM Resident
    s121_M=0.918*a1,      # Non breeding  survival juvenile and adult  IoM Migrant 
    s122_R=0.922*a1,      # Non breeding  survival juvenile and adult  BoB Resident
    s122_M=0.912*a1,      # Non breeding  survival juvenile and adult  BoB Migrant
    
    s211_R=0.984,         # Breeding  survival adult    IoM Resident
    s211_M=0.984,         # Breeding  survival adult    IoM Migrant 
    s212_R=0.993,         # Breeding  survival adult    BoB Resident
    s212_M=0.993,         # Breeding  survival adult    BoB Migrant
    
    s221_R=0.917,         # Non breeding  survival juvenile and adult   IoM Resident
    s221_M=0.918,         # Non breeding  survival juvenile and adult  IoM Migrant 
    s222_R=0.922,         # Non breeding  survival juvenile and adult  BoB Resident
    s222_M=0.912,         # Non breeding  survival juvenile and adult  BoB Migrant
    
    
    m01=m[1]*scn,           # Migration from IoM 
    m02=m[2]*scn,           # Migration from BoB
    
    w0=sqrt(0.99)*a0,     # Survival migration newborn
    w1=sqrt(0.99)*a1,     # Survival migration yearling
    w2=sqrt(0.99),        # Survival migration adults
    d=0.1                 # Dispersal probability
  )
  
  ##### GET LAMBDA #####
  
  A<-matrix(sapply(AA, eval,vr , NULL), nrow=12,ncol=12, byrow=TRUE)
  lambda<-max(Re(eigen(t(A))$values))        # long-term geometric rate of population growth
  SSS<-stable.stage(A)
  l[1,scn]<-lambda
  sss[,scn]<-SSS
  ##### GET ELASTICITY     #####
  elasticities<-vitalsens(AA, vr)
  elasticities$names<-  c("f1_R"  , "f1_M" ,  "f2_R"  , "f2_M" ,  "s011_R" ,"s011_M" ,"s012_R" ,"s012_M" ,"s021_R" ,"s021_M",
                          "s022_R", "s022_M" ,"s111_R", "s111_M" ,"s112_R" ,"s112_M", "s121_R", "s121_M" ,"s122_R", "s122_M",
                          "s211_R" ,"s211_M", "s212_R", "s212_M" ,"s221_R" ,"s221_M", "s222_R", "s222_M" ,"zz m01"   , "zz m02"  , 
                          "z w0"   ,  "z w1" ,    "z w2"    , "zzz d")   
  elasticities$t<-rep(c(1,2),17)
  elasticities[29:34,5]<-3
  elasticities$t<-as.factor(elasticities$t)
  elasticities$colony<-c(rep(c(1,1,2,2),7),3,3,3,3,3,3)
  
  elasticities$colony <- c("IoM"  , "IoM" ,  "BoB"  , "BoB" ,  "IoM" ,"IoM" ,"BoB" ,"BoB" ,"IoM" ,"IoM",
                           "BoB", "BoB" ,"IoM", "IoM" ,"BoB" ,"BoB", "IoM", "IoM" ,"BoB", "BoB",
                           "IoM" ,"IoM", "BoB", "BoB" ,"IoM" ,"IoM", "BoB", "BoB" ,"na"   , "na"  , 
                           "na"   ,  "na" ,    "na"    , "na")
  
  elasticities$strategy <- c("R"  , "M" ,  "R"  , "M" ,  "R" ,"M" ,"R" ,"M" ,"R" ,"M",
                             "R", "M" ,"R", "M" ,"R" ,"M", "R", "M" ,"R", "M",
                             "R" ,"M", "R", "M" ,"R" ,"M", "R", "M" ,"na"   , "na"  , 
                             "na"   ,  "na" ,    "na"    , "na")
  elasticities$para <- c("f"  , "f" ,  "f"  , "f" ,  "s" ,"s" ,"s" ,"s" ,"s" ,"s",
                         "s", "s" ,"s", "s" ,"s" ,"s", "s", "s" ,"s", "s",
                         "s" ,"s", "s", "s" ,"s" ,"s", "s", "s" ,"m"   , "m"  , 
                         "w"  ,  "w" ,    "w"    , "d")
  elasticities$age <- c("2"  , "2" ,  "2"  , "2" , 
                        "0" ,"0" ,"0" ,"0" ,"0" ,"0","0", "0" ,
                        "1", "1" ,"1" ,"1", "1", "1" ,"1", "1",
                        "2" ,"2", "2", "2" ,"2" ,"2", "2", "2" ,"na"   , "na"  , 
                        "na"   ,  "na" ,    "na"    , "na")
  elasticities$season<- c("B"  , "B" , "B"  , "B" , 
                          "B"  , "B" , "B"  , "B" ,"NB" ,"NB" ,"NB" ,"NB" ,
                          "B"  , "B" ,"B"  , "B" , "NB" ,"NB" ,"NB" ,"NB" ,
                          "B"  , "B" ,"B"  , "B" , "NB" ,"NB" ,"NB" ,"NB" ,
                          "na","na","na","na","na","na")
  axisLabels.x <- c("f1_R"  , "f1_M" ,  "f2_R"  , "f2_M" ,  "s011_R" ,"s011_M" ,"s012_R" ,"s012_M" ,"s021_R" ,"s021_M",
                    "s022_R", "s022_M" ,"s111_R", "s111_M" ,"s112_R" ,"s112_M", "s121_R", "s121_M" ,"s122_R", "s122_M",
                    "s211_R" ,"s211_M", "s212_R", "s212_M" ,"s221_R" ,"s221_M", "s222_R", "s222_M" ,"m01"   , "m02"  , 
                    "w0"   ,  "w1" ,    "w2"    , "d")
  labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
  elasticitiesmeta<-elasticities
  
  meta<-ggplot(data=elasticities, aes(x=names, y=elasticity, fill=t)) +
    geom_bar(stat="identity", colour="black") +
    #scale_fill_manual(values=c( "deepskyblue4", "cyan","skyblue4"))+ 
    scale_fill_manual(values = c("black", "darkgrey","white"),labels= c("R","M","Other"), guide = guide_legend(reverse = TRUE))+
    #  scale_x_discrete(labels=axisLabels.x)+
    ylim(-0.1, 0.6)+
    geom_hline(yintercept=0, linetype=1, color = "black")+
    theme(axis.text.x = element_text(angle = 90),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Change axis line
          axis.line = element_line(colour = "black")
    )
  data_meta<- elasticitiesmeta
  data_meta_IoM<-data_meta[which(elasticities$colony=="IoM"),]
  data_meta_BoB<-data_meta[which(elasticities$colony=="BoB"),]
  data_meta_other<-data_meta[which(elasticities$colony=="na"),]
  data_meta_IoM$what<-data_meta_BoB$what<-c("f",   
                                            "f",   
                                            "s011", 
                                            "s011",
                                            "s021",
                                            "s021", 
                                            "s111", 
                                            "s111", 
                                            "s121", 
                                            "s121", 
                                            "s211", 
                                            "s211",
                                            "s221",
                                            "s221")
  data_meta_IoM$order<-c(1,   
                         1,   
                         2, 
                         2,
                         3,
                         3, 
                         4, 
                         4, 
                         5, 
                         5, 
                         6, 
                         6,
                         7,
                         7)
  
  data_meta_BoB$order<-c(8,   
                         8,   
                         9,
                         9,
                         10,
                         10,
                         11,
                         11, 
                         12, 
                         12, 
                         13, 
                         13, 
                         14, 
                         14)
  data_meta_other$what<-c("m1",   
                          "m2", 
                          "w0",
                          "w1",
                          "w2",
                          "d")
  data_meta_other$order<-c(18,   
                           19, 
                           15,
                           16,
                           17,
                           20)
  data<-rbind(data_meta_IoM,data_meta_BoB,data_meta_other)
  data$order<-as.factor(data$order)
  e[,,scn]<-as.matrix(data)
  p<-ggplot(data=data, aes(x=order, y=elasticity, fill=t)) +
    geom_bar(stat="identity", colour="black",position="stack") +
    #scale_fill_manual(values=c( "deepskyblue4", "cyan","skyblue4"))+ 
    scale_fill_manual(values = c( "darkgrey","black","white"),labels= c("R","M","Other"), guide = guide_legend(reverse = TRUE))+
    #  scale_x_discrete(labels=axisLabels.x)+
    geom_hline(yintercept=0, linetype=1, color = "black")+
    geom_vline(xintercept=7.5, linetype=2, color = "grey")+
    geom_vline(xintercept=14.5, linetype=2, color = "grey")+
    ylim(-0.25, 1)+
    geom_text(x=2, y=1, label="IoM",size=8)+
    geom_text(x=9, y=1, label="BoB",size=8)+
    geom_text(x=17, y=1, label="Common",size=8)+
    geom_text(x=2.5, y=-0.08, label=paste(lab[scn]))+
    geom_text(x=2.7, y=-0.18, label=paste("m1=", scn*m[1],"  m2= ", scn*m[2]))+
    
    scale_x_discrete(labels= c("f1","s011","s021","s111","s121","s211","s221",
                               "f2","s012","s022","s112","s122","s212","s222",
                               "w0","w1", "w2","m1", "m2", "d"  )) +
    xlab(expression("Parameter ("*theta*")")) +
    ylab(expression("Elasticity of "*lambda*" to "*theta)) + 
    theme(axis.text.x = element_text(angle = 0),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Change axis line
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")
    )
  
assign(paste("plot",scn, sep="_"),p )

  }    
grid.arrange( plot_1,plot_2)
d<-as.data.frame(rbind(e[,,1],e[,,2]))
when<-c(rep(1,34),rep(2,34))
d<-cbind(d,when)
names(d)<-c(names(data),"when")
d$when<-as.factor(d$when)
d$elasticity<-as.numeric(paste(d$elasticity))
d$order<-as.numeric(paste(d$order))
facet_names<- c("1"="f1","2"="s011","3"="s021","4"="s111","5"="s121","6"="s211","7"="s221",
                "8"= "f2","9"="s012","10"="s022","11"="s112","12"="s122","13"="s212","14"="s222",
                "15"="w0","16"="w1","17"= "w2","18"="m1", "19"="m2", "20"="d"  )
ggplot(data=d, aes(x=when, y=elasticity, group=strategy)) +
  geom_line(aes(colour = factor(strategy)),linetype="dotted") +
  geom_point(aes(colour=factor(strategy),fill=factor(strategy)),shape=21,size=2) +
  scale_color_manual(values=c('black','black',"grey"))+
  scale_fill_manual(values=c("black", "white","grey")) + 
  geom_hline(yintercept=0, linetype=1, color = "grey")+
  scale_x_discrete(name = "Scenario", labels = c("m","2m"))+
  ylab(expression("Elasticity of "*lambda*" to "*theta)) + 
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'),
        axis.text.x = element_text(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        legend.position = "none"
  )+facet_wrap( ~ order,nrow=1, labeller = as_labeller(facet_names))


l
sss
seq(from=1, to=12,by=2)
sss[,seq(from=1, to=12,by=2)]
migrants<-matrix(0,ncol=2,nrow=6)
odd <- function(x) x%%2 != 0
even <- function(x) x%%2 == 0
migrants<-sss[which(even(c(1:12))==T),]
colSums(migrants)
colSums(migrants[1:3,])
colSums(migrants[4:6,])
e

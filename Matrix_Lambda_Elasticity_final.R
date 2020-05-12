
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
#### Define model characteristics #### 
nparameters<-18                        # number of parameters to estimate elasticities
scenarios<-7                           # number of spatio-temporal scenarios
species<-c("Shortlived",               # 1-SHORTLIVED  
           "Longlived")                # 2-LONGLIVED
nLHS=2                                 # number of species Shortlived/Longlived
a=c(1,3)                               # Age classes - age at first reproducion for each LHS
k=2                                    # Number of patches hosting the metapopulation.
                                       #P1-Patch 1 
                                       #P2-Patch 2 
                                       #P3- Metapopulation
nmodels<-3                             # Number of models
states<-c(1,1,2)                       # States dimension tag
ntag<-c(1,1,2)                         # Matrix dimension tag


#define the models in a matrix
# Define Model Matrices ####
# Matrices are defined for three different temporal extents of seasonal movement plasticity
    #Mseason
    #Mannual
    #Mlife
#These matrices are specified for P1,P2 and the whole metapopulation for each model and life-history strategy.

#####   Shortlived species  #####
#Mseason
      #Metapopulation
      MDSM =expression(                  s111*s121*(m1 - 1)^2 + f1*s011*s021*(m01 - 1)^2 + m1*m2*s111*s122*w1*w2 + f1*m01*m02*s011*s022*w01*w02, - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1) - f2*m02*s012*s021*w01*(m01 - 1) - f2*m02*s012*s022*w02*(m02 - 1),
                               - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1) - f1*m01*s011*s021*w02*(m01 - 1) - f1*m01*s011*s022*w01*(m02 - 1),                  s112*s122*(m2 - 1)^2 + f2*s012*s022*(m02 - 1)^2 + m1*m2*s112*s121*w1*w2 + f2*m01*m02*s012*s021*w01*w02)
      #Patch1  
      MDSP1 = expression(s111*s121*(m1 - 1)^2 + f1*s011*s021*(m01 - 1)^2 - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1) + m1*m2*s111*s122*w1*w2 - f2*m02*s012*s021*w01*(m01 - 1) - f2*m02*s012*s022*w02*(m02 - 1) + f1*m01*m02*s011*s022*w01*w02)
      #Patch2      
      MDSP2 =expression(s112*s122*(m2 - 1)^2 + f2*s012*s022*(m02 - 1)^2 - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1) + m1*m2*s112*s121*w1*w2 - f1*m01*s011*s021*w02*(m01 - 1) - f1*m01*s011*s022*w01*(m02 - 1) + f2*m01*m02*s012*s021*w01*w02)

#Myear
      #Metapopulation
      M1SM =expression( m1*s111*s122*w1*w2 - f1*s011*s021*(m01 - 1) - s111*s121*(m1 - 1) + f1*m01*s011*s022*w01*w02,0,
              0, m2*s112*s121*w1*w2 - f2*s012*s022*(m02 - 1) - s112*s122*(m2 - 1) + f2*m02*s012*s021*w01*w02)
      
      #Patch1  
      M1SP1 =expression(m1*s111*s122*w1*w2 - f1*s011*s021*(m01 - 1) - s111*s121*(m1 - 1) + f1*m01*s011*s022*w01*w02)
      #Patch2      
      M1SP2 =expression(m2*s112*s121*w1*w2 - f2*s012*s022*(m02 - 1) - s112*s122*(m2 - 1) + f2*m02*s012*s021*w01*w02)
#Mlife
      #Metapopulation
      M12SM=expression(s111*s121-f1*s011*s021*(m01-1),               -f1*s011*s021*(m01-1), 0,0,
                       f1*m01*s011*s022*w01*w02,s111*s122*w1*w2+f1*s011*s022*m01*w01*w02,    0,0,
                       
                       0,0,  s112*s122-f2*s012*s022*(m02-1),-f2*s012*s022*(m02-1),
                       0,0,            f2*s012*s021*m02*w01*w02,s112*s121*w1*w2+f2*m02*s012*s021*w01*w02)
      
      #Patch1
      M12SP1=expression(s111*s121-f1*s011*s021*(m01-1),                    -f1*s011*s021*(m01-1),
                        f1*m01*s011*s022*w01*w02,s111*s122*w1*w2+f1*m01*s011*s022*w01*w02)
      #Patch2      
      M12SP2=expression(s112*s122-f2*s012*s022*(m02-1),                    -f2*s012*s022*(m02-1),
                        f2*m02*s012*s021*w01*w02, s112*s121*w1*w2+f2*m02*s012*s021*w01*w02)
##### Long-live species   #####
#Mseason
      #Metapopulation
      MDLM=expression(0,0,f1*s011*s021*(m01 - 1)^2 + f1*m01*m02*s011*s022*w01*w02,0,0, - f2*m02*s012*s021*w01*(m01 - 1) - f2*m02*s012*s022*w02*(m02 - 1),          s111*s121*(m1 - 1)^2 + m1*m2*s111*s122*w1*w2,0,0, - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1),0,0,0,s111*s121*(m1 - 1)^2 + m1*m2*s111*s122*w1*w2,s111*s121*(m1 - 1)^2 + m1*m2*s111*s122*w1*w2,0, - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1),- m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1),                                                     0,                                                     0, - f1*m01*s011*s021*w02*(m01 - 1) - f1*m01*s011*s022*w01*(m02 - 1),                                                     0,                                                     0,           f2*s012*s022*(m02 - 1)^2 + f2*m01*m02*s012*s021*w01*w02, - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1),                                                     0,                                                                 0,          s112*s122*(m2 - 1)^2 + m1*m2*s112*s121*w1*w2,                                                     0,                                                                 0,                                                     0, - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1),             - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1),                                                     0,          s112*s122*(m2 - 1)^2 + m1*m2*s112*s121*w1*w2,                      s112*s122*(m2 - 1)^2 + m1*m2*s112*s121*w1*w2)
      #Patch1
      MDLP1 =expression(                                                                                                  0,                                                                                                  0, f1*s011*s021*(m01 - 1)^2 - f2*m02*s012*s021*w01*(m01 - 1) - f2*m02*s012*s022*w02*(m02 - 1) + f1*m01*m02*s011*s022*w01*w02, s111*s121*(m1 - 1)^2 - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1) + m1*m2*s111*s122*w1*w2,                                                                                                  0,                                                                                                                         0,                                                                                                  0, s111*s121*(m1 - 1)^2 - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1) + m1*m2*s111*s122*w1*w2,                        s111*s121*(m1 - 1)^2 - m2*s112*s121*w1*(m1 - 1) - m2*s112*s122*w2*(m2 - 1) + m1*m2*s111*s122*w1*w2)
      #Patch2      
      MDLP2= expression(                                                                                                  0,                                                                                                  0, f2*s012*s022*(m02 - 1)^2 - f1*m01*s011*s021*w02*(m01 - 1) - f1*m01*s011*s022*w01*(m02 - 1) + f2*m01*m02*s012*s021*w01*w02, s112*s122*(m2 - 1)^2 - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1) + m1*m2*s112*s121*w1*w2,                                                                                                  0,                                                                                                                         0,                                                                                                  0, s112*s122*(m2 - 1)^2 - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1) + m1*m2*s112*s121*w1*w2,                        s112*s122*(m2 - 1)^2 - m1*s111*s121*w2*(m1 - 1) - m1*s111*s122*w1*(m2 - 1) + m1*m2*s112*s121*w1*w2)
#Myear   
      #Metapopulation
      M1LM =expression(                                       0,                                       0, f1*m01*s011*s022*w01*w02 - f1*s011*s021*(m01 - 1),                                       0,                                       0,                                                 0, m1*s111*s122*w1*w2 - s111*s121*(m1 - 1),                                       0,                                                 0,                                       0,                                       0,                                                 0,                                       0, m1*s111*s122*w1*w2 - s111*s121*(m1 - 1),           m1*s111*s122*w1*w2 - s111*s121*(m1 - 1),                                       0,                                       0,                                                 0,                                       0,                                       0,                                                 0,                                       0,                                       0, f2*m02*s012*s021*w01*w02 - f2*s012*s022*(m02 - 1),                                       0,                                       0,                                                 0, m2*s112*s121*w1*w2 - s112*s122*(m2 - 1),                                       0,                                                 0,                                       0,                                       0,                                                 0,                                       0, m2*s112*s121*w1*w2 - s112*s122*(m2 - 1),           m2*s112*s121*w1*w2 - s112*s122*(m2 - 1))
      #Patch1
      M1LP1 =expression(                                       0,                                       0, f1*m01*s011*s022*w01*w02 - f1*s011*s021*(m01 - 1), m1*s111*s122*w1*w2 - s111*s121*(m1 - 1),                                       0,                                                 0,                                       0, m1*s111*s122*w1*w2 - s111*s121*(m1 - 1),           m1*s111*s122*w1*w2 - s111*s121*(m1 - 1))
      #Patch2      
      M1LP2= expression(                                       0,                                       0, f2*m02*s012*s021*w01*w02 - f2*s012*s022*(m02 - 1), m2*s112*s121*w1*w2 - s112*s122*(m2 - 1),                                       0,                                                 0,                                       0, m2*s112*s121*w1*w2 - s112*s122*(m2 - 1),           m2*s112*s121*w1*w2 - s112*s122*(m2 - 1))
#Mlife 
      #Metapopulation
      M12LM=expression(0,0,0,0,-f1*s011*s021*(m01-1),-f1*s011*s021*(m01-1),0,0,0,0,0,0,
                       0,0,0,0,f1*m01*s011*s022*w01*w02,f1*m01*s011*s022*w01*w02,0,0,0,0,0,0,
                       s111*s121,0,0,0,0,0,0,0,0,0,0,0,
                       0,s111*s122*w1*w2,0,0,0,0,0,0,0,0,0,0,
                       0,0,s111*s121,0,s111*s121,0,0,0,0,0,0,0,
                       0,0,0,s111*s122*w1*w2,0,s111*s122*w1*w2,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,-f2*s012*s022*(m02-1),-f2*s012*s022*(m02-1),
                       0,0,0,0,0,0,0,0,0,0,f2*m02*s012*s021*w01*w02,f2*m02*s012*s021*w01*w02,
                       0,0,0,0,0,0,s112*s122,0,0,0,0,0,
                       0,0,0,0,0,0,0,s112*s121*w1*w2,0,0,0,0,
                       0,0,0,0,0,0,0,0,s112*s122,0,s112*s122,0,
                       0,0,0,0,0,0,0,0,0,s112*s121*w1*w2,0,s112*s121*w1*w2)
      #Patch1
      M12LP1=expression(0,0,0,0,-f1*s011*s021*(m01-1),-f1*s011*s021*(m01-1),
                        0,0,0,0,f1*m01*s011*s022*w01*w02,f1*m01*s011*s022*w01*w02,
                        s111*s121,0,0,0,0,0,
                        0,s111*s122*w1*w2,0,0,0,0,
                        0,0,s111*s121,0,s111*s121,0,
                        0,0,0,s111*s122*w1*w2,0,s111*s122*w1*w2)
       #Patch2      
       M12LP2=expression(  0,0,0,0,-f2*s012*s022*(m02-1),-f2*s012*s022*(m02-1),
                          0,0,0,0,f2*m02*s012*s021*w01*w02,f2*m02*s012*s021*w01*w02,
                          s112*s122,0,0,0,0,0,
                          0,s112*s121*w1*w2,0,0,0,0,
                          0,0,s112*s122,0,s112*s122,0,
                          0,0,0,s112*s121*w1*w2,0,s112*s121*w1*w2)# Define mimatrix ####
# mimatrix allows my to call the corresponding matrix for  each model/patch/lhs in the coming fuctions
mimatrix<-array(c(
  "MDSP1",
  "MDSP2",
  "MDSM",
  
  "MDLP1",
  "MDLP2",
  "MDLM",
  
  "M1SP1",
  "M1SP2",
  "M1SM", 
  
  "M1LP1",
  "M1LP2",
  "M1LM",

  "M12SP1",
  "M12SP2",
  "M12SM", 
  
  "M12LP1",
  "M12LP2",
  "M12LM"
 
  
),dim=c(k+1,nLHS,nmodels))
#Get metamatrix with all parameters 
getmeta<-function(X){for (i in 1:dim(X)[1]){            # i denotes each possible combination of mobile fractions and mobile fraction survival
  if(epsilondir==0){ ## B>NB
    survivalmigration_0B   =as.numeric(X[i,2])  
    survivalmigration_0NB  =epsilon*survivalmigration_0B  
    cost0B<-1-survivalmigration_0B
    cost0NB<-1-survivalmigration_0NB
      
    survivalmigration_B=1-(cost0B*alphaam)
    survivalmigration_NB=1-(cost0NB*alphaam)
    
    tag<-paste("cw=",epsilon,"*cs",sep="")
  }
  else {## NB>B
    survivalmigration_0NB  =as.numeric(X[i,2])
    survivalmigration_0B   =epsilon*survivalmigration_0NB  
    cost0B<-1-survivalmigration_0B
    cost0NB<-1-survivalmigration_0NB
    survivalmigration_NB=1-(cost0NB*alphaam)
    survivalmigration_B=1-(cost0B*alphaam)
    tag<-paste("cs=",epsilon,"*cw",sep="")
  }
  
  if(gammadir==0){## P1>P2
    mig_ratio1=as.numeric(X[i,1])                                    #Extract mobile fraction
    mig_ratio01=alphamp*mig_ratio1
    mig_ratio2=gamma*mig_ratio1 
    mig_ratio02=alphamp*mig_ratio2
    tag2<-paste("m2=",gamma,"*m1",sep="")
  }
  else{## P1<P2
    mig_ratio2=as.numeric(X[i,1])
    mig_ratio02=alphamp*mig_ratio2
    mig_ratio1=gamma*mig_ratio2
    mig_ratio01=alphamp*mig_ratio1
    tag2<-paste("m1=",gamma,"*m2",sep="")
  }
  for( sp in 1:nLHS){                                                   #Loop through life-history strategies
    alphaa<-alphatag[sp]
    s<-sar11[sp]                                #Define baseline survival
    f<-far11[sp]                                #Define baseline fecundity
    ####defining the matrix ####
    m<-matrix(c( 
      f,                                         f,                     f,                       f,                      f,                                f,                                              f,     #1  Fecundity P1
      f,                                         f,                     f,                alphag*f,               alphag*f,                         alphag*f,                                        alphag*f,  #2  Fecundity P2
      
      s*alphaa,                           s*alphaa,              s*alphaa,                s*alphaa,               s*alphaa,                         s*alphaa,                                        s*alphaa,  #3  Newborn breedng survival  P1
      s*alphaa,                           s*alphaa,       alphag*s*alphaa,         alphag*s*alphaa, alphas*alphag*s*alphaa,    alphaR*alphas*alphag*s*alphaa,                                 alphag*s*alphaa,  #4  Newborn breedng survival P2
      s*alphaa,                    alphas*s*alphaa,              s*alphaa,         alphas*s*alphaa, alphas*alphag*s*alphaa,           alphas*alphag*s*alphaa,                alphaR*   alphas*alphag*s*alphaa,   #5  Newborn non breeding survival P1
      s*alphaa,                    alphas*s*alphaa,       alphag*s*alphaa,  alphas*alphag*s*alphaa,               s*alphaa,                  alphaR*s*alphaa,                          alphas*alphaR*s*alphaa,   #6  Newborn non breeding survival P2
      
      s,                                         s,                     s,                       s,                      s,                                s,                                               s,  #7  Adult breeding survival P1
      s,                                         s,              alphag*s,                alphag*s,        alphas*alphag*s,           alphaR*alphas*alphag*s,                                        alphag*s,  #8  Adult breeding survival P2
      s,                                  alphas*s,                     s,                alphas*s,        alphas*alphag*s,                  alphas*alphag*s,                          alphaR*alphas*alphag*s,   #9  Adult non-breeding survival P1
      s,                                  alphas*s,              alphag*s,         alphas*alphag*s,                      s,                         alphaR*s,                                 alphas*alphaR*s, #10  Adult non-breeding survival P2
      
      mig_ratio01,                       mig_ratio01,            mig_ratio01,              mig_ratio01,             mig_ratio01,                       mig_ratio01,                               mig_ratio01,    #11  Newborn Mobile fraction   P1
      mig_ratio02,                       mig_ratio02,            mig_ratio02,              mig_ratio02,             mig_ratio02,                       mig_ratio02,                               mig_ratio02,    #12  Newborn Mobile fraction   P2
      survivalmigration_0B ,     survivalmigration_0B,  survivalmigration_0B,     survivalmigration_0B,    survivalmigration_0B,              survivalmigration_0B,                      survivalmigration_0B,    #13  Newborn Breeding season Mobile fraction survival  P1
      survivalmigration_0NB,    survivalmigration_0NB, survivalmigration_0NB,    survivalmigration_0NB,   survivalmigration_0NB,             survivalmigration_0NB,                     survivalmigration_0NB,       #14 Newborn Non-Breeding season Mobile fraction survival  P1
      
      mig_ratio1,                       mig_ratio1,            mig_ratio1,              mig_ratio1,             mig_ratio1,                       mig_ratio1,                                      mig_ratio1,     #15  Mobile fraction   P1
      mig_ratio2,                       mig_ratio2,            mig_ratio2,              mig_ratio2,             mig_ratio2,                       mig_ratio2,                                      mig_ratio2,     #16  Mobile fraction   P2
      survivalmigration_B ,     survivalmigration_B,  survivalmigration_B,     survivalmigration_B,    survivalmigration_B,              survivalmigration_B,                             survivalmigration_B,   #17  Breeding season Mobile fraction survival  P1
      survivalmigration_NB,    survivalmigration_NB, survivalmigration_NB,    survivalmigration_NB,   survivalmigration_NB,             survivalmigration_NB,                            survivalmigration_NB   #18  Non-Breeding season Mobile fraction survival  P1
    )
    , nrow=nparameters,ncol=scenarios,byrow=TRUE)
    meta[,,i,sp]<-m    #Stores each combination in meta
  }                    #End for sp
}                      #End i
  return(meta)}        #End getmeta()





THEFUNCTION<-function(X,           # Vector containing combinations of mobile fraction and mobile fractionsurvival
                      k,            # Number of subpopulations
                      a,            # Age classes
                      meta,         # parameter values
                      scenario,
                      sim,
                      model,
                      strategy,      # LHS 
                      p
){
  mimat<-get(mimatrix[p,strategy,model])
  a=c(1,3)          # Age classes - age at first reproducion for each LHS
  states<-c(1,1,2)
  ntag<-c(1,1,2) 
  if (model==3){npara=nparameters-2}else{npara=nparameters}
  elas<-array(0,dim=c(npara,10,dim(X)[1]))
  LBD<-array(0,dim=c(1,9,dim(X)[1]))                     # array to store lambdas for each scenario and combinations of mobile fraction and mobile fraction survival
  mimat<-get(mimatrix[p,strategy,model])
  SSS<-array(0,dim=c(a[strategy]*ntag[p]*states[model],13,dim(X)[1]))           # array to store stable stage distribution for each scenario and combinations of mobile fraction and mobile fraction survival
  for (i in 1:dim(X)[1]){
       if (model==3){
        npara=nparameters-2
        tag<-c("f1","f2","s011","s012","s021","s022","s111","s112","s121","s122","m01","m02","w01","w02","w1","w2")
        vr<- list(
          f1=meta[1,scenario,i,strategy],          # Mean offspring
          f2=meta[2,scenario,i,strategy],          # Mean offspring
          s011=meta[3,scenario,i,strategy],        # Breeding  survival newborn  P1
          s012=meta[4,scenario,i,strategy],        # Breeding  survival newborn   P2 
          s021= meta[5,scenario,i,strategy],       # Non breeding  survival newborn   P1
          s022= meta[6,scenario,i,strategy],       # Non breeding  survival newborn    P2
          s111=meta[7,scenario,i,strategy],        # Breeding  survival adult    P1
          s112=meta[8,scenario,i,strategy],        # Breeding  survival adult   P2
          s121= meta[9,scenario,i,strategy],       # Non breeding  survival juvenile and adult   P1
          s122= meta[10,scenario,i,strategy],      # Non breeding  survival juvenile and adult  P2
          m01=meta[11,scenario,i,strategy],         # Breeding migration 
          m02=meta[12,scenario,i,strategy],         # Non Breeding migration
          w01=meta[13,scenario,i,strategy],        # Survival migration 
          w02=meta[14,scenario,i,strategy],         # Survival migration 
          w1=meta[17,scenario,i,strategy],        # Survival migration 
          w2=meta[18,scenario,i,strategy]         # Survival migration 

        ) }
       else{
        npara=nparameters
        tag<-c("f1","f2","s011","s012","s021","s022","s111","s112","s121","s122","m01","m02","w01","w02","m1","m2","w1","w2")
        vr<- list(
          f1=meta[1,scenario,i,strategy],          # Mean offspring
          f2=meta[2,scenario,i,strategy],          # Mean offspring
          s011=meta[3,scenario,i,strategy],        # Breeding  survival newborn  P1
          s012=meta[4,scenario,i,strategy],        # Breeding  survival newborn   P2 
          s021= meta[5,scenario,i,strategy],       # Non breeding  survival newborn   P1
          s022= meta[6,scenario,i,strategy],       # Non breeding  survival newborn    P2
          s111=meta[7,scenario,i,strategy],        # Breeding  survival adult    P1
          s112=meta[8,scenario,i,strategy],        # Breeding  survival adult   P2
          s121= meta[9,scenario,i,strategy],       # Non breeding  survival juvenile and adult   P1
          s122= meta[10,scenario,i,strategy],      # Non breeding  survival juvenile and adult  P2
          m01=meta[11,scenario,i,strategy],         # Breeding migration 
          m02=meta[12,scenario,i,strategy],         # Non Breeding migration
          w01=meta[13,scenario,i,strategy],        # Survival migration 
          w02=meta[14,scenario,i,strategy],         # Survival migration 
          m1=meta[15,scenario,i,strategy],         # Breeding migration 
          m2=meta[16,scenario,i,strategy],         # Non Breeding migration
          w1=meta[17,scenario,i,strategy],        # Survival migration 
          w2=meta[18,scenario,i,strategy]         # Survival migration 
        ) 
      }
       
##### GET A
  A<-matrix(sapply(mimat, eval,vr , NULL),   # A matrix: Evaluate the matrix  for each model (mimat) with the corresponding vital rates (vr) for each combination, scenario and life history strategy for P1,P2 and Metapopulation.
               nrow=a[strategy]*ntag[p]*states[model],
               ncol=a[strategy]*ntag[p]*states[model], 
               byrow=TRUE)
   
##### GET LAMBDA
     LBD[,1,i]<- max(Re(eigen(t(A))$values))        #Store lambda
     LBD[,3,i]<- meta[15,scenario,i,strategy]       #Store m
     LBD[,4,i]<- meta[17,scenario,i,strategy]       #Store w
     LBD[,5,i]<- scenario                           #Store scenario
     LBD[,6,i]<- sim                                #Store scenario
     LBD[,7,i]<- model                              #Store LHS
     LBD[,8,i]<- strategy                           #Store strategy
     LBD[,9,i]<- p                                  #Store simulation
     
##### GET Stable Stage Structure
 if(model%in%c(1,2)){SS<-rep(0,a[strategy]*states[p])}
 if(model%in%c(3)){SS<-rep(c(1,2),a[strategy]*states[p])}
 if(p%in%c(1,2)){sspatch<-p}
 if(p==3){       sspatch<-rep(c(1:k),each =a[strategy]*states[model] )}
     SSS[,1,i]<-stable.stage(A)            #Store stable stage structure
     SSS[,2,i]<-SS                                 #resident/migrant
     SSS[,3,i]<-rep(1:a[strategy],each =states[p]) #age                    
     SSS[,4,i]<-sspatch                            #patch
     if (model%in%c(1,2)){
       SSS[,5,i]<-meta[15,scenario,i,strategy]   #migrant fraction
     }else{
       SSS[,5,i]<-sum(SSS[which(SSS[,2,i]==2),1,i])  #migrant fraction
     }
     SSS[,6,i]<-max(Re(eigen(t(A))$values))        #Store lambda
     SSS[,7,i]<-meta[15,scenario,i,strategy]       #Store m
     SSS[,8,i]<-meta[17,scenario,i,strategy]       #Store w
     SSS[,9,i]<- scenario                    
     SSS[,10,i]<- sim                        #Input possible stage structure
     SSS[,11,i]<- model
     SSS[,12,i]<-strategy                       #Store p .
     SSS[,13,i]<-p

     if (model%in%c(1,2)){
       LBD[,2,i]<- meta[15,scenario,i,strategy]   #migrant fraction
     }else{
       LBD[,2,i]<- sum(SSS[which(SSS[,2,i]==2),1,i])  #Store real migrant fraction
     }
##### GET ELASTICITY    
    x<-vitalsens(mimat, vr)
    x[,4]<-meta[15,scenario,i,strategy]
    x[,5]<-meta[17,scenario,i,strategy]
    x[,6]<-scenario
    x[,7]<-sim
    x[,8]<-model
    x[,9]<-strategy
    x[,10]<-p
    
    names(x)<-c("value", "sensitivity","elasticity","m1","w1","scn","sim","model","strategy","patch")
    elas[,,i]<-as.matrix(x)
  }
  l<- LBD[,,1]
  ss<- SSS[,,1]
  elasticit<- elas[,,1]
  for (i in 2:dim(X)[1]){ #i loops through combinations
    l<-rbind(l,LBD[,,i])
    ss<-rbind(ss,SSS[,,i])
    elasticit<-rbind(elasticit,elas[,,i])
 }# end loop i
  lambd<-as.data.frame(l)
  ss<-as.data.frame(ss)
  elasticit<-as.data.frame(elasticit)
  elasticit[,11] <-rep(tag,dim(elasticit)[1]/length(tag))
  names(lambd)<-c("lambda","migrants", "m1","w1","scenario","sim","model","strategy","patch")
  names(ss)<-c("sss","migstrat","age","where","migratory fraction","lambda", "m1","w1","scenario","sim","model","strategy","patch")
  names(elasticit)<-c("value", "sensitivity","elasticity","m1","w1","scn","sim","model","strategy","patch","para")
  return(list(lambd,ss,elasticit))
  rm(lambd)
  rm(ss)
  rm(elasticit)
  closeAllConnections()
}


LAMBDA<-function(X,k,a,meta,scenario,sim){
  Model<-c("MD","M1","M12")          #Movement model labels
  Patch<-c("P1","P2","Meta")                #Movement model labels
  for (strategy in 1:nLHS){
    for(model in 1:nmodels){
      for (p in 1:(k+1)){
        assign(paste("l_",Model[model], "_", Patch[p], sep=""), as.data.frame(THEFUNCTION(X,k,a,meta,scenario,sim,model,strategy,p)[1]))
      }
    } 
    lambd<-rbind(l_MD_Meta,l_MD_P1,l_MD_P2,l_M1_Meta,l_M1_P1,l_M1_P2,l_M12_Meta,l_M12_P1,l_M12_P2)
    if(strategy==1){
      LHS="Shortlived"
    }
    if(strategy==2){
      LHS="Longlived"
    }
    assign(paste("lambda_",LHS, sep=""), lambd)
    
  }
  lambd<-rbind(lambda_Shortlived,lambda_Longlived)
  
  return(lambd) 
  rm(lambd) 
  
}

SSS<-function(X,meta,scenario,sim){
  for (strategy in 1:nLHS){
    for(model in 1:nmodels){
      for (p in 1:(k+1)){
        assign(paste("sssp_", p,"_m_", model,"_s_", strategy,"_sim_", sim, sep=""),THEFUNCTION(X,k,a,meta,scenario,sim,model,strategy,p)[2])
      }
    }
  }
  sssall<- do.call(rbind, mget(ls(pattern="sssp")))
  return(sssall)
  rm(list = ls(pattern="sssp"))
  rm(list = ls(pattern="sssall"))
}

ELASTICITY<-function(X,meta,scenario,sim){
  for (strategy in 1:nLHS){
    for(model in 1:nmodels){
      for (p in 1:(k+1)){
        assign(paste("elasticitiesp_", p,"_m_", model,"_s_", strategy,"_sim_", sim, sep=""),as.data.frame(THEFUNCTION(X,k,a,meta,scenario,sim,model,strategy,p)[3]))
      }
    }
  }
  elasticitiesall<- do.call(rbind, mget(ls(pattern="elasticitiesp")))
  return(elasticitiesall)
  rm(list = ls(pattern="elasticitiesp"))
  rm(list = ls(pattern="elasticitiesall"))
}


######################################
sar11=c(0.729,0.995)
far11=c(4.5,2)              #Fecundity      of shortlived and longlived in patch 1 during the breeding season. 
#Movement fraction and movement survival
migvalues <-c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999)#seq(0,1,0.1)      #Mobile fraction range of values
survvalues<-seq(0,1,0.1)#seq(0,1,0.5)       #Mobile fraction survival range of values

###Generate all possible combinations of mobile fraction and mobile fraction survival values. 
X<-expand.grid(migvalues,survvalues) 

#Dispersal
alphatag=c(0.63,0.775) #Age factor-           Newborn/adult survival ratio this is seasonal annual (alphatag^2) is 0.4 for shortlived species and 0.6 for long-lived.
aam=0.5          #Age factor-           Newborn/adult movement survival ratio   
aamp=0.95       #Age factor-           Newborn/adult movement survival ratio   
setwd("directory name")
data<-read.csv("sims.csv", sep=",", header=T)
scenario<-7
for (sim in 1:dim(data)[1]){#dim(data)[1]
  alphas=data[sim,1]                  #Seasonality gradient- Non-breeding/Breeding season ratio 
  alphag=data[sim,2]                 #Patch gradient-       P2/P1 ratio                       - alphag
  alphaR=data[sim,3]                  #Reciprocity-          Allows to define scenario 5.      - alphaR
  epsilon=data[sim,4]          #Ratio for migration survival between patches
  epsilondir=data[sim,5]       # 0 p1>p2, 1 p2>p1
  gamma=data[sim,6]            #Ratio for migration fraction between patches
  gammadir=data[sim,7]       # 0 p1>p2, 1 p2>p1
  activeagedifferenceinmigrationsurvival=data[sim,8]
  activepatchdifferenceinmigration=data[sim,9]
  d=data[sim,11]
  
  if(activeagedifferenceinmigrationsurvival==TRUE){
    alphaam=aam                  #Age factor-           Newborn/adult movement survival ratio   
  }else{    alphaam=1   }
  if(activepatchdifferenceinmigration==TRUE){
    alphamp=aamp                  #Age factor-           Newborn/adult movement survival ratio   
  }else{alphamp=1}
  meta<-array(0, dim=c(nparameters,scenarios,dim(X)[1],nLHS))
  meta<-getmeta(X)
  assign(paste("lambda_",sim, sep=""), LAMBDA(X,k,a,meta,scenario,sim))
  assign(paste("elasticity_",sim, sep=""), ELASTICITY(X,meta,scenario,sim))
}

lambda<-do.call(rbind, lapply( paste0("lambda_",1:dim(data)[1]) , get) )
elasticity<-do.call(rbind, lapply( paste0("elasticity_",1:dim(data)[1]) , get) )
getwd()
write.table(lambda,"lambda.csv",sep=",")
write.table(elasticity,"elasticity.csv",sep=",")



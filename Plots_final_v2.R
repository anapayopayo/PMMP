###  Seasonally mobile metapopulation lambda and elasticities plot funcitons       #### 
###  Author: Ana Payo Payo                                                         ####
###  Date:      07/08/2019                                                         ####

####    House keeping    #### 
#### Clean the workspace 

#this script produces plots from the article that are arranged differently that in the 
# manuscript. We exported this figures and modified their layout (but not their content)
# in a grapth design software. 
rm(list=ls())
ls()

setwd("directory name")

library(ggplot2)
library(viridis)
library(gridExtra)
library("RColorBrewer")
library(viridis)
e<-read.csv("elasticity.csv",sep=",")
l<-read.csv("lambda.csv",sep=",")
head(l)
l$dispersers<-0
l[which(l$model==1),10]<-l[which(l$model==1),3]*2*(1-l[which(l$model==1),3])
l[which(l$model==1),2]<-l[which(l$model==1),3]*l[which(l$model==1),3]
l$moving<-l$migrants+l$dispersers

#test that proves that simulations with different parameters provide different persistences see sims_art - Copy.csv for parameters
# e<-read.csv("elasticitytest.csv",sep=",")
# l<-read.csv("lambdatest.csv",sep=",")
l$modelname<-"NA"
names(l)
l[l$model==1,12]<-" Mseason"
l[l$model==2,12]<-"Mannual"
l[l$model==3,12]<-"Mlife"
l$LHSname<-"NA"

l[l$strategy==1,13]<-" Short"
l[l$strategy==2,13]<-"Long"
l$pname<-"NA"
l[l$patch==1,14]<-"P1"
l[l$patch==2,14]<-"P2"
l[l$patch==3,14]<-"Meta"


e$modelname<-"NA"
names(e)
e[e$model==1,12]<-" Mseason"
e[e$model==2,12]<-"Mannual"
e[e$model==3,12]<-"Mlife"
levels(e$model) <- c(" Mseason","Mannual","Mlife")

e$LHSname<-"NA"

e[e$strategy==1,13]<-"Short"
e[e$strategy==2,13]<-"Long"
levels(e$strategy) <- c(" Short","Long")

e$pname<-"NA"
e[e$patch==1,14]<-"P1"
e[e$patch==2,14]<-"P2"
e[e$patch==3,14]<-"Meta"

e<-e[!e$para=="d",]

e$plottag<-0
e$plottag[e$model==1&e$strategy==1] = 1
e$plottag[e$model==1&e$strategy==2] = 4
e$plottag[e$model==2&e$strategy==1] = 2
e$plottag[e$model==2&e$strategy==2] = 5
e$plottag[e$model==3&e$strategy==1] = 3
e$plottag[e$model==3&e$strategy==2] = 6
e$plottag =factor(e$plottag)
levels(e$plottag) <- c(" Mseason","Mannual","Mlife","  Mseason"," Mannual"," Mlife")

e$para =factor(e$para)
e$para = factor(e$para,levels(e$para)[c(1:2,7:18,3:6)])
levels(e$para) <- c("f[1]", " f[2]", expression("s"["011"]), expression(" s"["012"]), expression(" s"["021"] ), expression("s"["022"]), " s[111]", " s[112]", " s[121]", " s[122]", expression(" w"["01"]), expression(" w"["02"]), " w[1]", " w[2] ", expression("m"["01"]), expression(" m"["02"]), " m[1]", " m[2]")



tag <-c("=",">","=") #survival
tag2<-c("=","=",">")#migration
dv<-c(0,0.05,0.05,0.05,0.05,0.05,0,0)#migration


      
###########
      for(sim in 1:4){
          datas<-l
          datas$mtick<-datas$moving
          datas<-datas[datas$sim==sim, ]
          datas<-datas[datas$patch==2, ]
          #datas<-datas[datas$patch==3, ]

          s <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)
 
          sc <- c("deeppink4","violetred4","violetred4","white","turquoise4","deepskyblue4","deepskyblue4")
          
          p_l<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
            geom_raster(data=datas, aes(fill=lambda),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
            scale_fill_gradientn(limits = c(0,2),
                                 colours=  sc,
                                 labels = s,
                                 breaks=s,
                                 name = expression(lambda))+
            theme(panel.background = element_blank(),
                  panel.spacing = unit(0, "lines"), 
                  strip.text.x = element_text(colour = "black", face = "bold"),
                  strip.background = element_rect(colour = "white", fill = "NA"),
                  # legend.position = "none",
                  axis.text.x = element_text(angle = 90))+
            facet_grid(strategy~modelname)+
            labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
                 y = ~ paste("w"),
                 x = ~ paste("m")
            )
           assign(paste("plot_original_",sim, sep="_"),p_l )
          
          s <- seq(0,1,by=0.1)
          sc <-viridis(11)
          table(datas$d)
          p_l_mtick<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
            geom_raster(data=datas, aes(fill=mtick),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
            scale_fill_gradientn(limits = c(0,1),
                                 colours=  sc,
                                 labels = s,
                                 breaks=s,
                                 name = paste("m"))+
            theme(panel.background = element_blank(),
                  panel.spacing = unit(0, "lines"), 
                  strip.text.x = element_text(colour = "black", face = "bold"),
                  strip.background = element_rect(colour = "white", fill = "NA"),
                  # legend.position = "none",
                  axis.text.x = element_text(angle = 90))+
            facet_grid(strategy~modelname,labeller = label_parsed)+
            labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
                 y = ~ paste("w"),
                 x = ~ paste("m")
            )
          assign(paste("plot",sim, sep="_"),p_l_mtick )
          
          datas$mtick<-datas$migrants
          p_l_mtick2<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
            geom_raster(data=datas, aes(fill=mtick),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
            scale_fill_gradientn(limits = c(0,1),
                                 colours=  sc,
                                 labels = s,
                                 breaks=s,
                                 name = paste("m"))+
            theme(panel.background = element_blank(),
                  panel.spacing = unit(0, "lines"), 
                  strip.text.x = element_text(colour = "black", face = "bold"),
                  strip.background = element_rect(colour = "white", fill = "NA"),
                  # legend.position = "none",
                  axis.text.x = element_text(angle = 90))+
            facet_grid(strategy~modelname,labeller = label_parsed)+
            labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
                 y = ~ paste("w"),
                 x = ~ paste("m")
            )
          assign(paste("plot",sim, sep="__"),p_l_mtick2 )
          
          
          datas$mtick<-datas$dispersers
          p_l_mtick3<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
            geom_raster(data=datas, aes(fill=mtick),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
            scale_fill_gradientn(limits = c(0,1),
                                 colours=  sc,
                                 labels = s,
                                 breaks=s,
                                 name = paste("m"))+
            theme(panel.background = element_blank(),
                  panel.spacing = unit(0, "lines"), 
                  strip.text.x = element_text(colour = "black", face = "bold"),
                  strip.background = element_rect(colour = "white", fill = "NA"),
                  # legend.position = "none",
                  axis.text.x = element_text(angle = 90))+
            facet_grid(strategy~modelname,labeller = label_parsed)+
            labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
                 y = ~ paste("w"),
                 x = ~ paste("m")
            )
          assign(paste("plot",sim, sep="___"),p_l_mtick3 )
          
      }
      grid.arrange( plot_original__1,plot_1)

      grid.arrange(plot_1,plot__1,plot___1)
      
      grid.arrange( plot_original__2,
                    plot_2)
      grid.arrange( plot_original__3,
                    plot_3)
      grid.arrange( plot_original__4,
                    plot_4)

   rm(list=ls(pattern="plot_*"))
  
   
   
   
   ###########
   for(sim in 1:4){
     datas<-l
     datas$mtick<-datas$moving
     datas<-datas[datas$sim==sim, ]
     datas<-datas[datas$patch==3, ]
     s <- seq(0,1,by=0.1)
     sc <-viridis(11)
     table(datas$d)
     p_l_mtick<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
       geom_raster(data=datas, aes(fill=mtick),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
       scale_fill_gradientn(limits = c(0,1),
                            colours=  sc,
                            labels = s,
                            breaks=s,
                            name = paste("m"))+
       theme(panel.background = element_blank(),
             panel.spacing = unit(0, "lines"), 
             strip.text.x = element_text(colour = "black", face = "bold"),
             strip.background = element_rect(colour = "white", fill = "NA"),
             # legend.position = "none",
             axis.text.x = element_text(angle = 90))+
       facet_grid(strategy~modelname,labeller = label_parsed)+
       labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
            y = ~ paste("w"),
            x = ~ paste("m")
       )
     assign(paste("plot",sim, sep="_"),p_l_mtick )
     
     datas$mtick<-datas$migrants
     p_l_mtick2<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
       geom_raster(data=datas, aes(fill=mtick),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
       scale_fill_gradientn(limits = c(0,1),
                            colours=  sc,
                            labels = s,
                            breaks=s,
                            name = paste("m"))+
       theme(panel.background = element_blank(),
             panel.spacing = unit(0, "lines"), 
             strip.text.x = element_text(colour = "black", face = "bold"),
             strip.background = element_rect(colour = "white", fill = "NA"),
             # legend.position = "none",
             axis.text.x = element_text(angle = 90))+
       facet_grid(strategy~modelname,labeller = label_parsed)+
       labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
            y = ~ paste("w"),
            x = ~ paste("m")
       )
     assign(paste("plot",sim, sep="__"),p_l_mtick2 )
     
     
     datas$mtick<-datas$dispersers
     p_l_mtick3<-ggplot(aes(x=m1, y=w1, z=lambda), data = datas) +              #Plot  lambdas
       geom_raster(data=datas, aes(fill=mtick),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
       scale_fill_gradientn(limits = c(0,1),
                            colours=  sc,
                            labels = s,
                            breaks=s,
                            name = paste("m"))+
       theme(panel.background = element_blank(),
             panel.spacing = unit(0, "lines"), 
             strip.text.x = element_text(colour = "black", face = "bold"),
             strip.background = element_rect(colour = "white", fill = "NA"),
             # legend.position = "none",
             axis.text.x = element_text(angle = 90))+
       facet_grid(strategy~modelname,labeller = label_parsed)+
       labs(subtitle=  paste("wB",tag[sim],"wNB"," m",tag2[sim],"m2",sep=""),
            y = ~ paste("w"),
            x = ~ paste("m")
       )
     assign(paste("plot",sim, sep="___"),p_l_mtick3 )
     
   }

   grid.arrange(plot_1,plot__1,plot___1)

   
   rm(list=ls(pattern="plot_*"))
   
   
   
   
    for (sim in 1:4){
     d  <-e[e$sim==sim, ]
     d  <-d[d$patch==3, ] ## to plot patch one change to 1 to plot patch two change to 2
     s  <- seq(-2,1,by=0.4)
     sc <- brewer.pal(n = 7, name = 'RdBu')
     dd <-d
     pd <- position_dodge(0.1) # move them .1 to the left and right
     p  <-ggplot(aes(x=m1, y=w1, z=elasticity), data = dd) +              #Plot  lambdas
          geom_raster(data=dd, aes(fill=elasticity),interpolate = TRUE) +   #Interpolate between points to have a smooth representation
          scale_fill_gradientn(limits = c(-2,1),
                            colours=  sc,
                            labels = s,
                            breaks=s,
                            name = paste("m"),
                            na.value = "deeppink4")+
          theme(panel.background = element_blank(),
               panel.spacing = unit(0, "lines"), 
               strip.text.x = element_text(colour = "black", face = "bold"),
               strip.background = element_rect(colour = "white", fill = "NA"),
               axis.ticks=element_blank(),# legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               
               strip.text.y = element_text(angle=0))+
          facet_grid(plottag~para,labeller = label_parsed)+
          labs(y = paste("w"),
              x = ~ paste("m")
          )
        assign(paste("plot_sim",sim, sep="_"),p )
    }
  
   plot_sim_1
   plot_sim_2
   plot_sim_3
   plot_sim_4
   
 
  c<- rev(viridis(11))
  
   for (sim in 1:4){
     e$w1<- as.factor(e$w1)
     d<-e
     d<-d[d$sim==sim, ]
     d<-d[d$patch==3, ]## to plot patch one change to 1 to plot patch two change to 2
     
     dd<-d
     pd <- position_dodge(0.1) # move them .1 to the left and right
     p<-ggplot(d, aes(x=m1)) + 
       geom_hline(yintercept = 0,col="grey")+
       geom_line(  aes(y=elasticity,group = w1,colour=w1),position=pd) +
       scale_color_manual(name="w",values=c)+
      ylim(-1,  1)+
       theme( 
         axis.text.x = element_text(colour="darkgrey",angle = 90,size = 8, hjust = 1),
         axis.text.y = element_text(colour="darkgrey",size = 8, hjust = 1),
         strip.text.y = element_text(angle = 0, size = 12,colour = "white", face = "bold"),
         strip.text.x = element_text(colour = "white",size = 12, face = "bold"),
         strip.background = element_rect(colour = "deepskyblue4", fill = "deepskyblue4"),
         panel.background = element_rect(colour="white",fill="white"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
       facet_grid(plottag~para,labeller = label_parsed)+
       labs(y = paste("w"),
            x = ~ paste("m")
       )
     assign(paste("plot_sim",sim, sep="_"),p )

   }

   plot_sim_1
   plot_sim_2
   plot_sim_3
   plot_sim_4
   
     
#====================================================================================================#

if (!require(c("bipartite",'dplyr'))) install.packages(c("bipartite",'dplyr'))
library(bipartite);library(dplyr)

#====================================================================================================#

#=================#
# READING DATASET #
#=================#

data<-read.table("Plant-hummingbird_interaction_PNSC.txt", header=T)

#===============================#
# CREATING INTERACTION MATRICES #
#===============================#

#---General interaction matrix---#
general_matrix<-tapply(data$Interaction, data[,3:2],sum)
general_matrix[is.na(general_matrix)] <- 0

#---Rupestrian grassland interaction matrix---#
rupestrian_matrix<-frame2webs(data, varnames = c("Plant","Hummingbird", "Phyto",
                                         "Interaction"), type.out="list", 
                      emptylist = F)[["rupestrian_field"]]

#---Riparian forest interaction matrix---#
forest_matrix<-frame2webs(data, varnames = c("Plant","Hummingbird", "Phyto",
                                         "Interaction"), type.out="list", 
                      emptylist = F)[["riparian_forest"]]

#====================================================================================================#

#=========================#
# NETWORK LEVEL ANALYSIS* #
#=========================#
#*Due to null models, z-scores vary with each run


#---General hummingbird-plant network---#

#-- Modulariry
module_general_network<-computeModules(general_matrix, method = "Beckett") 
plotModuleWeb(module_general_network) #Q=0.3633; 4 modules

nulls_general<-nullmodel(general_matrix, N=1000, method = "vaznull")

#-Testing for significance
modules.nulls_general<-sapply(nullmodel(general_matrix, N=100, method = "vaznull")
                         , computeModules, method="Beckett")

like.nulls_general <- sapply(modules.nulls_general, function(x) x@likelihood)
(z_modules_general <- (module_general_network@likelihood - 
                         mean(like.nulls_general))/sd(like.nulls_general))
#z-score=17.58

#-Pvalue
2*pnorm(-abs(z_modules_general))

#--Network specialization H2 
networklevel(general_matrix, index = 'H2') #0.4396

#--Testing for significance
(z_h2_general<-(networklevel(general_matrix, index = "H2")- 
                      mean(unlist(sapply(nulls_general, networklevel, 
                      index="H2"))))/sd(unlist(sapply(nulls_general, 
                                                      networklevel,index="H2"))))
#z-score=17.12

#-Pvalue
2*pnorm(-abs(z_h2_general))


#---Network nestedness wNODF
networklevel(general_matrix, index = 'weighted NODF') #33.488

#--Testing for significance
(z_wNODF_general<-(networklevel(general_matrix, index = "weighted NODF")- 
                  mean(unlist(sapply(nulls_general, networklevel, 
                                     index="weighted NODF"))))
        /sd(unlist(sapply(nulls_general, networklevel,index="weighted NODF"))))

#z-score=-3.72

#--Pvalue
2*pnorm(-abs(z_wNODF_general))

#=============================================#  
#=Rupestrian fields hummingbird-plant network=#
#=============================================#  

#--Modularity

module_rupestrian_network<-computeModules(rupestrian_matrix, method = "Beckett")
plotModuleWeb(module_rupestrian_network)
#Q=0.38; 3 modules

#-Testing for significance

nulls_rupestrian<-nullmodel(rupestrian_matrix, N=1000, method = "vaznull")

modules.nulls_rupestrian<-sapply(nullmodel(rupestrian_matrix, N=100, 
                        method = "vaznull"), computeModules, method="Beckett")

like.nulls_rupestrian <- sapply(modules.nulls_rupestrian, function(x) x@likelihood)
(z_modules_rupestrian <- (module_rupestrian_network@likelihood - 
                         mean(like.nulls_rupestrian))/sd(like.nulls_rupestrian))
#z-score=9.85

#-Pvalue
2*pnorm(-abs(z_modules_rupestrian))


#--Network specialization H2 
networklevel(rupestrian_matrix, index = 'H2') #0.55

#--Testing for significance
(z_h2_rupestrian<-(networklevel(rupestrian_matrix, index = "H2")- 
                  mean(unlist(sapply(nulls_rupestrian, networklevel, 
                               index="H2"))))/sd(unlist(sapply(nulls_rupestrian, 
                                                    networklevel,index="H2"))))
#z-score=9.29

#-Pvalue
2*pnorm(-abs(z_h2_rupestrian))

#---Network nestedness wNODF
networklevel(rupestrian_matrix, index = 'weighted NODF') #44.636

#--Testing for significance
(z_wNODF_rupestrian<-(networklevel(rupestrian_matrix, index = "weighted NODF")- 
                     mean(unlist(sapply(nulls_rupestrian, networklevel, 
                                        index="weighted NODF"))))
  /sd(unlist(sapply(nulls_rupestrian, networklevel,index="weighted NODF"))))

#z-score=-1.23

#-Pvalue
2*pnorm(-abs(z_wNODF_rupestrian))

#===============================================#  
#===Riparian forest hummingbird-plant network===#
#===============================================# 

#--Modularity

module_forest_network<-computeModules(forest_matrix, method = "Beckett")
plotModuleWeb(module_forest_network)
#Q=0.346; 4 modules

#-Testing for significance

nulls_forest<-nullmodel(forest_matrix, N=1000, method = "vaznull")

modules.nulls_forest<-sapply(nullmodel(forest_matrix, N=100, 
                        method = "vaznull"), computeModules, method="Beckett")

like.nulls_forest<- sapply(modules.nulls_forest, function(x) x@likelihood)
(z_modules_forest <- (module_forest_network@likelihood - 
                         mean(like.nulls_forest))/sd(like.nulls_forest))

#z-score=15.86
#-Pvalue
2*pnorm(-abs(z_modules_forest))

#--Network specialization H2 
networklevel(forest_matrix, index = 'H2') #0.41

#--Testing for significance
(z_h2_forest<-(networklevel(forest_matrix, index = "H2")- 
                     mean(unlist(sapply(nulls_forest, networklevel, 
                                  index="H2"))))/sd(unlist(sapply(nulls_forest, 
                                                    networklevel,index="H2"))))
#z-score=15.01
#-Pvalue
2*pnorm(-abs(z_h2_forest))

#---Network nestedness wNODF
networklevel(forest_matrix, index = 'weighted NODF') #24.02

#--Testing for significance
(z_wNODF_forest<-(networklevel(forest_matrix, index = "weighted NODF")- 
                        mean(unlist(sapply(nulls_forest, networklevel, 
                                           index="weighted NODF"))))
  /sd(unlist(sapply(nulls_forest, networklevel,index="weighted NODF"))))
#z-score=-5.12

#-Pvalue
2*pnorm(-abs(z_wNODF_forest))


#====================================================================================================#

#============================#  
#===SPECIES LEVEL ANALYSIS===#
#============================#

specieslevel_general<-specieslevel(general_matrix, index = c("d","species strength"), level="higher")
specieslevel_general<-data.frame(names=row.names(specieslevel_general), specieslevel_general)

specieslevel_rupestrian<-specieslevel(rupestrian_matrix, index = c("d","species strength"), level="higher")
specieslevel_rupestrian<-data.frame(names=row.names(specieslevel_rupestrian), specieslevel_rupestrian)

specieslevel_forest<-specieslevel(forest_matrix, index = c("d","species strength"), level="higher")
specieslevel_forest<-data.frame(names=row.names(specieslevel_forest), specieslevel_forest)

#-Combining species level data from all networks

df_list<-list(specieslevel_general,specieslevel_rupestrian, specieslevel_forest)
specieslevel_all<-Reduce(function(x, y) merge(x, y, all=TRUE, by="names"), 
                         df_list, accumulate=FALSE)
colnames(specieslevel_all)<-c("Hummingbirds", "ss_general", "d_general", 
                              "ss_rupestrian","d_rupestrian", "ss_forest", 
                              "d_forest")

#--Testing whether specieslevel index differ between phytophysiognomies

normality_test<-do.call(rbind, lapply(specieslevel_all[,4:7], 
                      function(x) shapiro.test(x)[c("statistic", "p.value")]))
#d index has a normal distribution, species strength index does not. 

wilcox.test(specieslevel_all$ss_rupestrian, specieslevel_all$ss_forest)
#W=23, p-value=0.4079

t.test(specieslevel_all$d_rupestrian, specieslevel_all$d_forest)
#t= -0.64989, df=11.541, p-value=0.5285

#====================================================================================================#


#===================#
# PLOTTING NETWORKS #
#===================#

#---Adjusting taxa names
colnames(rupestrian_matrix)<-gsub('_', '. ', colnames(rupestrian_matrix))
row.names(rupestrian_matrix) <- gsub('_', '. ', row.names(rupestrian_matrix))
row.names(rupestrian_matrix)<-  gsub('Lessingianthus. sp.', 
                                     'Lessingianthus sp.', 
                                     row.names(rupestrian_matrix))
row.names(rupestrian_matrix)<-  gsub('Lychnophora. sp', 
                                     'Lychnophora sp.', 
                                     row.names(rupestrian_matrix))
row.names(rupestrian_matrix)<-  gsub('Eremanthus. sp', 
                                     'Eremanthus sp.', 
                                     row.names(rupestrian_matrix))
row.names(forest_matrix) <- gsub('Asteraceae. sp.', 'Asteraceae sp.', 
                                 row.names(forest_matrix))
colnames(forest_matrix)<-gsub('_', '. ', colnames(forest_matrix))


#---Rupestrian fields hummingbird-plant network---#
rupestrian_network<-plotweb(rupestrian_matrix, text.rot = 90, 
                col.interaction = "#e0dede",  col.high = "#757474", 
                col.low = "#000000", labsize = 1.1, low.y = 0.7, high.y = 1.4)


#---Riparian forest hummingbird-plant network---#
forest_network<-plotweb(forest_matrix, text.rot = 90, 
                 col.interaction = "#e0dede",  col.high = "#757474", 
                 col.low = "#000000", labsize = 1.1, low.y = 0.7, high.y = 1.4)


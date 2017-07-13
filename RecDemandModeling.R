#R code for aspatial and spatial Poisson models of freshwater aquatic species richness as reported in Table 1 of Davis & Darling (2017) 
##load packages
library(pander)
library(spdep)


##read in data
all_Aug8_resfits <- read.csv("//Aa.ad.epa.gov/ord/RTP/Users/A-D/adavis03/Net MyDocuments/Projects/EPA_Invasives/popVsRecDemd/data/ManuscriptFinalData/all_Aug8_resfits.csv")
attach(all_Aug8_resfits)
all<-all_Aug8_resfits 

##Standarize predictors
FFSUMarea<-FFSUM/AreaSqKm
FFSUMarealogStd<-scale(log(FFSUMarea+0.5),scale=TRUE,center=TRUE)
popdens<-DasyPop/AreaSqKm
popdenslogStd<-scale(log(popdens+0.5),scale=TRUE,center=TRUE)
FFSUMarealogStd1<-FFSUMarealogStd[,1]
popdenslogStd1<-popdenslogStd[,1]

##Derive aspatial models for full data and temporal subsets to examine effect of temporally  mismatched data

### models using full data set
m1<-glm(new_ar~FFSUMarealogStd1,offset=perAqHab,family=quasipoisson,data=all)
m1a<-glm(new_pr~FFSUMarealogStd1,offset=perAqHab,family=quasipoisson,data=all)

m2<-glm(new_ar~popdenslogStd1,offset=perAqHab,family=quasipoisson,data=all)
m2a<-glm(new_pr~popdenslogStd1,offset=perAqHab,family=quasipoisson,data=all)

### models using the 2005+ temporal subset
m3<-glm(new_ar05~FFSUMarealogStd1,offset=perAqHab,family=quasipoisson,data=all)
m3a<-glm(new_pr05~FFSUMarealogStd1,offset=perAqHab,family=quasipoisson,data=all)

m4<-glm(new_ar05~popdenslogStd1,offset=perAqHab,family=quasipoisson,data=all)
m4a<-glm(new_pr05~popdenslogStd1,offset=perAqHab,family=quasipoisson,data=all)

##Obtain 

aspatial_list<-list(m3,m3a,m4,m4a)
aspatial_CI_list1<-lapply(aspatial_list1, function(x) exp(confint(x,2)))
aspatial_CI_list<-lapply(aspatial_list, function(x) exp(confint(x,2)))
ci_mat<-do.call(rbind,aspatial_CI_list)
ci_mat2<-do.call(rbind,aspatial_CI_list1)

Predictor<-list("Freshwater fishing","Freshwater fishing","Population density","Population density")
Response<-list("Animal Richness","Plant Richness")
ci_mat1<-cbind(Response,Predictor,ci_mat)
ci_mat3<-cbind(Response,Predictor,ci_mat2)


panderOptions("digits",2)
pander(ci_mat3, caption="Table 1a: 95% confidence intervals for freshwater fishing and population density using aspatial GLMs with the full dataset")
pander(ci_mat1, caption="Table 1b: 95% confidence intervals for freshwater fishing and population density using aspatial GLMs with the 2005+ temporal subset")

####Not run
###1. create 2nd order queen contiguity matrix in the open source GeoDa available from: https://geodacenter.github.io/
###2. Create negibhor list in spdep
#queen2_listw<-nb2listw(huc8.nb2, glist=NULL, style="B", zero.policy=TRUE)
#3.assess spatial autocorrelation in aspatial and spatial models
#moran.mc(residuals(yourmodel,type="deviance"), listw=queen2_listw, 999, zero.policy=TRUE, alternative="greater")
#Apply moran eigenvector filtering function in spdep package to models using full data and 2005+ subset
##(Note these models can be lenghthy to run and the derived eigenvectors for each model are included with the data, see below)
# me1<-ME(new_ar~FFSUMarealogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me1a<-ME(new_pr~FFSUMarealogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me1b<-ME(new_tot~FFSUMarealogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)


#me2<-ME(new_ar~popdenslogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me2a<-ME(new_pr~popdenslogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me2b<-ME(new_tot~popdenslogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)


#me3<-ME(new_ar05~FFSUMarealogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me3a<-ME(new_pr05~FFSUMarealogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me3b<-ME(new_tot05~FFSUMarealogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)


#me4<-ME(new_ar05~popdenslogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me4a<-ME(new_pr05~popdenslogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)
#me4b<-ME(new_tot05~popdenslogStd1,offset=perAqHab,family="poisson",data=all,listw=queen2_listw,alpha=0.2,stdev=TRUE,verbose=FALSE)

###Run spatial models using pre-derived eigenvectors

#Extract eigenvectors corresponding to each model
me1afit<-all[,grep("me1afit",names(all),value=FALSE)]
me1fit<-all[,grep("me1fit",names(all),value=FALSE)]
me2fit<-all[,grep("me2fit",names(all),value=FALSE)]
me2afit<-all[,grep("me2afit",names(all),value=FALSE)]
me3fit<-all[,grep("me3fit",names(all),value=FALSE)]
me3afit<-all[,grep("me3afit",names(all),value=FALSE)]
me4fit<-all[,grep("me4fit",names(all),value=FALSE)]
me4afit<-all[,grep("me4afit",names(all),value=FALSE)]

#assign model parameters and eigenvectors for each  model to individual dataframes
sm1.df<-data.frame(cbind(new_ar,FFSUMarealogStd1,me1fit,perAqHab))
sm1a.df<-data.frame(cbind(new_pr,FFSUMarealogStd1,me1afit,perAqHab))
sm2.df<-data.frame(cbind(new_ar,popdenslogStd1,me2fit,perAqHab))
sm2a.df<-data.frame(cbind(new_pr,popdenslogStd1,me2afit,perAqHab))
sm3.df<-data.frame(cbind(new_ar05,FFSUMarealogStd1,me3fit,perAqHab))
sm3a.df<-data.frame(cbind(new_pr05,FFSUMarealogStd1,me3afit,perAqHab))
sm4.df<-data.frame(cbind(new_ar05,popdenslogStd1,me4fit,perAqHab))
sm4a.df<-data.frame(cbind(new_pr05,popdenslogStd1,me4afit,perAqHab))

#spatial models with pre-derived eigenvectors
sm1<-glm(new_ar~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm1.df)
sm1a<-glm(new_pr~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm1a.df)
sm2<-glm(new_ar~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm2.df)
sm2a<-glm(new_pr~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm2a.df)
sm3<-glm(new_ar05~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm3.df)
sm3a<-glm(new_pr05~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm3a.df)
sm4<-glm(new_ar05~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm4.df)
sm4a<-glm(new_pr05~.-perAqHab,offset=perAqHab,family=quasipoisson,data=sm4a.df)


spatial_listFull<-list(sm1,sm1a,sm2,sm2a)
spatial_listSub<-list(sm3,sm3a,sm4,sm4a)
spatial_CI_listFull<-lapply(spatial_listFull, function(x) exp(confint(x,2)))
spatial_CI_listSub<-lapply(spatial_listSub, function(x) exp(confint(x,2)))

ci_matFull<-do.call(rbind,spatial_CI_listFull)
ci_matSub<-do.call(rbind,spatial_CI_listSub)

Predictor<-list("Freshwater fishing","Freshwater fishing","Population density","Population density")
Response<-list("Animal Richness","Plant Richness")
ci_matFull1<-cbind(Response,Predictor,ci_matFull)
ci_matSub1<-cbind(Response,Predictor,ci_matSub)


panderOptions("digits",2)
pander(ci_matFull1, caption="Table 1c: 95% confidence intervals for freshwater fishing and population density using spatial Poisson models with the full dataset")
pander(ci_matSub1, caption="Table 1d: 95% confidence intervals for freshwater fishing and population density using spatial Poisson models with the 2005+ temporal subset")

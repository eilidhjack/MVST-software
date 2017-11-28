#####R Script for Paper######
#####Read in data, source files, run model######
#####Eilidh Jack, 21/11/2017######

###Load library's
library(Rcpp)
library(RcppArmadillo)
library(CARBayes)
library(shapefiles)
library(mvtnorm)
library(truncdist)
library(MCMCpack)




library(spdep)
library(coda)
library(MASS)
library(classInt)
library(plyr)
library(colorspace)



##Read in the Data##
ThreeDisease <- read.csv("Ten Years Three Disease - IG.csv", header=T)

##Load in W matrix which has been altered to connect islands to the mainland
load("W.RData")


##Create Health Board factor for mean centering by HB
CHD<-ThreeDisease[which(ThreeDisease$Disease=="CHD"),]
CHD.oneyear<-CHD[which(CHD$Year=="2012"),]
HB<-factor(CHD.oneyear$HealthBoard)
##Create X matrix to mean centre the phi's
X.phi<-model.matrix(~HB-1)

##Create Health Board factor for temporal trends
HB.2<-factor(CHD$HealthBoardYear)
X.HB<-model.matrix(~HB.2-1)

##Create HB factor to add to formula
HB.run<-factor(ThreeDisease$HealthBoardYear)

##source functions
source("MVLerouxAR.R")
sourceCpp("Leroux.phi.cpp")
sourceCpp("kronecker.cpp")
sourceCpp("AR1.HB.cpp")

ThreeDisease<-space.time(ThreeDisease$Yi~ThreeDisease$JSA.mean+ThreeDisease$Percent.Asian.mean.log+
                                 ThreeDisease$Percent.Black.mean.log+factor(ThreeDisease$Urban)+offset(log(ThreeDisease$Ei)) + HB.run,
                               samples=150000, burnin=50000, prior.mean.beta=NULL, prior.var.beta=NULL,  
                               progress.bar=T, thin=5, block=5, d=3,W=W, tp=10, S=14, X.HB=X.HB, X.phi=X.phi) 

##Arguments
##formula - Insert formula for the covariate part of the model
##samples - number of samples to be drawn
##burnin - number of samples to be removed from the chain
##prior.mean.beta - a vector for the prior mean of each beta (defaults to a vector of 0's)
##prior.var.beta - a vector for the prior variances for each beta (defaults to a vector of 100's)
##progress.bar - whether progress bar appears or not, default = T
##thin - value to thin by, default = 1
##block - number of parameters per block, default=5
##d - number of diseases in data
##W - W neighborhood matrix
##tp - number of time points
##S - number of Health Boards
##X.HB - Matrix of HB by year
##X.phi - Design matrix of factors for HB to allow for phi to be mean centered by HB
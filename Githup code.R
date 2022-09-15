#@mfsallam, 8SEP2022, mosquito data from Florida#
install.packages("devtools")
library(devtools)
devtools::install_github("nicholasjclark/MRFcov", force = TRUE)
library(MRFcov)
library(MRFcov)
library(devtools)
library(glmnet)
library(ggplot2)
library(igraph)

setwd("Drive name")
getwd()

######2020#####
########################
# load data
#5km
dat5k<-read.csv("file name")
dim(dat5k)
#10km
dat10k<-read.csv("file name")
dim(dat10k)

##load coordinates 
#5km
coords<-dat5k[,2:3]
View(coords)
#10km
coords10k<-dat10k[,2:3]
View(coords10k)

##preparing number of nodes 
#5km
mos<-dat5k[, c(4:33)]
mos<-mos[, colSums(mos != 0) > 0]##Select mosquito columns that their summations not equal to zero and greater than zero
dim(mos)
names(mos)
#10km
mos10k<-dat10k[, c(4:33)]
dim(mos10k)
names(mos10k)

##preparing n_covariates
#2020
#5km
n_covariates20<-dat5k[, c(34:42)]
mos<-cbind(mos,n_covariates20)
#10km
n_covariates10k<-dat10k[, c(34:42)]
moscov10k<-cbind(mos10k,n_covariates10k)

##n_nodes = number of mosquito species columns
##n_covariates = number of covariate columns
## prep_splines = T means it is including the spatial splines based on the coordinates
#5km
MRFcov20_5k_spatial<- MRFcov_spatial(data= mos, symmetrise = "mean", prep_covariates = TRUE, n_nodes = 30, n_cores = 23, n_covariates = 9,
                                     family = "poisson", coords = coords, prep_splines = TRUE, progress_bar = TRUE)
MRFcov20_5k_spatial
#10km
MRFcov20_10k_spatial<- MRFcov_spatial(data= moscov10k, symmetrise = "mean", prep_covariates = TRUE, n_nodes = 30, n_cores = 23, n_covariates = 9,
                                      family = "poisson", coords = coords10k, prep_splines = TRUE, progress_bar = TRUE)
MRFcov20_10k_spatial

#MRF without covariates
#5km
moswo<-subset(mos[, c(1:30)])
View(moswo)
MRFmod20wo<-MRFcov_spatial(data=moswo, symmetrise = "mean", prep_covariates = T, n_nodes = 30, family = "poisson", coords = coords, prep_splines = T, progress_bar = T)
MRFmod20wo
#10km
moswo10k<-subset(mos10k[, c(1:30)])
View(moswo10k)
MRFmod20wo10k<-MRFcov_spatial(data=moswo10k, symmetrise = "mean", prep_covariates = T, n_nodes = 30, family = "poisson", coords = coords10k, prep_splines = T, progress_bar = T)
MRFmod20wo10k



# We can estimate the regression coefficients
MRFmod20wo$key_coefs
MRFcov20_5k_spatial$Key_coefs

## Plot heatmap for MRFmod
#2020
plotMRF_hm(MRF_mod = MRFcov20_5k_spatial, main = 'Estimated node interaction in response to covariates', plot_observed_vals = FALSE, data = MRFcov20_5k_spatial$mrf_data) #FALSE becasue model family is not binomial, cannot plot observed occurrences and co-occurrences

plotMRF_hm(MRF_mod = MRFmod20wo, main = 'Estimated node interaction without covariates', plot_observed_vals = FALSE, data = MRFmod20wo$mrf_data) #FALSE becasue model family is not binomial, cannot plot observed occurrences and co-occurrences


#unpenalized MRF model withuot covariates
#2020
MRF.null.coefs <- plotMRF_hm(MRF_mod = MRFmod20wo, main = "Unregularized MRF")
MRF.null.coefs


## cross validation: As a check of the models estimation accuracy, we run the model using functions in David Harrisâ€™rosaliapackage, which directly solves small MRF networks by maximum likelihood. We specify a weakly informative logistic prior to ensure finite parameter estimates (following Harris, 2016
#2020
#5km
lambda.plot <- cv_MRF_diag(data = mos [, 1:30], n_nodes = 30, n_cores = 23, family = 'poisson', compare_null = TRUE) # working code

#10km
lambda.plot <- cv_MRF_diag(data = moswo10k, n_nodes = 30, n_cores = 23, family = 'poisson', compare_null = TRUE) # working code

##Bootsrap testing for uncertainity
#We can now explore interaction coeffients using bootsrap for uncertainity
#2020 bootstrap MRF without covariates


### THIS IS THE UPDATED BOOTSTRAP -- see the sample_prop=0.8
#Bootstrap without covariates ### Working code###
#5km
booted_MRFwoR <- bootstrap_MRF(data = moswo,n_nodes = 30, sample_seed = 5, sample_prop = 1, n_bootstraps = 500, n_cores = 23, family = 'poisson', spatial = T, coords = coords) 
booted_MRFwoR
plotMRF_hm(booted_MRFwoR)
#10km
booted_MRFwoR <- bootstrap_MRF(data = moswo10k,n_nodes = 30, sample_seed = 5, sample_prop = 1, n_bootstraps = 500, n_cores = 23, family = 'poisson', spatial = T, coords = coords10k) 
booted_MRFwoR
plotMRF_hm(booted_MRFwoR)

#booted_MRFwoR8 <- bootstrap_MRF(data = moswo,n_nodes = 30, sample_seed = 5, sample_prop = 0.8, n_bootstraps = 500, n_cores = 23, family = 'poisson', spatial = T, coords = coords) 
#booted_MRFwoR8

MRF.cov.coefs <- plotMRF_hm(MRF_mod = booted_MRFwoR, main = "Predicted interactions (95% CIs)")
grid::grid.draw(MRF.cov.coefs)


#2020 bootstrap CRF with covariates
#5km
booted_MRF_20_wcR <- bootstrap_MRF(data = mos, n_nodes = 30, sample_seed = 5, sample_prop = 1 , n_covariates = 9, n_bootstraps = 500, n_cores = 23, family = 'poisson', spatial = T, coords = coords)#working code
booted_MRF_20_wcR
MRF.cov.coefs <- plotMRF_hm(MRF_mod = booted_MRF_20_wcR, main = "95% CI of Predicted interactions in response to covariates")
grid::grid.draw(MRF.cov.coefs)

#10km
booted_10k_MRFwcR <- bootstrap_MRF(data = moscov10k,n_nodes = 30, sample_seed = 5, sample_prop = 1, n_covariates = 9, n_bootstraps = 500, n_cores = 23, family = 'poisson', spatial = T, coords = coords10k) 
booted_10k_MRFwcR
plotMRF_hm(booted_10k_MRFwcR)

## plot net of MRFmod
#2020 CRF with covariates
net20<- igraph::graph.adjacency(MRFcov20_5k_spatial$graph, weighted = T, mode = "undirected")
# NEW CODE FOR EDITED FIGURES 
CRF20 <- ggraph(net20, layout = 'igraph', algorithm = 'circle') +         
  # to use igraph network model
  geom_edge_link0(aes(width = weight, color= weight < 0), 
                  alpha = 0.5) +                                          
  # add weighted edges to the graph-does not work with theme_void()
  scale_edge_width(range = c(.4, 2)) +                                  
  # control range of edge size
  scale_edge_colour_manual(values=c('firebrick1', 'dodgerblue3')) +
  # for conditional fill, background color white and edge_link aes(width = weight, color= weight < 0) +
  ## scale_edge_colour_gradient2(low = "black", mid = "white", high = "red",
  ##                             midpoint = 0) + 
  ## for gradient fill, change background color and edge_link aes(width = weight, color= weight)
  geom_node_point(color='black', 
                  fill='grey40', 
                  size=4, 
                  shape = 21) +                                           
  # add nodes to the graph
  coord_fixed() +                                                         
  # ensures the network remains a round circle
  expand_limits(x = c(-2, 2), y = c(-2, 2.3)) +
  # to add space for labels
  ## geom_node_text(aes(angle = -((-node_angle(x, y)+90)%%180)+90, 
  ##                    label = name), 
  ##               size=3, 
  ##                hjust='outward') 
  ## geom_node_text(aes(label = name, hjust='outward'), 
  ##                nudge_x = net20$data$x * .1, 
  ##                nudge_y = net20$data$y * .1, repel = T)
  geom_node_text(aes(x = x*1.1, y=y*1.1, 
                     label=name, 
                     hjust='outward', 
                     angle = -((-node_angle(x, y)+90)%%180)+90), 
                 size=3.5) +                                              
  # species name labels at 90° angle to node, text flipped at 180°
  # to not have it stand on its head on the left side of the network
  theme_graph() +                                                         
  # theme has white background
  theme(legend.position = "none")                                         
# remove legend
## labs(title = 'Title', subtitle = 'subtitle')
# adds title and subtitle over plot
CRF20
# install the latest version of graphical extremes with
## install.packages("devtools")
#devtools::install_github("sebastian-engelke/graphicalExtremes")
# Note that only the gitHub version contains the full flights data
library(graphicalExtremes)
library(nloptr)
library(igraph)
library(matrixcalc)
library(MVN)
library(cluster)
library(xtable)

source("help_functions.R")
source("loc_metr_algorithm.R")

# Filter and select observations without NAs
us_cont=flights$airports$IATA[which(flights$airports$Latitude > 24& flights$airports$Latitude < 50  & flights$airports$Longitude < -60  & flights$airports$Longitude > -130)]
airports=which(rowSums(apply(flights$flightCounts[us_cont,us_cont,1:16], c(3),colSums)>1000)==16&rowSums(apply(flights$flightCounts[us_cont,us_cont,1:16], c(3),rowSums)>1000)==16)
airports=names(airports)
d=length(airports)



# Data and exploratory analysis

data=na.omit(flights$delays[,airports,"arrivals"]+flights$delays[,airports,"departures"])
nrow(data)

#Clustering

nClusters <- 4
diss <- emp_vario(data,p=0.95)

kmed <- pam(
  diss,
  diss = TRUE,
  k = nClusters,
  medoids = 'random',
  nstart = 10,
  variant = 'original'
)


# First cluster has largest average variogram

vario_cluster1=diss[which(kmed$clustering==1),which(kmed$clustering==1)]
mean(vario_cluster1[upper.tri(vario_cluster1)])

vario_cluster2=diss[which(kmed$clustering==2),which(kmed$clustering==2)]
mean(vario_cluster2[upper.tri(vario_cluster2)])

vario_cluster3=diss[which(kmed$clustering==3),which(kmed$clustering==3)]
mean(vario_cluster3[upper.tri(vario_cluster3)])

vario_cluster4=diss[which(kmed$clustering==4),which(kmed$clustering==4)]
mean(vario_cluster4[upper.tri(vario_cluster4)])

# Plot airports
flightsplot=plotFlights(airports[which(kmed$clustering==1)],plotConnections = FALSE,map="state",clipMap = 1.3,useAirportNFlights = TRUE,returnGGPlot = TRUE)
flightsplot=flightsplot + scale_size_continuous(name="Number of flights", range=c(2,10))+guides(size = guide_legend(override.aes = list(alpha = 1)))
plot(flightsplot)

pdf("airports_map.pdf",width = 5.5,height = 4.5)
plotFlights(airports[which(kmed$clustering==1)],plotConnections = FALSE,map="state",clipMap = 1.2,useAirportNFlights = TRUE)
dev.off()

# Check first cluster for Hüsler--Reiss

cluster=which(kmed$clustering==1)
airports_cluster = airports[cluster]
m=length(cluster)

plist=seq(0.94,0.97,length.out=7)
pvalues=rep(0,length(plist))
samples=rep(0,length(plist))

for (i in 1:length(plist)) {
  data_exp_margins=log(data2mpareto(data[,cluster],p=plist[i]))
  data_halfspace= data_exp_margins[which(rowSums(data_exp_margins)>0),]
  data_hyperplane= data_halfspace%*%(diag(m)-matrix(1/m,m,m))
  samples[[i]]=nrow(data_hyperplane)
  
  test=mvn(data_hyperplane[,-1],mvn_test = "hz")
  pvalues[i]=test$multivariate_normality$p.value
}

test_results = matrix(0,3,length(plist))
test_results[1,]=plist
test_results[2,]=samples
test_results[3,]=pvalues

table = xtable(test_results)
print(table,type = "latex")

# Exploratory analysis for positive dependence

p=0.95
vario_emp=emp_vario(data[,cluster],p=p)
Theta_emp=Gamma2Theta(vario_emp)

graph_full=make_full_graph(m)
lm = loc_metr(vario_emp,graph_full)
lm$lm
prop_emp=1-length(lm$failed_ineq)/(lm$number_valid_ineq+length(lm$failed_ineq))
prop_emp

mean(Theta2Qvec(Theta_emp)>=0)


##########################
## Estimation under LMP ##
##########################

# Separating data into training and validation data
dates <- as.Date(rownames(data[,cluster]))
years <- as.numeric(format(dates, "%Y"))
fold_id <- floor((years - min(years)) / 4) + 1
folds <- split(seq_len(nrow(data)), fold_id)
indices = combn(1:4, 2)

# Estimation

results1_list <- list()
results2_list <- list()
Gamma1_listlist <- list()
Gamma2_listlist <- list()
Gamma1_best_list <- list()
Gamma2_best_list <- list()
loglik1_training_listlist <-list()
loglik2_training_listlist <-list()
proportion_LMP_listlist <- list()
proportion_EMTP2_listlist <- list()
ecount_listlist <- list()
dual_dim_list <- list()
sample_size1 <-list()
sample_size2 <-list()


rholist=seq(0, 0.3, length.out = 11)


for (k in 1:ncol(indices)) {
  train_idx <- unlist(folds[indices[,k]], use.names = FALSE)
  test_idx <- unlist(folds[setdiff(1:4,indices[,k])],use.names = FALSE)
  
  training_data <- data[train_idx,cluster ]
  validation_data  <- data[test_idx,cluster ]
  
  
  #Estimate empirical variogram
  Y=data2mpareto(training_data,p=p)
  vario_emp=emp_vario(Y)
  Theta_emp=Gamma2Theta(vario_emp)
  sample_size1[[k]]=nrow(Y)
  sample_size2[[k]]=nrow(data2mpareto(validation_data,p=p))
  
  ##########################
  ### eglearn estimation ###
  ##########################
  
  lasso.est = eglearn(training_data,p=p,rholist=rholist)
  
  #First step estimate and HR likelihood evaluation on validation data
  
  Gamma1_list <- list()
  loglik_step1 <- list()
  proportion_LMP_list <- list()
  proportion_EMTP2_list <- list()
  ecount_list <- list()
  loglik1_training_list<-list()
  
  for (i in 1:length(rholist)){
    Gamma1_list[[i]] = complete_Gamma(vario_emp,graph=lasso.est$graph[[i]],final_tol=1e-6, N=75000)
    loglik_step1[[i]] <- loglik_HR(data=validation_data,p=p, graph = lasso.est$graph[[i]],
                                   Gamma = Gamma1_list[[i]], cens = FALSE)[1]
    loglik1_training_list[[i]] <- loglik_HR(data=training_data,p=p, graph = lasso.est$graph[[i]],
                                            Gamma = Gamma1_list[[i]], cens = FALSE)[1]
    lm = loc_metr(Gamma1_list[[i]],lasso.est$graph[[i]])
    dual_dim_list[[i]]=lm$number_valid_ineq+length(lm$failed_ineq)
    proportion_LMP_list[[i]]=(1-length(lm$failed_ineq)/dual_dim_list[[i]])*100
    proportion_EMTP2_list[[i]]=(mean(Theta2Qvec(Gamma2Theta(Gamma1_list[[i]]))>=-1e-5))*100
    ecount_list[[i]] = ecount(lasso.est$graph[[i]]) 
  }
  
  results1_list[[k]]=loglik_step1
  Gamma1_listlist[[k]] = Gamma1_list
  proportion_LMP_listlist[[k]]=proportion_LMP_list
  proportion_EMTP2_listlist[[k]]=proportion_EMTP2_list
  ecount_listlist[[k]]=ecount_list
  loglik1_training_listlist[[k]]=loglik1_training_list
  
  lasso.best=which(unlist(loglik_step1)==max(unlist(loglik_step1)))
  Gamma1_best_list[[k]]=Gamma1_list[[lasso.best]]
  
  
  
  
  #####################################
  ## Second step for lasso estimates ##
  #####################################
  
  #compute second step estimates and list of HR log-likelihoods on validation data.
  Gamma2_list <- list()
  loglik_step2 <- list()
  loglik2_training_list <- list()
  for (i in 1:length(rholist)){
    res = loc_metr_alg(Gamma1_list[[i]],lasso.est$graph[[i]])
    Gamma2_list[[i]]=res$Gamma_opt
    if(all(Gamma1_list[[i]]==Gamma2_list[[i]])){loglik_step2[[i]]=loglik_step1[[i]]}
    else{loglik_step2[[i]] <- loglik_HR(data=validation_data,p=p, graph = lasso.est$graph[[i]],
                                        Gamma = Gamma2_list[[i]], cens = FALSE)[1]}
    loglik2_training_list[[i]] <- loglik_HR(data=training_data,p=p, graph = lasso.est$graph[[i]],
                                            Gamma = Gamma2_list[[i]], cens = FALSE)[1]
  }
  Gamma2_listlist[[k]] =Gamma2_list 
  lasso.best2=which(unlist(loglik_step2)==max(unlist(loglik_step2)))
  Gamma2_best_list[[k]]=Gamma2_list[[lasso.best2]]
  results2_list[[k]]=loglik_step2
  loglik2_training_listlist[[k]]=loglik2_training_list
  print(k)
}


############

results1= matrix(unlist(results1_list), ncol = length(results1_list), byrow=FALSE)
results2= matrix(unlist(results2_list), ncol = length(results1_list), byrow=FALSE)

results2-results1

results1_max=apply(results1,2,max)
results2_max=apply(results2,2,max)
results2_max-results1_max



mean1 <- rowMeans(results1, na.rm = TRUE)
mean2 <- rowMeans(results2, na.rm = TRUE)

which(mean1==max(mean1))
which(mean2==max(mean2))


edges=rowMeans(matrix(unlist(ecount_listlist), ncol = length(ecount_listlist), byrow=FALSE))
proportions_LMP = rowMeans(matrix(unlist(proportion_LMP_listlist), ncol = length(proportion_LMP_listlist), byrow=FALSE))
proportions_EMTP2 = rowMeans(matrix(unlist(proportion_EMTP2_listlist), ncol = length(proportion_EMTP2_listlist), byrow=FALSE))


proportions = matrix(0, nrow = 7, ncol = length(rholist))
proportions[1,]=rholist
proportions[2,]= edges
proportions[3,]=unlist(dual_dim_list)
proportions[4,]= proportions_LMP
proportions[5,]= proportions_EMTP2
proportions[6,]= mean1
proportions[7,]= mean2

table = xtable(round(proportions,2))
print(table,type = "latex")

# Training data likelihood

loglik_train1= matrix(unlist(loglik1_training_listlist), ncol = length(loglik1_training_listlist), byrow=FALSE)
loglik_train2= matrix(unlist(loglik2_training_listlist), ncol = length(loglik2_training_listlist), byrow=FALSE)

mean1 <- rowMeans(loglik_train1, na.rm = TRUE)
mean2 <- rowMeans(loglik_train2, na.rm = TRUE)

matplot(rholist, cbind(mean1, mean2), type = "b",
        xlab = expression(paste("regularization parameter ", rho )),
        ylab = "Log-likelihood", main = expression(paste("Training log-likelihood vs. ", rho )),
        col=c('black','blue'), pch=19, lwd=2)

# Sample size

mean(unlist(sample_size1))
mean(unlist(sample_size2))

print(xtable(rbind(unlist(sample_size1),unlist(sample_size2))),type = "latex")



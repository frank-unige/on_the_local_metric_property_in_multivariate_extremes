# install the latest version of graphical extremes with
## install.packages("devtools")
#devtools::install_github("sebastian-engelke/graphicalExtremes")
# Note that only the gitHub version contains the full flights data
library(graphicalExtremes)
library(nloptr)
library(igraph)
library(matrixcalc)
source("help_functions.R")
source("loc_metr_algorithm.R")

# Filter and select observations without NAs

airports=which(rowSums(apply(flights$flightCounts[,,1:16], c(3),colSums)>2000)==16&rowSums(apply(flights$flightCounts[,,1:16], c(3),rowSums)>2000)==16)
d=length(airports)

# Plot airports

pdf("airports_map.pdf",width = 7,height = 4.5)
plotFlights(names(airports),plotConnections = FALSE,map="world",clipMap = 1.3,useAirportNFlights = TRUE)
dev.off()

# Separating data into training and validation data

training_data=na.omit(flights$delays[1:2191,airports,"arrivals"]+flights$delays[1:2191,airports,"departures"])
validation_data=na.omit(flights$delays[2192:5601,airports,"arrivals"]+flights$delays[2192:5601,airports,"departures"])

#Estimate empirical variogram
p=0.85
vario_emp=emp_vario(training_data,p=p)
Theta_emp=Gamma2Theta(vario_emp)

#Exploratory analysis for metric property and EMTP2

graph_full=make_full_graph(d)
lm = loc_metr(vario_emp,graph_full)
lm$lm
prop_emp=1-length(lm$failed_ineq)/(lm$number_valid_ineq+length(lm$failed_ineq))
prop_emp

mean(Theta2Qvec(Theta_emp)>=0)

##########################
### eglearn estimation ###
##########################

rholist=seq(0, 0.2, length.out = 21)
lasso.est = eglearn(training_data,p=p,rholist=rholist)

#First step estimate and HR likelihood evaluation on validation data

Gamma1_list <- list()
loglik_step1 <- list()
for (i in 1:length(rholist)){
  Gamma1_list[[i]] = complete_Gamma(vario_emp,graph=lasso.est$graph[[i]],final_tol=1e-6)
  loglik_step1[[i]] <- loglik_HR(data=validation_data,p=p, graph = lasso.est$graph[[i]],
                                 Gamma = Gamma1_list[[i]], cens = FALSE)[1]
}

lasso.best=which(unlist(loglik_step1)==max(unlist(loglik_step1)))
loglik_step1[lasso.best]
lasso.est$graph[[lasso.best]]

pdf("best_graph.pdf",width = 7,height = 4.5)
plotFlights(names(airports),graph=lasso.est$graph[[lasso.best]],map="world",clipMap = 1.3)
dev.off()

# First step estimates partially satisfy the local metric property, but not the EMTP2 constraints

for (i in 1:length(rholist)) {
  print(loc_metr(Gamma1_list[[i]],lasso.est$graph[[i]])$lm)
}

for (i in 1:length(rholist)) {
  print(mean(Theta2Qvec(Gamma2Theta(Gamma1_list[[i]]))>=-1e-5)
)
}

# For comparison: Performance of eglearn estimates wrt information criteria is lower than the best above
Gamma1_MBIC = complete_Gamma(vario_emp,graph=lasso.est$graph_ic$mbic,final_tol=1e-6)
Gamma1_BIC = complete_Gamma(vario_emp,graph=lasso.est$graph_ic$bic,final_tol=1e-6)
Gamma1_AIC = complete_Gamma(vario_emp,graph=lasso.est$graph_ic$aic,final_tol=1e-6)
loglik_step1_MBIC = loglik_HR(data=validation_data,p=p, graph = lasso.est$graph_ic$mbic,
                             Gamma = Gamma1_MBIC, cens = FALSE)[1]
loglik_step1_BIC = loglik_HR(data=validation_data,p=p, graph = lasso.est$graph_ic$bic,
                            Gamma = Gamma1_BIC, cens = FALSE)[1]
loglik_step1_AIC = loglik_HR(data=validation_data,p=p, graph = lasso.est$graph_ic$aic,
                             Gamma = Gamma1_AIC, cens = FALSE)[1]
loc_metr(Gamma1_MBIC,lasso.est$graph_ic$mbic)$lm
loc_metr(Gamma1_BIC,lasso.est$graph_ic$bic)$lm
loc_metr(Gamma1_AIC,lasso.est$graph_ic$aic)$lm


#####################################
## Second step for lasso estimates ##
#####################################

#compute second step estimates and list of HR log-likelihoods on validation data.
Gamma2_list <- list()
loglik_step2 <- list()
for (i in 1:length(rholist)){
  Gamma2_list[[i]] = loc_metr_alg(Gamma1_list[[i]],lasso.est$graph[[i]])$Gamma_opt
  if(all(Gamma1_list[[i]]==Gamma2_list[[i]])){loglik_step2[[i]]=loglik_step1[[i]]}
  else{loglik_step2[[i]] <- loglik_HR(data=validation_data,p=p, graph = lasso.est$graph[[i]],
                                 Gamma = Gamma2_list[[i]], cens = FALSE)[1]}
}

# Plotting log-likelihoods vs. Rho
pdf("plot_lik.pdf",width = 7,height = 4.5)
par(cex = 1.25, cex.lab = 1.3, cex.axis = 1, cex.main = 1.5,
    mar = c(4,4,3,2) +.1)
matplot(rholist, cbind(unlist(loglik_step1), unlist(loglik_step2)), type = "b",
        xlab = expression(paste("regularization parameter ", rho )),
        ylab = "log-likelihood", main = expression(paste("Log-likelihood vs. ", rho )),
        ylim = c(min(min(unlist(loglik_step1)), min(unlist(loglik_step2))),
                 max(max(unlist(loglik_step1)), max(unlist(loglik_step2)))),
        col=c('black', 'blue'), pch=19, lwd=2)
grid()
dev.off()



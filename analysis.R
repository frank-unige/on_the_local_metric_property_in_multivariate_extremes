library(tidyverse)
library(latex2exp)
library(xtable)


############################
### Performance analysis ###
############################


data_sim <- read_rds("output/sim_study_speed_high_dim.rds") %>% 
  unnest(cols = c("perf")) %>%
  group_by(type,d) %>% 
  summarise(mean_value = mean(value))%>%
  pivot_wider(names_from = d, values_from = mean_value)

results = matrix(nrow=4,ncol = 3)
colnames(results)=colnames(data_sim)[2:4]
rownames(results)=c("Dual dimension","Start proportion","Duality gap","Time")

results[1,]=as.matrix(data_sim[5,2:4])
results[2,]=as.matrix(data_sim[6,2:4])
results[3,]=as.matrix(data_sim[4,2:4])
results[4,]=as.matrix(data_sim[7,2:4])

results    

table = xtable(round(results,2))
print(table,type = "latex")  



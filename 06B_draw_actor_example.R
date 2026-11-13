source('06B_predict_many_tournaments.R')
z <- readRDS('06B_many_results/three_tournament_population.rds')
actor <- c(1,1,1) # edit these three talents
first <- many_allocate(actor,z$beta,z$phi,z$B,z$fixed$scaler)
responsive <- many_allocate(actor,z$beta,z$phi,z$B,tail(z$checks,1)[[1]]$responsive$scaler)
tab <- rbind(first$allocation,responsive$allocation)
rownames(tab) <- c('Fixed equal E_t reference','Responsive E_t');print(round(tab,3))
write.csv(data.frame(condition=rownames(tab),tab),
  '06B_many_results/actor_example.csv',row.names=FALSE)
draw <- function(){par(mfrow=c(1,2),mar=c(.5,.5,2,.5));many_pie(first,expression(Fixed~equal~E[t]~reference));many_pie(responsive,expression(Responsive~E[t]))}
png('Figure_many_asymmetric_actor_pies.png',width=8,height=4,units='in',res=600,type='cairo');draw();dev.off()
pdf('Figure_many_asymmetric_actor_pies.pdf',width=8,height=4);draw();dev.off()

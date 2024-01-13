

 
## data processing
require(readxl)
require(data.table)
require(tidycensus)
require(lme4)
require(splines)
setwd('/Users/chongliangluo/Dropbox/R/covax/politic/')

# Data sources
# https://vaers.hhs.gov/data/datasets.html 
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/42MVDX 
# https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Age-and-Sex-Trends-in-the-Uni/5i5k-6cmh/explore/
# https://www.cdc.gov/flu/fluvaxview/interactive-general-population.htm
# https://data.census.gov/
# https://cran.r-project.org/web/packages/tidycensus/index.html
 
load('Data_Covid_AE_politics.rda')

######### self-defined functions for data processing ########################
logit = function(x) log(x/(1-x))
## weighted LM: weights from sample size and heteroscadecity
lm_heteroscadecity <- function(x,y,weights){
  lm0 <- lm.wfit(x,y,weights)
  wt1 = weights/sum(weights) 
  wt2 = 1/abs(lm0$res)/sum(1/abs(lm0$res))   
  wt = wt1 * wt2 / sum(wt1 * wt2) * length(y)
  tt = lm(y~0+x,weights=wt,x=T,y=T) 
  r2 = 1 - sum((y - x%*%tt$coef)^2*weights) / sum((y-mean(y))^2*weights) #
  return(cbind(summary(tt)$coef[,c(1,4)], r2))
}  
 
## produce main effects from interaction term in regression
mylincom = function(fit, lc.coef, outcome=''){ 
  b=summary(fit)$coef[,1]
  S=summary(fit)$vcov
  B=lc.coef%*%b  
  SD = sqrt(diag(lc.coef%*%S%*%t(lc.coef)))
  pv = pnorm(-abs(B/SD))*2
  data.frame(Outcome=outcome,
             Factor=apply(lc.coef,1,function(a) paste0(a[a!=0],'*',names(b)[a!=0],collapse = '+')),
             OR=round(exp(B),3), 
             CI95=paste0('(', round(exp(B-1.96*SD),3), ', ', round(exp(B+1.96*SD),3), ')'),
             # CI95_lb=round(exp(B-1.96*SD),3), 
             # CI95_ub=round(exp(B+1.96*SD),3), 
             pval=round(pv,3))
}
## produce OR and 95% CI from regression output
tabfit <- function(fit=NULL, fit.coef=NULL, outcome=NULL, digits = 3){  # table of fit coef: OR, CI, p-val
  if(is.null(fit.coef)) fit.coef = summary(fit)$coef
  tt = fit.coef
  tab = data.frame(Outcome=outcome, Factor=row.names(tt),
                   OR = round(exp(tt[,1]),digits),
                   CI95 = paste0('(', round(exp(tt[,1]-1.96*tt[,2]),digits), ', ', round(exp(tt[,1]+1.96*tt[,2]),digits), ')'),
                   pval = round(tt[,4],3))
  row.names(tab)=NULL
  return(tab)
}
#############################################################################


############################## main analysis ################################
## logistic regression adjusting for FLU reporting rate (logit scale with continuity correction)
res = data.table()
aa = '>=18 Years'
sta = st[age==aa,] # &state_po!='DC'
fit0 = sta[vaccine=='COVID19',glm(cbind(vaers_report,total_admin-vaers_report)~  
                                    logit(sta[vaccine=='FLU',report_admin/10000+0.5/total_admin])  
                                  +sex_MFratio+MedianAge_st+ 
                                    percent_rep_st,family='binomial')]  
res = rbind(res, tabfit(fit.coef=summary(fit0)$coef,outcome=paste0('report_admin, ', aa)))

fit0 = sta[vaccine=='COVID19',glm(cbind(vaers_severe,total_admin-vaers_severe)~
                                    logit(sta[vaccine=='FLU',severe_admin/10000+0.5/total_admin])
                                  +sex_MFratio+ MedianAge_st+ 
                                    percent_rep_st,family='binomial')]  
res = rbind(res, tabfit(fit.coef=summary(fit0)$coef,outcome=paste0('severe_admin, ', aa)))

fit0 = sta[vaccine=='COVID19',glm(cbind(vaers_severe,vaers_report-vaers_severe)~
                                    logit(sta[vaccine=='FLU',severe_report/100+0.5/vaers_report])
                                  +sex_MFratio+ MedianAge_st+ 
                                    percent_rep_st,family='binomial')] 
res = rbind(res, tabfit(fit.coef=summary(fit0)$coef,outcome=paste0('severe_report, ', aa)))

## age-stratified
for(aa in c('18-49 Years','50-64 Years','>=65 Years')){
  sta = st[age==aa,] 
  fit0 = sta[vaccine=='COVID19',glm(cbind(vaers_report,total_admin-vaers_report)~
                                      logit(sta[vaccine=='FLU',report_admin/10000+0.5/total_admin])
                                    +sex_MFratio+  
                                      percent_rep_st,family='binomial')]  
  res = rbind(res, tabfit(fit.coef=summary(fit0)$coef,outcome=paste0('report_admin, ', aa)))
  
  fit0 = sta[vaccine=='COVID19',glm(cbind(vaers_severe,total_admin-vaers_severe)~
                                      logit(sta[vaccine=='FLU',severe_admin/10000+0.5/total_admin])
                                    +sex_MFratio+  
                                      percent_rep_st,family='binomial')]  
  res = rbind(res, tabfit(fit.coef=summary(fit0)$coef,outcome=paste0('severe_admin, ', aa)))
  
  fit0 = sta[vaccine=='COVID19',glm(cbind(vaers_severe,vaers_report-vaers_severe)~
                                      logit(sta[vaccine=='FLU',severe_report/100+0.5/vaers_report])
                                    +sex_MFratio+  
                                      percent_rep_st,family='binomial')] 
  res = rbind(res, tabfit(fit.coef=summary(fit0)$coef,outcome=paste0('severe_report, ', aa)))
}
# write.csv(res, 'tmp.csv')
#############################################################################



####################### sensitivity analysis #1 ############################# 
require(lme4)
# state-level GLMM, COVID vs FLU
res = data.table()
aa = '>=18 Years' 
fit.st = st[age==aa,glmer(cbind(vaers_report,total_admin-vaers_report)~MedianAge_st+sex_MFratio+  
                            vaccine*percent_rep_st+(1|state_po),family='binomial')]
b.covid = mylincom(fit.st, matrix(c(0,0,0,0,1,1),nrow=1), outcome=paste0('report_admin, ',aa)) 
res = rbind(res, tabfit(fit.st, outcome=paste0('report_admin, ',aa), digits = 3), b.covid)

fit.st = st[age==aa,glmer(cbind(vaers_severe,total_admin-vaers_severe)~MedianAge_st+sex_MFratio+  
                            vaccine*percent_rep_st+(1|state_po),family='binomial')]
b.covid = mylincom(fit.st, matrix(c(0,0,0,0,1,1),nrow=1), outcome=paste0('severe_admin, ',aa))
res = rbind(res, tabfit(fit.st, outcome=paste0('severe_admin, ',aa), digits = 3), b.covid)

fit.st = st[age==aa,glmer(cbind(vaers_severe,vaers_report-vaers_severe)~MedianAge_st+sex_MFratio+  
                            vaccine*percent_rep_st+(1|state_po),family='binomial')]
b.covid = mylincom(fit.st, matrix(c(0,0,0,0,1,1),nrow=1), outcome=paste0('severe_report, ',aa))
res = rbind(res, tabfit(fit.st, outcome=paste0('severe_report, ',aa), digits = 3), b.covid)
for(aa in c('18-49 Years','50-64 Years','>=65 Years')){ 
  fit.st = st[age==aa,glmer(cbind(vaers_report,total_admin-vaers_report)~sex_MFratio+ 
                              vaccine*percent_rep_st+(1|state_po),family='binomial')]
  b.covid = mylincom(fit.st, matrix(c(0,0,0, 1,1),nrow=1), outcome=paste0('report_admin, ',aa)) 
  res = rbind(res, tabfit(fit.st, outcome=paste0('report_admin, ',aa), digits = 3), b.covid)
  
  fit.st = st[age==aa,glmer(cbind(vaers_severe,total_admin-vaers_severe)~sex_MFratio+ 
                              vaccine*percent_rep_st+(1|state_po),family='binomial')]
  b.covid = mylincom(fit.st, matrix(c(0,0,0,1,1),nrow=1), outcome=paste0('severe_admin, ',aa))
  res = rbind(res, tabfit(fit.st, outcome=paste0('severe_admin, ',aa), digits = 3), b.covid)
  
  fit.st = st[age==aa,glmer(cbind(vaers_severe,vaers_report-vaers_severe)~sex_MFratio+ 
                              vaccine*percent_rep_st+(1|state_po),family='binomial')]
  b.covid = mylincom(fit.st, matrix(c(0,0,0,1,1),nrow=1), outcome=paste0('severe_report, ',aa))
  res = rbind(res, tabfit(fit.st, outcome=paste0('severe_report, ',aa), digits = 3), b.covid)
}
# write.csv(res[Factor%in%c('vaccineCOVID19:percent_rep_st')], file='tmp2.csv')
# write.csv(res[Factor%in%c('1*percent_rep_st+1*vaccineCOVID19:percent_rep_st')], file='tmp3.csv')
############################################################################# 


####################### sensitivity analysis #2 #############################
## VAERS individual-level, outcome 3, COVID vs Others
covax[,vax0:=ifelse(vax=='COVID19','COVID19','Other')]
covax[,vax0:=factor(vax0,levels = c('Other','COVID19'))]
covax[,table(vax0)]
fit1t = covax[,glmer(severe~vax0*percent_rep_st+
                       bs(age_st,df=4,degree=3)*percent_rep_st+
                       SEX*percent_rep_st+
                       history*percent_rep_st+(1|state_po), family='binomial')]
round(summary(fit1t)$coef[,-3],3)
tabfit(fit1t, outcome='')

## politic effect-modifiers (interactions)
x = covax$age_st
xbs = bs(x, df=4, degree=3)  
cbs = predict(xbs,seq(-3,3,1)) # age 20 - 80  
mylincom = function(fit, lc.coef){ 
  b=summary(fit)$coef[,1]
  S=summary(fit)$vcov
  B=lc.coef%*%b  
  SD = sqrt(diag(lc.coef%*%S%*%t(lc.coef)))
  pv = pnorm(-abs(B/SD))*2
  data.table(LC=apply(lc.coef,1,function(a) paste0(a[a!=0],'*',names(b)[a!=0],collapse = '+')),
             Estimate=round(B,3), 
             SD=round(SD,3),
             OR=round(exp(B),2), 
             CI95_lb=round(exp(B-1.96*SD),2), 
             CI95_ub=round(exp(B+1.96*SD),2), 
             pvalue=round(pv,3))
} 
OR_politic_age = mylincom(fit1t, cbind(0,0,1,0,0,0,0,0,0,1,cbs,0,0))
OR_politic_vaxOthers = mylincom(fit1t, cbind(0,0,1,0,0,0,0,0,0,0,cbs[4,,drop=F],0,0))
OR_politic_sexM = mylincom(fit1t, cbind(0,0,1,0,0,0,0,0,0,1,cbs[4,,drop=F],1,0))
OR_politic_history1 = mylincom(fit1t, cbind(0,0,1,0,0,0,0,0,0,1,cbs[4,,drop=F],0,1))
OR_politic = rbind(OR_politic_vaxOthers,OR_politic_history1,OR_politic_sexM,OR_politic_age)
OR_politic = OR_politic[c(7,1:6,8:10),] # age=50 is baseline
OR_politic[,label:=c('Baseline', 'Vax = Others', 'History = Yes', 'Sex = M',
                     'Age = 20', 'Age = 30', 'Age = 40', 'Age = 60', 'Age = 70', 'Age = 80')]
OR_politic[,OR:=sprintf('%.2f', OR.V1),
][,pvalue:=as.character(pvalue.V1)
][pvalue=='0',pvalue:='<0.001'
][,CI95:=sprintf('(%.2f, %.2f)', CI95_lb.V1, CI95_ub.V1)]
# write.csv(OR_politic[,-1], file='tmp.csv')
#############################################################################


#########################  main Figure  ####################################
st1 = st[age=='>=18 Years'] # overall
names(st1)[7] = 'percent_REP' # c('report_admin', 'severe_admin', 'severe_report')
st1 = st1[order(percent_REP)]

#### LOESS curves
# pdf('plot_politic_3outcomes_noage.pdf', height=7, width=15)
# pdf('plot_politic_3outcomes_Fig1_loess.pdf', height=7, width=15)
par(mfrow=c(1,3))
st1[vaccine=='COVID19',plot(report_admin~percent_REP, ylim=c(0,30),
                            main='(1) Reporting rate of VAERS among vaccinated', 
                            xlab='% Republican in 2020 presidential election', # yaxt='n',
                            ylab = "Reporting rate (1/10,000)",
                            cex.axis=1.6, cex.lab=1.6, cex.main=1.5,
                            col= 'red', cex=total_admin/max(total_admin)*8)] 
fit_loess = st1[vaccine=='COVID19', loess(report_admin~percent_REP, weights=total_admin,span = 1)]
plx<-predict(fit_loess, se=T)
lines(plx$fit~fit_loess$x)
lines(fit_loess$x,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(fit_loess$x,plx$fit + qt(0.975,plx$df)*plx$se, lty=2) 
st1[vaccine=='COVID19',][rank(report_admin)>=42|rank(report_admin)<2|percent_REP<30|total_admin>4e7,
                         text(percent_REP,report_admin-0.4,state_po,cex=1,col='black')]

st1[vaccine=='COVID19',plot(severe_admin~percent_REP, ylim=c(0,5.5),
                            main='(2) Reporting rate of severe among vaccinated',
                            xlab='% Republican in 2020 presidential election',  
                            ylab = "Reporting rate (1/10,000)", 
                            cex.axis=1.6, cex.lab=1.6,cex.main=1.5,
                            col= 'red', cex=vaers_report/max(vaers_report)*8)]
fit_loess = st1[vaccine=='COVID19', loess(severe_admin~percent_REP, weights=total_admin,span = 1)]
plx<-predict(fit_loess, se=T)
lines(plx$fit~fit_loess$x)
lines(fit_loess$x,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(fit_loess$x,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
st1[vaccine=='COVID19',][rank(severe_admin)>=42|rank(severe_admin)<2|percent_REP<30|total_admin>4e7,
                         text(percent_REP,severe_admin-0.1,state_po,cex=1,col='black')]

st1[vaccine=='COVID19',plot(severe_report~percent_REP, ylim=c(0,32),
                            main='(3) Reporting rate of severe among VAERS',
                            xlab='% Republican in 2020 presidential election',  
                            ylab = "Reporting rate (1/100)", 
                            cex.axis=1.6, cex.lab=1.6,cex.main=1.5,
                            col= 'red', cex=vaers_report/max(vaers_report)*8)]
fit_loess = st1[vaccine=='COVID19', loess(severe_report~percent_REP, weights=vaers_report,span = 1)]
plx<-predict(fit_loess, se=T)
lines(plx$fit~fit_loess$x)
lines(fit_loess$x,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(fit_loess$x,plx$fit + qt(0.975,plx$df)*plx$se, lty=2) 
st1[vaccine=='COVID19',][rank(severe_report)>=42|rank(severe_report)<2|percent_REP<30|total_admin>4e7,
                         text(percent_REP,severe_report-1,state_po,cex=1,col='black')]

# pdf('plot_politic_3outcomes_FigS1_loess_FLU.pdf', height=7, width=15)
par(mfrow=c(1,3))
st1[vaccine=='FLU',plot(report_admin~percent_REP, ylim=c(0,0.8),
                        main='(1) Reporting rate of VAERS among vaccinated', 
                        xlab='% Republican in 2020 presidential election', # yaxt='n',
                        ylab = "Reporting rate (1/10,000)",
                        cex.axis=1.6, cex.lab=1.6, cex.main=1.5,
                        col= 'red', cex=total_admin/max(total_admin)*8)] 
fit_loess = st1[vaccine=='FLU', loess(report_admin~percent_REP, weights=total_admin,span = 1)]
plx<-predict(fit_loess, se=T)
lines(plx$fit~fit_loess$x)
lines(fit_loess$x,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(fit_loess$x,plx$fit + qt(0.975,plx$df)*plx$se, lty=2) 
st1[vaccine=='FLU',][rank(report_admin)>=42|rank(report_admin)<2|percent_REP<30|total_admin>4e7,
                     text(percent_REP,report_admin-0.4,state_po,cex=1,col='black')]

st1[vaccine=='FLU',plot(severe_admin~percent_REP, ylim=c(0,0.08),
                        main='(2) Reporting rate of severe among vaccinated',
                        xlab='% Republican in 2020 presidential election',  
                        ylab = "Reporting rate (1/10,000)", 
                        cex.axis=1.6, cex.lab=1.6,cex.main=1.5,
                        col= 'red', cex=vaers_report/max(vaers_report)*8)]
fit_loess = st1[vaccine=='FLU', loess(severe_admin~percent_REP, weights=total_admin,span = 1)]
plx<-predict(fit_loess, se=T)
lines(plx$fit~fit_loess$x)
lines(fit_loess$x,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(fit_loess$x,plx$fit + qt(0.975,plx$df)*plx$se, lty=2) 
st1[vaccine=='FLU',][rank(severe_admin)>=42|rank(severe_admin)<2|percent_REP<30|total_admin>4e7,
                     text(percent_REP,severe_admin-0.1,state_po,cex=1,col='black')]

st1[vaccine=='FLU',plot(severe_report~percent_REP, ylim=c(0,15),
                        main='(3) Reporting rate of severe among VAERS',
                        xlab='% Republican in 2020 presidential election',  
                        ylab = "Reporting rate (1/100)", 
                        cex.axis=1.6, cex.lab=1.6,cex.main=1.5,
                        col= 'red', cex=vaers_report/max(vaers_report)*8)]
fit_loess = st1[vaccine=='FLU', loess(severe_report~percent_REP, weights=vaers_report,span = 1)]
plx<-predict(fit_loess, se=T)
lines(plx$fit~fit_loess$x)
lines(fit_loess$x,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(fit_loess$x,plx$fit + qt(0.975,plx$df)*plx$se, lty=2) 
st1[vaccine=='FLU',][rank(severe_report)>=42|rank(severe_report)<2|percent_REP<30|total_admin>4e7,
                     text(percent_REP,severe_report-1,state_po,cex=1,col='black')]
dev.off()



### weighted LM line 
bb = matrix(NA, 3*4, 3) 
# pdf('plot_politic_3outcomes_noage_rev.pdf', height=7, width=15)
par(mfrow=c(1,3))
st1[vaccine=='COVID19',plot(report_admin~percent_REP, ylim=c(0,26),
                            main='(1) Reporting rate of VAERS among vaccinated', 
                            xlab='% Republican in 2020 presidential election', # yaxt='n',
                            ylab = "Reporting rate (1/10,000)",
                            cex.axis=1.6, cex.lab=1.6, cex.main=1.5,
                            col= 'red', cex=total_admin/max(total_admin)*8)] 
bb[1:2,] = st1[vaccine=='COVID19', lm_heteroscadecity(cbind(1,percent_REP), report_admin, weights=total_admin)]
abline(bb[1:2,1],col='red')  
st1[vaccine=='COVID19',][rank(report_admin)>=42|rank(report_admin)<2|percent_REP<30|total_admin>4e7,
                         text(percent_REP,report_admin-0.4,state_po,cex=1,col='black')]
bb[3:4,] = st1[vaccine=='FLU', lm_heteroscadecity(cbind(1,percent_REP), report_admin, weights=total_admin)]
abline(bb[3:4,1],col='grey')
legend('topleft', col=c('grey','red'),lty=1, legend=c('FLU','COVID19'),title='Vaccine category',cex = 1.5)

st1[vaccine=='COVID19',plot(severe_admin~percent_REP, ylim=c(0,5.5),
                            main='(2) Reporting rate of severe among vaccinated',
                            xlab='% Republican in 2020 presidential election',  
                            ylab = "Reporting rate (1/10,000)", 
                            cex.axis=1.6, cex.lab=1.6,cex.main=1.5,
                            col= 'red', cex=vaers_report/max(vaers_report)*8)]
bb[5:6,] = tt=st1[vaccine=='COVID19', lm_heteroscadecity(cbind(1,percent_REP), severe_admin, weights=total_admin)]
abline(bb[5:6,1],col='red')
st1[vaccine=='COVID19',][rank(severe_admin)>=42|rank(severe_admin)<2|percent_REP<30|total_admin>4e7,
                         text(percent_REP,severe_admin-0.1,state_po,cex=1,col='black')]
bb[7:8,] = st1[vaccine=='FLU', lm_heteroscadecity(cbind(1,percent_REP), severe_admin, weights=total_admin)]
abline(bb[7:8,1],col='grey')
legend('topleft', col=c('grey','red'),lty=1, legend =c('FLU','COVID19'),title='Vaccine category',cex = 1.5)

st1[vaccine=='COVID19',plot(severe_report~percent_REP, ylim=c(0,32),
                            main='(3) Reporting rate of severe among VAERS',
                            xlab='% Republican in 2020 presidential election',  
                            ylab = "Reporting rate (1/100)", 
                            cex.axis=1.6, cex.lab=1.6,cex.main=1.5,
                            col= 'red', cex=vaers_report/max(vaers_report)*8)]
bb[9:10,] = st1[vaccine=='COVID19', lm_heteroscadecity(cbind(1,percent_REP), severe_report, weights=vaers_report)]
abline(bb[9:10,1],col='red')
st1[vaccine=='COVID19',][rank(severe_report)>=42|rank(severe_report)<2|percent_REP<30|total_admin>4e7,
                         text(percent_REP,severe_report-0.5,state_po,cex=1,col='black')]
bb[11:12,] = st1[vaccine=='FLU', lm_heteroscadecity(cbind(1,percent_REP), severe_report, weights=vaers_report)]
abline(bb[11:12,1],col='grey')
legend('topleft', col=c('grey','red'),lty=1, legend =c('FLU','COVID19'),title='Vaccine category',cex = 1.5)
dev.off()
# round(bb[seq(2,12,2),],4)
# [1,]  0.0653 0.0013 0.1080
# [2,]  0.0015 0.0701 0.2295
# [3,]  0.0218 0.0000 0.0892
# [4,] -0.0001 0.0641 0.0804
# [5,]  0.1895 0.0000 0.0668
# [6,] -0.0537 0.0004 0.0666
##############################################################################

### save data 
save(st, covax, file='Data_Covid_AE_politics.rda')


### Table S1
tmp = copy(st)
tmp[,vaers_report:=paste0(vaers_report, ' (', round(report_admin,2),')')
][,vaers_severe:=paste0(vaers_severe, ' (', round(severe_admin,2),')')]
tmp = cbind(tmp[age=='>=18 Years'&vaccine=='COVID19',][order(percent_rep), c(6,8,7,12,10,11)],
            tmp[age=='>=18 Years'&vaccine=='FLU',][order(percent_rep), c(12,10,11)])
# write.csv(tmp, 'tmp1.csv')

### Table S2, fill out
covax = politic[,c('state_po','percent_rep','elected')][covax,,on='state_po']
covax[,table(vax) ]
covax[,table(vax,severe) ]
covax[vax=='COVID19',table(elected,severe) ]
covax[vax=='COVID19',table(SEX,severe) ]
covax[vax=='COVID19',prop.table(table(elected)),by=severe ]
covax[,table(vax,elected) ]
covax[vax=='COVID19',prop.table(table(elected)) ]
covax[, paste0(round(mean(percent_rep),1),' (', round(sd(percent_rep),1), ')'), by=.(vax,severe)]
covax[, paste0(round(mean(AGE_YRS),1),' (', round(sd(AGE_YRS),1), ')'), by=.(vax,severe)]
covax[, prop.table(table(SEX)), by=.(vax,severe)]
covax[, prop.table(table(SEX)), by=.(vax)]
covax[, table(SEX), by=.(vax)]

covax[vax=='FLU',table(elected,severe) ]
covax[vax=='FLU',prop.table(table(elected)),by=severe ]
covax[vax=='FLU',prop.table(table(elected)) ]
covax[vax=='FLU',t.test(percent_rep~severe) ]
covax[vax=='FLU',t.test(AGE_YRS~severe) ]
covax[vax=='FLU',chisq.test(table(elected,severe)) ]
covax[vax=='FLU',chisq.test(table(SEX,severe)) ]
covax[vax=='FLU',table(SEX,severe) ]


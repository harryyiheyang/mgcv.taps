library(mgcv)
library(devtools)
library(jmcm)
document()
data("aids")

fit0 = gam(cd4~s(time,bs="Fixed",fx=T,xt=list(getA=linearity_piecewise_discontinuity,para=list(para1=1.5,para=0)))+s(age,bs="gp")+s(cesd,bs="gp")+packs+drugs+sex, data=aids, family=quasipoisson(), method="REML", gamma=100)
plot(fit0)

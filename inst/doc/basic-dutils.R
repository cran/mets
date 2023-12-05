## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets)
data(melanoma)

## -----------------------------------------------------------------------------
is.data.frame(melanoma)

## -----------------------------------------------------------------------------
dmean(melanoma,~thick+I(log(thick)))

## -----------------------------------------------------------------------------
dmean(melanoma,~thick+I(log(thick))|I(days>500))

## -----------------------------------------------------------------------------
dmean(melanoma,thick+I(log(thick))~sex|I(days>500))

## -----------------------------------------------------------------------------
dmean(melanoma,thick+I(log(thick))~I(dcut(days)))

## -----------------------------------------------------------------------------
dmean(melanoma,"s*"+"*a*"~sex|I(days>500))

## -----------------------------------------------------------------------------
melanoma=drename(melanoma,tykkelse~thick)
names(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
melanoma=drm(melanoma,~thick+sex)
names(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
melanoma=ddrop(melanoma,~thick+sex)
names(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
melanoma=dkeep(melanoma,~thick+sex+status+days)
names(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
ddrop(melanoma) <- ~thick+sex
names(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
names(melanoma)
melanoma=dkeep(melanoma,~days+status+.)
names(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
dstr(melanoma)

## -----------------------------------------------------------------------------
dlist(melanoma)

## -----------------------------------------------------------------------------
dlist(melanoma, ~.|sex==1)

## -----------------------------------------------------------------------------
dlist(melanoma, ~ulc+days+thick+sex|sex==1)

## -----------------------------------------------------------------------------
dsummary(melanoma)

## -----------------------------------------------------------------------------
dsummary(melanoma,~thick+status+sex)

## -----------------------------------------------------------------------------
dsummary(melanoma,thick+days+status~sex)

## -----------------------------------------------------------------------------
dsummary(melanoma,thick+days+status~sex|thick<97)

## -----------------------------------------------------------------------------
dsummary(melanoma,thick+status~+1|sex==1)

## -----------------------------------------------------------------------------
dsummary(melanoma,~thick+status|sex==1)

## -----------------------------------------------------------------------------
dsummary(melanoma,thick+days+status~sex|I(thick<97 & sex==1))

## -----------------------------------------------------------------------------
dtable(melanoma,~status+sex)

## -----------------------------------------------------------------------------
dtable(melanoma,~status+sex+ulc,level=2)

## -----------------------------------------------------------------------------
dtable(melanoma,~status+sex+ulc,level=1)

## -----------------------------------------------------------------------------
dtable(melanoma,~status+sex+ulc+dcut(days)+I(days>300),level=1)

## -----------------------------------------------------------------------------
data(melanoma)
mel= dsort(melanoma,~days)
dsort(melanoma) <- ~days
head(mel)

## -----------------------------------------------------------------------------
dsort(melanoma) <- ~days-status
head(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
melanoma= transform(melanoma, thick2=thick^2, lthick=log(thick) ) 
dhead(melanoma)

## -----------------------------------------------------------------------------
 melanoma=dtransform(melanoma,ll=thick*1.05^ulc,sex==1)  
 melanoma=dtransform(melanoma,ll=thick,sex!=1)  
 dmean(melanoma,ll~sex+ulc)

## -----------------------------------------------------------------------------
melanoma=dcut(melanoma,~thick,breaks=c(0,200,500,800,2000))

## -----------------------------------------------------------------------------
dlevels(melanoma)

## -----------------------------------------------------------------------------
dtable(melanoma,~thickcat.0)

## -----------------------------------------------------------------------------
dcut(melanoma,breaks=c(0,200,500,800,2000)) <- gr.thick1~thick
dlevels(melanoma)

## -----------------------------------------------------------------------------
dcut(melanoma) <- ~ thick  # new variable is thickcat.4
dlevels(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
dcut(melanoma,breaks=2) <- ~ thick  # new variable is thick.2
dlevels(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
mela= dcut(melanoma,thickcat4+dayscat4~thick+days,breaks=4)
dlevels(mela)

## -----------------------------------------------------------------------------
data(melanoma)
dcut(melanoma,breaks=4) <- thickcat4+dayscat4~thick+days
dlevels(melanoma)

## -----------------------------------------------------------------------------
melanoma$gthick = cut(melanoma$thick,breaks=c(0,200,500,800,2000))
melanoma$gthick = cut(melanoma$thick,breaks=quantile(melanoma$thick),include.lowest=TRUE)

## -----------------------------------------------------------------------------
data(melanoma)
dcut(melanoma,breaks=4) <- thickcat4~thick
dlevels(melanoma) 

## -----------------------------------------------------------------------------
dtable(melanoma,~thickcat4)
melanoma = drelevel(melanoma,~thickcat4,ref="(194,356]")
dlevels(melanoma)

## -----------------------------------------------------------------------------
melanoma = drelevel(melanoma,~thickcat4,ref=2)
dlevels(melanoma)

## -----------------------------------------------------------------------------
melanoma = drelevel(melanoma,~thickcat4,newlevels=1:3)
dlevels(melanoma)

## -----------------------------------------------------------------------------
dkeep(melanoma) <- ~thick+thickcat4
melanoma = drelevel(melanoma,gthick2~thickcat4,newlevels=list(1:2,3:4))
dlevels(melanoma)

## -----------------------------------------------------------------------------
dfactor(melanoma,levels=c(3,1,2,4)) <-  thickcat4.2~thickcat4
dlevel(melanoma,~ "thickcat4*")
dtable(melanoma,~thickcat4+thickcat4.2)

## -----------------------------------------------------------------------------
melanoma=drelevel(melanoma,gthick3~thickcat4,newlevels=list(group1.2=1:2,group3.4=3:4))
dlevels(melanoma)

## -----------------------------------------------------------------------------
data(melanoma)
melanoma = dfactor(melanoma,~status, labels=c("malignant-melanoma","censoring","dead-other"))
melanoma = dfactor(melanoma,sexl~sex,labels=c("females","males"))
dtable(melanoma,~sexl+status.f)

## -----------------------------------------------------------------------------
melanoma = dnumeric(melanoma,~sexl)
dstr(melanoma,"sex*")
dtable(melanoma,~'sex*',level=2)

## -----------------------------------------------------------------------------
sessionInfo()


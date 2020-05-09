Carnivore Coexistence Analysis
================
Clayton T. Lamb, Liber Ero Postdoctoral Fellow, University of British Columbia
09 May, 2020

### Load Data, Functions and Cleanup Data

``` r
##Load Packages
lapply(c("gtsummary", "gridExtra", "grid", "sjPlot", "huxtable", "lwgeom","ggpubr","ggspatial","ggmap", "ggrepel","rgdal","raster","sf","velox","here","survival","lubridate","MuMIn","lme4","ggeffects", "adehabitatHR", "rnaturalearth","rnaturalearthdata","viridisLite", "units", "knitr", "suncalc", "tidyverse"), require, character.only = TRUE)

options(scipen=999)

##ggplot color/fillscales
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
##custom ggplot theme
theme_Publication <- function(...){
  theme_bw()+
  theme(plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold",size = rel(1.3)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(size = rel(1.1)), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.5, "cm"),
        legend.spacing = unit(0.05, "cm"),
        legend.text = element_text(size = rel(1)),
        legend.title = element_blank(),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))
}
##################################################################
##load data
##################################################################
cap <- read.csv(here::here("data", "cap.csv"))
reloc <- read.csv(here::here("data", "reloc.csv"))
mort <- read.csv(here::here("data", "mort.csv"))
repro <- read.csv(here::here("data", "repro.csv"))


##################################################################
##Cleanup
##################################################################
##remove bears not collared
cap <- cap%>%filter(!EndReason%in%c("DidntCollar"))

##Get all into Date Format
cap$CollarStartDate <- ymd(cap$CollarStartDate)
cap$CollarEndDate <- ymd(cap$CollarEndDate)


##Calculate time monitored
cap$time <- as.numeric(cap$CollarEndDate-cap$CollarStartDate)

##add mort info into reloc
mort <- mort%>%
  mutate(Cause_CL=case_when(Cause_CL%in%"dolp" ~ "attract./conflict",
                            TRUE~as.character(Cause_CL)))%>%
  filter(CollarOn_CL%in%"Y")


reloc <- reloc%>%
  select("BearID","Date", "Time", "X", "Y")%>%
  mutate(mort=0,
         Cause_CL=NA,
         HumanSuspected_CL=NA)%>%
  rbind(mort%>%mutate(Time=NA)%>%select("BearID","Date", "Time", "X", "Y", "Cause_CL", "HumanSuspected_CL")%>%mutate(mort=1))


##add sex,age, ageclass,and study area into reloc
cap <- cap%>%group_by(BearID)%>%mutate(b_y=year(CollarStartDate)-Age)%>%
  summarize(birthyear=min(b_y, na.rm=TRUE))%>%
  left_join(cap,by="BearID")%>%
  ungroup()

reloc <- reloc%>%
  left_join(cap%>%select(BearID,Sex,birthyear, StudyArea)%>%
  distinct(BearID,Sex,birthyear, StudyArea),by="BearID")%>%
  mutate(Age=year(Date)-birthyear)%>%
  mutate(ageclass=case_when(
    Age %in% c(0,1,2) ~ "Dependent_Young",
    Age %in% c(3:6) ~ "Subadult",
    Age >6 ~ "Adult"
  ))

##add month day year
reloc <- reloc%>%
  mutate(yr=year(Date),
         m=month(Date),
         d=day(Date))

##only keep active months
reloc <- reloc %>% filter(m%in%c(4:11))

##keep only adult and subadult data
reloc <- reloc%>%filter(ageclass%in%c("Adult","Subadult"))
```

Add Spatial Attributes
----------------------

``` r
##make bear data spatial
reloc <- st_as_sf(reloc, 
                  coords = c("X","Y"),
                  crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

##reproject for MOVEBANK
# reloc%>%
#   filter(StudyArea%in%"Elk Valley-Lamb")%>%
#   st_transform("+init=epsg:26911")%>%
#   cbind(st_coordinates(.))%>%
#   as_tibble()%>%
#   filter(!is.na(Time))%>%
#   mutate(dt=paste(Date, Time, sep=" "))%>%
#   select(-geometry)%>%
#   write_csv("reloc_mb.csv")
# 
# cap%>%
#   filter(StudyArea%in%"Elk Valley-Lamb")%>%
#   mutate(sp="Ursus arctos")%>%
#   write_csv("cap_mb.csv")

## HUMAN INFLUENCE INDEX
##make smaller
hii<- raster(here::here("data", "spatial", "hii.tif"))

##make velox (faster)
hiiv <- velox(hii)

##extract hii to reloc
reloc <- reloc%>%mutate(hii=hiiv$extract_points(sp = reloc))

rm(hii)

## NDVI
##make smaller
ndvi<- raster(here::here("data", "spatial", "ndvi.tif"))

##make velox (faster)
ndviv <- velox(ndvi)

##extract hii to reloc
reloc <- reloc%>%mutate(ndvi=ndviv$extract_points(sp = reloc))

rm(ndvi)

###make Tibble
reloc <- reloc%>%
  as_tibble()
```

Explore Mortality Data
----------------------

``` r
cod <-reloc%>%
  filter(mort==1 & ageclass%in%c("Adult","Subadult"))%>%
  rename(Cause=Cause_CL)%>%
  mutate(Cause=case_when(Cause%in%"unknown"~"unk.",
                         Cause%in%"attract./conflict"~"confl.",
                         Cause%in%"natural"~"nat.",
                         Cause%in%"road/train"~"collis.",
                         Cause%in%"poached"~"nat.",
                         Cause%in%"harvest"~"harv.",
                         TRUE~Cause))%>%
  mutate(Cause=fct_relevel(Cause, "confl.", "collis.", "ill.", "harv.", "nat."))%>%
  ggplot(aes(hii, stat(count), fill = Cause)) +
  geom_density(position = "fill", n=10)+
  ylab("Mortalities (prop.)")+
  xlab("Human Influence Index")+
  xlim(0,40)+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_Publication()+
  guides(fill=guide_legend(nrow=1))+
  theme(legend.key.size= unit(0.25, "cm"),
        legend.spacing = unit(0.01, "cm"),
        legend.text = element_text(size = rel(0.9)),
        legend.position = "none")
cod
```

![](coexistence_sunrise_files/figure-markdown_github/Explore%20Mortality%20Data-1.png)

``` r
##N by type
reloc%>%
  filter(mort==1 & ageclass%in%c("Adult","Subadult"))%>%
  rename(Cause=Cause_CL)%>%
  mutate(Cause=case_when(Cause%in%"unknown"~"unk.",
                         Cause%in%"attract./conflict"~"confl.",
                         Cause%in%"natural"~"nat.",
                         Cause%in%"road/train"~"collis.",
                         Cause%in%"poached"~"nat.",
                         Cause%in%"harvest"~"harv.",
                         TRUE~Cause))%>%
  mutate(Cause=fct_relevel(Cause, "confl.", "collis.", "ill.", "harv.", "nat."))%>%
  group_by(Cause)%>%
  summarise(n=n())

##N human caused? 76%
hc <-mort%>%
  group_by(HumanSuspected_CL)%>%
  tally()
```

### 76 % of grizzly bear mortalities are human caused

Explore Reproduction Data
-------------------------

``` r
##do subadults reproduce much?
repro%>%
  mutate(yr=year(Date))%>%
  left_join(cap%>%mutate(yr=year(CollarStartDate))%>%distinct(BearID, birthyear))%>%
  mutate(age=yr-birthyear)%>%
  mutate(Count=ifelse(is.na(Count),0,Count))%>%
  #mutate(Count=ifelse(Count>=1,1,0))%>%
  filter(age>=0)%>%
  group_by(age)%>%
  summarise(mean=mean(Count))%>%
  ggplot()+
  geom_bar(aes(x=age, y=mean), stat="identity")+
  ylab("# offspring")+
  theme_Publication()
```

![](coexistence_sunrise_files/figure-markdown_github/Explore%20Reproduction%20Data-1.png)

``` r
cub.prop <-repro%>%
  mutate(yr=year(Date))%>%
  left_join(cap%>%mutate(yr=year(CollarStartDate))%>%distinct(BearID, birthyear))%>%
  mutate(age=yr-birthyear)%>%
  mutate(Count=ifelse(is.na(Count),0,Count))%>%
  mutate(Count=ifelse(Count>=1,1,0))%>%
  filter(age>0 & age<=6)%>%
  group_by(BearID)%>%
  summarize(mean=mean(as.numeric(Count), na.rm=TRUE))%>%
  mutate(Count=ifelse(mean>0,1,0))%>%
  summarize(prop=round(mean(Count)*100,0))
```

### 16 % of subadults(&lt;=6) reproduce

Prep Monthly Survival Data
--------------------------

``` r
##pull out exact mort locs
kills <- reloc%>%filter(mort%in%1)

##summarize live locs by month
reloc.m <- reloc%>%  
  filter(mort%in%0)%>%
  group_by(BearID,yr,m,Sex,Age, ageclass,StudyArea, mort)%>%
  summarize_at(vars(hii,ndvi), mean, na.rm = TRUE)%>%
  as_tibble()%>%
  ungroup()

##assign 30 days to each, as they are each a month
reloc.m$time <- 30

##length of survival on last month of life= day of month
kills$time<-kills$d
##match columns
kills <- kills%>%select(colnames(reloc.m))

##bind and join last month with live and dead if needed
reloc.m <- rbind(reloc.m,kills)%>%
  group_by(BearID,Sex,Age, ageclass, StudyArea,yr,m)%>%
  summarize(mort=max(mort, na.rm=TRUE),
            hii=mean(hii, na.rm=TRUE),
            ndvi=mean(ndvi, na.rm=TRUE),
            time=min(time, na.rm=TRUE))%>%
  ungroup()
```

Run Cox Models
--------------

``` r
##################################################################
##Fit Models
##################################################################

fit1 <- coxph(Surv(time, mort)~ Sex +
                cluster(BearID), data = reloc.m)
fit1.1 <- coxph(Surv(time, mort)~ ageclass+
                  cluster(BearID), data = reloc.m)
fit1.2 <- coxph(Surv(time, mort)~ Sex + ageclass +
                  cluster(BearID), data = reloc.m)
fit2 <- coxph(Surv(time, mort)~  Sex * ageclass + 
                cluster(BearID), data = reloc.m)
fit3 <- coxph(Surv(time, mort)~ Sex + ageclass + hii +
                cluster(BearID), data = reloc.m)
fit4 <- coxph(Surv(time, mort)~ (Sex+ageclass*hii) +
                cluster(BearID), data = reloc.m)
fit5 <- coxph(Surv(time, mort)~ (Sex*hii + ageclass) +
                cluster(BearID), data = reloc.m)
fit6 <- coxph(Surv(time, mort)~ (Sex*ageclass*hii) +
                cluster(BearID), data = reloc.m)
fit7 <- coxph(Surv(time, mort)~ Sex + ageclass + hii +  + ndvi + 
                cluster(BearID), data = reloc.m)
fit8 <- coxph(Surv(time, mort)~ (Sex+ageclass*hii) + ndvi + 
                cluster(BearID), data = reloc.m)
fit9 <- coxph(Surv(time, mort)~ (Sex*ageclass*hii) + ndvi + 
                cluster(BearID), data = reloc.m)

##prep model selection table
mod.sel <- model.sel(fit1.1, fit1.2, fit1,fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9,  rank="AICc")

# replace Model name with formulas
for(i in 1:nrow(mod.sel)) mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] 

mod.sel <- mod.sel%>%
  select(Model,df,AICc,delta, weight)%>%
  rename(dAICc=delta)%>%
  mutate(Model=str_replace(Model,"\\|","]"))%>%
  mutate_if(is.numeric, function(x) round(x, 2))%>%
  as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_col_width(c(4.2,0.5,1,0.7,0.7))%>%
  set_width(0.6)

huxtable::number_format(mod.sel)[, 3] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)

quick_docx(mod.sel,file=here::here("tables", "MortHazard_AIC.docx"), open=FALSE)


##################################################################
##Test Assumptions
##################################################################
print(cox.zph(fit8))


##################################################################
##Plot Effects, weighted by top models
##################################################################
##make prediction data
surv.pred <- expand.grid(Sex=c("M","F"), ageclass=c("Subadult","Adult"), hii=c(0:40), ndvi=mean(reloc.m$ndvi, na.rm=TRUE))

###predict each of 4 top mods
surv.pred$pred6 <- predict(fit6, newdata=surv.pred, type="risk")
surv.pred$pred.se6 <- predict(fit6, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

surv.pred$pred7 <- predict(fit7, newdata=surv.pred, type="risk")
surv.pred$pred.se7 <- predict(fit7, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

surv.pred$pred8 <- predict(fit8, newdata=surv.pred, type="risk")
surv.pred$pred.se8 <- predict(fit8, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

surv.pred$pred9 <- predict(fit9, newdata=surv.pred, type="risk")
surv.pred$pred.se9 <- predict(fit9, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit


##average each model prediction by model weight
wts <- Weights(AIC(fit6, fit7, fit8, fit9))
surv.pred$pred <- (surv.pred$pred6*wts[[1]])+(surv.pred$pred7*wts[[2]])+(surv.pred$pred8*wts[[3]])+(surv.pred$pred9*wts[[4]])
surv.pred$pred.se <- (surv.pred$pred.se6*wts[[1]])+(surv.pred$pred.se7*wts[[2]])+(surv.pred$pred.se8*wts[[3]])+(surv.pred$pred.se9*wts[[4]])


##relevel factors for plotting
surv.pred <- surv.pred%>%
  mutate(ageclass=fct_relevel(ageclass,levels=c("Subadult", "Adult")),
         Sex=fct_relevel(Sex,levels=c("M", "F")))

##summarize across sexes to produce simpler plot
surv.pred<- surv.pred%>%
  group_by(ageclass, hii)%>%
  summarise(pred=mean(pred,na.rm=TRUE),
            pred.se=mean(pred.se, na.rm=TRUE))
```

``` r
kable(mod.sel[-1,], escape=FALSE)
```

|     | Model                                                                                  | df  | AICc    | dAICc | weight |
|-----|:---------------------------------------------------------------------------------------|:----|:--------|:------|:-------|
| 2   | Sex + ageclass + hii + ndvi + ageclass:hii                                             | 5   | 1085.82 | 0     | 0.57   |
| 3   | Sex + ageclass + hii + ndvi + Sex:ageclass + Sex:hii + ageclass:hii + Sex:ageclass:hii | 8   | 1088.06 | 2.23  | 0.19   |
| 4   | Sex + ageclass + hii + ndvi                                                            | 4   | 1089.18 | 3.36  | 0.11   |
| 5   | Sex + ageclass + hii + ageclass:hii                                                    | 4   | 1089.66 | 3.83  | 0.08   |
| 6   | Sex + ageclass + hii + Sex:ageclass + Sex:hii + ageclass:hii + Sex:ageclass:hii        | 7   | 1091.17 | 5.34  | 0.04   |
| 7   | Sex + ageclass + hii                                                                   | 3   | 1093.64 | 7.81  | 0.01   |
| 8   | Sex + hii + ageclass + Sex:hii                                                         | 4   | 1094.95 | 9.12  | 0.01   |
| 9   | Sex + ageclass                                                                         | 2   | 1108.98 | 23.16 | 0      |
| 10  | Sex + ageclass + Sex:ageclass                                                          | 3   | 1109    | 23.18 | 0      |
| 11  | ageclass                                                                               | 1   | 1111.67 | 25.84 | 0      |
| 12  | Sex                                                                                    | 1   | 1116.68 | 30.86 | 0      |

Cox Models Results
------------------

``` r
###plot standardized response
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

plot_model( 
  coxph(
    Surv(time, mort)~ (Sex+ageclass*hii) + ndvi + cluster(BearID),
    data = reloc.m%>%ungroup()%>%mutate(hii=range01(hii, na.rm=TRUE))
  ),
  sort.est = TRUE,
  title = "",
  axis.title=c("Mortality Hazard"))+
  geom_hline(yintercept = 1, linetype="dotted")+
  theme_Publication()
```

![](coexistence_sunrise_files/figure-markdown_github/Cox%20Models%201%20Results-1.png)

``` r
ggsave(here::here("plots", "modelplots","MortHazard_HII_Coeff.png"), width=6, height=5, units="in")

##coef table
tbl_regression(fit8)


surv <- ggplot(surv.pred,aes(x=hii,y=pred,fill=ageclass, color=ageclass))+
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin = pred-pred.se, ymax = pred+pred.se), colour = NA,alpha = 0.3) +
  ylab("Mortality Risk")+
  xlab("Human Influence Index")+
  xlim(0,40)+
  #ylim(0,23)+
  scale_colour_manual(values=c("#386cb0", "#fdb462"))+
  scale_fill_manual(values=c("#386cb0", "#fdb462"))+
  theme_Publication()+
  theme(legend.position = c(0.4,0.9))

#ggsave(here::here("plots","hazard.png"), width=4, height=5, units="in")


surv
```

![](coexistence_sunrise_files/figure-markdown_github/Cox%20Models%201%20Results-2.png)

``` r
###summary stats
#reloc%>%
#  summarise(n=n_distinct(BearID))

##ratio of hazard between adults and subads
sa_risk <-surv.pred%>%spread(ageclass,pred)%>%
  group_by(hii)%>%
  summarize(s=mean(Subadult, na.rm=TRUE),
            a=mean(Adult, na.rm=TRUE))%>%
  mutate(dif=(s/a))%>%
  filter(hii==40)%>%
  pull(dif)
```

### Compared to adults, subadults faced an 17.1x higher mortality risk in the highest human-influenced areas

Individual Spatial Avoidance of HII
-----------------------------------

``` r
##assign median risk individual experienced overall
risk.bears <- reloc.m%>%
  filter(ageclass%in%c("Subadult", "Adult"))%>%
  group_by(BearID)%>%
  summarize(med.hii=median(hii, na.rm=TRUE))%>%
  dplyr::select(BearID,med.hii)

##assign median risk individual experienced as subadult
risk.bears2 <- reloc.m%>%
  filter(ageclass%in%c("Subadult"))%>%
  group_by(BearID)%>%
  summarize(med.hii2=median(hii, na.rm=TRUE))%>%
  dplyr::select(BearID,med.hii2)

##prep model data
a <- reloc.m%>%
  filter(ageclass%in%c("Adult","Subadult"))%>%
  mutate(Age2=as.numeric(paste(Age,m-3,sep=".")))%>%
  left_join(risk.bears, by=c("BearID"))%>%
  left_join(risk.bears2, by=c("BearID"))%>%
  mutate(hii2=hii-med.hii)%>%
  ungroup()%>%
  dplyr::select(BearID, ageclass, hii2, med.hii, med.hii2, ndvi, m, Age2)%>%
  drop_na()

a%>%
  group_by(BearID)%>%
  summarise(n=n_distinct(Age2))%>%
  filter(n>=10)%>%
  pull(BearID)->keep

a <- a%>%
  filter(BearID%in% keep)

unique(a$BearID)%>%length()

##model
spat.learn <- lmer(hii2 ~ Age2 +(1 | BearID) + (1|m), data = a,
                   REML = FALSE)

spat.learn2 <- lmer(hii2 ~ 1 + (1 | BearID), data = a,
                    REML = FALSE)

spat.learn3 <- lmer(hii2 ~ 1 + (1 | BearID) + (1|m), data = a,
                    REML = FALSE)

spat.learn4 <- lmer(hii2 ~ Age2 + med.hii2 + (1 | BearID) + (1|m), data = a,
                    REML = FALSE)

spat.learn5 <- lmer(hii2 ~ Age2*med.hii2 + (1 | BearID) + (1|m), data = a,
                    REML = FALSE)

spat.learn6 <- lmer(hii2 ~ Age2*med.hii2 + ndvi + (1 | BearID) + (1|m), data = a,
                    REML = FALSE)


##model table
mod.sel <- model.sel(spat.learn, spat.learn2, spat.learn3, spat.learn4, spat.learn5, spat.learn6)

for(i in 1:nrow(mod.sel)) {mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] }

mod.sel <- mod.sel%>%
  select(Model,df,AICc,delta, weight)%>%
  rename(dAICc=delta)%>%
  mutate(Model=str_replace_all(Model,"\\|","]"))%>%
  mutate_if(is.numeric, function(x) round(x, 2))%>%
  as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_col_width(c(4.2,0.5,1,0.7,0.7))%>%
  set_width(0.6)

huxtable::number_format(mod.sel)[, 3] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)

quick_docx(mod.sel, file=here::here("tables", "spat.learn_AIC.docx"), open=FALSE)




###PREDICT

pred.spat.learn <- expand.grid(Age2=unique(a$Age2),  BearID="Bailey_L_EV", m=7, med.hii2=c(0:40), ndvi=mean(a$ndvi, na.rm=TRUE))


pred.spat.learn$pred <- predict(spat.learn6, newdata = pred.spat.learn)


pred.spat.learn2<- pred.spat.learn%>%
  mutate(ageclass=case_when(Age2<=6~"Subadult",
                            Age2>6~"Adult"))%>%
  dplyr::select(ageclass,med.hii2,pred)%>%
  group_by(ageclass,med.hii2)%>%
  summarise(pred=mean(pred,na.rm=TRUE))%>%
  spread(ageclass,pred)%>%
  mutate(change=Adult-Subadult,
         ageclass="Adult")%>%
  dplyr::select(med.hii2,change, ageclass)%>%
  as.tibble()%>%
  rbind(data.frame(med.hii2=c(0:40), change=0, ageclass="Subadult"))%>%
  mutate(ageclass=fct_relevel(ageclass, levels=c("Subadult","Adult")),
         change2=(change/med.hii2)*100)%>%
  filter(med.hii2%in%seq(5,40,by=5))
```

``` r
kable(mod.sel[-1,])
```

|     | Model                                              | df  | AICc    | dAICc | weight |
|-----|:---------------------------------------------------|:----|:--------|:------|:-------|
| 2   | Age2 \* med.hii2 + ndvi + (1 \] BearID) + (1 \] m) | 8   | 7166.97 | 0     | 0.96   |
| 3   | Age2 + med.hii2 + (1 \] BearID) + (1 \] m)         | 6   | 7174.33 | 7.36  | 0.02   |
| 4   | Age2 \* med.hii2 + (1 \] BearID) + (1 \] m)        | 7   | 7175.34 | 8.37  | 0.01   |
| 5   | Age2 + (1 \] BearID) + (1 \] m)                    | 5   | 7180.25 | 13.28 | 0      |
| 6   | 1 + (1 \] BearID) + (1 \] m)                       | 4   | 7180.72 | 13.75 | 0      |
| 7   | 1 + (1 \] BearID)                                  | 3   | 7189.37 | 22.4  | 0      |

Results: Individual Spatial Avoidance of HII
--------------------------------------------

``` r
# replace Model name with formulas

###plot standardized response
  
plot_model( 
  lmer(hii2 ~ Age*`HII (median, subadult)` + ndvi + (1 | BearID) + (1|m), data = a%>%ungroup()%>%mutate(med.hii2=range01(med.hii2, na.rm=TRUE),
                                                                                                        ndvi=range01(ndvi, na.rm=TRUE),
                                                                                                        Age2=range01(Age2, na.rm=TRUE))%>%rename(`HII (median, subadult)`=med.hii2, Age=Age2),
       REML = FALSE),
  sort.est = TRUE,
  title = "",
  axis.title=c("Coefficient (std.)")
)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_Publication()
```

![](coexistence_sunrise_files/figure-markdown_github/Spatial%20Hab%20Use%20Results-1.png)

``` r
ggsave(here::here("plots", "modelplots","spat.learn_Coeff.png"), width=5, height=4, units="in")

##coefficient table
tbl_regression(spat.learn6)

hab.change <- ggplot(pred.spat.learn2%>%mutate(change2=case_when(med.hii2<10~NA_real_, TRUE~change2)), aes(x=change2, y=med.hii2,  group=med.hii2))+
  geom_line(aes(color=ageclass),size=0.5)+
  geom_point(aes(fill=ageclass),size=2, colour="black",pch=21) +
  xlim(-13, 57)+
  ylim(0,43)+
  xlab("Habitat Use (% Change)")+
  ylab("")+
  theme_bw()+
  #ggtitle("Change in space use with age: within individuals")+
  scale_fill_Publication()+
  scale_colour_Publication()+
  theme_Publication()+
  theme(legend.position = c(0.7,0.08))
hab.change
```

![](coexistence_sunrise_files/figure-markdown_github/Spatial%20Hab%20Use%20Results-2.png)

``` r
#ggsave(here::here("plots","Space_Within2.png"), width=3, height=3.5, units="in")
```

Individual Temporal Avoidance of HII
------------------------------------

``` r
# ###ADD STEP LENGTHS (SAVE AND RELOAD AS IT TAKES A WHILE)
# ##make bear data spatial
# reloc <- st_as_sf(reloc,
#                   crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
# 
# ##prep loop components
# ind <-unique(reloc$BearID)
# reloc.step <- list()
# ###Calculate Step Length
# for(i in 1:length(ind)){
#   a <- reloc%>%
#     filter(BearID%in%ind[i])%>%
#     mutate(dt=ymd_hms(paste(Date,Time, sep=" ")),
#            dist=NA,
#            timedif=NA)%>%
#     arrange(dt)
# 
#   if(nrow(a)>1){
#   for(j in 2:nrow(a)){
#     a[j,"timedif"] <- difftime(a[j,"dt"][[1]],a[j-1,"dt"][[1]],units='hours')%>%as.numeric()
#     a[j,"dist"] <- st_distance(a[j,], a[j-1,])[[1]]%>%as.numeric()
#   }
#     reloc.step[[i]] <- as_tibble(a)
#   print(i)
#   }
# }
# reloc.step <- bind_rows(reloc.step)
# 
# reloc.step <- reloc.step%>%
#    st_as_sf(crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")%>%
#   cbind(st_coordinates(.))
# 
# 
# write_csv(reloc.step%>%select(-geometry), "cap_fix_dist.csv")
# 
# move<- reloc.step


##load
move <- read_csv("cap_fix_dist.csv")

##get into lat long
move <- st_as_sf(move,
                 coords=c("X","Y"),
                  crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")%>%
  st_transform("+init=epsg:4326")%>%
  cbind(st_coordinates(.))%>%
  as_tibble()%>%
  select(-geometry)


##convert all times to PST
move <- move%>%
  mutate(dt2=case_when(StudyArea%in%c("Westslopes-McLellan", "Flathead-McLellan","Elk Valley-Apps","Elk Valley-Lamb")~ force_tz(move$dt,"America/Edmonton")%>%with_tz("America/Vancouver"),
        TRUE~force_tz(move$dt,"America/Vancouver"))) 

#mapview(move%>%filter(StudyArea%in%"Elk Valley-Lamb"))

##assign day and night
daym1 <-getSunlightTimes(data=move%>%mutate(Date=Date-1)%>%rename(lat=Y, lon=X, date=Date)%>%select(lat,lon,date),
  keep = c("dusk"), tz = "America/Vancouver")

day <-getSunlightTimes(data=move%>%rename(lat=Y, lon=X, date=Date)%>%select(lat,lon,date),
  keep = c("dawn", "dusk"), tz = "America/Vancouver")

dayp1 <-getSunlightTimes(data=move%>%mutate(Date=Date+1)%>%rename(lat=Y, lon=X, date=Date)%>%select(lat,lon,date),
  keep = c("dawn"), tz = "America/Vancouver")


##make day/night intervals
move$night1= daym1$dusk%--%day$dawn
move$day= day$dawn%--%day$dusk
move$night2= day$dusk%--%dayp1$dawn


###make relocation intervals
move<-move%>%
  filter(timedif<=8 & timedif>0.5)%>%drop_na(timedif)

move$reloc.interval <- (move$dt2-hours(move$timedif%>%as.integer()))%--%move$dt2

move.select <- move%>%
              mutate(day.night=case_when(reloc.interval %within% night1~ as.character("Night"),
                             reloc.interval %within% night2 ~ as.character("Night"),
                             reloc.interval %within% day ~ as.character("Day")), 
         Age2=as.numeric(paste(Age,(month(Date)-3),sep=".")),
         m=as.factor(month(month(Date))))%>%
  filter(!is.na(day.night))%>%
  filter(StudyArea %in% c( "Flathead-McLellan", "Parsnip-Ciarniello", "SelkirkPurcell-Proctor", "Elk Valley-Lamb", "Lilloet-Iredale", "BabineKispiox-MacHutchonMahon", "Elk Valley-Apps", "Tulsequah-Diemart"))

ggplot(move.select, aes(dist/1000))+
  geom_histogram()+
  theme_Publication()+
  xlab("Distance (km)")
```

![](coexistence_sunrise_files/figure-markdown_github/Temporal%20Avoidance:prep%20data-1.png)

``` r
##remove a couple outliers
move.select<- move.select%>% filter(dist>0 & dist<15000) ##can't really be 0 between successive collar relocations, as even if in same spot collar error should mean >0..0.001..
ggplot(move.select, aes(dist/1000))+
  geom_histogram()+
  theme_Publication()+
  xlab("Distance (km)")
```

![](coexistence_sunrise_files/figure-markdown_github/Temporal%20Avoidance:prep%20data-2.png)

``` r
##save
write_csv(move.select, "relocs_nocturnality.csv")

length(unique(move.select$BearID))
##plot

move.select%>%
  select(day, StudyArea)%>%
  mutate(dawn=int_start(day)%>%with_tz(tz="MST")%>%hour()+(int_start(day)%>%with_tz(tz="MST")%>%minute())/60,
         dusk=int_end(day)%>%with_tz(tz="MST")%>%hour()+(int_end(day)%>%with_tz(tz="MST")%>%minute())/60,
         day=int_start(day)%>%with_tz(tz="MST")%>%ymd_hms()%>%yday())%>%
  mutate(dusk=case_when(dusk<2~dusk+24,TRUE~dusk))%>%
  gather(class,hour,-day,-StudyArea)%>%
  arrange(StudyArea, day,class)%>%
  ggplot()+
  geom_path(aes(x=day,y=hour, group=StudyArea, color=StudyArea))+
  facet_grid(.~class)+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%d-%b")) +
  scale_y_reverse()+
  theme_bw()
```

![](coexistence_sunrise_files/figure-markdown_github/Temporal%20Avoidance:prep%20data-3.png)

``` r
ggsave(here::here("plots", "rev1", "dawn_dusk.png"), width=9, height=5)
```

Individual Temporal Avoidance of HII
------------------------------------

``` r
##################
###Change in time use WITHIN individuals
##################


dat <- move.select%>%
  mutate(m=as.factor(month(Date)),
         y=year(Date),
         rate=dist/timedif)%>%
  select(StudyArea,BearID, Sex, Age2,y,m,d, hii, ndvi, rate, day.night)%>%
  group_by(StudyArea,BearID, Sex, Age2,y,m,d, day.night)%>%
  summarise(rate=mean(rate, na.rm=TRUE),
            hii=mean(hii, na.rm=TRUE),
            ndvi=mean(ndvi, na.rm=TRUE))%>%
  ungroup()%>%
  spread(day.night, rate)%>%
  group_by(StudyArea,BearID, Sex, Age2,y,m,d)%>%
  summarise(Day=mean(Day, na.rm=TRUE),
            Night=mean(Night, na.rm=TRUE),
            hii=mean(hii, na.rm=TRUE),
            ndvi=mean(ndvi, na.rm=TRUE))%>%
  ungroup()%>%
  mutate(noc=(Night/(Night + Day))*100)%>%
  drop_na(noc)

# unique(dat$StudyArea)



##save
write_csv(dat, "relocs_nocturnality_monthly.csv")


library(viridisLite)
dat%>%
  mutate(Age2=round(Age2,0))%>%
  group_by(BearID,Age2)%>%
  summarise(noc=mean(noc),
            hii=mean(hii))%>%
    filter(hii<40)%>%
  ggplot(aes(x=Age2,y=hii, color=noc))+
  geom_jitter(alpha=0.7, size=4)+
  scale_colour_viridis_c(option="C")+
  theme_bw()
```

![](coexistence_sunrise_files/figure-markdown_github/Temporal%20Avoidance-1.png)

``` r
keep <- dat%>%
  group_by(BearID)%>%
  summarise(n=n_distinct(Age2))%>%
  filter(n>=10)%>%
  pull(BearID)

dat <- dat%>%
  filter(BearID%in% keep)

#hist(dat$Age2)

temp.learn.all <- lmer(noc~Age2*hii +(1| BearID) +(1|m) , data =dat,
                       REML = FALSE)
temp.learn.all2 <- lmer(noc~Age2+hii + (1| BearID) +(1|m) , data =dat,
                        REML = FALSE)
temp.learn.all3 <- lmer(noc~Age2+(1| BearID) +(1|m) , data =dat,
                        REML = FALSE)
temp.learn.all4 <- lmer(noc~hii+ (1| BearID) +(1|m) , data =dat,
                        REML = FALSE)
temp.learn.all5 <- lmer(noc~1+(1| BearID) +(1|m) , data =dat,
                        REML = FALSE)
temp.learn.all6 <- lmer(noc~ (1| BearID) +(1|m) , data =dat,
                        REML = FALSE)
temp.learn.all7 <- lmer(noc~Age2*hii +  ndvi +(1| BearID) +(1|m), data =dat,
                        REML = FALSE)



mod.sel <- model.sel(temp.learn.all,temp.learn.all2, temp.learn.all3, temp.learn.all4, 
                     temp.learn.all5, temp.learn.all6, temp.learn.all7,  rank="AICc")


# replace Model name with formulas 
for(i in 1:nrow(mod.sel)) {mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] }

mod.sel <- mod.sel%>%
  select(Model,df,AICc,delta, weight)%>%
  rename(dAICc=delta)%>%
    mutate(Model=str_replace_all(Model,"\\|","]"))%>%
  mutate_if(is.numeric, function(x) round(x, 2))%>%
  as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_col_width(c(4.2,0.5,1,0.7,0.7))%>%
  set_width(0.6)

huxtable::number_format(mod.sel)[, 3] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)

quick_docx(mod.sel, file=here::here("tables", "temp.learn_AIC.docx"), open=FALSE)


##PREDICT

pred.temp.learn<-expand.grid(hii=as.numeric(c(0:64)), m="10",Age2=unique(dat$Age2),
                              BearID="Dr. Mike_A_EV", ndvi=mean(dat$ndvi, na.rm=TRUE))




pred.temp.learn$pred <-predict(temp.learn.all7,newdata=
                                 pred.temp.learn, type="response")


# pred.temp.learn2 <- pred.temp.learn%>%
#   mutate(Age2=case_when(Age2<6~"Subadult",
#                         Age2>10~"Adult"))%>%
#   group_by(hii,Age2)%>%
#   summarise(pred=mean(pred, na.rm=TRUE))%>%
#   ungroup()%>%
#   spread(Age2,pred)%>%
#   mutate(change=(Adult-Subadult)/Subadult*100,
#          Age2="Adult")%>%
#   as.tibble()%>%
#   select(hii,Age2,change)%>%
#   rbind(data.frame(hii=c(0:40),Age2="Subadult",change=0))%>%
#   mutate(Age2=fct_relevel(Age2, levels=c("Subadult", "Adult")))

pred.temp.learn2 <- pred.temp.learn%>%
  mutate(Age2=case_when(Age2<=6~"Subadult",
                        Age2>6~"Adult"))%>%
  group_by(hii,Age2)%>%
  summarise(pred=mean(pred, na.rm=TRUE))%>%
  as.tibble()%>%
    drop_na(Age2)%>%
  mutate(Age2=fct_relevel(Age2, levels=c("Subadult", "Adult")))
  
  
  pred.temp.learn3 <- pred.temp.learn%>%
  mutate(Age2=case_when(Age2<=6~"Subadult",
                        Age2>6~"Adult"))%>%
  group_by(hii,Age2)%>%
  summarise(pred=mean(pred, na.rm=TRUE))%>%
  as.tibble()%>%
    drop_na(Age2)%>%
  spread(Age2,pred)%>%
    mutate(pred.dif=(Adult-Subadult)/Subadult*100,
           Age2="Adult")%>%
    select(hii,pred.dif, Age2)%>%
    rbind(tibble(hii=c(0:64),pred.dif=0, Age2="Subadult"))%>%
  mutate(Age2=fct_relevel(Age2, levels=c("Subadult", "Adult")))
```

``` r
kable(mod.sel[-1,])
```

|     | Model                                         | df  | AICc     | dAICc  | weight |
|-----|:----------------------------------------------|:----|:---------|:-------|:-------|
| 2   | Age2 \* hii + ndvi + (1 \] BearID) + (1 \] m) | 8   | 80006.6  | 0      | 0.75   |
| 3   | Age2 \* hii + (1 \] BearID) + (1 \] m)        | 7   | 80008.82 | 2.21   | 0.25   |
| 4   | hii + (1 \] BearID) + (1 \] m)                | 5   | 80045.23 | 38.62  | 0      |
| 5   | Age2 + hii + (1 \] BearID) + (1 \] m)         | 6   | 80046.4  | 39.8   | 0      |
| 6   | 1 + (1 \] BearID) + (1 \] m)                  | 4   | 80541.92 | 535.31 | 0      |
| 7   | (1 \] BearID) + (1 \] m)                      | 4   | 80541.92 | 535.31 | 0      |
| 8   | Age2 + (1 \] BearID) + (1 \] m)               | 5   | 80543.12 | 536.52 | 0      |

Results: Individual Temporal Avoidance of HII
---------------------------------------------

``` r
plot_model( 
  lmer(noc~Age*hii + ndvi +(1| BearID) +(1|m) , data =dat%>%ungroup()%>%mutate(hii=range01(hii, na.rm=TRUE),
                                                                                         ndvi=range01(ndvi, na.rm=TRUE),
                                                                                         Age2=range01(Age2, na.rm=TRUE))%>%rename(Age=Age2),
       REML = FALSE),
  sort.est = TRUE,
  title = "",
  axis.title=c("Coefficient (std.)")
)+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_Publication()
```

![](coexistence_sunrise_files/figure-markdown_github/Temporal%20Avoidance%20Results-1.png)

``` r
ggsave(here::here("plots", "modelplots","temp.learn_Coeff.png"), width=5, height=4, units="in")

##coefficient table
tbl_regression(temp.learn.all7)

noc.change2 <- ggplot(pred.temp.learn2%>%filter(hii%in%seq(0,40,by=5))%>%mutate(pred=case_when(hii<10~NA_real_, TRUE~pred)),aes(x=pred,y=hii, fill=Age2, color=Age2, group=hii))+
  geom_path(size=0.5) +
  geom_point(colour="black",pch=21,size=2)+
  xlab("Nocturnality (%)")+
  ylab("Human Influence Index")+
  xlim(30,100)+
  ylim(0,43)+
  theme_bw()+
  #ggtitle("Change in time use with age: within individuals")+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_Publication()+
  theme(legend.position = c(0.7,0.1))
noc.change2
```

![](coexistence_sunrise_files/figure-markdown_github/Temporal%20Avoidance%20Results-2.png)

``` r
noc.change <- ggplot(pred.temp.learn3%>%filter(hii%in%seq(0,40,by=5)),aes(x=pred.dif,y=hii, fill=Age2, color=Age2, group=hii))+
  geom_path(size=0.5) +
  geom_point(colour="black",pch=21,size=2)+
  xlab("Nocturnality (% Change)")+
  ylab("Human Influence Index")+
  xlim(-13,80)+
  ylim(0,43)+
  theme_bw()+
  #ggtitle("Change in time use with age: within individuals")+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_Publication()+
  theme(legend.position = c(0.7,0.08))
#noc.change
#ggsave(here::here("plots","Time_Within2.png"), width=3, height=3.5, units="in")
```

Immigration Required
--------------------

``` r
##load results from source-sink analysis
imi <- read_csv(here::here("immi_rqrd.csv"))

imi.rqrd <- imi%>%
  filter(hii%in%seq(0,40,by=5))%>%
  ggplot(aes(x=hii, y=im*100))+
  geom_col(col="#386cb0", fill="#386cb0")+
  coord_flip()+
  ylim(-13,57)+
  xlim(0,43)+
  ylab("Immigration Required (%)")+
  xlab("")+
  scale_fill_Publication()+
  theme_Publication()
imi.rqrd
```

![](coexistence_sunrise_files/figure-markdown_github/imi%20rqrd-1.png)

``` r
ggarrange(noc.change2, hab.change,imi.rqrd,
          ncol = 3, nrow = 1,
          labels=c("D", "E", "F"))
```

![](coexistence_sunrise_files/figure-markdown_github/imi%20rqrd-2.png)

``` r
ggsave(here::here("plots","Fig2b.png"), width=10, height=4, units="in")
```

Nocturnality and Mortality Analysis
-----------------------------------

``` r
dat <- move.select%>%
  mutate(m=as.factor(month(Date)),
         y=year(Date),
         rate=dist/timedif)%>%
  select(StudyArea,BearID, Sex, Age2,m,d, hii, ndvi, rate, day.night)%>%
  group_by(StudyArea,BearID, Sex, Age2,m,d, day.night)%>%
  summarise(rate=mean(rate, na.rm=TRUE),
            hii=mean(hii, na.rm=TRUE),
            ndvi=mean(ndvi, na.rm=TRUE))%>%
  ungroup()%>%
  spread(day.night, rate)%>%
  group_by(StudyArea,BearID, Sex, Age2,m,d)%>%
  summarise(Day=mean(Day, na.rm=TRUE),
            Night=mean(Night, na.rm=TRUE),
            hii=mean(hii, na.rm=TRUE),
            ndvi=mean(ndvi, na.rm=TRUE))%>%
  ungroup()%>%
  mutate(noc=(Night/(Night + Day))*100)%>%
  drop_na(noc)



length(unique(dat$BearID))
##nocturnality
bear.mean <- dat%>%mutate(Age=round(Age2,0))%>%group_by(BearID, Age)%>%summarise(noc2=mean(noc,na.rm=TRUE))
reloc.m.noc <- reloc.m%>%
  filter(BearID%in%unique(dat$BearID))%>%
  mutate(Age2=paste0(Age,".",m-3)%>%as.numeric(),
         timedif=1,
         m=as.factor(m),
         Sex=Sex)%>%
  left_join(dat%>%group_by(BearID,Age2)%>%summarise(noc=mean(noc,na.rm=TRUE)), by=c("BearID", "Age2"))%>%
   left_join(bear.mean, by=c("BearID", "Age"))%>%
   mutate(noc=case_when(is.na(noc)~noc2, 
                        TRUE~noc))%>%
  group_by(BearID,Sex, Age)%>%
  summarise(noc=mean(noc,na.rm=TRUE),
            mort=max(mort,na.rm=TRUE),
            hii=mean(hii,na.rm=TRUE))%>%
  filter(hii>=10)%>%
  mutate(time=365)%>%
  drop_na(noc)

length(unique(reloc.m.noc$BearID))

##fit survival models
nocfit4 <- coxph(Surv(time, mort)~ noc + hii+ 
                   cluster(BearID), data = reloc.m.noc)
nocfit5 <- coxph(Surv(time, mort)~ Age+ noc + hii+ 
                   cluster(BearID), data = reloc.m.noc)
nocfit6 <- coxph(Surv(time, mort)~ noc * hii+ 
                   cluster(BearID), data = reloc.m.noc)
nocfit7 <- coxph(Surv(time, mort)~ Sex+ noc + hii+ 
                   cluster(BearID), data = reloc.m.noc)


###model selection
mod.sel <- model.sel(nocfit4,nocfit5,nocfit6,nocfit7,  rank="AICc")


# replace Model name with formulas 
for(i in 1:nrow(mod.sel)) {mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] }

mod.sel <- mod.sel%>%
  select(Model,df,AICc,delta, weight)%>%
  rename(dAICc=delta)%>%
  mutate_if(is.numeric, function(x) round(x, 2))%>%
  as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_col_width(c(4.2,0.5,1,0.7,0.7))%>%
  set_width(0.6)

huxtable::number_format(mod.sel)[, 3] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)

quick_docx(mod.sel, file=here::here("tables", "MortHazard_Nocturnal_AIC.docx"), open=FALSE)




##make prediction data
nocsurv.pred <- expand.grid(hii=c(10,20, 30, 40),noc=c(0:100), time=as.numeric(365), mort=0, Age=10)

nocsurv.pred$pred <- predict(nocfit4, newdata=nocsurv.pred, type="risk")
nocsurv.pred$pred.se <- predict(nocfit4, newdata=nocsurv.pred, type="risk", se.fit=TRUE)$se.fit

nocsurv.pred$hii <- as.factor(nocsurv.pred$hii)
```

``` r
kable(mod.sel[-1,])
```

|     | Model               | df  | AICc  | dAICc | weight |
|-----|:--------------------|:----|:------|:------|:-------|
| 2   | noc + hii           | 2   | 57.99 | 0     | 0.97   |
| 3   | Sex + noc + hii     | 3   | 66.55 | 8.56  | 0.01   |
| 4   | Age + noc + hii     | 3   | 66.81 | 8.82  | 0.01   |
| 5   | noc + hii + noc:hii | 3   | 67.95 | 9.96  | 0.01   |

Results: Nocturnality and Mortality Analysis
--------------------------------------------

``` r
###plot standardized response

plot_model( 
  coxph(
    Surv(time, mort)~  noc + hii + 
                   cluster(BearID), 
    data = reloc.m.noc%>%ungroup()%>%
      mutate(hii=range01(hii, na.rm=TRUE),
             noc=range01(noc, na.rm=TRUE))),
  sort.est = TRUE,
  title = "",
  axis.title=c("Mortality Hazard")
)+
  geom_hline(yintercept = 1, linetype="dotted")+
  theme_Publication()
```

![](coexistence_sunrise_files/figure-markdown_github/noc+survival%20results-1.png)

``` r
ggsave(here::here("plots", "modelplots","MortHazard_Nocturnal_Coeff.png"), width=6, height=4, units="in")

##coefficient table
tbl_regression(nocfit4)

noc.surv <- ggplot(nocsurv.pred,aes(x=noc,y=pred,fill=hii, color=hii, group=hii))+
  geom_path(size=1.5) +
  geom_ribbon(aes(ymin = pred-pred.se, ymax = pred+pred.se), colour = NA,alpha = 0.3) +
  ylab("Mortality Risk")+
  xlab("Nocturnal (%)")+
  scale_colour_manual(values=c("orange","#a6cee3","#fb9a99","#7fc97f"))+
  scale_fill_manual(values=c("orange","#a6cee3","#fb9a99","#7fc97f"))+
  theme_Publication()+
  theme(legend.position = c(0.65,0.9))
noc.surv
```

![](coexistence_sunrise_files/figure-markdown_github/noc+survival%20results-2.png)

``` r
###Survival increase per unit increase in noc
nocsurv.pred <- expand.grid(hii=c(30), noc=c(0:100), time=as.numeric(365), mort=0)
nocsurv.pred$pred.surv4 <- exp(-predict(nocfit4, newdata=nocsurv.pred, type="expected"))

##slope
noc.slope <- ((max(nocsurv.pred$pred.surv4)-min(nocsurv.pred$pred.surv4))/100)*100




###prep for basin plot
###predict each of 2 top mods
nocsurv.pred <- expand.grid(hii=c(0:50), noc=c(50,100), time=as.numeric(365), mort=0)
nocsurv.pred$pred.surv4 <- exp(-predict(nocfit4, newdata=nocsurv.pred, type="expected"))




ggplot(nocsurv.pred%>%group_by(hii,noc)%>%summarise(pred.surv=mean(pred.surv4)), aes(x=hii, y=pred.surv, color=as.factor(noc), group=as.factor(noc)))+geom_line()+theme_bw()+xlim(0,40)
```

![](coexistence_sunrise_files/figure-markdown_github/noc+survival%20results-3.png)

``` r
nocsurv.pred%>%group_by(hii,noc)%>%summarise(pred.surv=mean(pred.surv4))%>%write_csv(here::here("data","surv_basin.csv"))
```

### For every 1% change in nocturnality, survival increases by 1%

``` r
ggarrange(surv, noc.surv,cod,
          ncol = 3, nrow = 1,
          labels="AUTO")
```

![](coexistence_sunrise_files/figure-markdown_github/plot%20fig2-1.png)

``` r
ggsave(here::here("plots","Fig2a.png"), width=10, height=4, units="in")
```

################################################ 

################################################ 

MAPS
----

################################################ 

################################################ 

Predict Immigration and Nocturnality across space
-------------------------------------------------

Summary stats for maps
----------------------

``` r
##Load Data
range <- st_read(here::here("data","spatial","formap","For_NYT","Current_LambEdits_poly.shp"))%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
ndvi <-raster(here::here("data","spatial","formap","NDVIspring_small.tif"))%>%
  projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
hii.plot <- raster(here::here("data","spatial","formap","hii_small.tif"))%>%
  projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")



#save and load completed version, is slow
ndvi.range <- raster(here::here("data","spatial","formap","NDVIspring_NA.tif"))%>%
  crop(extent(range%>%st_transform("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))%>%
  disaggregate(fact=5)%>%
  projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# saveRDS(ndvi.range, here::here("data","spatial","ndvi.range.rds"))
# ndvi.range <- readRDS(here::here("data","spatial","ndvi.range.rds"))
  
  
hii.range <- raster(here::here("data","spatial","formap","hii_n_amer"))%>%
  crop(extent(range%>%st_transform("+proj=longlat +ellps=clrk66 +no_defs")))%>%
  projectRaster(ndvi.range)
#saveRDS(hii.range, here::here("data","spatial","hii.range.rds"))
# hii.range <- readRDS(here::here("data","spatial","hii.range.rds"))

##create habitat mask, clip by bear range
hab <- raster::mask(ndvi.range,as(st_zm(range), "Spatial"))
hab[hab<=0.2]<-NA


hii.range <- raster::mask(hii.range%>%crop(hab),hab)

###Calculate % area that is sink
sink.perc <-(table(values(hii.range)>=13)[2]/ (table(values(hii.range)>=13)[2]+table(values(hii.range)>=13)[1]))*100
##1.298731


###denote sinks
hii.sink <- hii.range
hii.sink[hii.sink<13]<-NA
hii.sink[hii.sink>=13]<-1

##make poly
hii.sink <- rasterToPolygons(hii.sink, n=4, fun=NULL, na.rm=TRUE)%>%
  st_as_sf()%>%
  summarise(d=mean(layer))

##add 20km buffer
hii.sink20 <- hii.sink%>%
  st_buffer(dist=20000, endCapStyle="ROUND")%>%
  st_intersection(range%>%st_make_valid())

###calculate area affected by sinks w/ 20km buffer.
sinkaffected.perc<-(st_area(hii.sink20)%>%set_units( km^2)/
    st_area(range)%>%set_units( km^2))*100 ##15.28608%


# st_area(range)%>%set_units( km^2)

#clip to southern range
range.south <- st_crop(range%>%st_transform("+proj=longlat +ellps=clrk66 +no_defs"), xmin =-140 , ymin = 43, xmax = -108, ymax = 60)%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
#mapview(range.south)

south.sinkaffected.perc <- (st_area(hii.sink20%>%st_crop(st_bbox(range.south)))%>%set_units( km^2)/
    st_area(range.south)%>%set_units( km^2))*100 
# st_area(range.south)%>%set_units( km^2)

###Calculate % area that is sink
#plot(hii.range>=15)
s.vals <- values(hii.range%>%crop(extent(range.south)))
south.sink.perc <- (table(s.vals>=13)[2]/ (table(s.vals>=13)[2]+table(s.vals>=13)[1]))*100
##3.272988


##calculate # of people inside southern basin of coexistence
hd <- raster(here::here("data","spatial","gpw-v4-population-density_2015.tif"))%>%
  crop(extent(range%>%st_transform("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))%>%
  projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


hdin <- mask(hd,as(hii.sink20%>%st_zm(), "Spatial"))
# 
# plot(hdin, col=rainbow(10))
# plot(as(hii.sink20%>%st_zm(), "Spatial"), add=TRUE)

sum(values(hdin), na.rm=TRUE)* (0.466 * 0.899)/1E6

##below 60 degrees
hdinsouth <- hdin%>%
  crop(extent(c(-2.4E6,-0.95E6,0.5E6,2.4E6)))

# plot(hdinsouth)
# plot(as(hii.sink20%>%st_zm(), "Spatial"), add=TRUE)

sum(values(hdinsouth), na.rm=TRUE)* (0.466 * 0.899)/1E6
```

### 1.7% of the grizzly bear range examined here is sink habitat (population growth&lt;1). This rises to 4.3% in the southern portion of the range.

### 19% of the grizzly bear range examined here is affected by sink habitat through the pull of immigrants from sources. This rises to 46.4% in the southern portion of the range.

### 1.22 million people live within this 19% of the grizzly bear range that is affected by sink habitat through the pull of immigrants from sources

Map it
------

``` r
##load data

world <- ne_countries(scale='medium',returnclass = 'sf')%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
range.exp <- st_read(here::here("data","spatial","formap","For_NYT","Expansion_LambEdits_poly.shp"))%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
  st_make_valid()
hist <- st_read(here::here("data","spatial","formap","For_NYT","Historic_poly.shp"))%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
  st_difference(st_buffer(range,10))
  

range.neg <-st_difference( st_make_valid(world), range.exp)

#range.exp <-st_difference(range), st_make_valid(range.exp))


##prep extents
big = matrix(c(-2.4E6,0.5E6,-2.4E6,2.4E6,-0.95E6,2.4E6,-0.95E6,0.5E6,-2.4E6,0.5E6),ncol=2, byrow=TRUE)
mtn = matrix(c(-1.62E6,1E6,-1.62E6,1.55E6,-1.13E6,1.55E6,-1.13E6,1E6,-1.62E6,1E6),ncol=2, byrow=TRUE)
wild = matrix(c(-1.88E6,2.32E6,-1.88E6,2.5E6,-1.68E6,2.5E6,-1.68E6,2.32E6,-1.88E6,2.32E6),ncol=2, byrow=TRUE)
pts = list(mtn)
pl1 = st_polygon(pts)
pl1 <- st_sfc(pl1, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

pts = list(wild)
pl2 = st_polygon(pts)
pl2 <- st_sfc(pl2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


##standardize to max and min within bear range
hii.plot[hii.plot>46] <- 46
hii.plot[hii.plot<0] <- 0


###NOCTURNAL
noc.mod <- lm(noc~hii*Age2 + ndvi, data =dat)
noc.mod2 <- lm(noc~hii*Age2 + ndvi+I(ndvi^2), data =dat)

model.sel(noc.mod,noc.mod2)

##prep age
Age2 <- hii.plot
Age2[!is.na(Age2)] <- 15

##fix NA's in ndvi-very cloudy areas up in high mtn areas, not bear hab generally, only a few spots
ndvi[is.na(ndvi)] <- 0.3

##assign anything <0.3 as 0.3 (just for mapping purposes, generally mtns, etc.)
values(ndvi)[values(ndvi)<0.3] <- 0.3

##stack
stack <- stack(hii.plot,Age2, ndvi)
names(stack) <- c("hii","Age2","ndvi")


###predict

noc.rast <- predict(stack,noc.mod)

##change to tibble
noc.rast<-  noc.rast%>%
  as.data.frame(xy = TRUE)%>%
  mutate(noc=round(layer,0),
         noc=case_when(noc<10~10,noc>90~90, TRUE~noc))%>%
  drop_na(noc)


noc.map <- ggplot() +
  geom_raster(data = noc.rast,
              aes(x = x, y = y, fill =noc))+
  geom_sf(data = world, fill =NA, color="grey") +
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range.neg,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data=range.exp,fill=NA, color = "purple", size=0.8)+
  geom_sf(data=range,fill =NA, color = "black", size=0.6)+
  geom_sf(data=pl1,fill =NA, color = "grey", size=0.5)+
  geom_sf(data=pl2,fill =NA, color = "grey", size=0.5)+
  #annotation_scale(location = "bl", width_hint = 0.25)+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma", direction=1)+
  #scale_fill_viridis_c(option="cividis")+
  scale_fill_viridis_c(option="viridis")+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]), name = "% Nocturnal")+
  coord_sf(xlim =c(-2.4E6,-0.95E6), ylim = c(0.48E6, 2.55E6), expand=FALSE,
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(linetype = "solid", color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")

noc.map.mtns <- ggplot() +
  geom_raster(data = noc.rast,
              aes(x = x, y = y, fill =noc))+
  geom_sf(data = world, fill =NA, color="grey") +
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range.neg,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data=range.exp,fill=NA, color = "purple", size=0.8)+
  geom_sf(data=range,fill =NA, color = "black", size=0.6)+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma", direction=-1)+
  #scale_fill_viridis_c(option="cividis")+
  #scale_fill_viridis_c(option="viridis", direction=-1)+
  #scale_fill_viridis_c(option="viridis", direction=1)+
  scale_fill_viridis_c(option="viridis")+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]), name = "% Nocturnal")+
  coord_sf(xlim =c(-1.62E6,-1.13E6), ylim = c(1E6, 1.55E6),expand=FALSE,
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(linetype = "solid", color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")


noc.map.wild <- ggplot() +
  geom_raster(data = noc.rast,
              aes(x = x, y = y, fill =noc))+
  geom_sf(data = world, fill =NA, color="grey") +
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range.neg,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data=range.exp,fill=NA, color = "purple")+
  geom_sf(data=range,fill =NA, color = "black")+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]), name = "% Nocturnal")+
  scale_fill_viridis_c(option="viridis")+
  coord_sf(xlim =c(-1.88E6,-1.68E6), ylim = c(2.32E6, 2.5E6), expand=FALSE,
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_rect(fill = "aliceblue"), 
         panel.border = element_rect(linetype = "solid", color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")





###IMI
lmod <- readRDS(here::here("lambdamod", "m1.rds"))

###predict
l.rast <- predict(stack,lmod)


##change to tibble
imi.rast<-  l.rast%>%
  as.data.frame(xy = TRUE)%>%
  mutate(l=round(layer,2),
         imi=case_when((1-l)>0&(1-l)<=0.6~(1-l),
                       (1-l)>0.6~0.6,
                       (1-l)<=0~0),
         l=case_when(l<0.5~0.5,
                     TRUE~l),
         imi=round(imi*100,0))%>%
  drop_na(imi)



imi.map <- ggplot() +
  geom_raster(data = imi.rast,
              aes(x = x, y = y, fill =imi))+
  geom_sf(data = world, fill =NA, color="grey") +
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range.neg,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data=range.exp,fill=NA, color = "purple", size=0.8)+
  geom_sf(data=hii.sink20,fill =NA, color = "tan",size=0.15)+
  geom_sf(data=range,fill =NA, color = "black", size=0.6)+
  geom_sf(data=pl1,fill =NA, color = "grey", size=0.5)+
  geom_sf(data=pl2,fill =NA, color = "grey", size=0.5)+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  scale_fill_viridis_c(option="viridis", direction=1)+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[c(2:7)]), name = "% Immigrants")+
  coord_sf(xlim =c(-2.4E6,-0.95E6), ylim = c(0.48E6, 2.55E6), expand=FALSE,
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(linetype = "solid", color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")

imi.map.mtns <- ggplot() +
  geom_raster(data = imi.rast,
              aes(x = x, y = y, fill =imi))+
  geom_sf(data = world, fill =NA, color="grey") +
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range.neg,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data=range.exp,fill=NA, color = "purple", size=0.8)+
  geom_sf(data=hii.sink20,fill =NA, color ="tan",size=0.4)+
  geom_sf(data=range,fill =NA, color = "black", size=0.9)+
 #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma", direction=)+
  #scale_fill_viridis_c(option="cividis")+
  #scale_fill_viridis_c(option="viridis", direction=-1)+
  scale_fill_viridis_c(option="viridis", direction=1)+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]), name = "% Immigrants")+
  coord_sf(xlim =c(-1.62E6,-1.13E6), ylim = c(1E6, 1.55E6),expand=FALSE,
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_rect(fill = "aliceblue"), 
         panel.border = element_rect(linetype = "solid", color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")

imi.map.wild <- ggplot() +
  geom_raster(data = imi.rast,
              aes(x = x, y = y, fill =imi))+
  geom_sf(data = world, fill =NA, color="grey") +
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range.neg,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data=range.exp,fill=NA, color = "purple", size=0.5)+
  geom_sf(data=hii.sink20,fill =NA, color = "tan",size=0.1)+
  geom_sf(data=range,fill =NA, color = "black")+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  scale_fill_viridis_c(option="viridis", direction=1)+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]), name = "% Immigrants")+
  coord_sf(xlim =c(-1.88E6,-1.68E6), ylim = c(2.32E6, 2.5E6), expand=FALSE,
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25,  text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_rect(fill = "aliceblue"), 
         panel.border = element_rect(linetype = "solid", color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")





# Extract the legend. Returns a gtable
noc.leg <- get_legend(ggplot() +
                        geom_raster(data = noc.rast,
                                    aes(x = x, y = y, fill =noc))+
                        scale_fill_viridis_c(option="viridis", direction=1)+
                        theme(legend.position = "top",
                              legend.title = element_blank())
)%>%as_ggplot()

imi.leg <- get_legend(ggplot() +
                        geom_raster(data = imi.rast,
                                    aes(x = x, y = y, fill =imi))+
                   scale_fill_viridis_c(option="viridis", direction=1)+
                        theme(legend.position = "top",
                              legend.title = element_blank())
)%>%as_ggplot()



##save individidually
ggsave(here::here("plots","Fig3", "imi.map.png"), width=2.5, height=3, units="in", bg="transparent", imi.map)
ggsave(here::here("plots","Fig3", "imi.map.mtns.png"), width=3, height=3, units="in", bg="transparent", imi.map.mtns)
ggsave(here::here("plots","Fig3", "imi.map.wild.png"), width=1.5, height=1.5, units="in", bg="transparent", imi.map.wild)
ggsave(here::here("plots","Fig3", "imi.leg.png"), width=1.5, height=0.4, units="in", bg="transparent", imi.leg)

ggsave(here::here("plots","Fig3", "noc.map.png"), width=2.5, height=3, units="in", bg="transparent", noc.map)
ggsave(here::here("plots","Fig3", "noc.map.mtns.png"), width=3, height=3, units="in", bg="transparent", noc.map.mtns)
ggsave(here::here("plots","Fig3", "noc.map.wild.png"), width=1.5, height=1.5, units="in", bg="transparent", noc.map.wild)
ggsave(here::here("plots","Fig3", "noc.leg.png"), width=1.5, height=0.6, units="in", bg="transparent", noc.leg)


ggarrange(imi.map,noc.map)
```

![](coexistence_sunrise_files/figure-markdown_github/Map:%20Immigration%20and%20Nocturnality-1.png)

Study Area Map
--------------

``` r
##Load Data and cleanup
sa_mcp <- st_read(here::here("data", "spatial", "telem_dna_sa_merge.shp"))
  

st_area(sa_mcp%>%st_buffer(dist=30E3))%>%  ##SCR mask of integration
  set_units("km^2")

# hii.big<- raster(here::here("data","spatial","formap","hii_n_amer"))%>%
#   aggregate(fact=10)%>%
#   projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# 
# #-140, -107, 41, 61
# #writeRaster(hii.big, '/Users/clayton.lamb/Google Drive/Documents/University/Geographic_Data/Human/hii_n_america_grid-2/hii_clip.tif', overwrite=TRUE)
# 
# hii.big<- as.data.frame(hii.big, xy = TRUE)%>%
#   filter(!is.na(hii_n_amer))%>%
#   mutate(hii_n_amer=case_when(hii_n_amer>40~40,
#                               TRUE~hii_n_amer))
# saveRDS(hii.big, here::here("data","spatial","hii.big.rds"))
hii.big <- readRDS(here::here("data","spatial","hii.big.rds"))

##prep extents
koot.extent = matrix(c(-1.5E6,1.12E6,-1.5E6,1.37E6,-1.2E6,1.37E6,-1.2E6,1.12E6,-1.5E6,1.12E6),ncol=2, byrow=TRUE)
coast.extent = matrix(c(-2.1E6,1.35E6,-2.1E6,2.4E6,-1.51E6,2.4E6,-1.51E6,1.35E6,-2.1E6,1.35E6),ncol=2, byrow=TRUE)
pts = list(koot.extent)
pl1 = st_polygon(pts)
pl1 <- st_sfc(pl1, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

pts = list(coast.extent)
pl2 = st_polygon(pts)
pl2 <- st_sfc(pl2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

mask <- st_make_valid(world)%>%
  st_difference(st_union(st_make_valid(range),st_make_valid(hist)))%>%
  select(continent)%>%
  filter(continent%in%"North America")%>%
  group_by(continent)%>%
  summarise()
  

##Prep maps
big <-ggplot() +
  geom_raster(data = hii.big,
              aes(x = x, y = y, fill =hii_n_amer))+
  geom_sf(data = world, fill =NA, color="grey90") +
   geom_sf(data=mask,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed", size=0.7) +
  geom_sf(data=range,fill =NA, color = "black", size=1)+
  geom_sf(data=pl1,fill =NA, color = "grey", size=0.9)+
  geom_sf(data=pl2,fill =NA, color = "grey", size=0.9)+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  scale_fill_viridis_c(option="viridis", direction=1)+
  geom_sf(data=sa_mcp, fill=NA,col="white", size=0.3)+
  #geom_point(data=reloc.df, aes(x = X, y = Y), fill=NA,col=viridis(5)[5], size=0.2)+
  #geom_polygon(data=bear_mcp , aes(x = long, y = lat, group = group), fill=NA,col="red", lwd=0.2)+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]))+
  coord_sf(xlim = c(-3.8E6, 2E6), ylim = c(-1.5E6, 4E6), crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  xlab("Longitude")+ ylab("Latitude")+
  annotation_scale(location = "br", width_hint = 0.25)+
  theme(axis.ticks=element_blank(),
        panel.background = element_blank(), 
        panel.border =  element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")

# hii.koot<- raster(here::here("data","spatial","formap","hii_n_amer"))%>%
#   crop(extent(c(-122, -110, 47, 54)))%>%
#   projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
#   as.data.frame(xy = TRUE)%>%
#   filter(!is.na(hii_n_amer))%>%
#   mutate(hii_n_amer=case_when(hii_n_amer>40~40,
#                               TRUE~hii_n_amer))
# saveRDS(hii.koot, here::here("data","spatial","hii.koot.rds"))
hii.koot <- readRDS(here::here("data","spatial","hii.koot.rds"))

koot <- ggplot() +
  geom_raster(data = hii.koot,
              aes(x = x, y = y, fill =hii_n_amer))+
  geom_sf(data = world, fill =NA, color="grey90") +
  geom_sf(data=mask,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range,fill =NA, color = "black", size=1)+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  scale_fill_viridis_c(option="viridis", direction=1)+
  geom_sf(data=sa_mcp, fill=NA,col="white", size=0.5)+
  #geom_polygon(data=bear_mcp , aes(x = long, y = lat, group = group), fill=NA,col="red", lwd=0.2)+
  #geom_point(data=reloc.df, aes(x = X, y = Y), fill=NA,col=viridis(5)[5], size=0.3)+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]))+
  coord_sf(xlim = c(-1.5E6,-1.2E6), ylim = c(1.12E6, 1.45E6),
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")


# hii.coast<- raster(here::here("data","spatial","formap","hii_n_amer"))%>%
#   crop(extent(c(-133, -110, 47, 60)))%>%
#   projectRaster(crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
#   as.data.frame(xy = TRUE)%>%
#   filter(!is.na(layer))%>%
#   filter(!is.na(hii_n_amer))%>%
#   mutate(hii_n_amer=case_when(hii_n_amer>40~40,
#                               TRUE~hii_n_amer))
#saveRDS(hii.coast, here::here("data","spatial","hii.coast.rds"))
hii.coast <- readRDS(here::here("data","spatial","hii.coast.rds"))


coast <- ggplot() +
  geom_raster(data = hii.coast,
              aes(x = x, y = y, fill =hii_n_amer))+
  geom_sf(data = world, fill =NA, color="grey30") +
  #geom_sf(data=mask,fill="grey", color = NA, alpha=0.3)+
  geom_sf(data = hist, fill =NA, color="black",linetype="dashed") +
  geom_sf(data=range,fill =NA, color = "black", size=1)+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  scale_fill_viridis_c(option="viridis", direction=1)+
  geom_sf(data=sa_mcp, fill=NA,col="white", size=0.5)+
  #geom_point(data=reloc.df, aes(x = X, y = Y), fill=NA,col=viridis(5)[5], size=0.3)+
  #geom_polygon(data=bear_mcp , aes(x = long, y = lat, group = group), fill=NA,col="red", lwd=0.2)+
  #scale_fill_gradientn(colours = rev(magma(8, alpha = 0.8)[2:7]))+
  coord_sf(xlim = c(-2.1E6, -1.51E6), ylim = c(1.35E6, 2.4E6),
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  annotation_scale(location = "bl", width_hint = 0.25, text_col="white")+
  theme(axis.ticks=element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(color="black",fill = NA,size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")



ggsave(here::here("plots","fig1","Fig1a.png"), width=3, height=3, units="in", koot, bg="transparent")
ggsave(here::here("plots","fig1","Fig1b.png"), width=6, height=6, units="in", big, bg="transparent")
ggsave(here::here("plots","fig1","Fig1c.png"), width=2.5, height=6, units="in", coast, bg="transparent")


ggarrange(coast, big, NULL, koot, 
          ncol = 2, nrow = 2, 
          widths = c(2, 3), heights = c(3, 1.55))
```

![](coexistence_sunrise_files/figure-markdown_github/Study%20Area%20Map-1.png)

WARP human wildlife conflict
----------------------------

``` r
hwc <- read_csv(here::here("data", "conflict", "encounter-265616-177058.csv"))%>%
  st_as_sf(coords = c("encounter_lng","encounter_lat"),
           crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")

#mapview(hwc)


clip <-st_intersection(st_make_valid(range%>%st_buffer(dist=5000)),
                       st_read("/Users/clayton.lamb/Google Drive/Documents/University/Geographic_Data/Administrative_Boundaries/Canada_Boundaries/BC_Boundary.shp")%>%
                         st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 "))
clip$group <- 1
clip <- clip%>%
  group_by(group)%>%
  summarise(Shape_Area=mean(Shape_Area, na.rm=TRUE))


##prep
ndvi2 <- ndvi%>%mask(clip)%>%aggregate(11)
hii2 <- hii.plot%>%mask(clip)%>%aggregate(11)

###
blank <- ndvi2
values(blank) <-0
hwc$count <- 1
r=rasterize(as(hwc%>%filter(species_name=="GRIZZLY BEAR"),"Spatial"),blank,field="count",fun=sum, update=TRUE)%>%
  rasterToPoints(spatial=TRUE)%>%
  st_as_sf()

##make velox (faster)
hiiv <- velox(hii2)

##extract hii to reloc
r <- r%>%mutate(hii=hiiv$extract_points(sp = as(r, "Spatial")))%>%
  drop_na(hii)

rm(hii2)

## NDVI
##make velox (faster)
ndviv <- velox(ndvi2)

##extract hii to reloc
r <- r%>%mutate(ndvi=ndviv$extract_points(sp = as(r, "Spatial")))%>%
  drop_na(ndvi)

rm(ndvi2)
r <- r%>%mutate(layer=layer/3,ndvi=as.numeric(ndvi), hii=as.numeric(hii))


##mods
hwc.mod <- glm(layer~hii + I(hii^2), data=r)
hwc.mod1 <- glm(layer~hii, data=r)
hwc.mod2 <- glm(layer~hii + I(hii^2) + ndvi, data=r)
hwc.mod3 <- glm(layer~hii +ndvi, data=r)
hwc.mod4 <- glm(layer~1, data=r)

mod.sel <-model.sel(hwc.mod,hwc.mod1,hwc.mod2,hwc.mod3,hwc.mod4)

##coefficient table
tbl_regression(hwc.mod)

# replace Model name with formulas 
for(i in 1:nrow(mod.sel)) {mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] }

mod.sel <- mod.sel%>%
  select(Model,df,AICc,delta, weight)%>%
  rename(dAICc=delta)%>%
  mutate_if(is.numeric, function(x) round(x, 2))%>%
  as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_col_width(c(4,0.5,1,0.9,0.7))%>%
  set_width(0.6)

huxtable::number_format(mod.sel)[, 3] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)
huxtable::number_format(mod.sel)[, 4] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)

quick_docx(mod.sel, file=here::here("tables", "warp.conflict_AIC.docx"), open=FALSE)
```

``` r
kable(mod.sel[-1,], espace=FALSE)
```

|     | Model                 | df  | AICc     | dAICc   | weight |
|-----|:----------------------|:----|:---------|:--------|:-------|
| 2   | hii + I(hii^2) + ndvi | 5   | 18015.94 | 0       | 0.58   |
| 3   | hii + I(hii^2)        | 4   | 18016.56 | 0.62    | 0.42   |
| 4   | hii + ndvi            | 4   | 19342.96 | 1327.01 | 0      |
| 5   | hii                   | 3   | 19358.67 | 1342.73 | 0      |
| 6   | 1                     | 2   | 20188.3  | 2172.36 | 0      |

``` r
##make prediction data
conf.pred <- expand.grid(ndvi=0.6, hii=c(0:40))
wts <- Weights(AIC(hwc.mod, hwc.mod2))
conf.pred$pred1 <-  predict(hwc.mod,conf.pred)
conf.pred$pred2 <-  predict(hwc.mod2,conf.pred)
 conf.pred$pred <- (conf.pred$pred1*wts[1]) + (conf.pred$pred2*wts[2])
 conf.pred$se <- (predict(hwc.mod,conf.pred,se=TRUE)$se.fit*wts[1]) + (predict(hwc.mod2,conf.pred,se=TRUE)$se.fit*wts[2])
 
 conf.pred%>%
  ggplot( aes(x=hii, y=pred)) + 
  geom_line()+
  geom_ribbon(aes(x=hii, ymin=pred-se, ymax=pred+se), alpha=0.5)+
  xlab("Human Influence Index")+
  ylab("Conflict Calls (/yr/100km2)")+
  theme_Publication()
```

![](coexistence_sunrise_files/figure-markdown_github/human%20wildlife%20conflict%20plot-1.png)

``` r
 ggsave(here::here("plots","conflictcalls_hii.png"), width=5, height=3, units="in")
```

human wildlife conflict: bear perspective
-----------------------------------------

``` r
###bear view
lambbears <- move.select%>%
  filter(StudyArea%in%c("Elk Valley-Lamb"))%>%
  distinct(BearID)%>%
  pull(BearID)

noc.conflict <-dat %>%
  filter(BearID%in%lambbears)%>%
  mutate(conflict=case_when(BearID%in%c("Donnie_L_EV", "EVGM88_L_EV", "Matt_L_EV", "Sid_L_EV", "Kathy_L_EV")~1, ##bears who have been in conflict
                            TRUE~0))%>%
  group_by(conflict, hii, Age2)%>%
  summarise(noc=mean(noc), ndvi=mean(ndvi))


  noc.conflict1 <- glm(conflict~1, data=noc.conflict, family="binomial")
  noc.conflict1.1 <- glm(conflict~hii, data=noc.conflict, family="binomial")
   noc.conflict1.2 <- glm(conflict~noc, data=noc.conflict, family="binomial")
  noc.conflict2 <- glm(conflict~hii+noc, data=noc.conflict, family="binomial")
  noc.conflict3 <- glm(conflict~hii+noc + hii*noc, data=noc.conflict, family="binomial")
  noc.conflict4 <- glm(conflict~hii+noc + hii*noc + ndvi, data=noc.conflict, family="binomial")
    noc.conflict5 <- glm(conflict~hii+noc + ndvi, data=noc.conflict, family="binomial")

mod.sel <- model.sel(noc.conflict1,noc.conflict1.1,noc.conflict1.2,noc.conflict2,noc.conflict3,noc.conflict4,noc.conflict5)

##coefficient table
tbl_regression(noc.conflict4,estimate_fun = function(x) style_ratio(x, digits = 3))
summary(noc.conflict4)
# replace Model name with formulas 
for(i in 1:nrow(mod.sel)) {mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] }

mod.sel <- mod.sel%>%
  select(Model,df,AICc,delta, weight)%>%
  rename(dAICc=delta)%>%
  mutate_if(is.numeric, function(x) round(x, 2))%>%
  as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_col_width(c(4.2,0.5,1,0.7,0.7))%>%
  set_width(0.6)

huxtable::number_format(mod.sel)[, 3] <- list(
  function(x)
    prettyNum(x, big.mark = ",",
              scientific = FALSE)
)

quick_docx(mod.sel, file=here::here("tables", "bear.conflict_AIC.docx"), open=FALSE)
```

``` r
kable(mod.sel[-1,], espace=FALSE)
```

|     | Model                         | df  | AICc   | dAICc | weight |
|-----|:------------------------------|:----|:-------|:------|:-------|
| 2   | hii + noc + hii \* noc + ndvi | 5   | 707.36 | 0     | 1      |
| 3   | hii + noc + ndvi              | 4   | 720.74 | 13.39 | 0      |
| 4   | hii + noc + hii \* noc        | 4   | 735.56 | 28.2  | 0      |
| 5   | hii + noc                     | 3   | 745.56 | 38.2  | 0      |
| 6   | hii                           | 2   | 756.55 | 49.19 | 0      |
| 7   | 1                             | 1   | 793.6  | 86.24 | 0      |
| 8   | noc                           | 2   | 794.69 | 87.33 | 0      |

``` r
##make prediction data
conf.pred <- expand.grid(ndvi=0.6, hii=c(10,20,30,40),noc=c(0:100))
wts <- Weights(AIC(noc.conflict4, noc.conflict5))
conf.pred$pred1 <-  predict(noc.conflict4,conf.pred, type="response")
conf.pred$pred2 <- predict(noc.conflict5,conf.pred, type="response")
 conf.pred$pred <- (conf.pred$pred1*wts[1]) + (conf.pred$pred2*wts[2])
conf.pred$se <- (predict(noc.conflict4,conf.pred, type="response", se=TRUE)$se.fit*wts[1]) + (predict(noc.conflict5,conf.pred, type="response", se=TRUE)$se.fit*wts[2])
 

  
conf.pred%>%
  mutate(HII=as.factor(hii))%>%
  ggplot( aes(x=noc, y=pred, group=as.factor(hii))) +
  geom_ribbon(aes(x=noc, ymin=pred-se, ymax=pred+se, fill=HII), alpha=0.5)+
  geom_line()+
  xlab("Nocturnal (%)")+
  ylab("p(conflict)")+
  theme_Publication()+
  theme(legend.title = element_text(),
        legend.position = "right",
         legend.direction = "vertical")+
    scale_colour_manual(values=c("orange","#a6cee3","#fb9a99","#7fc97f"))+
  scale_fill_manual(values=c("orange","#a6cee3","#fb9a99","#7fc97f"))
```

![](coexistence_sunrise_files/figure-markdown_github/human%20wildlife%20conflict:%20bear%20perspective%20plot-1.png)

``` r
ggsave(here::here("plots","conflict_nocturnality.png"), width=5, height=3.5, units="in")

##how nuch conflict is adult nocturnality reducing? predict using 57% change in nocturnality between subdadult and adult, as predicted from nocturnality analysis above (HII=40).
pred <- conf.pred%>%filter(hii==40 & noc%in%c(41,93))%>%pull(pred)
  


reduce <- round(((pred[1]-pred[2])/pred[1])*100,0)
```

### Using a subset of 40 GPS-collared bears for which we have maintained records of conflict incidents with people (as reported by the BC govt), we show, in the highest human influenced areas, adults’ increases in nocturnality reduces their chance of conflict by 71% compared to subadults

additional analyses following revision 1: n bears inside sink
-------------------------------------------------------------

``` r
dens <- raster(here::here("secr_hii_oscr","outputs","dens_surf.tif"))
#plot(dens)

##clip to southern range (<60 degrees latitude-most appropriate extrapolation area given the data used)
dens <- mask(dens, range.south%>%st_transform("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")%>%as("Spatial"))
plot(dens)
```

![](coexistence_sunrise_files/figure-markdown_github/additional%20analyses%20following%20revision%201:%20n%20bears%20inside%20sink-1.png)

``` r
##how many bears in area?
sum(values(dens), na.rm=TRUE)

sink.bears <- mask(dens, 
                   hii.sink%>%
                       st_transform("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")%>%
                     st_buffer(6000)%>% ##capture those edge animals
                     as("Spatial"))

plot(sink.bears)
```

![](coexistence_sunrise_files/figure-markdown_github/additional%20analyses%20following%20revision%201:%20n%20bears%20inside%20sink-2.png)

``` r
###what % of bears experience sink?
sum(values(sink.bears),na.rm=TRUE)/sum(values(dens), na.rm=TRUE)

###n bears outside of sink
sum(values(dens), na.rm=TRUE)-sum(values(sink.bears),na.rm=TRUE)

###n bears inside
sum(values(sink.bears),na.rm=TRUE)
```

### 4158 (or 16.1%) of brown bears in the southern range live in, or are exposed to, sink habitats (home range center &lt; 6 km from sink). There are approximately 21633 bears outside the sink areas.

Check sensitivity of results to alternate paramaterizations
===========================================================

``` r
################
################
###Do results change with mixed effects model?? no
################
################
library(coxme)
surv.me <- reloc.m%>%filter(ageclass%in%c("Adult","Subadult"))%>%drop_na()

fitme1 <- coxme(Surv(time, mort)~ Sex +
                (1|StudyArea), data = surv.me)
fitme1.1 <- coxme(Surv(time, mort)~ ageclass+
                  (1|StudyArea), data = surv.me)
fitme1.2 <- coxme(Surv(time, mort)~ Sex + ageclass +
                  (1|StudyArea), data =surv.me)
fitme2 <- coxme(Surv(time, mort)~  Sex * ageclass + 
                (1|StudyArea), data = surv.me)
fitme3 <- coxme(Surv(time, mort)~ Sex + ageclass + hii +
                (1|StudyArea), data = surv.me)
fitme4 <- coxme(Surv(time, mort)~ (Sex+ageclass*hii) +
                (1|StudyArea), data = surv.me)
fitme5 <- coxme(Surv(time, mort)~ (Sex*hii + ageclass) +
                (1|StudyArea), data = surv.me)
fitme6 <- coxme(Surv(time, mort)~ (Sex*ageclass*hii) +
                (1|StudyArea), data = surv.me)
fitme7 <- coxme(Surv(time, mort)~ (Sex + Age*hii) +
                (1|StudyArea), data = surv.me)


model.sel(fitme1, fitme1.1, fitme1.2,fitme2,fitme3,fitme4,fitme5,fitme6,fitme7)



###predict each of 3 top mods
surv.me$pred3 <- predict(fitme3, type="risk")

surv.me$pred4 <- predict(fitme4, type="risk")

surv.me$pred6 <- predict(fitme6, type="risk")



##average each model prediction by model weight
wts <- Weights(AIC(fitme3, fitme4, fitme6))
surv.me$pred <- (surv.me$pred3*wts[[1]])+(surv.me$pred4*wts[[2]])+(surv.me$pred6*wts[[3]])



##relevel factors for plotting
surv.me <- surv.me%>%
  ungroup()%>%
  mutate(ageclass=fct_relevel(ageclass,levels=c("Subadult", "Adult")),
         Sex=fct_relevel(Sex,levels=c("M", "F")))

##summarize across sexes to produce simpler plot
surv.me<- surv.me%>%
  group_by(ageclass, hii)%>%
  summarise(pred=mean(pred,na.rm=TRUE))



#ranef(fitme6)

ggplot(surv.me,aes(x=hii,y=pred,fill=ageclass, color=ageclass))+
  geom_point()+
  geom_smooth(size=2, method = "loess") + 
  ylab("Risk Factor (x above baseline)")+
  xlab("Human Influence Index")+
  xlim(0,33)+
  ylim(0,30)+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_bw() +theme(plot.title = element_text(face = "bold",
                                              size = rel(1.2), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1.3)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(size = rel(1.1)), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.key.size= unit(0.6, "cm"),
                    legend.margin = unit(0.1, "cm"),
                    legend.text = element_text(size = rel(1)),
                    legend.title = element_blank(),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold"))
```

![](coexistence_sunrise_files/figure-markdown_github/check%20sensitivity%20of%20results-1.png)

``` r
ggsave(here::here("plots","ranef_hazard.png"), width=4, height=5, units="in")



################
################
###Do results change with road density? no
################
################

##################################################################
##Add Spatial Data
##################################################################
##make bear data spatial
reloc.rd <- reloc%>%st_as_sf()

## HUMAN INFLUENCE INDEX
##make smaller
rd<- raster(here::here("data","spatial","rd_USEdens_8km.tif"))%>%
  projectRaster(crs="+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") 


##make velox (faster)
rdv <- velox(rd)


##extract hii to reloc
reloc.rd <- reloc.rd%>%mutate(rd=rdv$extract_points(sp = reloc.rd))%>%
  as_tibble()

rm(rd)

##################################################################
##Smooth over monthly periods
##################################################################

##pull out exact mort locs
kills.rd <- reloc.rd%>%filter(mort%in%1)

##summarize live locs by month
reloc.rd.m <- reloc.rd%>%  
  filter(mort%in%0)%>%
  group_by(BearID,yr,m,Sex,Age, ageclass,StudyArea, mort)%>%
  summarize(rd=mean(rd,na.rm=TRUE))%>%
  as_tibble()

##assign 3 days to each, as they are each a month
reloc.rd.m$time <- 30

##length of survival on last month of life= day of month
kills.rd$time<-kills.rd$d

##match columns
kills.rd <- kills.rd%>%select(colnames(reloc.rd.m))

##bind and join last month with live and dead if needed
reloc.rd.m <- rbind(reloc.rd.m,kills.rd)%>%
  group_by(BearID,Sex,Age, ageclass, StudyArea,yr,m)%>%
  summarize(mort=max(mort, na.rm=TRUE),
            rd=mean(rd, na.rm=TRUE),
            time=min(time, na.rm=TRUE))



##################################################################
##Fit Models
##################################################################

fit1 <- coxph(Surv(time, mort)~ Sex +
                cluster(BearID), data = reloc.rd.m)
fit1.1 <- coxph(Surv(time, mort)~ ageclass+
                  cluster(BearID), data = reloc.rd.m)
fit1.2 <- coxph(Surv(time, mort)~ Sex + ageclass +
                  cluster(BearID), data = reloc.rd.m)
fit2 <- coxph(Surv(time, mort)~  Sex * ageclass + 
                cluster(BearID), data = reloc.rd.m)
fit3 <- coxph(Surv(time, mort)~ Sex + ageclass + rd +
                cluster(BearID), data = reloc.rd.m)
fit4 <- coxph(Surv(time, mort)~ (Sex+ageclass*rd) +
                cluster(BearID), data = reloc.rd.m)
fit5 <- coxph(Surv(time, mort)~ (Sex*rd + ageclass) +
                cluster(BearID), data = reloc.rd.m)
fit6 <- coxph(Surv(time, mort)~ (Sex*ageclass*rd) +
                cluster(BearID), data = reloc.rd.m)
fit7 <- coxph(Surv(time, mort)~ (Sex + Age*rd) +
                cluster(BearID), data = reloc.rd.m)

model.sel(fit1.1, fit1.2, fit1,fit2, fit3, fit4, fit5, fit6, fit7)


##################################################################
##Test Assumptions
##################################################################
cox.zph(fit3)
#ggcoxzph(cox.zph(fit3))

##################################################################
##Plot Effects, weighted by top models
##################################################################

##make prediction data
surv.pred <- expand.grid(Sex=c("M","F"), ageclass=c("Subadult","Adult"), rd=c(0:6))

###predict each of 3 top mods
surv.pred$pred3 <- predict(fit3, newdata=surv.pred, type="risk")
surv.pred$pred.se3 <- predict(fit3, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

surv.pred$pred4 <- predict(fit4, newdata=surv.pred, type="risk")
surv.pred$pred.se4 <- predict(fit4, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

surv.pred$pred6 <- predict(fit6, newdata=surv.pred, type="risk")
surv.pred$pred.se6 <- predict(fit6, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

surv.pred$pred5 <- predict(fit5, newdata=surv.pred, type="risk")
surv.pred$pred.se5 <- predict(fit5, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit

##average each model prediction by model weight
wts <- Weights(AIC(fit3, fit4, fit6, fit5))
surv.pred$pred <- (surv.pred$pred3*wts[[1]])+(surv.pred$pred4*wts[[2]])+(surv.pred$pred6*wts[[3]])+(surv.pred$pred5*wts[[4]])
surv.pred$pred.se <- (surv.pred$pred.se3*wts[[1]])+(surv.pred$pred.se4*wts[[2]])+(surv.pred$pred.se6*wts[[3]])+(surv.pred$pred.se5*wts[[4]])


##relevel factors for plotting
surv.pred <- surv.pred%>%
  mutate(ageclass=fct_relevel(ageclass,levels=c("Subadult", "Adult")),
         Sex=fct_relevel(Sex,levels=c("M", "F")))

##summarize across sexes to produce simpler plot
surv.pred<- surv.pred%>%
  group_by(ageclass, rd)%>%
  summarise(pred=mean(pred,na.rm=TRUE),
            pred.se=mean(pred.se, na.rm=TRUE))

# ##ratio of hazard between adults and subads
# surv.pred%>%spread(ageclass,pred)%>%
#   group_by(hii)%>%
#   summarize(s=mean(Subadult, na.rm=TRUE),
#             a=mean(Adult, na.rm=TRUE))%>%
#   mutate(dif=s/a)%>%
#   ggplot(aes(x=hii,y=dif))+geom_line()


ggplot(surv.pred,aes(x=rd,y=pred,fill=ageclass, color=ageclass))+
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin = pred-pred.se, ymax = pred+pred.se), colour = NA,alpha = 0.3) +
  ylab("Mortality Risk (x above baseline)")+
  xlab("Road Density + Use Index")+
  #ylim(0,23)+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_bw() +theme(plot.title = element_text(face = "bold",
                                              size = rel(1.2), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1.3)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(size = rel(1.1)), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.key.size= unit(0.6, "cm"),
                    legend.margin = unit(0.1, "cm"),
                    legend.text = element_text(size = rel(1)),
                    legend.title = element_blank(),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold"))
```

![](coexistence_sunrise_files/figure-markdown_github/check%20sensitivity%20of%20results-2.png)

``` r
ggsave(here::here("plots","rd_hazard.png"), width=4, height=5, units="in")


library(ggridges)
ggplot(reloc.rd.m, aes(x = rd, y = StudyArea)) + 
  geom_density_ridges( jittered_points = TRUE,
                       position = position_points_jitter(width = 0.05, height = 0),
                       point_shape = '|', point_size = 1.5, point_alpha = 1, alpha = 0.7)+
  xlab("Road Density + Use Index")+
  theme_ridges()
```

![](coexistence_sunrise_files/figure-markdown_github/check%20sensitivity%20of%20results-3.png)

``` r
ggsave(here::here("plots","rd_use_byproj.png"), width=10, height=6.5, units="in")










################
################
###Do results change with weekly analyses? no, and also violate model assumptions, so proceed with a month
################
################

##################################################################
##Smooth over WEEKLY periods
##################################################################

##pull out exact mort locs
kills <- reloc%>%  
  mutate(w=week(Date))%>%filter(mort%in%1)

##summarize live locs by week
reloc.m <- reloc%>%  
  mutate(w=week(Date))%>%
  filter(mort%in%0)%>%
  group_by(BearID,yr,w,Sex,Age, ageclass,StudyArea, mort)%>%
  summarize(hii=mean(hii,na.rm=TRUE))%>%
  as_tibble()

##assign 7 days to each, as they are each a week
reloc.m$time <- 7

##length of survival on last month of life= day of month
kills$time<-kills$d

##match columns
kills <- kills%>%select(colnames(reloc.m))

##bind and join last month with live and dead if needed
reloc.m <- rbind(reloc.m,kills)%>%
  group_by(BearID,Sex,Age, ageclass, StudyArea,yr,w)%>%
  summarize(mort=max(mort, na.rm=TRUE),
            hii=mean(hii, na.rm=TRUE),
            time=min(time, na.rm=TRUE))



##################################################################
##Fit Models
##################################################################

fit1 <- coxph(Surv(time, mort)~ Sex +
                cluster(BearID), data = reloc.m)
fit1.1 <- coxph(Surv(time, mort)~ ageclass+
                  cluster(BearID), data = reloc.m)
fit1.2 <- coxph(Surv(time, mort)~ Sex + ageclass +
                  cluster(BearID), data = reloc.m)
fit2 <- coxph(Surv(time, mort)~  Sex * ageclass + 
                cluster(BearID), data = reloc.m)
fit3 <- coxph(Surv(time, mort)~ Sex + ageclass + hii +
                cluster(BearID), data = reloc.m)
fit4 <- coxph(Surv(time, mort)~ (Sex+ageclass*hii) +
                cluster(BearID), data = reloc.m)
fit5 <- coxph(Surv(time, mort)~ (Sex*hii + ageclass) +
                cluster(BearID), data = reloc.m)
fit6 <- coxph(Surv(time, mort)~ (Sex*ageclass*hii) +
                cluster(BearID), data = reloc.m)


model.sel(fit1.1, fit1.2, fit1,fit2, fit3, fit4, fit5, fit6)


##################################################################
##Test Assumptions
##################################################################
cox.zph(fit3)



##make prediction data
surv.pred <- expand.grid(Sex=c("M","F"), ageclass=c("Subadult","Adult"), hii=c(0:40), ndvi=mean(reloc.m$ndvi, na.rm=TRUE))

surv.pred$pred <- predict(fit4, newdata=surv.pred, type="risk")
surv.pred$pred.se <- predict(fit4, newdata=surv.pred, type="risk", se.fit=TRUE)$se.fit


##relevel factors for plotting
surv.pred <- surv.pred%>%
  mutate(ageclass=fct_relevel(ageclass,levels=c("Subadult", "Adult")),
         Sex=fct_relevel(Sex,levels=c("M", "F")))

##summarize across sexes to produce simpler plot
surv.pred<- surv.pred%>%
  group_by(ageclass, hii)%>%
  summarise(pred=mean(pred,na.rm=TRUE),
            pred.se=mean(pred.se, na.rm=TRUE))

ggplot(surv.pred,aes(x=hii,y=pred,fill=ageclass, color=ageclass))+
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin = pred-pred.se, ymax = pred+pred.se), colour = NA,alpha = 0.3) +
  ylab("Mortality Risk")+
  xlab("Human Influence Index")+
  xlim(0,40)+
  #ylim(0,23)+
  scale_colour_manual(values=c("#386cb0", "#fdb462"))+
  scale_fill_manual(values=c("#386cb0", "#fdb462"))+
  theme_Publication()+
  theme(legend.position = c(0.4,0.9))
```

![](coexistence_sunrise_files/figure-markdown_github/check%20sensitivity%20of%20results-4.png)

``` r
ggsave(here::here("plots","weekly_hazard.png"), width=4, height=5, units="in")
```

additional analyses following revision 1: nocturnality and movement rates across space and hii
----------------------------------------------------------------------------------------------

``` r
noc.dat<-read_csv("relocs_nocturnality_monthly.csv")

# 
# ###model
# m1 <-glm(noc~ndvi+Age2+hii + hii*Age2 + I(hii^2) +Sex, data=noc.dat)
# m2 <-glm(noc~ndvi+Age2+hii, data=noc.dat)
# m3 <-glm(noc~Age2+hii + hii*Age2, data=noc.dat)
# m4 <-glm(noc~ndvi+Age2+hii + hii*Age2 , data=noc.dat)
# 
# 
# 
# model.sel(m1,m2,m3,m4,  rank="AICc")
# 
# pred.noc.r1<-expand.grid(hii=as.numeric(c(0:40)),
#                          Age2=3:25,
#                          ndvi=mean(noc.dat$ndvi, na.rm=TRUE),
#                          Sex="F")
# 
# 
# 
# 
# pred.noc.r1$nocturnality <-predict(m1,newdata=
#                                  pred.noc.r1, type="response")
# 
# 
# ggplot(pred.noc.r1,
#   aes(x=hii, y=Age2, z=nocturnality, fill=nocturnality))+
#   ylab("Age")+
#   xlab("Human Influence Index")+
#   geom_tile()+
#   geom_contour()


##noc raw data
noc.dat%>%
  mutate(Age2=round(Age2,0))%>%
  group_by(BearID,Age2)%>%
  summarise(nocturnality=mean(noc),
            hii=mean(hii))%>%
    filter(hii<40)%>%
  ggplot(aes(x=Age2,y=hii, color=nocturnality))+
  geom_jitter(alpha=0.7, size=4)+
  scale_colour_viridis_c(option="C")+
    xlab("Age")+
  ylab("Human Influence Index")+
  theme_bw()
```

![](coexistence_sunrise_files/figure-markdown_github/additional%20analyses%20following%20revision%201:%20nocturnality%20and%20movement%20rates%20across%20space%20and%20hii-1.png)

``` r
ggsave(here::here("plots", "rev1", "nocturnal_raw.png"), width=6, height=5)

##movement rate raw data
noc.dat%>%
  mutate(Age2=round(Age2,0))%>%
  group_by(BearID,Age2)%>%
  summarise(day=mean(Day),
            night=mean(Night),
            hii=mean(hii))%>%
    filter(hii<40)%>%
  ungroup()%>%
  select(Age2,day, night, hii)%>%
  gather(time, speed,-Age2,-hii)%>%
  ggplot(aes(color=Age2,x=hii, y=speed))+
  geom_jitter(alpha=0.5, size=4)+
  scale_colour_viridis_c(option="C")+
  ylab("Speed (m/hr)")+
  xlab("Human Influence Index")+
facet_grid(.~time)+
  theme_bw()
```

![](coexistence_sunrise_files/figure-markdown_github/additional%20analyses%20following%20revision%201:%20nocturnality%20and%20movement%20rates%20across%20space%20and%20hii-2.png)

``` r
ggsave(here::here("plots", "rev1", "movement_raw.png"), width=9, height=5)

##habitat use raw data
reloc%>%
  filter(hii<=40)%>%
  ggplot(aes(x=hii, y=Age))+
  geom_point(alpha=0.01, size=1)+
  ylab("Age")+
  xlab("Human Influence Index")+
  theme_bw()
```

![](coexistence_sunrise_files/figure-markdown_github/additional%20analyses%20following%20revision%201:%20nocturnality%20and%20movement%20rates%20across%20space%20and%20hii-3.png)

``` r
ggsave(here::here("plots", "rev1", "habuse_raw.png"), width=6, height=5)
```

plot known dispersers
---------------------

``` r
#load bear data

#mort data
mort.disp <-mort%>%
  filter(BearID%in% c("Tina_FH", "Matt_L_EV", "Donnie_L_EV", "Sid_L_EV"))%>%
  st_as_sf(coords=c("X","Y"),
                   crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")%>%
    select(BearID)%>%
  rbind(
  data.frame(BearID=c("Tina_FH", "Matt_L_EV"),
                   X=c(-113.987,-115.465252),
                   Y=c(49.425,49.516892))%>%
  st_as_sf(coords=c("X","Y"),
                   crs = 4326)%>%
    st_transform("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))%>%
  st_transform(crs = 4326)%>%
  cbind(st_coordinates(.))


###known far dispersers  
disp <-st_as_sf(reloc,
                   crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")%>%
  filter(BearID%in% c("Tina_FH", "Matt_L_EV", "Donnie_L_EV", "Sid_L_EV"))%>%
  st_transform(crs = 4326)%>%
  mutate(Age2=paste0(Age,".",m)%>%as.numeric(),
         type="dispersed",
         Date=ymd(Date))%>%
  select(BearID, Date, Age, type)

##Donnie AB data
disp <- disp%>%
  rbind(
    st_read(here::here("data","spatial","dispersbears","Bear163", "locs_2016.shp"))%>%
  mutate(BearID="Donnie_L_EV",
         Date=ymd_hms(DateObs),
         Age=year(Date)-2014,
         type="natal")%>%
    st_transform(crs = 4326)%>%
    select(BearID, Date, Age, type))%>%
  rbind(
  st_read(here::here("data","spatial","dispersbears","Bear163", "locs_2017.shp"))%>%
  mutate(BearID="Donnie_L_EV",
         Date=ymd_hms(DateObs),
         Age=year(Date)-2014,
         type="natal")%>%
    st_transform(crs = 4326)%>%
    select(BearID, Date, Age, type))
```

    ## Reading layer `locs_2016' from data source `/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/data/spatial/dispersbears/Bear163/locs_2016.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 252 features and 48 fields
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: -115.1735 ymin: 50.59842 xmax: -115.0693 ymax: 50.92693
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## Reading layer `locs_2017' from data source `/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/data/spatial/dispersbears/Bear163/locs_2017.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 62 features and 62 fields
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: -115.1842 ymin: 50.62482 xmax: -115.109 ymax: 50.71597
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs

``` r
#Add Sid AB data
disp <- disp%>%
  rbind(
    st_read(here::here("data","spatial","dispersbears","Bear162", "Bear162.shp"))%>%
  mutate(BearID="Sid_L_EV",
         Date=ymd(LMT_Date),
         Age=year(Date)-2012,
         type="natal")%>%
    st_transform(crs = 4326)%>%
    select(BearID, Date, Age, type))
```

    ## Reading layer `Bear162' from data source `/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/data/spatial/dispersbears/Bear162/Bear162.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 896 features and 32 fields
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: -179769300000000064855122377195140780868567068050500468753735806562255881568525585520900552438967460767694165825208694872191919032324418090261363859732219097655018930368758270583513610708936040396203771265380480757050279273579187930844689739574764971740729296678532397561030338063719570761854957992502563438592 ymin: -179769300000000064855122377195140780868567068050500468753735806562255881568525585520900552438967460767694165825208694872191919032324418090261363859732219097655018930368758270583513610708936040396203771265380480757050279273579187930844689739574764971740729296678532397561030338063719570761854957992502563438592 xmax: -114.8253 ymax: 50.78919
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs

``` r
#Add Matt USA data
disp <- disp%>%
  rbind(
    st_read(here::here("data","spatial","dispersbears","GB810_MattsMom_GPSVHFLocs", "GB810_MattsMom_GPSVHFLocs.shp"))%>%
  mutate(BearID="Matt_L_EV",
         Date=ymd(Date),
         Age=year(Date)-2014,
         type="natal")%>%
    st_transform(crs = 4326)%>%
    select(BearID, Date, Age, type))%>%
    rbind(
    st_read(here::here("data","spatial","dispersbears","GB18986_Matt_DNAHits", "GB18986_Matt_DNAHits.shp"))%>%
  mutate(BearID="Matt_L_EV",
         Date=ymd(Date),
         Age=year(Date)-2014,
         type="natal")%>%
    st_transform(crs = 4326)%>%
    select(BearID, Date, Age, type))
```

    ## Reading layer `GB810_MattsMom_GPSVHFLocs' from data source `/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/data/spatial/dispersbears/GB810_MattsMom_GPSVHFLocs/GB810_MattsMom_GPSVHFLocs.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 1408 features and 18 fields
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: -116.1643 ymin: 48.78457 xmax: -115.9272 ymax: 48.95825
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs
    ## Reading layer `GB18986_Matt_DNAHits' from data source `/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/data/spatial/dispersbears/GB18986_Matt_DNAHits/GB18986_Matt_DNAHits.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 16 features and 16 fields
    ## geometry type:  POINT
    ## dimension:      XY
    ## bbox:           xmin: -116.1262 ymin: 48.81362 xmax: -116.02 ymax: 48.935
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs

``` r
##add Age2
disp<-disp%>%
  cbind(st_coordinates(.))%>%
  mutate(Age2=paste0(Age,".", month(Date)-3)%>%as.numeric(),
         type=case_when((BearID%in%"Tina_FH" & Age2<5.5)~"natal",
                        (BearID%in%"Tina_FH" & Age2>=5.5)~"dispersed",
                        BearID%in%"Sid_L_EV" & Y<50.5~"dispersed",
                        TRUE~type),
         sex=case_when(BearID%in%"Tina_FH"~"Female",
                       TRUE~"Male"))%>%
  filter(Y>0)


##load other spatial data
border <- st_read(here::here("data", "spatial", "borders", "North_America.shp"))%>%
  mutate(lab=case_when(FID_canada%in%9~"British Columbia, Canada",
                         FID_canada%in%8~"Alberta, Canada",
                         FID_usa%in%2~"Montana, USA"))%>%
  drop_na(lab)%>%
    st_transform(crs = 4326)
```

    ## Reading layer `North_America' from data source `/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/data/spatial/borders/North_America.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 70 features and 2 fields
    ## geometry type:  MULTIPOLYGON
    ## dimension:      XY
    ## bbox:           xmin: -178.2176 ymin: 18.92179 xmax: -52.62783 ymax: 83.11506
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs

``` r
hii.koot2<- raster(here::here("data","spatial","formap","hii_n_amer"))%>%
  crop(extent(c(-117.5, -112, 46, 51.5)))%>%
  projectRaster(crs="+proj=longlat +datum=WGS84 +no_defs")%>%
  as.data.frame(xy = TRUE)%>%
  filter(!is.na(hii_n_amer))%>%
  mutate(hii_n_amer=case_when(hii_n_amer>40~40,
                              TRUE~hii_n_amer))



register_google("AIzaSyCOwGx2D77XOqRgGhKmcb5F4Kt_S61tCLI")


plot.mvmt <- function(bear.name, zoomed, map.coords, label.df, displ){
bear <- disp%>%filter(BearID%in% c(bear.name))
bear.line <- bear%>%arrange(Date)%>%summarize(do_union=FALSE) %>% st_cast("LINESTRING")
a <- qmap(map.coords, zoom = zoomed, color = "bw")+
    geom_raster(data = hii.koot2,
              aes(x = x, y = y, fill=hii_n_amer),  alpha=0.4)+
    geom_sf(data=border,inherit.aes = FALSE, color="black",size=1.4, fill=NA)+
  geom_sf(data=bear.line,inherit.aes = FALSE, color="white", alpha=0.6)+
  geom_point(data=bear%>%as_tibble(),aes(x=X, y=Y), color="white", size=1.2,  alpha=0.4)+
  geom_point(data=mort.disp%>%filter(BearID%in% c(bear.name))%>%as_tibble(),aes(x=X, y=Y), color="black", size=3)+
  geom_text(data=label.df,
  aes(x, y, label = text), color="white")+
  labs(title=paste0("Brown Bear ", bear$sex[1], " - ", displ, " km shift"))+
  scale_fill_viridis_c()+
  theme(legend.position = "none")+
    ggrepel::geom_label_repel(
    data = bear%>%filter(type=="natal")%>%as_tibble()%>%select(X,Y)%>%summarise_all(mean)%>%mutate(lab="Natal Range")%>%
      rbind(mort.disp%>%filter(BearID%in% c(bear.name))%>%as_tibble()%>%select(X,Y)%>%mutate(lab="Mortality")),
    aes(label = lab, x=X, y=Y),
    segment.size  = 0.4,
    nudge_x=0.4,
    size=3
  )+
    annotation_scale(location = "br", width_hint = 0.25, text_cex=1, text_col="white")

return(a)
}



zoomed=9
a <- plot.mvmt(bear.name="Donnie_L_EV",
          map.coords=c(-114.9391, 50.31506),
          label.df=data.frame(
                               x = c(-114.4, -115.5),
                               y = c(50.5, 50),
                               text = c("AB, Can.", "BC, Can.")),
          displ=77
)

zoomed=9

zoomed=9
b <- plot.mvmt(bear.name="Sid_L_EV",
          map.coords=c(-114.9391, 50.25506),
          label.df=data.frame(
                               x = c(-114.4, -115.5),
                               y = c(50.5, 50),
                               text = c("AB, Can.", "BC, Can.")),
          displ=84
)


zoomed=8
c <- plot.mvmt(bear.name="Matt_L_EV",
          map.coords=c(-115.5076, 49.62778),
          label.df=data.frame(
                               x = c(-114.4, -116.5, -114.5, -116.7),
                               y = c(50.5, 50, 48.8, 48.7),
                               text = c("AB, Can.", "BC, Can.", "MT, USA", "ID, USA")),
          displ=95
)

zoomed=9
d <- plot.mvmt(bear.name="Tina_FH",
          map.coords=c(-114.1729, 49.20561),
          label.df=data.frame(
                               x = c(-114, -114.8, -114.4),
                               y = c(49.7, 49.2, 48.8),
                               text = c("AB, Can.", "BC, Can.", "MT, USA")),
          displ=53
)

hii.leg <- get_legend( ggplot()+geom_raster(data = hii.koot2%>%rename(HII=hii_n_amer),
              aes(x = x, y = y, fill=HII),  alpha=0.4)+
              scale_fill_viridis_c()+
                        theme(legend.position = "top"))%>%as_ggplot()


ggarrange(ggarrange(a,b,c,d, labels="AUTO"),ggarrange(hii.leg), nrow=2, heights=c(1,0.1))
```

![](coexistence_sunrise_files/figure-markdown_github/plot%20known%20dispersers-1.png)

``` r
ggsave(here::here("plots", "dispersal","Selected.png"), width=8, height=8, units="in")


###all known FH cap and dead locs
fh.disp <- read_csv(here::here("data", "spatial","dispersbears", "CAPTURE TO DEATH LOCATIONS APRIL 2020.csv"))%>%
  st_as_sf(coords=c("LONG","LAT"),
                   crs = 4326)

fh.disp.line <- fh.disp%>%group_by(NAME, SEX)%>%summarize(do_union=FALSE) %>% st_cast("LINESTRING")%>%ungroup()%>%mutate(l=st_length(.)%>%as.numeric())

hist <- ggplot(fh.disp.line%>%as_tibble())+
  geom_histogram(aes(x=l/1000))+
   geom_vline(data = fh.disp.line%>%as_tibble()%>%group_by(SEX)%>%summarise(l=mean(l/1000)), aes(xintercept=l),linetype="dashed")+
  theme_bw()+
  ylab("Count")+
  xlab("Dispersal distance (km)")+
  facet_wrap(vars(SEX), scales="free")+
  theme( axis.title = element_text(face = "bold",size = rel(1.3)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(size = rel(1.1)), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line())


qmap(c(-114.1729, 48.6), zoom = 8, color = "bw")+
    geom_raster(data = hii.koot2%>%rename(HII=hii_n_amer),
              aes(x = x, y = y, fill=HII),  alpha=0.4)+
    geom_sf(data=border,inherit.aes = FALSE, color="black",size=1.4, fill=NA)+
  geom_sf(data=fh.disp.line,inherit.aes = FALSE)+
  geom_sf(data=fh.disp,inherit.aes = FALSE,aes(color=TYPE), size=1.2)+
  scale_fill_viridis_c()+
    geom_text(data=data.frame(x = c(-112.9, -115.6, -113),
                               y = c(49.2, 49.2, 48.7),
                               text = c("AB, Can.", "BC, Can.", "MT, USA")),
  aes(x, y, label = text), color="white")+
      annotation_scale(location = "tr", width_hint = 0.25, text_cex=1, text_col="white")+
  annotation_custom(ggplotGrob(hist), xmin =-114.6, xmax = -112.4, ymin = 47.4, ymax = 48.3)
```

![](coexistence_sunrise_files/figure-markdown_github/plot%20known%20dispersers-2.png)

``` r
ggsave(here::here("plots", "dispersal","FH_all.png"), width=6, height=5, units="in")
```

dispersal distance (km) summaries
---------------------------------

``` r
kable(fh.disp.line%>%
  as_tibble()%>%
  group_by(SEX)%>%
  summarise(n=n(),mean=mean(l/1000), max=max(l/1000)))
```

| SEX |    n|      mean|        max|
|:----|----:|---------:|----------:|
| F   |   10|  12.08326|   52.72722|
| M   |   21|  31.52357|  159.93389|

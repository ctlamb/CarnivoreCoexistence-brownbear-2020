Source-Sink/Immigration Analysis
================
Clayton T. Lamb, Liber Ero Postdoctoral Fellow, University of British Columbia
09 May, 2020

Load Data
---------

``` r
set.seed(2019)
##Load Packages
##Load Packages
lapply(c("gtsummary", "gridExtra", "grid", "sjPlot", "huxtable", "demogR","ggpubr","ggspatial","ggmap", "rgdal","raster","sf","velox","here","survival","lubridate","MuMIn","lme4",
         "ggeffects", "secr", "knitr","forcats", "ggridges","fuzzyjoin", "tidyverse"), require, character.only = TRUE)

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
        legend.title = element_blank(),
        legend.key.size= unit(0.6, "cm"),
        legend.spacing = unit(0.1, "cm"),
        legend.text = element_text(size = rel(1)),
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
repro <- read.csv(here::here("data", "repro_CLedits.csv"))

####Parse repro to annual measures
repro <- repro%>%
  mutate(Year=year(Date))%>%
  group_by(BearID,AgeClass, StudyArea, Year)%>%
  summarise(Count=max(Count, na.rm=TRUE),
            litter=min(Litter, na.rm=TRUE),
            S=mean(S, na.rm=TRUE),
            ndead=mean(ndead, na.rm=TRUE),
            Yrl_dispersUnknown=mean(as.numeric(Yrl_dispersUnknown),na.rm=TRUE))
```

Prep Survival Data
------------------

``` r
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
  mutate(Cause_CL=case_when(Cause_CL%in%"dolp" ~ "Attractant/conflict",
                            TRUE~as.character(Cause_CL)))%>%
  filter(CollarOn_CL%in%"Y")


reloc <- reloc%>%
  select("BearID","Date", "Time", "X", "Y")%>%
  mutate(mort=0,
         Cause_CL=NA)%>%
  rbind(mort%>%mutate(Time=NA)%>%select("BearID","Date", "Time", "X", "Y", "Cause_CL")%>%mutate(mort=1))


##add sex,age, ageclass,and study area into reloc
cap <- cap%>%group_by(BearID)%>%mutate(b_y=year(CollarStartDate)-Age)%>%
  summarize(birthyear=min(b_y, na.rm=TRUE))%>%
  left_join(cap,by="BearID")

reloc <- reloc%>% left_join(cap%>%select(BearID,Sex,birthyear, StudyArea)%>%distinct(BearID,Sex,birthyear, StudyArea),by="BearID")%>%
  mutate(Age=year(Date)-birthyear)%>%
  mutate(ageclass=case_when(
    Age %in% c(0) ~ "Age (0)",
    Age %in% c(1) ~ "Age (1)",
    Age %in% c(2:4) ~ "Age (2:4)",
    Age %in% c(5:9) ~ "Age (5:9)",
    Age %in%c(10:35) ~ "Age (>=10)"
  ))

##add month day year
reloc <- reloc%>%
  mutate(yr=year(Date),
         m=month(Date),
         d=day(Date))

##only keep active months
reloc <- reloc %>% filter(m%in%c(4:11))

##drop bears with unknown ages
reloc <- reloc%>%
  filter(!Age<0)


##################################################################
##Add Spatial Data
##################################################################
##make bear data spatial
reloc <- st_as_sf(reloc, 
                  coords = c("X","Y"),
                  crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

########
## HUMAN INFLUENCE INDEX
########

hii<- raster(here::here("data", "spatial", "hii.tif"))

##make velox (faster)
hiiv <- velox(hii)


##extract hii to reloc
reloc <- reloc%>%mutate(hii=hiiv$extract_points(sp = reloc))

########
## VEGETATION
########

##stright summer NDVI
r <- raster(here::here("data", "spatial", "NDVIspring.tif"))

##make velox (faster)
dndviv <- velox(r)


##extract to reloc
reloc <- reloc%>%mutate(ndvi=dndviv$extract_points(sp = reloc))%>%
  mutate(ndvi=case_when(is.na(ndvi)~0,
                         TRUE~ndvi))



##################################################################
##Smooth over monthly periods
##################################################################

##pull out exact mort locs
kills <- reloc%>%filter(mort%in%1)

##summarize live locs by month
reloc.m <- reloc%>%  
  filter(mort%in%0)%>%
  group_by(BearID,yr,m,Sex,Age, ageclass,StudyArea, mort)%>%
summarise_at(vars(hii, ndvi), mean, na.rm = TRUE)%>%
  as_tibble()

##assign 30 days to each, as they are each a month
reloc.m$time <- 30

##length of survival on last month of life= day of month
kills$time<-kills$d
##match columns
kills <- kills%>%select(colnames(reloc.m))

##bind and join last month with live and dead if needed
reloc.m <- rbind(as_tibble(reloc.m),as_tibble(kills))%>%
  group_by(BearID,Sex,Age, ageclass, StudyArea,yr,m)%>%
  summarise_at(vars(hii, ndvi), mean, na.rm = TRUE)%>%
  left_join(rbind(as_tibble(reloc.m),as_tibble(kills))%>%group_by(BearID,Sex,Age, ageclass, StudyArea,yr,m)%>%summarise(mort=max(mort, na.rm=TRUE),
            time=min(time, na.rm=TRUE)))
```

Summary Stats
-------------

``` r
###############
###Summary Tables
###############
inds <- unique(reloc$BearID)
move.sum <- data.frame()
for(i in 1:length(inds)){
a <- reloc%>%filter(BearID%in%inds[i])%>%
  select(BearID, StudyArea, Date, Age, yr, hii,mort)%>%
  as_tibble()%>%
  mutate(Date=ymd(Date))%>%
  select(-geometry)%>%
  fuzzy_left_join(
   cap%>%filter(BearID%in%inds[i])%>%select(BearID,CollarStartDate,CollarEndDate,CollarType)%>%mutate(BearID=as.character(BearID)),
  by = c(
    "BearID" = "BearID",
    "Date" = "CollarStartDate",
    "Date" = "CollarEndDate"
  ),
  match_fun = list(`==`, `>=`, `<=`)
) 

move.sum <- rbind(move.sum, a)
print(i)
}


gps.ind <- move.sum%>%gather(key, value, CollarType) %>%
  distinct(StudyArea, BearID.x,value,key)%>%
  group_by(StudyArea, key, value) %>%
  tally %>% 
  spread(value, n, fill = 0)%>%
  mutate(GPSindividuals_perc=round((GPS/(VHF+GPS))*100),0)

gps.loc <-move.sum%>%gather(key, value, CollarType) %>%
  group_by(StudyArea, key, value) %>%
  tally %>% 
  spread(value, n, fill = 0)%>%
  mutate(GPSrelocations_perc=round((GPS/(VHF+GPS))*100),0)


  
join.gps <- gps.ind%>%select(StudyArea, GPSindividuals_perc)%>%
  left_join(gps.loc%>%select(StudyArea, GPSrelocations_perc))%>%
  ungroup()%>%
  select(-key)

sum <- reloc%>%
  as.tibble()%>%
  ungroup()%>%
  group_by(StudyArea)%>%
  summarize(animals=n_distinct(BearID),
            #age_mean=mean(Age)%>%round(0),
            ages=paste0(min(Age, na.rm=TRUE)%>%round(0),"-",max(Age, na.rm=TRUE)%>%round(0)),
            #years_n=(max(yr)-min(yr))+1,
            years=paste0(min(yr, na.rm=TRUE),"-",max(yr, na.rm=TRUE)),
            #hii_mean=mean(hii, na.rm=TRUE)%>%round(1),
            hii=paste0(min(hii, na.rm=TRUE)%>%round(1),"-",max(hii, na.rm=TRUE)%>%round(1)),
            morts=sum(mort),
            relocations=(n()-morts)
            )%>%
  left_join(join.gps%>%select(StudyArea, GPSindividuals_perc))%>%
  left_join(repro%>%group_by(StudyArea)%>%summarise(repro=n()))%>%
  ungroup()%>%
  select(StudyArea, animals, ages,  years, hii, morts, repro, relocations, GPSindividuals_perc)
  

 sum <- add_row(sum, StudyArea="SUMMARY",
            animals = sum(sum$animals),
            ages=paste0(min(reloc%>%as.tibble()%>%ungroup()%>%pull(Age), na.rm=TRUE)%>%round(0),
                        "-",
                        max(reloc%>%as.tibble()%>%ungroup()%>%pull(Age), na.rm=TRUE)%>%round(0)),
            years=paste0(min(reloc%>%as.tibble()%>%ungroup()%>%pull(yr), na.rm=TRUE)%>%round(0),
                         "-",
                         max(reloc%>%as.tibble()%>%ungroup()%>%pull(yr), na.rm=TRUE)%>%round(0)),
            hii=paste0(min(reloc%>%as.tibble()%>%ungroup()%>%pull(hii), na.rm=TRUE)%>%round(0),
                         "-",
                         max(reloc%>%as.tibble()%>%ungroup()%>%pull(hii), na.rm=TRUE)%>%round(0)),
            morts=sum(sum$morts),
            repro=sum(sum$repro),
            relocations=sum(sum$relocations),
          GPSindividuals_perc=(sum(gps.ind$GPS)/sum(c(gps.ind$GPS, gps.ind$VHF))*100)%>%round(0))%>%
   as_hux(add_colnames = TRUE,
         scientific=FALSE)%>%
  theme_article()%>%
  set_width(0.9)

huxtable::number_format(sum)[, c(8)] <- list(
  function(x)
    prettyNum(x,
              scientific = FALSE))

quick_docx(sum,file=here::here("tables", "summary.docx"), open=FALSE)
write_csv(sum,here::here("tables", "summary.csv"))

###ADDITIONAL STATS
##how many bear-years of monitoring?
sum(cap$time, na.rm=TRUE)/365  ##808

##average monitoring time
cap%>%group_by(BearID)%>%summarise(sum=sum(time,na.rm=TRUE)/365)%>%pull(sum)%>%mean(na.rm=TRUE) ##1.7
cap%>%group_by(BearID)%>%summarise(sum=sum(time,na.rm=TRUE)/365)%>%pull(sum)%>%max(na.rm=TRUE) ##21.02

##average project time
mean((sum[2:13,]$yr_max%>%as.numeric()-sum[2:13,]$yr_min%>%as.numeric())+1)

ggplot(reloc.m%>%drop_na(StudyArea), aes(x = hii, y = StudyArea)) + 
  geom_density_ridges( jittered_points = TRUE,
                       position = position_points_jitter(width = 0.05, height = 0),
                       point_shape = '|', point_size = 1.5, point_alpha = 1, alpha = 0.7)+
  xlab("Human Influence Index")+
  theme_ridges()
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/Summary%20Stats-1.png)

``` r
ggsave(here::here("plots","hab_use_byproj.png"), width=10, height=6.5, units="in")
```

Fit Survival Models
-------------------

``` r
##################################################################
##Fit Survival Models
##################################################################
reloc.m.ad <- reloc.m %>% filter(!ageclass %in% c("Age (0)", "Age (1)"))%>%drop_na(ageclass)

fit1 <- coxph(Surv(time, mort)~ Age+ 
               cluster(BearID), data = reloc.m.ad)
fit2 <- coxph(Surv(time, mort)~ Sex + Age+ 
               cluster(BearID), data = reloc.m.ad)
fit3 <- coxph(Surv(time, mort)~ Sex + Age + hii +
               cluster(BearID), data = reloc.m.ad)
fit4 <- coxph(Surv(time, mort)~ (Sex+Age*hii) +
               cluster(BearID), data = reloc.m.ad)
fit5 <- coxph(Surv(time, mort)~ (Sex*hii + Age) +
                cluster(BearID), data = reloc.m.ad)
fit8 <- coxph(Surv(time, mort)~ (Sex*ageclass*hii + ndvi) +
                cluster(BearID), data = reloc.m.ad)
fit9 <- coxph(Surv(time, mort)~ (Sex+Age*hii + ndvi) +
                cluster(BearID), data = reloc.m.ad)
fit10 <- coxph(Surv(time, mort)~ (Sex+Age*hii + I(Age^2) + ndvi) +
                cluster(BearID), data = reloc.m.ad)
fit11 <- coxph(Surv(time, mort)~ (Sex+Age*hii + I(Age^2)  + I(Age^3) + ndvi) +
                cluster(BearID), data = reloc.m.ad)
fit11.5 <- coxph(Surv(time, mort)~ (Sex+Age*hii + I(Age^2)  + I(Age^3) + I(Age^4) + ndvi) +
                cluster(BearID), data = reloc.m.ad)
fit12 <- coxph(Surv(time, mort)~ (Sex*Age*hii + I(Age^2) + ndvi) +
                cluster(BearID), data = reloc.m.ad)
fit13 <- coxph(Surv(time, mort)~ (Sex*Age*hii + I(Age^2)  + I(Age^3) + ndvi) +
                 cluster(BearID), data = reloc.m.ad)
```

Export table
------------

``` r
##prep model selection table
mod.sel <- model.sel(fit1,fit2,fit3,fit4,fit5,fit8, fit9, fit10, fit11, fit11.5, fit12, fit13,  rank="AICc")

# replace Model name with formulas little tricky so be careful 
for(i in 1:nrow(mod.sel)) mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] 

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

##model coefficient table
tbl_regression(fit10,estimate_fun = function(x) style_ratio(x, digits = 3))

quick_docx(mod.sel,file=here::here("tables", "Survival_AIC.docx"), open=FALSE)
```

``` r
kable(mod.sel[-1,])
```

|     | Model                                                                                    | df  | AICc    | dAICc | weight |
|-----|:-----------------------------------------------------------------------------------------|:----|:--------|:------|:-------|
| 2   | Sex + Age + hii + I(Age^2) + ndvi + Age:hii                                              | 6   | 1125.82 | 0     | 0.51   |
| 3   | Sex + Age + hii + I(Age^2) + I(Age^3) + ndvi + Age:hii                                   | 7   | 1127.45 | 1.63  | 0.23   |
| 4   | Sex + Age + hii + I(Age^2) + I(Age^3) + I(Age^4) + ndvi + Age:hii                        | 8   | 1127.96 | 2.14  | 0.18   |
| 5   | Sex + Age + hii + I(Age^2) + ndvi + Sex:Age + Sex:hii + Age:hii + Sex:Age:hii            | 9   | 1130.42 | 4.61  | 0.05   |
| 6   | Sex + Age + hii + I(Age^2) + I(Age^3) + ndvi + Sex:Age + Sex:hii + Age:hii + Sex:Age:hii | 10  | 1131.39 | 5.57  | 0.03   |
| 7   | Sex + Age + hii + ndvi + Age:hii                                                         | 5   | 1138.61 | 12.8  | 0      |
| 8   | Sex + Age + hii + Age:hii                                                                | 4   | 1141.41 | 15.59 | 0      |
| 9   | Sex + Age + hii                                                                          | 3   | 1142.6  | 16.78 | 0      |
| 10  | Sex + hii + Age + Sex:hii                                                                | 4   | 1144.56 | 18.74 | 0      |
| 11  | Sex + ageclass + hii + ndvi + Sex:ageclass + Sex:hii + ageclass:hii + Sex:ageclass:hii   | 12  | 1147.11 | 21.29 | 0      |
| 12  | Sex + Age                                                                                | 2   | 1159.93 | 34.11 | 0      |
| 13  | Age                                                                                      | 1   | 1164.21 | 38.39 | 0      |

Boot Survival Models
--------------------

``` r
##################################################################
##Fit Survival Models
##################################################################


pred.boot <- data.frame()
for(i in 1:500){
##bootstrap
reloc.m.i <- reloc.m.ad%>%sample_frac(1, replace=TRUE)


fit10 <- coxph(Surv(time, mort)~ (Sex+Age*hii + I(Age^2) + ndvi) +
                cluster(BearID), data = reloc.m.i)
fit11 <- coxph(Surv(time, mort)~ (Sex+Age*hii + I(Age^2)  + I(Age^3) + ndvi) +
                cluster(BearID), data = reloc.m.i)
fit11.5 <- coxph(Surv(time, mort)~ (Sex+Age*hii + I(Age^2)  + I(Age^3) + I(Age^4) + ndvi) +
                cluster(BearID), data = reloc.m.i)


##################################################################
##Plot Effects, weighted by top models
##################################################################

##make prediction data
surv.pred <- expand.grid(Sex=c("F"), Age=c(2:30), hii=c(0:40), time=30, mort=0, ndvi=seq(0.5,0.7, by=0.1))
  

###predict each of 2 top mods


###predict each of 3 top mods
surv.pred$pred10 <- exp(-predict(fit10, newdata=surv.pred, type="expected"))

surv.pred$pred11 <- exp(-predict(fit11, newdata=surv.pred, type="expected"))

surv.pred$pred115 <- exp(-predict(fit11.5, newdata=surv.pred, type="expected"))


##average each model prediction by model weight
wts <- Weights(AIC(fit10, fit11, fit11.5))
surv.pred$pred <- ((surv.pred$pred10*wts[[1]])+(surv.pred$pred11*wts[[2]])+(surv.pred$pred115*wts[[3]]))^6

surv.pred$pred <- surv.pred$pred * ((1-0.003)^6)



surv.pred$iter<-i


pred.boot <- rbind(pred.boot,surv.pred%>%select(Sex,Age,hii,ndvi,time,mort,pred,iter))
}
```

CUB + YRL SURVIVAL
------------------

``` r
reloc.cub <- reloc%>%
  filter(Sex=="F")%>%
  group_by(BearID,yr,Sex,Age)%>%
  summarise_at(vars(hii, ndvi), mean, na.rm = TRUE)%>%
  rename(Year=yr)%>%
  as_tibble()

##Dead F's with COY, cubs dead too.
mort.cubs <- mort%>%
  mutate(Year=year(Date), mort=1)%>%
  select(BearID, Year, mort)


##bind repro data together
repro.surv <- repro%>%
  ungroup()%>%
  filter(AgeClass%in%c("COY", "1YO"))%>%
  mutate(AgeClass=droplevels(AgeClass))%>%
  left_join(mort.cubs, by=c("BearID","Year"))%>%
  mutate(S=case_when((mort%in%1 & AgeClass%in%c("COY")) ~ 0,
                     TRUE~as.numeric(S)),
         ndead=case_when((mort%in%1 & AgeClass%in%c("COY"))  ~ Count,
                     TRUE~as.integer(ndead)))%>%
  drop_na(S)
  

##add in spatial data from mom
repro.surv <- repro.surv%>%
  left_join(reloc.cub, by=c("BearID","Year"))

##cleanup
repro.surv <- repro.surv%>%
  drop_na(hii)%>%
  filter(Age>=0)


##expand by individual cub

repro.expand <- data.frame()
for(i in 1:nrow(repro.surv)){
  a <- repro.surv[i,]
  b <- bind_rows(replicate(n=a[1,"Count"][[1]], a, simplify = FALSE))
  b$fate <- c(rep(1,a[1,"ndead"]), rep(0,(a[1,"Count"]-a[1,"ndead"])))
  repro.expand <- rbind(repro.expand,b)
}


##add time
repro.expand$time<-365
```

Fit Cub + Yrl Survival Models
-----------------------------

``` r
##################################################################
##Fit Survival Models
##################################################################
repro.expand <- repro.expand%>%filter(AgeClass%in%c("COY", "1YO"))

fit1 <- coxph(Surv(time, fate)~ hii +
               cluster(BearID), data = repro.expand)
fit2 <- coxph(Surv(time, fate)~ ndvi +
               cluster(BearID), data = repro.expand)
fit3 <- coxph(Surv(time, fate)~ hii + ndvi +
               cluster(BearID), data = repro.expand)
fit4 <- coxph(Surv(time, fate)~ hii + AgeClass +
               cluster(BearID), data = repro.expand)
fit5 <- coxph(Surv(time, fate)~ ndvi + AgeClass +
               cluster(BearID), data = repro.expand)
fit6 <- coxph(Surv(time, fate)~ hii + ndvi + AgeClass +
               cluster(BearID), data = repro.expand)
fit7 <- coxph(Surv(time, fate)~ hii + (ndvi*AgeClass) +
               cluster(BearID), data = repro.expand)
```

Export table
------------

``` r
##prep model selection table
mod.sel <- model.sel(fit1, fit2, fit3, fit4, fit5, fit6, fit7,  rank="AICc")

# replace Model name with formulas little tricky so be careful 
for(i in 1:nrow(mod.sel)) mod.sel$Model[i]<- as.character(formula(paste(rownames(mod.sel)[i])))[3] 

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


quick_docx(mod.sel,file=here::here("tables", "CubSurvival_AIC.docx"), open=FALSE)

##model coefficient table
tbl_regression(fit5,estimate_fun = function(x) style_ratio(x, digits = 3))
```

``` r
kable(mod.sel[-1,])
```

|     | Model                                 | df  | AICc   | dAICc | weight |
|-----|:--------------------------------------|:----|:-------|:------|:-------|
| 2   | ndvi                                  | 1   | 946.29 | 0     | 0.28   |
| 3   | hii                                   | 1   | 946.44 | 0.15  | 0.26   |
| 4   | ndvi + AgeClass                       | 2   | 947.52 | 1.24  | 0.15   |
| 5   | hii + AgeClass                        | 2   | 947.62 | 1.33  | 0.14   |
| 6   | hii + ndvi                            | 2   | 948.36 | 2.07  | 0.1    |
| 7   | hii + ndvi + AgeClass                 | 3   | 949.66 | 3.37  | 0.05   |
| 8   | hii + ndvi + AgeClass + ndvi:AgeClass | 4   | 951.22 | 4.93  | 0.02   |

Boot Cub Survival Models
------------------------

``` r
##################################################################
##Fit Survival Models
##################################################################


for(i in 1: 500){
#bootstrap
repro.expand.i <- repro.expand%>%sample_frac(1, replace=TRUE)


fit4 <- coxph(Surv(time, fate)~ hii + AgeClass +
               cluster(BearID), data = repro.expand.i)
fit5 <- coxph(Surv(time, fate)~ ndvi + AgeClass +
               cluster(BearID), data = repro.expand.i)
fit6 <- coxph(Surv(time, fate)~ hii + ndvi + AgeClass +
               cluster(BearID), data = repro.expand.i)


##################################################################
##Plot Effects, weighted by top models
##################################################################

##make prediction data
surv.pred <- expand.grid(Sex="pooled",AgeClass=c("COY", "1YO"), hii=c(0:40), time=365, fate=0, ndvi=seq(0.5,0.7, by=0.1))%>%
  mutate(Age=case_when(AgeClass==c("COY")~0,AgeClass==c("1YO")~1))
  

###predict each of 3 top mods
surv.pred$pred4 <- exp(-predict(fit4, newdata=surv.pred, type="expected"))

surv.pred$pred5 <- exp(-predict(fit5, newdata=surv.pred, type="expected"))

surv.pred$pred6 <- exp(-predict(fit6, newdata=surv.pred, type="expected"))



##average each model prediction by model weight
wts <- Weights(AIC(fit4, fit5, fit6))
surv.pred$pred <- (surv.pred$pred4*wts[[1]])+(surv.pred$pred5*wts[[2]])+(surv.pred$pred6*wts[[3]])


surv.pred$iter<-i

pred.boot <- rbind(pred.boot,surv.pred%>%select(Sex,Age,hii,ndvi,time,fate,pred,iter)%>%rename("mort"="fate"))

}
```

### scale YRL survival to published values

``` r
pred.boot%>%as_tibble()%>%ungroup()%>%filter(Age%in%1 & hii==9)%>%summarise(pred=mean(pred))

#(0.86-0.632)/0.632

pred.boot<- pred.boot%>%
  mutate(pred=case_when(Age%in%1~pred+(pred*0.36),TRUE~pred))%>%
  mutate(pred=case_when(Age%in%1 & pred>1~1,TRUE~pred))

pred.boot%>%as_tibble()%>%ungroup()%>%filter(Age%in%0 & hii==9)%>%summarise(pred=mean(pred))
#(0.755-0.689)/0.689

pred.boot<- pred.boot%>%
  mutate(pred=case_when(Age%in%0~pred+(pred*0.1),TRUE~pred))%>%
  mutate(pred=case_when(Age%in%0 & pred>1~1,TRUE~pred))
```

PLOT
----

``` r
ggplot()+
  geom_path(data=pred.boot%>%filter(Sex%in%c("F", "pooled"))%>%group_by(Age,hii,iter)%>%
              summarise(pred=mean(pred))%>%
              filter(hii%in%c(0,10,20,30,40)),
            aes(x=Age,y=pred,group=iter),
            size=0.1, alpha=0.05, col="red") +
    geom_line(data=pred.boot%>%filter(Sex%in%c("F", "pooled"))%>%group_by(Age,hii)%>%
              summarise(pred=mean(pred))%>%
              filter(hii%in%c(0,10,20,30,40)),
            aes(x=Age,y=pred),
            size=1, col="black") +
  ylab("Survival")+
  xlab("Age")+
  xlim(0,25)+
  ylim(0,1)+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_Publication()+
  facet_grid(hii~., labeller = label_both)
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/Plot%20Cub%20Survival%20Models-1.png)

``` r
ggsave(here::here("plots","Fsurvival.png"), width=5, height=9, units="in")
```

REPRODUCTION
------------

``` r
repro.rate <- repro%>%
  ungroup()%>%
  mutate(Count=case_when(!AgeClass%in%("COY")~as.integer(0),AgeClass%in%("COY")~Count))%>%
  left_join(reloc.cub, by=c("BearID","Year"))%>%
  drop_na(hii)%>%
  mutate(class=case_when(Age<9~"<9",Age>=9~"9+"))%>%
  filter(Age>=4)


m1 <- lm(Count~1, data=repro.rate)
m2 <- lm(Count~Age  , data=repro.rate)
m3 <- lm(Count~Age + I(Age^2), data=repro.rate)
m4 <- lm(Count~Age + I(Age^2)+ hii , data=repro.rate)
m5 <- lm(Count~Age + I(Age^2) + ndvi, data=repro.rate)
m6 <- lm(Count~Age + I(Age^2) + ndvi + hii, data=repro.rate)
m7 <- lm(Count~class, data=repro.rate)
```

Export table
------------

``` r
##prep model selection table
mod.sel <- model.sel(m1,m2,m3,m4,m5,m6,m7,  rank="AICc")

# replace Model name with formulas little tricky so be careful 
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

quick_docx(mod.sel,file=here::here("tables", "Repro_AIC.docx"), open=FALSE)

##model coefficient table
tbl_regression(m3,estimate_fun = function(x) style_ratio(x, digits = 3))
```

``` r
kable(mod.sel[-1,])
```

|     | Model                       | df  | AICc    | dAICc | weight |
|-----|:----------------------------|:----|:--------|:------|:-------|
| 2   | Age + I(Age^2)              | 4   | 1422.98 | 0     | 0.28   |
| 3   | Age + I(Age^2) + hii        | 5   | 1423.27 | 0.28  | 0.24   |
| 4   | Age + I(Age^2) + ndvi + hii | 6   | 1423.72 | 0.74  | 0.19   |
| 5   | Age + I(Age^2) + ndvi       | 5   | 1423.82 | 0.84  | 0.18   |
| 6   | class                       | 3   | 1425.87 | 2.89  | 0.07   |
| 7   | 1                           | 2   | 1427.92 | 4.93  | 0.02   |
| 8   | Age                         | 3   | 1429.12 | 6.13  | 0.01   |

Boot Repro Models
-----------------

``` r
pred.boot.repro <- data.frame()
for(i in 1: 500){
#bootstrap
repro.rate.i <- repro.rate%>%ungroup()%>%sample_frac(1, replace=TRUE)

m3 <- lm(Count~Age + I(Age^2), data=repro.rate.i)
m4 <- lm(Count~Age + I(Age^2)+ hii , data=repro.rate.i)
m6 <- lm(Count~Age + I(Age^2) + ndvi + hii, data=repro.rate.i)

##################################################################
##Plot Effects, weighted by top models
##################################################################

##make prediction data
repro.pred <- expand.grid(Age=c(4:28), 
                          hii=c(0:40), 
                          ndvi=seq(0.5,0.7, by=0.1))
  

###predict each of 2 top mods
repro.pred$pred3 <- predict(m3, newdata=repro.pred, type = "response")

repro.pred$pred4 <- predict(m4, newdata=repro.pred, type = "response")

repro.pred$pred6 <- predict(m6, newdata=repro.pred, type = "response")


##average each model prediction by model weight
wts <- Weights(AIC(m3, m4, m6))
repro.pred$pred <- (repro.pred$pred3*wts[[1]])+(repro.pred$pred4*wts[[2]])+(repro.pred$pred6*wts[[3]])


repro.pred$iter<-i

pred.boot.repro <- rbind(pred.boot.repro,repro.pred)
}


##replace<0 w/ 0

pred.boot.repro <- pred.boot.repro%>%
  mutate(pred=case_when(pred<0~0,
                        TRUE~as.numeric(pred)))
```

PLOT
----

``` r
ggplot()+
  geom_path(data=pred.boot.repro%>%filter(hii==20 & ndvi==0.5),aes(x=Age,y=pred/2, group=iter),size=0.1, alpha=0.1, col="red") +
   geom_path(data=pred.boot.repro%>%filter(hii==20 & ndvi==0.5)%>%group_by(Age)%>%summarise(pred=mean(pred)),aes(x=Age,y=pred/2),size=1, col="black")+
  ylab("Reproduction Rate (F cubs)")+
  xlab("Age")+
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

![](Source_sink_boot_contsAge_files/figure-markdown_github/Plot%20Repro%20Models-1.png)

``` r
ggsave(here::here("plots","Frepro.png"), width=7, height=5, units="in")
```

Make lifetables
---------------

``` r
lambda.loop <- data.frame()
ndvi.seq=seq(0.5,0.7, by=0.1)
for(k in 1:length(ndvi.seq)){  ###NDVI
  a <-pred.boot%>%mutate(ndvi2=as.integer(ndvi*10))%>%filter(ndvi2==ndvi.seq[k]*10 & Sex%in% c("F", "pooled"))
  b <- pred.boot.repro%>%mutate(ndvi2=as.integer(ndvi*10))%>%filter(ndvi2==ndvi.seq[k]*10)

  for(j in 0:40){              ###HII
  
  a.j <-a%>%filter(hii==j)
  b.j <- b%>%filter(hii==j)
  
  
    for(i in 1:500){             ###BOOT
    a.i <-a.j%>%filter(iter==i)%>%
      arrange(Age)

  
    b.i <- b.j%>%filter(iter==i)
    
Px <- round(a.i%>%pull(pred),4)
#Px[3] <-(Px[2]+Px[3]) /2
Fx <- round(c(0,0, 0, 0, b.i%>%mutate(pred=pred/2)%>%pull(pred),0,0, 0),4)
L <- odiag(Px, -1)
L[1, ] <- Fx


lambda.loop <- rbind(lambda.loop,
                     data.frame(ndvi=ndvi.seq[k],
                                hii=j,
                                iter=i,
                                # cub=Px[1],
                                # one=Px[2],
                                # sa=mean(Px[3:5]),
                                # ad=mean(Px[6:31]),
                                # m5=Fx[6],
                                # m6_8=mean(Fx[7:9]),
                                # m9_12=mean(Fx[10:13]),
                                # m13_7=mean(Fx[14:18]),
                                # m17p=mean(Fx[19:31]),
                                l=eigen(L)$values[1]))
  }
}
}
```

PLOT
----

``` r
ggarrange(
ggplot()+
  geom_path(data=lambda.loop%>%
              group_by(hii,iter)%>%
              summarise(l=as.numeric(mean(l))),aes(x=as.numeric(hii),y=as.numeric(l), group=iter),
            size=0.1, alpha=0.08, col="red") +
  geom_line(data=lambda.loop%>%
              group_by(hii)%>%
              summarise(l=mean(l)),aes(x=as.numeric(hii),y=as.numeric(l)),
            size=1,  col="black")+
  #geom_point(data=data.frame(study="McLellan FH all years", hii=9.3, l=1.046), aes(x=hii, y=l,label=study))+
  #geom_text(data=data.frame(hii=9.4,l=1.05), aes(hii, l), label="McLellan FH all years", vjust=-1)+
  ylab("Population Growth")+
  xlab("Human Influence Index")+
  xlim(0,40)+
  geom_hline(yintercept = 1, linetype="dashed")+
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
                    strip.text = element_text(face="bold")),



ggplot()+
  geom_path(data=lambda.loop%>%
              group_by(hii,iter)%>%
              summarise(l=as.numeric(mean(l)))%>%
              mutate(im=case_when((1-l)>0 ~ 1-l,
                      (1-l)<=0 ~ 0)*100),aes(x=as.numeric(hii),y=as.numeric(im), group=iter),
            size=0.1, alpha=0.08, col="red") +
  geom_line(data=lambda.loop%>%
              group_by(hii)%>%
              summarise(l=as.numeric(mean(l)))%>%
              mutate(im=case_when((1-l)>0 ~ 1-l,
                      (1-l)<=0 ~ 0)*100),
            aes(x=as.numeric(hii),y=as.numeric(im)),
            size=1,  col="black")+
  #geom_point(data=data.frame(study="McLellan FH all years", hii=9.3, l=1.046), aes(x=hii, y=l,label=study))+
  #geom_text(data=data.frame(hii=9.4,l=1.05), aes(hii, l), label="McLellan FH all years", vjust=-1)+
  ylab("Immigration Required (%)")+
  xlab("Human Influence Index")+
  xlim(0,40)+
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
                    strip.text = element_text(face="bold")),
          ncol = 1, nrow = 2,
          labels="AUTO")
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/Plot%20Lambda-1.png)

``` r
ggsave(here::here("plots","Lambda_immi.png"), width=7, height=10, units="in")
```

Plot Together
-------------

``` r
s.plot <- ggplot(pred.boot%>%filter(Sex%in%c("F", "pooled"))%>%group_by(Age,hii)%>%
              summarise(pred=mean(pred))%>%rename(Survival=pred),
       aes(x=Age, y=hii, fill=Survival))+
  ylab("Human Influence Index")+
  xlab("")+
  geom_tile()+
  ggtitle("Survival")+
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme_Publication()+
    theme(legend.position = "right",
        legend.direction = "vertical")

r.plot <- ggplot(pred.boot.repro%>%
                   group_by(Age,hii)%>%
                   summarise(pred=mean(pred)/2)%>%
                   ungroup()%>%
                   rbind(expand.grid(Age=c(0:3,29:31), hii=c(0:40),pred=0)%>%as_tibble())%>%
                   rename(Reproduction=pred),
       aes(x=Age, y=hii, fill=Reproduction))+
  ylab("Human Influence Index")+
  geom_tile()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme_Publication()+
  ggtitle("Reproduction")+
  theme(legend.position = "right",
        legend.direction = "vertical")


F_VitalRates<- ggarrange(s.plot, r.plot,
          ncol=1, nrow=2)

ggsave(here::here("plots","F_VitalRates.png"), width=6, height=8, units="in")
```

Plot Together-error
-------------------

``` r
s.plot <- ggplot(pred.boot%>%filter(Sex%in%c("F", "pooled"))%>%group_by(Age,hii)%>%
              summarise(pred=sd(pred))%>%rename(Survival=pred),
       aes(x=Age, y=hii, fill=Survival))+
  ylab("Human Influence Index")+
  xlab("")+
  geom_tile()+
  labs(title="Survival Standard Error")+
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme_Publication()+
    theme(legend.position = "right",
        legend.direction = "vertical")

r.plot <- ggplot(pred.boot.repro%>%
                   group_by(Age,hii)%>%
                   summarise(pred=sd(pred/2))%>%
                   ungroup()%>%
                   rbind(expand.grid(Age=c(0:3,29:31), hii=c(0:40),pred=0)%>%as_tibble())%>%
                   rename(Reproduction=pred),
       aes(x=Age, y=hii, fill=Reproduction))+
  ylab("Human Influence Index")+
  geom_tile()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme_Publication()+
  labs(title="Reproduction Standard Error")+
  theme(legend.position = "right",
        legend.direction = "vertical")


F_VitalRates_error <- ggarrange(s.plot, r.plot,
          ncol=1, nrow=2)

ggsave(here::here("plots","F_VitalRates_error.png"), width=6, height=8, units="in")

ggarrange(F_VitalRates,F_VitalRates_error,
          ncol=2, nrow=1,
          labels="AUTO")
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/plot%20repro%20and%20surv%20together-error-1.png)

``` r
ggsave(here::here("plots","F_VitalRates_combined.png"), width=12, height=8, units="in")
```

Create lambda model to plot map, used in SurvivalMechanisms analysis
--------------------------------------------------------------------

``` r
lambda.dat <- lambda.loop%>%
  group_by(ndvi,hii)%>%
  summarise(l=mean(l))

m1 <- lm(l~ndvi + hii, data=lambda.dat)
summary(m1)

saveRDS(m1, "/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/lambdamod/m1.rds")
```

immigration required, used in SurvivalMechanisms analysis
---------------------------------------------------------

``` r
lambda.loop%>%
  group_by(hii)%>%
  summarise(l=mean(l)%>%as.numeric())%>%
  mutate(im=case_when((1-l)>0 ~ 1-l,
                      (1-l)<=0 ~ 0))%>%
  write_csv("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/immi_rqrd.csv")
```

How long to make a survivor?
----------------------------

``` r
# pred.boot%>%filter(Sex%in%c("F", "pooled"))%>%group_by(Age,hii)%>%
#               summarise(pred=mean(pred))%>%rename(Survival=pred)%>%filter(hii==14 & Age %in% (7:25))%>%pull(Survival)%>%mean()
##aim for >90% survival

ggplot()+
    geom_path(data=pred.boot%>%filter(Age%in%c(3:25)&
                     hii%in%c(0,15,25,40)) , aes(x=Age, y=pred, color=as.factor(hii), group=iter),size=0.1, alpha=0.05) +
    geom_line(data=pred.boot%>%filter(Age%in%c(3:25)&
                     hii%in%c(0,15,25,40))%>%
                group_by(Age,hii)%>%
                summarise(pred=mean(pred)), aes(x=Age, y=pred, color=as.factor(hii)),size=2) +
  #geom_point(data=data.frame(study="McLellan FH all years", hii=9.3, l=1.046), aes(x=hii, y=l,label=study))+
  #geom_text(data=data.frame(hii=9.4,l=1.05), aes(hii, l), label="McLellan FH all years", vjust=-1)+
  ylab("Survival")+
  xlab("Age")+
  xlim(3,23)+
  geom_hline(yintercept = 0.90, linetype="dashed")+
  scale_colour_Publication()+
  scale_fill_Publication()+
  theme_Publication()+
  theme(legend.title = element_text())
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/survivor%20time-1.png)

``` r
ggsave(here::here("plots","survtime.png"), width=5, height=4, units="in")
 
 

 ###For 5 14yo
s.rates <- pred.boot%>%filter(Age%in%c(0:14)&
                     hii%in%c(5))%>%
                group_by(Age,hii)%>%
                summarise(pred=mean(pred))%>%
  pull(pred)
  
n=100 
for(i in 1:length(s.rates)){
  n<-n*s.rates[i]
  print(n)
}                       
                          
n.bears10 <- 100/n

###For 25  
s.rates <- pred.boot%>%filter(Age%in%c(0:14)&
                     hii%in%c(25))%>%
                group_by(Age,hii)%>%
                summarise(pred=mean(pred))%>%
  pull(pred)
  
n=100 
for(i in 1:length(s.rates)){
  n<-n*s.rates[i]
  print(n)
}                       
                          
n.bears <- 100/n
```

Immigrants supply the raw material for coexistence- for every successful 14-year-old coexister (defined as an adult bear that reaches an adult survival rate is similar to lambda&gt;=1 backcountry areas~ survival of 0.90), there will be approximately 29 dead conflict bears (HII=25). In comparison, only 4 bears would die before reaching 14 in wilder areas (HII=5)
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Immigration analysis
--------------------

``` r
##list to populate
ch.list <- list()

#projects
projects <- list.files(here::here("data","CapHists_SECR", "CapData"))



sex.dat <- data.frame()
for(i in 1:length(projects)){
  ##Read capture history files
sex.dat  <- rbind(sex.dat ,read_csv(here::here("data","CapHists_SECR", "CapData", paste0(projects[[i]]))))
                     
                     }

for(i in 1:length(projects)){
  ##Read capture history files
  a <- read.capthist(here::here("data","CapHists_SECR", "CapData", paste0(projects[[i]])), 
                     here::here("data","CapHists_SECR", "TrapData", paste0(projects[[i]])),
                     detector ="proximity",
                     trapcovnames="Trap_Type",
                     covnames="Sex",
                     sep = ",")
  
  ch.list[[i]] <- a
  names(ch.list)[i]<- str_replace(projects[[i]],".txt","")
}


grizzCH <- MS.capthist(ch.list)


##fix up sex covariate factor levels
for (i in 1:length(grizzCH))
    levels(covariates(grizzCH[[i]])$Sex) <- c('F','M','U')

grizzCH <- shareFactorLevels(grizzCH, columns=1)

for (i in 1:length(grizzCH)){
  levels(covariates(traps(grizzCH[[i]]))$Trap_Type) <- c('BS','RT')
}
rm(ch.list)


##make into DF
caps <- plyr::ldply(as.data.frame(grizzCH), data.frame)
traps <- data.frame()
for(i in 1:length(session(grizzCH))){
  traps <-  rbind(traps, 
                  data.frame(Sess=session(grizzCH)[i],TrapID=dimnames(traps(grizzCH)[[i]])[[1]], traps(grizzCH)[[i]]))
}
traps <- distinct(traps)


###join
df <- left_join(caps, traps, by="TrapID")

##clean
df <- df%>%mutate(year=str_sub(Session, -4, -1)%>%as.numeric())%>%as_tibble()

keep <- df%>%group_by(ID)%>%summarise(len=max(year)-min(year))%>%filter(len>0)%>%pull(ID)

df2 <- df%>%
  filter(ID %in% keep)


df2 <- st_as_sf(df2, 
                coords = c("x","y"),
                crs = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")%>%
  st_transform("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#mapview(df2)


##extract hii to imi
df2 <- df2%>%mutate(hii=hiiv$extract_points(sp = df2))



df3 <- df2%>%group_by(ID,year)%>%summarize(hii=mean(hii,na.rm=TRUE))%>%
  filter(year%in%c(max(year), min(year)))%>%
  mutate(cap=case_when(year==min(year)~"first",
                       year==max(year)~"last"))%>%
  spread(cap, hii)%>%
  group_by(ID)%>%
  summarise(first=mean(first,na.rm=TRUE),
            last=mean(last,na.rm=TRUE))%>%
  mutate(dif=last-first)



names <- df3%>%filter(last>13 & dif>5)%>%pull(ID)
names <- df3%>%pull(ID)
df4 <- df2%>%filter(ID%in%names)%>%group_by(ID,year)%>%summarize(hii=mean(hii,na.rm=TRUE))%>%
  filter(year%in%c(max(year), min(year)))%>%
  mutate(cap=case_when(year==min(year)~"first",
                       year==max(year)~"last"))



inds <- unique(df4$ID)
dist.df <- data.frame()
for( i in 1:length(inds)){
  a <- df4%>%filter(ID%in%inds[i])
  dist.df <- rbind(dist.df,
                   data.frame(ID=inds[i],
                              dist=st_distance(a[1,], a[2,])))
}

dist.df <- dist.df%>%left_join(df3)


dist.df<- dist.df%>%
  mutate(Type=case_when(last<13~"Source",
                        last>=13~"Sink"))


# ggplot(dist.df%>%left_join(sex.dat), aes(y=last,x=as.numeric(dist)/1000, color=Type))+
#   geom_point(aes(shape=Sex))+
#   stat_ellipse(type="norm", level=0.95)+
#   xlim(-1,84)+
#   xlab("Distance (km)")+
#   ylab("Settled HII")+
#   theme_Publication()

#ggsave(here::here("plots","source_sink_movement_DNA.png"), width=5, height=4, units="in")

dist.summary <- dist.df%>%
  group_by(Type)%>%
  mutate(dist=dist/1000)%>%
  summarise(mean=mean(dist),
            max=max(dist),
            min=min(dist),
            quantile10=quantile(dist, probs = c(0.1)),
                         quantile90=quantile(dist, probs = c(0.9)))%>%
  mutate(Type=fct_relevel(Type, "Source",  "Sink"))

dist.summary%>%
  ggplot(aes(x=Type, y=mean, color=Type))+
    ylab("Distance (km)")+
  geom_point(size=3)+
  geom_errorbar(aes(x=Type, ymin=quantile10, ymax=quantile90), width=0.2, size=1) + 
  theme_Publication()+
    theme(legend.position = "none")+
    coord_flip()
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/Immigration%20analysis-1.png)

``` r
  ggsave(here::here("plots","source_sink_movement_DNA2.png"), width=5, height=4, units="in")
```

### Summary of distance between first and last annual home range center for 430 animals who were monitored via genetic tags for &gt;1 year. Sink means their home range center for their last year of capture was in an areas where population growth was &lt;1.

``` r
kable(dist.summary)
```

| Type   |       mean|       max|  min|  quantile10|  quantile90|
|:-------|----------:|---------:|----:|-----------:|-----------:|
| Sink   |  16.143385|  81.27057|    0|   0.6974934|    39.73104|
| Source |   8.028924|  69.71263|    0|   0.0000000|    21.43988|

Prep data for basin of coexistence figure
-----------------------------------------

``` r
surv <- read_csv("surv_basin.csv")

basin.loop <- data.frame()
for(i in 1:nrow(surv)){
    a <- surv[i,]

if(a$noc==100){
Px <- c(0.7, 0.87, rep(a$pred.surv, times=29))
}
    
 if(a$noc==50){
Px <- c(0.7*.89, 0.87*.89, rep(a$pred.surv, times=29))
 }

Fx <- round(pred.boot.repro%>%
                   group_by(Age)%>%
                   summarise(pred=mean(pred)/2)%>%
                   rbind(data.frame(Age=c(0:3,29:31),pred=0))%>%
                   rename(repro=pred)%>%
              arrange(Age)%>%
              pull(repro),4)
L <- odiag(Px, -1)
L[1, ] <- Fx



basin.loop <- rbind(basin.loop,
                     data.frame(hii=a$hii,
                                noc=a$noc,
                                surv=a$pred.surv,
                                l=eigen(L)$values[1]))
}

m2 <- lm(l~ hii*noc, data=basin.loop )
summary(m2)



write_csv(basin.loop%>%mutate(l=as.numeric(l)), "/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/coexist_basin/lambda/lambda_noc.csv")

saveRDS(m2, "/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/coexist_basin/lambda/lambda_noc.rds")

basin.loop%>%
  mutate(l=as.numeric(l))%>%
ggplot(aes(x=hii, y=l, color=as.factor(noc), group=as.factor(noc)))+geom_hline(yintercept = 1, type="dashed")+geom_line()+theme_bw()+xlim(0,40)
```

additional analyses following revision 1
========================================

offspring observed across HII and nocturnal gradient
----------------------------------------------------

``` r
##calculate Ro for high HII and nocturnality animals
dat <- read_csv("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Telemetry_Prov_Grizz/Analyses/coexistence/relocs_nocturnality_monthly.csv")

repro.rate.noc <- repro%>%
  ungroup()%>%
  filter(AgeClass%in%c("COY", "1YO", "2YO", NA))%>%
  mutate(AgeClass=droplevels(AgeClass))%>%
  left_join(dat%>%rename(Year=y)%>%group_by(BearID,Year)%>%summarise(noc.mean=mean(noc,na.rm=TRUE),
                                                                     hii.mean=mean(hii,na.rm=TRUE),
                                                                     ndvi.mean=mean(ndvi,na.rm=TRUE),
                                                                     Age=round(mean(Age2),0)), by=c("BearID","Year"))%>%
  filter(Age>5)%>%
  drop_na(noc.mean)


##Dead F's
##bind repro data together
repro.rate.noc<- repro.rate.noc%>%
  ungroup()%>%
  left_join(mort.cubs, by=c("BearID","Year"))%>%
  mutate(Status=case_when(is.na(mort)~"lived", mort==1~"died"))

repro.rate.noc%>%
  filter(Age>5)%>%
  ggplot( aes(y=noc.mean, x=hii.mean, size=Count, color=Status)) +
  geom_jitter(fill=NA,  shape=21)+
  ylab("Nocturnal (%)")+
  xlab("Human Influence Index")+
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
        legend.key.size= unit(0.6, "cm"),
        legend.spacing = unit(0.1, "cm"),
        legend.text = element_text(size = rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))+
    guides(size=guide_legend(title="Offspring (#)"))+
  scale_colour_manual(values=c("red","forestgreen"))+
  scale_fill_manual(values=c("red","forestgreen"))
```

![](Source_sink_boot_contsAge_files/figure-markdown_github/calc%20Ro-1.png)

``` r
ggsave(here::here("plots", "noc_hii_cubs.png"), width=5, height=3.5, units="in")



  glm(Count~hii.mean,
      data=repro.rate.noc%>%
  filter(Age>5),
  family="poisson")%>%summary()%>%print()
```

    ## 
    ## Call:
    ## glm(formula = Count ~ hii.mean, family = "poisson", data = repro.rate.noc %>% 
    ##     filter(Age > 5))
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.340  -1.256  -1.147   1.001   1.921  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept) -0.10734    0.17233  -0.623    0.533
    ## hii.mean    -0.01793    0.01707  -1.050    0.294
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 151.71  on 110  degrees of freedom
    ## Residual deviance: 150.55  on 109  degrees of freedom
    ## AIC: 274.08
    ## 
    ## Number of Fisher Scoring iterations: 5

Offspring per female
--------------------

``` r
  glm(Count~group, family="poisson", data=repro.rate.noc%>%
              ungroup()%>%
              mutate(group=case_when(hii.mean>=20 & noc.mean>65~"Noc.Risky",
              hii.mean<=10~"Wild"))%>%
              drop_na(group))%>%summary()
```

    ## 
    ## Call:
    ## glm(formula = Count ~ group, family = "poisson", data = repro.rate.noc %>% 
    ##     ungroup() %>% mutate(group = case_when(hii.mean >= 20 & noc.mean > 
    ##     65 ~ "Noc.Risky", hii.mean <= 10 ~ "Wild")) %>% drop_na(group))
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.633  -1.293  -1.293   1.078   1.828  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)   0.2877     0.3536   0.814    0.416
    ## groupWild    -0.4675     0.3744  -1.249    0.212
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 120.01  on 84  degrees of freedom
    ## Residual deviance: 118.63  on 83  degrees of freedom
    ## AIC: 222.32
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
repro.rate.noc%>%
              ungroup()%>%
              mutate(group=case_when(hii.mean>=20 & noc.mean>65~"Noc.Risky",
              hii.mean<=10~"Wild"))%>%
              drop_na(group)%>%
  group_by(group)%>%
  summarise(count=mean(Count))%>%print()
```

    ## # A tibble: 2 x 2
    ##   group     count
    ##   <chr>     <dbl>
    ## 1 Noc.Risky 1.33 
    ## 2 Wild      0.835

Additional analyses following revision 2
========================================

S,R and lambda table
--------------------

``` r
pred.boot%>%
  filter(Sex%in%c("F", "pooled"))%>%
  group_by(Age,hii)%>%
  summarise(predict=mean(pred),
                        Surv.lower95=quantile(pred,0.05),
                        Surv.upper95=quantile(pred,0.95))%>%
  rename(Survival=predict)%>%
  mutate_all(round, 3)%>%
  left_join(
pred.boot.repro%>%
                   group_by(Age,hii)%>%
                   summarise(predict=mean(pred/2),
                        Repro.lower95=quantile(pred/2,0.05),
                        Repro.upper95=quantile(pred/2,0.95))%>%
                    ungroup()%>%
                   rbind(expand.grid(Age=c(0:3,29:31), hii=c(0:40),predict=0,Repro.lower95=0,Repro.upper95=0)%>%as_tibble())%>%
                  rename(Reproduction=predict)%>%
                  arrange(Age)%>%
                  mutate_all(round, 3),
by=c("Age","hii"))%>%
  arrange(hii,Age)%>%
  write_csv(here::here("tables", "vital_rates.csv"))




lambda.loop%>%
              group_by(hii)%>%
              summarise(lambda=mean(l)%>%as.numeric(),
                        lower95=quantile(l,0.05)%>%as.numeric(),
                        upper95=quantile(l,0.95)%>%as.numeric())%>%
              mutate_all(round, 3)%>%
  write_csv(here::here("tables", "population_growth.csv"))
```

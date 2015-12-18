# CHaMP Extrapolation Code
# Matt Nahorniak

##Read GAA files.  These are BIG and take a long time to read.  Read them once per session,
## then comment them out so as to no re-read them every time through the code.
#GAA.data = read.csv("IC_Sites_withMetrics_20151016.csv")
#Saved.GAA.data = GAA.data
#GAA_sheds = read.csv("GAA_Sites_20151006.csv", header=T)
#Saved.GAA_sheds = GAA_sheds


nrow(saved.GAA.data)
windows(record=T)
detach(data)


w=1


# Specify all watersheds that have CHaMP data and that we want to include
# in the model fitting process.
All.ModeledWatersheds = c("Tucannon", "Entiat", "Yankee Fork", 
                       "Upper Grande Ronde","Wenatchee", 
                      "Methow","South Fork Salmon", "Lemhi","John Day")
#All.ModeledWatersheds = c("Tucannon", "Entiat", "Yankee Fork", 
#                       "Upper Grande Ronde")

ModeledWatersheds = All.ModeledWatersheds
PredictedWatersheds = c("South Fork Clearwater", "Lower Clearwater","Lolo","Lochsa","Upper Salmon River Tributaries Above Redfish Lake")


# Comment out next four lines unless we're running through cross-watershed validation
for (w in 1:length(All.ModeledWatersheds)){
print(All.ModeledWatersheds[w])
ModeledWatersheds = All.ModeledWatersheds[1:length(All.ModeledWatersheds)!=w]
PredictedWatersheds = All.ModeledWatersheds[w]


cat.vars = NULL
cat.vars

GAA.data = Saved.GAA.data

# Some code to clean up " " values in "ChanlType"
levels(GAA.data$ChanlType)=c("Unknown", levels(GAA.data$ChanlType))
GAA.data$ChanlType[GAA.data$ChanlType==" "] = "Unknown"
GAA.data$ChanlType = factor(GAA.data$ChanlType)
levels(GAA.data$ChanlType)

# Some code to clean up " " values in "NatClass"
levels(GAA.data$NatClass)=c("Unknown", levels(GAA.data$NatClass))
GAA.data$NatClass[GAA.data$NatClass==" "] = "Unknown"
GAA.data$NatClass = factor(GAA.data$NatClass)
levels(GAA.data$NatClass)


# Some code to clean up " " values in "DistClass"
levels(GAA.data$DistClass)=c("Unknown", levels(GAA.data$DistClass))
GAA.data$DistClass[GAA.data$DistClass==" "] = "Unknown"
GAA.data$DistClass = factor(GAA.data$DistClass)
levels(GAA.data$DistClass)



GAA_sheds=Saved.GAA_sheds
nrow(GAA_sheds)
names(GAA_sheds)

names(GAA_sheds)
levels(GAA_sheds$EP_SHED)

#temp=GAA_sheds[GAA_sheds$EP_SHED %in% PredictedWatersheds,]
#plot(temp$X_ALBERS, temp$Y_ALBERS,
#col=match(temp$EP_SHED, levels(factor(temp$EP_SHED))))
## To Do Friday: Model for SubLT6, LWD (?)
## fit vs predicted plots and R-squared for Sin, SubLT2, SubLT6, LWD(?)
## extrapolations for Sin, SubLT2, SubLT6, LWD(?) to Jean

# Note: USalmon most like Yankee Fork according to coverage plots



#Load required libraries
#library(gstat)
library(survey)
#library(spsurvey)
library(gplots)







########3

#Pick the watershed(s) in which to develop model
#ModeledWatersheds = c("Tucannon", "Entiat", "Yankee Fork", 
#                       "Upper Grande Ronde","Minam","Wenatchee", 
#                      "Methow","South Fork Salmon", "Lemhi","John Day")
#ModeledWatersheds = c("Tucannon",  
#                       "Upper Grande Ronde","Minam", "Wenatchee"
#                      "Methow","South Fork Salmon", "Lemhi","John Day")


#ModeledWatersheds = c("Upper Grande Ronde")
#ModeledWatersheds = c("Yankee Fork")

#PredictedWatersheds = c("Tucannon", "Entiat", "Yankee Fork", 
#                       "Upper Grande Ronde", "South Fork Salmon", 
#                       "Upper Salmon River Tributaries Above Redfish Lake",
#                       "South Fork Clearwater", "Lower Clearwater","Lolo","Lochsa")


#PredictedWatersheds = c("South Fork Clearwater")
#PredictedWatersheds = c("Lower Clearwater")
#PredictedWatersheds = c("Lolo")
#PredictedWatersheds = c("Lochsa")
#PredictedWatersheds = c("Yankee Fork")


#PredictedWatersheds = c("Upper Salmon River Tributaries Above Redfish Lake")
#PredictedWatersheds = c("Upper Grande Ronde")
#PredictedWatersheds = c("Asotin", "Entiat", "Yankee Fork")
#PredictedWatersheds = c("South Fork Clearwater", "Lower Clearwater","Lolo","Lochsa","Upper Salmon River Tributaries Above Redfish Lake")


#PredictedWatersheds = c("Upper Salmon River Tributaries Above Redfish Lake",
#                       "South Fork Clearwater", "Lower Clearwater","Lolo","Lochsa")


#ModeledWatersheds = c("Entiat")
#ModeledWatersheds = c("Tucannon")



#Specify which metrics need to be log transformed
logmetrics = c(
"Sin",
"SubLT2",
"SubLT6",
"HSI_CH_Juv_WUA_per_m",
"HSI_CH_Spawn_WUA_per_m",
"HSI_St_Juv_WUA_per_m",
"HSI_St_Spawn_WUA_per_m",
"LWFreq_Bf",
"PoolResidDpth",
"DpthThlwg_UF_CV")
#,"SlowWater_Pct")


#####for (metric.name in metric.list) {


#Metric Names Here
#metric.name = "Sin" 
#metric.name = "SubLT2"
#metric.name = "SubLT6"
#metric.name = "HSI_CH_Juv_WUA_per_m"
#metric.name = "HSI_CH_Spawn_WUA_per_m"
#metric.name = "HSI_St_Juv_WUA_per_m"
#metric.name = "HSI_St_Spawn_WUA_per_m"
#metric.name = "NREI.dens.est.fpm"

# Additional metrics - non EPP but important for extrapolation
metric.name = "LWFreq_Bf" # Sucky Model
#metric.name = "PoolResidDpth"  # 
#metric.name = "DpthThlwg_UF_CV"  # This One needs a model
#metric.name = "SlowWater_Freq" ##  This One needs a model or the next one
#metric.name = "SlowWater_Pct" ##  This One needs a model
#metric.name = "SubD84"
#metric.name = "SubD50"


# Note:
# Need to first run grts analysis in order to generate 
# "METRICS_Ave_of_All_Years.csv" and  "CHaMP_Data_Mean_AdjWgt_by_Metric.csv" files

# Read weights and data files - these files are produced by the GRTS and Var Decomp process
data.wgt = read.csv("CHaMP_Data_Mean_AdjWgt_by_Metric.csv", header=T)
data = read.csv("METRICS_Ave_of_All_Years.csv", header=T)

nrow(data)
nrow(data.wgt)

# Assign appropriate data to "metric" and link weight to data data.frame
metric.idx = match(metric.name, names(data))
data$metric = data[,metric.idx]
nrow(data)

# Assign the weight fromt the weights file for the specific metric
metric.wgt.idx = match(metric.name, names(data.wgt))
data$wgt = data.wgt[,metric.wgt.idx]



# For some metrics, it's best to use the natural log of the values.  List is defined above.
if (metric.name %in% logmetrics){
data$metric = log(1+data$metric)}


# Read the GAA data (supplied by Carol) and assign "GAA" to the start of each
# name to easiliy identify these variables as potential GAA's
#
# Note: There have been lots of GAA files!  Effort required to continually ensure we have the
# latest, greatest version

#	GAA.data = read.csv("GAA_data_4_30_2014.csv", header=T)
#	GAA.data = read.csv("GAA_Site_WithMetrics_20150810.csv", header=T)
#      GAA.data = read.dbf("IC_PriorityWatersheds_GAAMetrics_20150924/IC_PriorityWatersheds_GAAMetrics_20150924.dbf")

#GAA.data = read.dbf("IC_Tucannon_Sites_20150917_Metrics/IC_Tucannon_Sites_20150917_Metrics.dbf")#$dbf
#GAA.data = read.dbf("Tucannon_GAA_Metrics_20150925/Tucannon_GAA_Metrics_20150925.dbf")
#GAA.data = read.csv("IC_Sites_withMetrics_20151016.csv")

#############
if (1==2) {
LWD_Sin_Conf.data = read.csv("GAA_EntTucYFGRUS_LWDSinConf_attributes_20151029.csv",header=T)
names(LWD_Sin_Conf.data)

# Deal w/ -9999 values - need to call them NA
GAA.data[GAA.data == -999]=NA
LWD_Sin_Conf.data[LWD_Sin_Conf.data == -999]=NA
GAA.data[GAA.data == -9999]=NA
LWD_Sin_Conf.data[LWD_Sin_Conf.data == -9999]=NA


temp=merge(GAA.data, LWD_Sin_Conf.data, by= "Site_ID")
nrow(GAA.data)

# temporarily comment this out
GAA.data = temp
}
nrow(GAA.data)
#############################

nrow(GAA_sheds)
nrow(GAA.data)
names(GAA_sheds)

# Read GAA_sheds at the top now, use either name "EP_shed" or "CHaMPsheds" to define watershed name
#GAA_sheds = read.csv("GAA_Sites_20151006.csv", header=T)

# Add "John Day" to "EP Sheds" since I'm using that to name my watershed.  May have to do this with
# Other watersheds too
levels(GAA_sheds$EP_SHED) = c(levels(GAA_sheds$EP_SHED), "John Day")
GAA_sheds$EP_SHED[GAA_sheds$CHaMPsheds=="John Day"] = "John Day"
levels(factor(GAA_sheds$EP_SHED))

shed.idx = match(GAA.data$Site_ID, GAA_sheds$Site_ID)
shed.idx
GAA.data$EP_SHED = GAA_sheds$EP_SHED[shed.idx] 



# Tucannon only as of 10/14/2015
# Tucannon - merge LWD attributes w/ other attributes.  
if (1==2) {
GAA.data2 = read.csv("Tucannon_GAA_Metrics_20150925/GAA_Tuc_LWDattributes_20150929.csv")
GAA.data=merge(GAA.data, GAA.data2, by.x="Site_ID", by.y="Site_ID")
}


# Tucannon only as of 10/14/2015
if (1==2) {
# Read the temperature data from Kris McNyset
GAA.Temp.Data = read.csv("GAA_Ave_Temp_Metrics.csv",header=T)
idx= match(GAA.data$Site_ID,GAA.Temp.Data$Site_ID)
idx
GAA.data$Pct12 = GAA.Temp.Data$Pct12[idx]
GAA.data$Pct13 = GAA.Temp.Data$Pct13[idx]
GAA.data$Pct16 = GAA.Temp.Data$Pct16[idx]
GAA.data$Pct18 = GAA.Temp.Data$Pct18[idx]
GAA.data$Pct20 = GAA.Temp.Data$Pct20[idx]
GAA.data$Pct22 = GAA.Temp.Data$Pct22[idx]
GAA.data$MxMx = GAA.Temp.Data$MxMx[idx]
GAA.data$sdMn = GAA.Temp.Data$sdMn[idx]
GAA.data$MnMx = GAA.Temp.Data$MnMx[idx]
}


names(GAA.data) = paste("GAA.", names(GAA.data), sep="")
names(GAA.data)
GAA.data$GAA.Site_ID

"GAA.Sin" %in% colnames(GAA.data)
# Remove sites w/ sinuosity = 0 if we're using sinuosity.  Hopefully a small subset of site - or else
# we're throwing out a lot of good data with the bad.
if ("GAA.Sinuosity" %in% colnames(GAA.data)){
GAA.data = GAA.data[GAA.data$GAA.Sin != 0,]}


# Vestiges of prior code - used terminology "Frame Data" instead of "GAA.data"
# because we use to use the pop frame as the source of the original GAA's.

#HERE!!!
names(GAA.data)
frame.data = GAA.data

names(GAA.data)
GAA.data$GAA.Site_ID
names(GAA.data)
data$SiteName
GAA.data$GAA.Site_ID
#########################################################
#Merge the GAA.data into the data
	idx = match(data$SiteName,GAA.data$GAA.Site_ID)
	data = data.frame(data, GAA.data[idx,])
data$SiteName
idx
idx[is.na(idx)==F]
# Reduce data to only data for which we have GAA info (or at least a matching site in the GAA file - we'll deal
# with N/A values later
data = data[is.na(idx)==F,]
levels(data$WatershedName)

##################################
# Reduce data to Watershed(s) selected at the top of script
	data = data[data$WatershedName %in% ModeledWatersheds,]
	data.saved=data

      val.predict = rep(NA, nrow(data))
length(val.predict)
nrow(data.saved)

val.loop=1
# For Cross Site Validation we'll loop through every CHaMP site for leave one out cross 
# validation.  If not doing leave one out cross validation, we'll just leave the first
# site out and go through the "loop" once.  Also do that if we're doing leave one watershed
# at a time out of the analysis.

## Use this when cross-validatint
#for (val.loop in 1:nrow(data.saved)){
## Use this when developing models
for (val.loop in 1:1){

v.data.idx= val.loop
val.data = data.saved[val.loop,]
nrow(val.data)
data = data.saved[(1:nrow(data)) %in% val.loop == FALSE,]
data = data.saved
print(paste("val loop", val.loop, "of", nrow(data.saved)))

#################################################################
####################################################################

# GAA's are all identified by "GAA" in the name.  Create an index of these.
	GAA.idx = (1:ncol(data))[grepl("GAA", names(data))]
	

# my.design provides design info (weights) to svyglm functions used for regression analysis
	my.design = svydesign(id= data$SiteName, weights=data$wgt, data=data, replacement=F)


#############################################3333
# Workspace for buiilding models

metric.name
#hist(data$metric)

attach(data)
#hist(metric)
mean(data$metric)
# Starter mod


mod = svyglm(metric ~ 1, design = my.design)
summary(mod)

idx=GAA.idx[1]


# GAA.idx is used in the "add1" function to list which GAA attributes we should consider.
# We have to weed out variables that have only a single value, or attributes with missing data, or the
# add1 function crashes.  That's what this mess, below, does.
for (idx in GAA.idx){
 #  print(length(levels(factor(na.omit(data[,idx])))))
    if (length(levels(factor(na.omit(data[,idx])))) < 2)
     {GAA.idx = GAA.idx[GAA.idx != idx]} 
       else {
          # print(length(na.omit(data[,idx])))
          if  (length(na.omit(data[,idx])) < (.8*length(metric))) {
            GAA.idx = GAA.idx[GAA.idx != idx]
           } 
          }
        }


# Don't include "Site_ID" as a potential GAA, cuz that's just stoopid.
GAA.idx = GAA.idx[GAA.idx != match("GAA.Site_ID", names(data))]

# Temp, while I don't have confinement or sinuosity across all basins
#GAA.idx = GAA.idx[GAA.idx != match("GAA.Confinem_3", names(data))]
#GAA.idx = GAA.idx[GAA.idx != match("GAA.Sinuosity", names(data))]


## plots are useful to, in addition to the "add1" function.
#for (a in GAA.idx){
#plot(data[,a],main= paste(a, colnames(data)[a]),metric,ylab=metric.name, xlab=colnames(data)[a],pch=19)
#}

# Use add1 and drop1, w/ AIC to see if it appears to make sense to add/subtract variables to model.

# don't do these lines while not building models...
if (1==2) {
add=add1(mod, scope = formula(data[,GAA.idx],2), data=data)
add.idx = order(add$AIC,decreasing=T)
add[add.idx,]
add[1,]
}

metric.name


if (metric.name == "SlowWater_Pct") {
mod = svyglm(metric ~1+GAA.StrmPwr+GAA.MeanU_v1+GAA.NatClass+GAA.ValleyClas+GAA.GradICB200, design=my.design)
cat.vars = c("GAA.NatClass", "GAA.ValleyClas")
#mod = svyglm(metric ~ 1+GAA.StrmPwr+GAA.MeanU_v1+GAA.ValleyClas+GAA.NLCD_24, design = my.design)
#mod = svyglm(metric ~1+GAA.StrmPwr+GAA.MeanU_v1+GAA.NatClass+GAA.Slope_GEO+GAA.NLCD_21, design=my.design)
#cat.vars = c("GAA.NatClass")
}


metric.name



if (metric.name=="DpthThlwg_UF_CV") {
#mod = svyglm(metric ~ 1+GAA.ValleyClas+GAA.Winter95+GAA.SO_v1+GAA.Winter2yr, design = my.design)
#cat.vars = c("GAA.ValleyClas")
mod = svyglm(metric ~ 1+GAA.ChanlType+GAA.Winter95+GAA.StrmPwr+GAA.NLCD_52+GAA.BFW_M, design = my.design)
cat.vars = c("GAA.ChanlType")

}





#D84, global model 12/10/15
if (metric.name=="SubD84"){
mod = svyglm(metric ~ 1+GAA.StrmPwr+GAA.NatPrin1+GAA.TRange+
GAA.ChanlType+GAA.MeanU_v1+GAA.WIDE_WW, design=my.design)
cat.vars = c("GAA.ChanlType")
}

#D50, using sampe covariates as D84 global model 12/10/15
if (metric.name=="SubD50"){
mod = svyglm(metric ~ 1+GAA.StrmPwr+GAA.NatPrin1+GAA.TRange+
GAA.ChanlType+GAA.MeanU_v1+GAA.WIDE_WW, design=my.design)
cat.vars = c("GAA.ChanlType")
}

##########################################
##########################################
# Best Fit Models 9/10/15 Tucannon Only


###########################################
#FishCovLW
if (metric.name == "FishCovLW"){
mod = svyglm(metric ~ 1+GAA.CURRSUSH+GAA.MAVELV, design = my.design)
}



########################
# PoolResidDpth - Lousy Model!!!! - Re-Do
if (metric.name == "PoolResidDpth"){
mod=svyglm(metric ~ 1+GAA.WIDE_WW+GAA.StrmPwr+GAA.DistClass+GAA.WIDE_BF, design = my.design)
cat.vars = c("GAA.DistClass")
}
summary(mod)
levels(factor(data$GAA.DistClass))

########################
# LWFreq_BF - Lousy Model!!! - ReDo!
if (metric.name == "LWFreq_Bf"){
#global LWFreq_BF MOdel
mod = svyglm(metric ~ 1+GAA.Summer20p+GAA.ChanlType+GAA.TRange+GAA.Winter95+GAA.MAFLOWV, design=my.design)
cat.vars = c("GAA.ChanlType")
}



########################
# SubLT6
if (metric.name == "SubLT6"){ 
# Tucannon Specific Model
#mod = svyglm(metric ~ 1+GAA.Summer95+GAA.NatClass, design = my.design)
#cat.vars = c("GAA.NatClass")

summary(mod)
#Global Model for SubLT2
mod = svyglm(metric ~ 1+GAA.NatClass+GAA.NLCD_11+GAA.NLCD_82+GAA.CURRSUSH+GAA.SO_v1, design = my.design)
cat.vars = c("GAA.NatClass")

#Upper Grande Ronde SubLT2
#mod = svyglm(metric ~ 1+GAA.NLCD_82+GAA.WIDE_WW+GAA.CumDran_v1+GAA.COMID_V1, design = my.design)

}

###################
# SubLT2
if (metric.name == "SubLT2"){
#Global Model for SubLT2
mod = svyglm(metric ~ 1+GAA.NatClass+GAA.NLCD_11+GAA.NLCD_82+GAA.CURRSUSH+GAA.SO_v1, design = my.design)
cat.vars = c("GAA.NatClass")

#Tucannon
#mod = svyglm(metric ~ 1+GAA.Summer95+GAA.NatClass, design = my.design)
#cat.vars = c("GAA.NatClass")

#Upper Grande Ronde SubLT2
#mod = svyglm(metric ~ 1+GAA.NLCD_82+GAA.WIDE_WW+GAA.CumDran_v1+GAA.COMID_V1, design = my.design)


}



###################
# "Sin"
if (metric.name == "Sin"){
# Global(ish) model for Sin
mod = svyglm(metric ~ 1+GAA.Sin+GAA.StrmPwr+GAA.ForHt5_10_VgHtLF_propPts65b+GAA.Conf, design = my.design)
#mod = svyglm(metric ~ 1+GAA.Sin+GAA.NewErode+GAA.NatClass+GAA.CumDran_v1, design = my.design)
#cat.vars = c("GAA.NatClass")
}

summary(mod)
################
## "HSI_St_Spawn_WUA_per_m"
if (metric.name == "HSI_St_Spawn_WUA_per_m"){
#HSI_St_Spawn_WUA_per_m all watersheds
mod = svyglm(metric ~ 1+GAA.MAFLOWU+GAA.NLCD_90+GAA.BFW_M+GAA.NLCD_23+GAA.NLCD_52, design = my.design)
#mod = svyglm(metric ~ 1+GAA.Channelflo+GAA.NatClass, design = my.design)
#cat.vars = c("GAA.NatClass")

}


####################3
## "HSI_St_Juv_WUA_per_m"
if (metric.name == "HSI_St_Juv_WUA_per_m"){

# Global Model
mod = svyglm(metric ~ 1+GAA.Channelflo+GAA.SO_v1+GAA.MAXELEVSMO+GAA.Winter2yr+GAA.BFW_M, design = my.design)

}


#############################
### "HSI_CH_Spawn_WUA_per_m"
if (metric.name == "HSI_CH_Spawn_WUA_per_m"){
#mod = svyglm(metric ~ 1+GAA.MAFLOWU+GAA.ChanlType+GAA.NLCD_23+GAA.NLCD_90+GAA.StrmPwr , design=my.design)
#cat.vars  = c("GAA.ChanlType")
#mod = svyglm(metric ~ 1 + GAA.Channelflo + GAA.NatClass, design = my.design)
#cat.vars = c("GAA.NatClass")
mod = svyglm(metric ~1+ GAA.Channelflo+GAA.MAFLOWU+GAA.CumDran_v1+GAA.WIDE_WW+GAA.StrmPwr+GAA.MeanU_v1+GAA.WIDE_BF, design = my.design)

}

#######################
#HSI_CH_Juv_WUA_per_m
if (metric.name == "HSI_CH_Juv_WUA_per_m"){
#mod = svyglm(metric ~ 1 + GAA.WIDE_WW+GAA.CHINRATE , design = my.design)
mod = svyglm(metric ~ 1 + GAA.WIDE_WW+GAA.CumDran_v1+GAA.MAVELU+GAA.DistClass+GAA.flow7q10+GAA.NLCD_41, design=my.design)
#mod = svyglm(metric ~ 1 + GAA.WIDE_WW+GAA.CumDran_v1+GAA.MAVELU+GAA.DistClass+GAA.flow7q10, design=my.design)
cat.vars = c("GAA.DistClass")

}
##########################

# NREI Model 8_28_2015
if (metric.name == "NREI.dens.est.fpm"){
#mod = svyglm(metric ~ 1+GAA.Channelflo+GAA.MINELEVSMO, design = my.design)
### NREI for all watersheds
mod = svyglm(metric ~ 1+GAA.BFW_M+GAA.Channelflo+GAA.CumDran_v1+GAA.WIDE_BF, design = my.design)
}

summary(mod)
########################
######################################################################3


## Consider this for count data or percent data
#mod = svyglm(metric/100 ~ 1+GAA_log_SPower+GAA_Riv_Style + GAA.ValleyClas , design=my.design, family = binomial(link='logit'))





###################################################
# Plot measured vs Fitted
#fitted(mod)
Measured=(mod$fitted.values+mod$residuals)
Fit = fitted(mod)
min(Fit,na.rm=T)
Fit[Fit < -50.6] = NA

## Set the plot color for plotting data
#data$plot.color = rainbow(length(levels(data$WatershedName)))[match(data$WatershedName, levels(data$WatershedName))]
#   data=data[data$wgt > 0,]



if (val.loop ==1) {
# Set the plot color for plotting data
levels(factor(data$WatershedName))

ModeledWatersheds
PredictedWatersheds
data$plot.color = rich.colors(length(levels(data$WatershedName)))[match(data$WatershedName, levels(data$WatershedName))]
   data=data[data$wgt > 0,]

length(levels(data$WatershedName))

data$plot.color[data$WatershedName=="Upper Grande Ronde"]="dark gray"

xlab = "Measured"
if (metric.name %in% logmetrics) {xlab = "log(1+Measured)"}
ylab = "Predicted"
if (metric.name %in% logmetrics) {ylab = "log(1+Predicted)"}



plot(Measured,Fit, main=paste("Extrapolation Model for",metric.name),cex=sqrt(.2+data$wgt/max(data$wgt)*3),
#plot(Fit, metric, main=paste("Extrapolation Model for",metric.name),cex=sqrt(.2+data$wgt/max(data$wgt)*3),
pch=19, col=data$plot.color,
xlab=xlab, ylab=ylab)

Measured

fit.mod = lm(Fit~Measured)

mod.rsquared = summary(fit.mod)$r.squared
mod.rsquared

fits=predict(fit.mod, newdata = data.frame(Measured=seq(0,1000, by=.1)), se=T)
fits
summary(fit.mod)

lines(seq(0,1000, by=.1),fits$fit)
lines(seq(0,1000, by=.1), fits$fit + 1.96* fits$se.fit, lt=2)
lines(seq(0,1000, by=.1), fits$fit - 1.96* fits$se.fit, lt=2)

levels(data$WatershedName)

#legend("topleft", levels(data$WatershedName), pch=19, col=rainbow(length(levels(data$WatershedName))),cex=.7)
leg.cols = rich.colors(length(levels(data$WatershedName)))
leg.cols[leg.cols=="#FFE201FF"] = "dark gray"
legend("topleft", levels(data$WatershedName), pch=19, 
#col=rich.colors(length(levels(data$WatershedName))),
col = leg.cols,
cex=.7)
leg.cols

}


### V here ####
# Predict at validation points
val.data
val.preds = predict.lm(mod, newdata = val.data)
val.predict[val.loop] = val.preds
val.predict[val.loop]

val.loop
detach(data)


}

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#
# Predict at new points
# Assign New data and clean it up for factors such that only sites having factors limited to
# those in the original regression are sites for which we attempt to make predictions
# Will have to include all factor variables in the fixes below.  May blow away a lot of sites
# if we don't have good coverage in our data for which models have been fit.
#
#
# Will want new data to be made from all GAA's (or everything within a given watershed or whatever), not just data.saved
# Match GAA sites IDs to master frame data (TBD)

# copy "frame.data" into newdata for manipulation
newdata = frame.data
newdata = frame.data[frame.data$GAA.EP_SHED %in% PredictedWatersheds,]

data$GAA.NatClass = factor(data$GAA.NatClass)
newdata$GAA.NatClass=factor(newdata$GAA.NatClass)


##########################################################################################################
## Need to find factor variables and reduce "newdata" to only those rows that have factors that
## are in the model.
###########################################################################################################

# Remove any data from "newdata" for which a categorical variable level exists for which
# we don't have a site in CHaMP data at that level - because we can't predict at that level
if (length(cat.vars)>0) {
for (i in 1:length(cat.vars)){
if (TRUE %in% (grepl(cat.vars[i], names(mod$coefficients))))
{
col.nd=match(cat.vars[i], colnames(newdata))
col.d =match(cat.vars[i], colnames(data))

idx = (newdata[,col.nd] %in% levels(factor(data[,col.d]))) *1
newdata = newdata[(idx==1),]
idx = newdata[,col.nd] != " " #Take out blanks as well, if there are any.
newdata = newdata[idx,]
newdata[,col.nd] = factor(newdata[,col.nd])
#actor(newdata[,col.nd], levels= levels(data[,col.d]))
}
}}



#"GAA.NatClass" %in% names(mod$coefficients)
#if ("GAA.NatClass" %in% names(mod$coefficients)){
#idx = (newdata$GAA.NatClass %in% levels(factor(data$GAA.NatClass))) *1
#newdata = newdata[(idx==1),]
#}
###################################

nrow(newdata)

## For models that use GAA.NatClass.  Screen to this (88% of sites covered.  Good enough!)
## Clean up data here: re-assign levels to match those used in regression model
#newdata$GAA.NatClass = factor(newdata$GAA.NatClass, levels=levels(factor(data$GAA.NatClass)))


levels(newdata$GAA.ChanlType)=levels(data$GAA.ChanlType)
levels(newdata$GAA.NatClass) =levels(data$GAA.NatClass)
levels(newdata$GAA.DistClass) =levels(data$GAA.DistClass)
levels(newdata$GAA.ValleyClas) =levels(data$GAA.ValleyClas)

# Make the Predictions!
Preds = predict(mod, newdata=newdata)


# Back-transform logged metrics
if (metric.name %in% logmetrics){
Preds = exp(Preds)-1
data$metric = exp(data$metric)-1
}
length(Preds)
nrow(newdata)

#
# Write Preds to a file....
#
# clean up, weird data format
Preds = Preds[1:length(Preds)]
#length(Preds)
#hist(data$metric)

length(Preds)
# Hokey bit to keep extrapolations within the limits of what's been observed (on the high side)
Preds[Preds > max(data$metric)] = max(data$metric)
Preds[Preds < 0] = 0
#Preds
#hist(Preds)

idx = match(names(Preds), rownames(newdata))

# Generate predictions dataframe
Predictions=data.frame("Site_ID"=newdata$GAA.Site_ID[idx],Preds,
"Chinook"=newdata$GAA.Chinook[idx],"Steelhead"=newdata$GAA.Steelhead[idx],
"EP_SHED" = newdata$GAA.EP_SHED[idx],
"X_ALBERS" = newdata$GAA.X_ALBERS[idx],
"Y_ALBERS" = newdata$GAA.Y_ALBERS[idx]
)


Test = Predictions
Test$GAA.BFW_M=newdata$GAA.BFW_M[idx]
Test$GAA.Channelflo = newdata$GAA.Channelflo[idx]
Test$GAA.CumDran_v1 =newdata$GAA.CumDran_v1[idx]
Test$GAA.WIDE_BF = newdata$GAA.WIDE_BF[idx]

summary(mod)
#Test$GAA.Channelflow
#plot(Test$Preds, Test$GAA.BFW_M)
#plot(Test$Preds, Test$GAA.Channelflo)
#plot(Test$Preds, Test$GAA.BFW_M)
#plot(Test$Preds, Test$GAA.BFW_M)



#Predictions
#nrow(Predictions)
#Predictions=data.frame("Site_ID"=newdata$GAA.Site_ID[idx],Preds,
#"Chinook"=newdata$GAA.Chinook[idx],"Steelhead"=newdata$GAA.Steelhead[idx])





# Add whether a site was a measured site to this output
Predictions$Measured = rep("FALSE", nrow(Predictions))
Predictions$Measured[match(data$SiteName, Predictions$Site_ID)]="TRUE"

# Change colname to metric name for metric
colnames(Predictions)[2]= paste("Predicted_", metric.name, sep="")
#colnames(Predictions)

# Write the file
nrow(Predictions)
names(Predictions)
head(Predictions)

WS.string = ""
for (k in 1:length(PredictedWatersheds)) {
WS.string = paste(WS.string, PredictedWatersheds[k],"_",sep="") 
}


# Only write predictions when we're really ready!
wd=getwd()
setwd("Predictions")
write.csv(Predictions, paste(metric.name,"_",WS.string,"Predicted.csv", sep=""))
setwd(wd)

if (w==1){All.Predictions = Predictions} else  {
All.Predictions = rbind(All.Predictions, Predictions)
}




nrow(Predictions)
names(Predictions)

Preds = Predictions[,2]
col.idx = 1+round(9*(Preds-min(Preds))/(max(Preds)-min(Preds)))
colors = terrain.colors(10)[col.idx]



#colors = match(newdata$GAA.EP_SHED,levels(factor(newdata$GAA.EP_SHED)))
#colors = 1+1*(newdata$GAA.Steelhead[idx]=="Chinook")+
#1*(newdata$GAA.Steelhead[idx]=="Steelhead")

Preds[Preds==0]
colors[Preds==0]="red"

nrow(newdata)
plot(newdata$GAA.X_ALBERS, newdata$GAA.Y_ALBERS)



#plot(newdata$GAA.X_ALBERS[idx],newdata$GAA.Y_ALBERS[idx],pch=19, cex=1)
points(newdata$GAA.X_ALBERS[idx][Predictions$Steelhead=="Steelhead"]
,
newdata$GAA.Y_ALBERS[idx][Predictions$Steelhead == "Steelhead"]
,pch=19, cex=1, 
col=colors[Predictions$Steelhead=="Steelhead"]
)

idx
points(newdata$GAA.X_ALBERS[idx],
newdata$GAA.Y_ALBERS[idx]
,pch=19, cex=1, 
col=colors
)



##idx
##Predictions$Chinook
#points(newdata$GAA.X_ALBERS[idx][Predictions$Chinook=="Chinook"]
#,newdata$GAA.Y_ALBERS[idx][Predictions$Chinook == "Chinook"]
#,pch=19, cex=1,# col="green")
#col=colors[Predictions$Chinook=="Chinook"]+1)



legend.text = round(seq(min(Preds), max(Preds), by=((max(Preds)-min(Preds))/9)),2)
legend("topright", legend= legend.text, col=terrain.colors(10),pch=19)



###############################





#####}
##################################################################################################


#PredictedWatersheds

GAA.data$GAA.X_Albers
p.data=GAA.data[GAA.data$GAA.EP_SHED == PredictedWatersheds,]
#plot(p.data$GAA.X_ALBERS, p.data$GAA.Y_ALBERS, pch=19, col="black")
st.p.data = p.data[p.data$GAA.Steelhead == "Steelhead",]
points(st.p.data$GAA.X_ALBERS, st.p.data$GAA.Y_ALBERS,pch=19, col="red")


ch.p.data = p.data[p.data$GAA.Chinook == "Chinook",]
points(ch.p.data$GAA.X_ALBERS, ch.p.data$GAA.Y_ALBERS,pch=19, col="blue")

legend("topleft", c(PredictedWatersheds, "St Only", "Chin+St"), pch=19, col=c("black","red","blue"))

print("HERE")
###############################


head(All.Predictions)

######
######
######



val.predict[val.predict < 0]= 0
hist(val.predict)
ylim=range(data.saved$metric, na.rm=T)

a=sum((val.predict-data.saved$metric)^2,na.rm=T)
b=sum((data.saved$metric-mean(data.saved$metric, na.rm=T))^2, na.rm=T)
1-a/b

print(metric.name)
print(paste("Cross Validation R-squared=",round((1-a/b),2)))


col.idx=match(data.saved$Watershed,levels(data.saved$Watershed))



data.saved$plot.color = rich.colors(length(levels(data.saved$WatershedName)))[col.idx]
data.saved$plot.color[data.saved$WatershedName=="Upper Grande Ronde"] = "dark gray"

plot(data.saved$metric,val.predict,main=paste("Leave one out cross validation of ",metric.name, sep=""),
xlab="Measured", ylab="Predicted",pch=19,col=data.saved$plot.color, ylim=ylim)
lines(c(0,100), c(0,100), lt=2)

color=rich.colors(length(levels(data.saved$WatershedName)))
levels(factor(data.saved$WatershedName))

color[levels(data.saved$WatershedName) == "Upper Grande Ronde"] = "dark gray"
leg.text=levels(data.saved$Watershed)
legend("topleft", leg.text, pch=19, ,cex=.7,bg="white",col=color)

#########################3
} # end of looping through watersheds for leave-one-watershed out cross-validation

metric.name

names(data)
names(All.Predictions)
idx= match(data$SiteName, All.Predictions$Site_ID)
WS.Val_Preds = All.Predictions[,2][idx]

col.idx=match(data$Watershed,levels(data$Watershed))
color = rich.colors(length(levels(data$WatershedName)))[col.idx]

color[data$WatershedName == "Upper Grande Ronde"] = "dark gray"

#Don't predict anything outside of or max range.
WS.Val_Preds[WS.Val_Preds > max(data$metric,na.rm=T)] = max(data$metric) 

plot(data$metric, WS.Val_Preds, main=paste(metric.name,"Leave one watershed out
 at a time cross-validation"),
col=color, pch=19,
ylim=range(data$metric,na.rm=T))

color.leg =rich.colors(length(levels(data.saved$WatershedName)))
color.leg[levels(data.saved$WatershedName) == "Upper Grande Ronde"] = "dark gray"
legend("topleft", leg.text, pch=19, ,cex=.7,bg="white",col=color.leg)

data$metric-WS.Val_Preds


use = (is.na(data$metric)==F)&(is.na(WS.Val_Preds)==F)
WS.Val_Preds[use]

data$metric-WS.Val_Preds
a =sum( (data$metric[use] - WS.Val_Preds[use])^2)
b = sum((data$metric[use] - mean(data$metric[use], na.rm=T))^2)
a
b
print(paste("Leave-One-Watershed-Out Rsquared = ",round(1-a/b,2)))

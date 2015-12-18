library(ggplot2)
library(spsurvey)
windows(record=T)

data = read.csv("MetricVisitInformation.csv", header=T)


levels(data$WatershedName)

#GAA = read.csv("GAA_Sites_PhilSubset_20150708.csv", header=T)
#GAA = read.csv("GAA_Site_WithMetrics_20150810.csv", header=T)
#GAA = read.dbf("IC_PriorityWatersheds_GAAMetrics_20150924/IC_PriorityWatersheds_GAAMetrics_20150924.dbf")

### File from Phil 12/9/2015
GAA = read.csv("PERMANOVA_Sites_20150708_original.csv", header=T)


GAA_sheds = read.csv("GAA_Sites_20151006.csv", header=T)

nrow(GAA)
nrow(GAA_sheds)
idx = match(GAA$Site_ID, GAA_sheds$Site_ID)

GAA$EP_SHED = GAA_sheds$EP_SHED[idx]
GAA$EP_SHED = factor(GAA$EP_SHED)

levels(GAA$EP_SHED)
#identify which sites in GAA dataset are champsites

data$SiteName
GAA$CHaMPSite = GAA$Site_ID %in% data$SiteName

#GAA_CHaMP = GAA[GAA$CHaMPSite == "TRUE",]

names(GAA)

GAA$Watershed=as.character(GAA$CHaMPsheds)
GAA$Watershed[is.na(GAA$Watershed == T)] = ""
GAA$EP_SHED = as.character(GAA$EP_SHED)

i=1
for (i in 1:nrow(GAA)) {
if (GAA$Watershed[i]==" ") {GAA$Watershed[i]= GAA$EP_SHED[i]}
if (GAA$Watershed[i]=="") {GAA$Watershed[i]= GAA$EP_SHED[i]}
}

GAA$Watershed = factor(GAA$Watershed)


metric.name = "StrmPwr"
metric.name = "Precip"
metric.name = "MAVELV"
metric.name = "TRange"
metric.name ="CumDran_v1"
metric.name ="Winter95"
metric.name ="Sinuosity"
metric.name ="Channelflo"
metric.name = "WIDE_WW"


metric.list = 
c(
#"NatClass", # Categorical
"Channelflo",
#"NLCD_82",
#"WIDE_WW",
"BFW_M",
"SO_v1",
#"CHINRATE",
#"CURRSUSH",
#"NLCD_11",
#"Conf",
"CumDran_v1",
"StrmPwr",
"WIDE_BF",
"MAFLOWU",
"NewErode",
"Trange",
"GDD",
"DistPrin1"
)

##HSI_St_Juv_WUA_per_m Continuous Variables

#metric.list = c("Channelflo","SO_v1",
#"MAXELEVSMO",
#"Winter2yr",
#"BFW_M")

metric.list %in% names(GAA)

##SubD50 Continuous variables
#metric.list = c("StrmPwr","NatPrin1","Trange","MeanU","WIDE_WW")

metric.list = c("StrmPwr", "Channelflo", MeanU_v1")
names(GAA)
metric.list %in% names(GAA)

Coverage_Limits = data.frame(metric.list, LL= rep(NA, length(metric.list)), UL = rep(NA, length(metric.list)))

metric.name = metric.list[3]
for (metric.name in metric.list) {
print(metric.name)
GAA$metric = GAA[,names(GAA)==metric.name]
GAA$metric[GAA$metric < -999] = NA
min(GAA$metric,na.rm=T)
GAA$metric



  # density plot comparing CHaMP and non-CHaMP sites
jpg.name = paste("Coverage_for_", metric.name,".jpg",sep="")
jpeg(jpg.name, 12,10, units='in', res=600)


length(GAA$metric[is.na(GAA$metric)==T])

theme_set(theme_gray(base_size = 18))
theme(text = element_text(size=20))

p1=ggplot(GAA, aes(x = metric, y = ..density.., color = CHaMPSite, fill = CHaMPSite)) +
    geom_density(alpha = 0.3) +
    theme(legend.position = 'bottom') +
    labs(
     # title = paste("CHaMP Coverages for", metric.name), 
  color = 'In CHaMP?', fill = 'In CHaMP?')
print(p1)
dev.off()


bp.name = paste("Coverage_for_", metric.name,"_boxplot.jpg",sep="")
jpeg(bp.name, 12,10, units='in', res=600)


GAA$metric
GAA$Watershed
par(mar=c(12,7,2,2))
boxplot((GAA$metric[GAA$CHaMPSite=="TRUE"]), #(GAA$metric) ,
(GAA$metric[GAA$Watershed == "Entiat"]),
(GAA$metric[GAA$Watershed == "Upper Grande Ronde"]),
(GAA$metric[GAA$Watershed == "Tucannon"]),
(GAA$metric[GAA$Watershed == "South Fork Salmon"]),
(GAA$metric[GAA$Watershed == "Yankee Fork"]),
(GAA$metric[GAA$Watershed == "South Fork Clearwater"]),
(GAA$metric[GAA$Watershed == "Lochsa"]),
(GAA$metric[GAA$Watershed == "Lolo"]),
(GAA$metric[GAA$Watershed == "Lower Clearwater"]),
(GAA$metric[GAA$Watershed == "Okanogan"]),
(GAA$metric[GAA$Watershed == "Imnaha"]),
(GAA$metric[GAA$Watershed == "Toppenish"]),
(GAA$metric[GAA$Watershed == "Walla Walla"]),
(GAA$metric[GAA$Watershed == "John Day"]),
(GAA$metric[GAA$Watershed == "Methow"]),
(GAA$metric[GAA$Watershed == "Wenatchee"]),
(GAA$metric[GAA$Watershed == "Lemhi"]),
(GAA$metric[GAA$Watershed == "Upper Salmon River Tributaries Above Redfish Lake"]),
names=c("CHaMP Sites",# "GAA File",
       "Entiat", "UGR", "Tucannon",
      "SF Salmon", "Yankee Fork", "SF Clearwater","Lochsa", "Lolo",
       "Lower Clearwater","Okanogan", "Imnaha", "Toppenish", "Walla Walla","John Day",
       "Methow", "Wenatchee", "Lemhi", "US Tribs above RF"),
#main=paste("Coverage for ",metric.name,sep=""),
col=c("magenta", 
#"blue", 
rep("light blue",18),
ylab=metric.name),
las=2, cex.axis=1.5, cex.main=2)



limits=quantile(GAA$metric[GAA$CHaMPSite=="TRUE"],c(.025, .975), na.rm=T)
lines(c(0,100), c(limits[1],limits[1]), lt=2)
lines(c(0,100), c(limits[2],limits[2]), lt=2)

dev.off()
Coverage_Limits[match(metric.name, metric.list),c(2,3)] = limits

col = rep("gray", nrow(GAA))
col[GAA$metric < limits[1]] = "blue"
col[GAA$metric > limits[2]] = "red"


map.name = paste("Coverage_for_", metric.name,"_map.jpg",sep="")
jpeg(map.name, 12,10, units='in', res=600)


plot(GAA$X_ALBERS, GAA$Y_ALBERS, col=col,pch=19, cex=.1,
#main=paste("Coverage for ",metric.name,sep=""),
xlab = "Easting", ylab="Northing", cex.axis= 1.2, cex.lab = 1.5)

legend("topright", c("Within CHaMP Range", "Below CHaMP Range","Above CHaMP Range"),
col=c("gray", "blue","red"), lt=1, cex=1.5)
dev.off()


#cols=heat.colors(10)
#cols=terrain.colors(10)
cols=topo.colors(10)
col.idx = 1+round(9*(GAA$metric-limits[1])/(limits[2]-limits[1]))
col.idx[col.idx < 1]= 1
col.idx[col.idx > 10] = 10
color= cols[col.idx]

col[GAA$metric < limits[1]] = "gray"
col[GAA$metric > limits[2]] = "gray"


Dmap.name = paste("Distribution_for_", metric.name,"_map.jpg",sep="")
jpeg(Dmap.name, 12,10, units='in', res=600)

#plot(GAA$X_ALBERS,GAA$Y_ALBERS)

plot(GAA$X_ALBERS, GAA$Y_ALBERS, col=color,pch=19, cex=.1,
#main=paste("Distribution of ",metric.name,sep=""),
xlab = "Easting", ylab="Northing")

legend("topright", legend= round(seq(limits[1], limits[2], 
by=(limits[2]-limits[1])/9),2),col=cols,pch=19,
cex=.8)

dev.off()

}


write.csv(Coverage_Limits, "coverage_limits.csv")





#################################################3
# Categorical Variable Coverage Maps

#cat.metric.list = c("DistClass","ChanlType","NatClass", "ValleyClas")
# SubD50 categorical variable list
cat.metric.list = c("ChanlType")

for (cat.metric.name in cat.metric.list) {
print(cat.metric.name)
GAA$cat.metric= GAA[,match(cat.metric.name,colnames(GAA))]
GAA$cat.metric


cat.metric.levels = levels(factor(GAA$cat.metric[GAA$CHaMPSite]))
counts.cat.metric.by.level = rep(0, length(cat.metric.levels))
k=1
nrow(GAA[GAA$CHaMPSite,])

for (k in 1:length(cat.metric.levels)){
counts.cat.metric.by.level[k] = sum(rep(1,nrow(GAA[GAA$CHaMPSite,]))[GAA[GAA$CHaMPSite,]$cat.metric == cat.metric.levels[k]])
}

data.frame(cat.metric.levels,counts.cat.metric.by.level)
covered.cat.metrics = factor(cat.metric.levels[counts.cat.metric.by.level > 4])
covered.cat.metrics

# set plot colors, gray if we have sites at a given categorical variable
# level, red if we don't
color = rep("gray", nrow(GAA))
color[(GAA$cat.metric %in% covered.cat.metrics)==FALSE] = "red"

map.name = paste("Coverage_for_", cat.metric.name,"_map.jpg",sep="")
jpeg(map.name, 12,10, units='in', res=600)


plot(GAA$X_ALBERS, GAA$Y_ALBERS, col=color,pch=19, cex=.1,
#main=paste("Distribution of Catogorical Variable ",cat.metric.name,sep=""),
xlab = "Easting", ylab="Northing", cex.axis= 1.2, cex.lab = 1.5)

legend("topright", c("Covered by CHaMP","Covereage <5 CHaMP Sites"),
col=c("gray", "red"), lt=1, cex=1.5)

dev.off()

}

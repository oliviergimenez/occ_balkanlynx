#-- load packages
library(rgdal)  
library(RColorBrewer)  
library(sp)
library(unmarked)
library(MuMIn)
library(rworldmap)

#-- load lynx det/non-det and interviewees' level of expertise
load('lynx.RData')
lynx_obs
dim(lynx_obs) # 392 grid cells, 1375 interviewees
# note: there are more cells than in Fig 1 of the paper due the extra cells 
#       that have been added with the area of knowledge
level_expertise
cell_labels

#-- read in grid 10km2 x 10km2 and site covariates
grid <- readOGR(".", "Grid_10km_Covariates")
head(grid@data)
print(proj4string(grid))

# get forest cover, sampling effort (number of interviewees), altitude and rugosity
# for each cell of the surveyed grid 
forest_cover = NULL
samp_eff = NULL
alti = NULL
rug = NULL
for (i in cell_labels){
	mask = (grid@data$Label == i)
	forest_cover = c(forest_cover,grid@data$Forest_Cov[mask])
	samp_eff = c(samp_eff,grid@data$NumAnswer[mask])
	alti = c(alti,grid@data$Altitude[mask])
	rug = c(rug,grid@data$Ruggedness[mask])
}

# standardize and store mean/sd
mean_fc = mean(forest_cover)
sd_fc = sd(forest_cover)
forest_cover = (forest_cover-mean_fc)/sd_fc

mean_alti = mean(alti)
sd_alti = sd(alti)
alti = (alti-mean_alti)/sd_alti

mean_rug = mean(rug)
sd_rug = sd(rug)
rug = (rug-mean_rug)/sd_rug

samp_eff = as.numeric(scale(samp_eff))
level_expertise = as.numeric(scale(level_expertise))

#-- fit occupancy models

# prepare data for occupancy analysis
lynx_occ = unmarkedFrameOccu(lynx_obs)
head(lynx_occ)
siteCovs(lynx_occ) <- data.frame(forest_cover = forest_cover, samp_eff = samp_eff, rug= rug,alti=alti)
obsCovs(lynx_occ) <- data.frame(quality_estim = level_expertise)

numSites(lynx_occ)
obsNum(lynx_occ)

# general model has:
# linear effect of forest cover and ruggedness, quadratic effect of altitude on occupancy
# linear effect of forest cover, ruggedness, sampling effort and level of expertise on detection
fm <- occu(~ samp_eff + quality_estim + forest_cover + rug ~ forest_cover + rug + alti + I(alti^2), lynx_occ)

# fit all possible combinations
# TAKES A WHILE: 
# dd <- dredge(fm)
load('dd.RData')

# drop models in which the quadratic term is included without the linear term
ind = (is.na(dd[,'psi(alti)']) + !is.na(dd[,'psi(I(alti^2))']))
mask = (ind==2)
dd2 = dd[!mask,]

# display first ten top-ranked models
dd2[1:10,]

# nb of models we run
length(rownames(dd2)) 

# coeff estimates
round(coef(dd2),2)

# model-averaging
res = model.avg(dd2, subset = delta < 2)

# model-averaged param estimates and ses
summary(res)

# confidence intervals
confint(res)

#-- predict occupancy everywhere
#   and get lower/upper bounds of 95% conf int
occup_pred = NULL
occup_lower = NULL
occup_upper = NULL
for (i in 1:nrow(grid@data)){
	l_psi = 1.508555 -0.8573843 * ((grid@data$Forest_Cov[i]-mean_fc)/sd_fc) + 1.205667 * ((grid@data$Ruggedness[i]-mean_rug)/sd_rug) -0.3388781 * ((grid@data$Altitude[i]-mean_alti)/sd_alti) -0.1818641 * ((grid@data$Altitude[i]-mean_alti)/sd_alti)^2
	se_lpsi = sqrt(0.26336^2 + ((grid@data$Forest_Cov[i]-mean_fc)/sd_fc)^2 * 0.29626^2 + ((grid@data$Ruggedness[i]-mean_rug)/sd_rug)^2 * 0.38331^2 + ((grid@data$Altitude[i]-mean_alti)/sd_alti)^2 * 0.27161^2 + (((grid@data$Altitude[i]-mean_alti)/sd_alti)^2)^2 * 0.10435^2)
	temp = 1/(1+exp(-l_psi))
	templ = 1/(1+exp(-(l_psi - 1.96 * se_lpsi)))
	tempu = 1/(1+exp(-(l_psi + 1.96 * se_lpsi)))
	occup_pred = c(occup_pred,temp)
	occup_lower = c(occup_lower,templ)
	occup_upper = c(occup_upper,tempu)
}

# show estimated occ with conf int
cbind(occup_lower, occup_pred, occup_upper)

# add predicted occ to shp file
grid@data$occ_new = occup_pred
grid@data$occ_lb = occup_lower
grid@data$occ_ub = occup_upper

#-- mapping 
#   (http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS3_MakingMaps_part1_mappingVectorData.html)
N = 9
at_regular = seq(0,1,length=N+1)
my.palette <- brewer.pal(n = N, name = "OrRd")
p1 = spplot(grid, c("occ_lb"), col.regions = my.palette, at = at_regular, main = c("CI lower limit"), col = "transparent")
p2 = spplot(grid, c("occ_new"), col.regions = my.palette, at = at_regular, main = c('estim occ'), col = "transparent")
p3 = spplot(grid, c("occ_ub"), col.regions = my.palette, at = at_regular, main = c("CI upper limit"), col = "transparent")
print(p1, split=c(1,1,3,1), more=TRUE)
print(p2, split=c(2,1,3,1), more=TRUE)
print(p3, split=c(3,1,3,1))

# note: to apply the Jenks method like in the paper
# check out https://stackoverflow.com/questions/5304057/partition-into-classes-jenks-vs-kmeans

#-- to get Fig 1 (or so), go through the following steps

# first, get map of surveyed countries (http://egallic.fr/european-map-using-r/)
worldMap <- getMap()
# worldMap@data$NAME
# proj4string(worldMap)

# countries involved
countries <- c("Greece","Albania","Serbia","Kosovo","Bosnia and Herz.","Montenegro","Macedonia")
# select only that countries
ind <- which(worldMap$NAME%in%countries)

balkans = worldMap[ind,]
# convert longlat -> utm
balkans = spTransform(balkans,CRS("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# proj4string(balkans)

p1 = spplot(grid[grid$filter==1,], c("occ_lb"), col.regions = my.palette, at = at_regular, main = c("CI lower limit"), col = "transparent",sp.layout = balkans)
p2 = spplot(grid[grid$filter==1,], c("occ_new"), col.regions = my.palette, at = at_regular, main = c('estim occ'), col = "transparent",sp.layout = balkans)
p3 = spplot(grid[grid$filter==1,], c("occ_ub"), col.regions = my.palette, at = at_regular, main = c("CI upper limit"), col = "transparent",sp.layout = balkans)
print(p1, split=c(1,1,3,1), more=TRUE)
print(p2, split=c(2,1,3,1), more=TRUE)
print(p3, split=c(3,1,3,1))


## 3D GRP model for Acton, Burr Oak, and Pleasant Hill
## predW should also be adjusted based on age and predator being modeled (Line 363)

##### Load Libraries #####
library(ggplot2)
library(gstat)
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(PBSmapping)
library(maps)
library(data.table)
library(grid)
library(tcltk)

rm(list = ls())                             ## clear environment

#####  Menu and read files code  #####
## More reservoirs can be added to this list - make sure folder names are the same
lake <- tk_select.list(c("Acton", "Burr Oak", "Pleasant Hill"), title = "Pick a lake:")

## Read in shoreline polygons saved in reservoir data files folder
shoreline <- readOGR(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/",
                                 lake, "/", lake, "_lake.shp", sep = ""))

## Read in depth data
depthdata.rg <- readOGR(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/",
                             lake, sep = ""), paste(lake, "_depths", sep = ""))

## Read in lake acoustic data
allHAD <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/", lake,
                         "/", lake, "_acoustics.csv", sep = ""))
allHAD$Date <- as.Date(allHAD$Date, "%m/%d/%Y")  ## changes date format


##  Read in water quality data
allWQ <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/", lake,
                     "/", lake, "_WQ.csv", sep = ""))
allWQ$Trip_Date <- as.Date(allWQ$Trip_Date, "%m/%d/%Y")  ## changes date format

##  Read in trawl data
alltrawl <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/", lake,
                           "/", lake, "_trawl.csv", sep = ""))
alltrawl$Trip_Date <- as.Date(alltrawl$Trip_Date, "%m/%d/%Y")  ## changes date format

#Read in transect data
xsectline <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/",
                             lake, "/", lake, "_xsect.csv", sep = ""))

## Input date for simulation and create 5-d window to include sampling not on exact day
sdates <- unique(allHAD$Date)
dayoichar <- tk_select.list(as.character(sdates), title = "Pick a date:")
dayoi <- as.Date(dayoichar)
day1 <- dayoi - 2
day3 <- dayoi + 2

## Subset data based on the date of interest, omits NA values from water quality data
dHAD <- subset(allHAD, Date == dayoi)
dHAD$Easting <- as.numeric(as.character(dHAD$Easting))
dWQ <- subset(allWQ, Trip_Date >= day1 & Trip_Date <= day3,
                 select = c(UTMEasting, UTMNorthing, ReadDepth, Temp, DO))
dWQDO <- na.omit(dWQ)
dtrawl <- subset(alltrawl, Trip_Date >= day1 & Trip_Date <= day3)

## Process trawl data  -- RIGHT NOW TRAWL DATA NOT SPATIALLY EXAMINED
gape <- tk_select.list(c("110.9 (LMB3)", "122.9 (SAE1)", "73.3 (WHC2)"), title = "Select gape size (mm):")
if(gape == "110.9 (LMB3)") yoy <- sum(dtrawl$Length < 110.9) ## sets gape limit for avg size of Age-2 LMB from Michaletz 1997
if(gape == "122.9 (SAE1)") yoy <- sum(dtrawl$Length < 122.9) ## sets gape limit for avg size of Age-3 LMB from Michaletz 1997
if(gape == "73.3 (WHC2)") yoy <- sum(dtrawl$Length < 73.3) ## sets gape limit for avg size of Age-1 SAE from Knight et al. 1984 (based on WAE gape)

## only prey <## mm are considered (use interactive predator model to determine gape of predator being modeled Gaeta et al. 2018)
preyprop <- yoy / nrow(dtrawl)  ## calculates the proportion of prey that are below gape


## Process acoustic data
dHAD$cDensity <- preyprop * dHAD$Density    ## adjusts for proportion YOY
dHAD$cBiomass <- preyprop * dHAD$Biomass

dHAD$cDensity <- dHAD$cDensity / 10000      ## converts #/ha to #/m3
dHAD$cBiomass <- dHAD$cBiomass * 0.1        ## converts kg/ha to g/m3

for(j in 0:max(dHAD$Layer_depth_min)) {     ## fills in vertical profile of abundance, 0-2 m depth filled with side-scan data
    if(!exists("preydata")) {
        preydata <- subset(dHAD, Layer_depth_min == j, select = c(Layer_depth_min, 
                            Easting, Northing, cDensity, cBiomass))
    } else {
        if(j == 1) {
            temp_dataset <- preydata
            temp_dataset$Layer_depth_min = 1
            preydata <- rbind(preydata, temp_dataset)
            rm(temp_dataset)
        } else {
            temp_dataset <- subset(dHAD, Layer_depth_min == j, select = 
                                       c(Layer_depth_min, Easting, Northing, 
                                         cDensity, cBiomass))
            preydata <- rbind(preydata, temp_dataset)
            rm(temp_dataset)
        }
    }
}

preydata$lnBiomass <- log(preydata$cBiomass + 1)    ## compute ln(biomass)
colnames(preydata)[1] <- "Depth"

##  Define min and max values for plots
maxTEMP = 28.1931199003548
medTEMP = 20.6965703694715
minTEMP = 13.2000208385881


maxDO = 13.3120167428226
medDO = 6.6560083714113
minDO = 0

maxPREY = 2.65427105720343
medPREY = 1.32713552860171
minPREY = 0

nxsect = max(dHAD$Interval)
maxzWQ = max(dWQ$ReadDepth)
maxzHAD = max(dHAD$Layer_depth_min)
maxz = max(maxzWQ, maxzHAD)  ## Find maximum sample depths from WQ and HAD

depthsTEMP = paste("mT", 0:maxz, sep = "") ## make character names for areal plots (one for every m depth)
depthsDO = paste("mD", 0:maxz, sep = "")
depthsBIOM = paste("mB", 0:maxz, sep = "")


##### Reads bioenergetic database #####
bioparam <- read.csv("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Reservoir data files/Bioenergetics.csv")


##### Subset and process depth data  #####
depthdata <- as.data.frame(depthdata.rg)
depthd <- as.data.frame(matrix(0, ncol = 3, nrow = nrow(depthdata)))
colnames(depthd) <- c("Easting", "Northing", "Depth")
depthd$Easting <- depthdata$Easting
depthd$Northing <- depthdata$Northing
if(lake != "Findlay1" & lake != "Findlay2" & lake != "Pleasant Hill") depthd$Depth <- -depthdata$Depth_m
if(lake == "Findlay1" | lake == "Findlay2" | lake == "Pleasant Hill") depthd$Depth <- -depthdata$depth_m

## determine if value = 0 and subset rows not equal to all 0's
row_sub = apply(depthd, 1, function(row) all(row !=0 ))
depthd <- depthd[row_sub,]

## subset for really large datasets then make depths positive, speeds up model, but less detail for depth
if(nrow(depthdata) >= 15000 & nrow(depthdata) < 50000) depthd <- depthd[seq(1, nrow(depthd), 10), ]
if(nrow(depthdata) >= 50000 & nrow(depthdata) < 200000) depthd <- depthd[seq(1, nrow(depthd), 20), ]
if(nrow(depthdata) >= 200000) depthd <- depthd[seq(1, nrow(depthd), 100), ]


##### Process Shoreline #####
## Change format to dataframe and define projection, if needed
if(lake == "Findlay1" | lake == "Findlay2" | lake == "Delaware") {
  proj4string(shoreline) <- CRS("+init=epsg:3728")
  shoreline <- spTransform(shoreline, CRS("+proj=utm +zone=17"))
}

if(lake != "Acton") proj4string(shoreline) <- CRS("+proj=utm +zone=17")
if(lake == "Acton") proj4string(shoreline) <- CRS("+proj=utm +zone=16") 
shoreoutline <-fortify(shoreline)


## Create zero depth data (i.e., shoreline) and add to depth data
zerodepth <- shoreoutline[,1:2]
zerodepth$Depth <- 0
colnames(zerodepth) <- c("Easting", "Northing", "Depth")
depthd <- rbind(depthd, zerodepth)
depthd3D <- depthd


## Set spatial coordinates to create a Spatial object
coordinates(depthd) <- ~Easting + Northing
if(lake != "Acton") proj4string(depthd) <- CRS("+proj=utm +zone=17")
if(lake == "Acton") proj4string(depthd) <- CRS("+proj=utm +zone=16")

## generates polygons from shoreline data, for making the areal maps and determining which cells are not land
if(lake == "Tappan") {  
    poly <- subset(shoreoutline, shoreoutline$group == 7.1, select = c(long, lat))
} else {
    poly <- subset(shoreoutline, shoreoutline$piece == 1, select = c(long, lat))
}
polyP <- Polygon(cbind(poly$long, poly$lat)) ## Add lat & lon to object
polySP <- SpatialPolygons(list(Polygons(list(polyP), ID=1)))


## Define grid extent for areal graphing and interpolating (grid is larger than reservoir, but land subtracted later)
xmax = round(max(depthd$Easting) + 100, -2)
xmin = round(min(depthd$Easting) - 100, -2)
ymax = round(max(depthd$Northing) + 100, -2)
ymin = round(min(depthd$Northing) - 100, -2)
x.range <- as.numeric(c(xmin, xmax))
y.range <- as.numeric(c(ymin, ymax))

## Create matrix of outside grid bounds and change to SpatialPolygons Dataframe
outside <- matrix(c(xmin,ymin, xmax,ymin, xmax,ymax, xmin,ymax, xmin,ymin),
                  ncol = 2, byrow = TRUE)
outside <- SpatialPolygons(list(Polygons(list(Polygon(outside)), ID = 1)))
outside <- SpatialPolygonsDataFrame(Sr=outside, data=data.frame(polyvar=357))

## Make an overlay grid for range from min to max by units of 50
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 50),
                   y = seq(from = y.range[1], to = y.range[2], by = 50)) 
coordinates(grd) <- ~x + y
if(lake != "Acton") proj4string(grd) <- CRS("+proj=utm +zone=17")
if(lake == "Acton") proj4string(grd) <- CRS("+proj=utm +zone=16")
gridded(grd) <- TRUE

## Define grid for transect data
grdxsect <- expand.grid(x = seq(from = 1, to = nrow(xsectline), by = 1), 
                        y = seq(from = 0, to = maxz, by = 1))

##### Interpolate depth using inverse distance weighted model  #####
idwd <- idw(formula = Depth ~ 1, locations = depthd, newdata = grd)
idwd.output = as.data.frame(idwd)     ## Output defined as dataframe with names
names(idwd.output)[1:3] <- c("Easting", "Northing", "Depth")
idwd.output$Depth <- round(idwd.output$Depth)                                   ## Round depth data

## Generate interpolated (IDW method) temp, DO and biomass data for each depth
for(doi in 0:maxz) {
    ## subsets data based on UTM
    idwTEMP <- subset(dWQ, dWQ$ReadDepth == doi, select = c(UTMEasting, UTMNorthing, Temp))
    idwDO <- subset(dWQDO, dWQDO$ReadDepth == doi, select = c(UTMEasting, UTMNorthing, DO))
    idwHAD <- subset(preydata, preydata$Depth == doi, select = c(Easting, Northing, cDensity))
    idwHAB <- subset(preydata, preydata$Depth == doi, select = c(Easting, Northing, lnBiomass))
    
    if(doi > maxzWQ) {       ## duplicates previous depth layer for deeper depths
        idwTEMP <- subset(dWQ, dWQ$ReadDepth == maxzWQ, select = c(UTMEasting, UTMNorthing, Temp))
        idwDO <- subset(dWQDO, dWQDO$ReadDepth == maxzWQ, select = c(UTMEasting, UTMNorthing, DO))
    }

    if(doi > maxzHAD) {     ## duplicates previous depth layer for deeper depths
        idwHAD <- subset(preydata, preydata$Depth == maxzHAD, select = c(Easting, Northing, cDensity))
        idwHAB <- subset(preydata, preydata$Depth == maxzHAD, select = c(Easting, Northing, lnBiomass))
    }
    
    ## Set spatial coordinates to create Spatial objects and interpolate
    coordinates(idwTEMP) <- ~UTMEasting + UTMNorthing
    if(lake != "Acton") proj4string(idwTEMP) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwTEMP) <- CRS("+proj=utm +zone=16")
    idwt <- idw(formula = Temp ~ 1, locations = idwTEMP, newdata = grd)
    idwt.output = as.data.frame(idwt)   ## Output defined as dataframe with names
    names(idwt.output)[1:3] <- c("Easting", "Northing", "Temp")
    
    coordinates(idwDO) <- ~UTMEasting + UTMNorthing
    if(lake != "Acton") proj4string(idwDO) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwDO) <- CRS("+proj=utm +zone=16")
    idwo <- idw(formula = DO ~ 1, locations = idwDO, newdata = grd)
    idwo.output = as.data.frame(idwo)   ## Output defined as dataframe with names
    names(idwo.output)[1:3] <- c("Easting", "Northing", "DO")
    
    coordinates(idwHAD) <- ~Easting + Northing
    if(lake != "Acton") proj4string(idwHAD) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwHAD) <- CRS("+proj=utm +zone=16")
    idwDEN <- idw(formula = cDensity ~ 1, locations = idwHAD, newdata = grd)
    idwDEN.output = as.data.frame(idwDEN)     ## Output defined as dataframe with names
    names(idwDEN.output)[1:3] <- c("Easting", "Northing", "Density")
    
    coordinates(idwHAB) <- ~Easting + Northing
    if(lake != "Acton") proj4string(idwHAB) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwHAB) <- CRS("+proj=utm +zone=16")
    idwBIO <- idw(formula = lnBiomass ~ 1, locations = idwHAB, newdata = grd)
    idwBIO.output = as.data.frame(idwBIO)     ## Output defined as dataframe with names
    names(idwBIO.output)[1:3] <- c("Easting", "Northing", "lnBiomass")
    
    if(!exists("WQdataset")) {  ## creates dataset with interpolated data
        WQdataset <- idwt.output
        WQdataset$DO <- idwo.output$DO
        WQdataset$Density <- idwDEN.output$Density
        WQdataset$lnBiomass <- idwBIO.output$lnBiomass
        WQdataset$Depth <- rep(doi,nrow(WQdataset))
    } else {                    ## appends data to dataset once original created
        temp_dataset <- idwt.output
        temp_dataset$DO <- idwo.output$DO
        temp_dataset$Density <- idwDEN.output$Density
        temp_dataset$lnBiomass <- idwBIO.output$lnBiomass
        temp_dataset$Depth <- rep(doi,nrow(idwt.output))
        WQdataset <-rbind(WQdataset, temp_dataset)
        rm(temp_dataset)
    }
}

WQdataset <- WQdataset[c(-4)]     ## removes added var1.var column

## Generate transect temp, DO and prey biomass data for each depth
xsectdata <- data.frame(0,0,0,0,0,0,0,0)     ## empty dataframe for xsect data

for(i in 1:nrow(xsectline)) {                ## loop through points on xsect
  easting = as.numeric(xsectline[i, 1])    ## easting point
  northing = as.numeric(xsectline[i, 2])   ## northing point
  
  ## find closest data for WQ
  distance <- sqrt((WQdataset$Easting - easting)^2 + (WQdataset$Northing - northing)^2)
  nearestdata <- WQdataset[which(distance == min(distance)), ]
  
  ## Finds closest bottom depth for ploting
  dnearest <- sqrt((depthd3D[, 1] - easting)^2 + (depthd3D[, 2] - northing)^2)
  maxdepth3D <- round(max(depthd3D[which(dnearest == min(dnearest)), 3]))
  
  ## Subset commands just in case there are multiple nearest data points
  if(length(unique(nearestdata$Easting)) > 1) {
    nearestdata <- nearestdata[which(nearestdata$Easting == min(nearestdata$Easting)), ]
  }
  if(length(unique(nearestdata$Northing)) > 1) {
    nearestdata <- nearestdata[which(nearestdata$Northing == min(nearestdata$Northing)), ]
  }
  
  for(iii in 0:maxz) {                   ## loop through each xsect depth  
    xsectdata[i + nrow(xsectline) * iii, ] <-              ## define data
      cbind(grdxsect[i + nrow(xsectline) * iii, 1:2], nearestdata[iii + 1, 1:6])        
    
    if(iii > maxdepth3D) xsectdata[i + nrow(xsectline) * iii, 5:8] = NA
  }
}
colnames(xsectdata) <- c("x", "y", "Easting", "Northing", "Temp", "DO",
                         "Density", "lnBiomass")
grdxsect <- na.omit(xsectdata)

#####  SHORE EXTRACTION AND PRINTING  #####
## Loops through number of shore polygons and creates map
for (i in 1:nlevels(shoreoutline$piece)) {
    if(lake == "Tappan" && i == 1) {
        poly <- subset(shoreoutline, shoreoutline$group == 7.1, select = c(long, lat))
    } else {
        poly <- subset(shoreoutline, shoreoutline$piece == i, select = c(long, lat))
    }
    polyP <- Polygon(cbind(poly$long, poly$lat)) ## Add lat & lon to object
    polySP <- SpatialPolygons(list(Polygons(list(polyP), ID=i)))
    outside <- gDifference(outside, polySP)
    nam <- paste("polySP", i, sep="")             ## Create sequential names
    assign(nam, polySP)                           ## Assign name
}

## Creates list of ggplot commands to plot all island polygons
islands <- list()
if(nlevels(shoreoutline$piece) >= 2) {
    for (ii in 2:nlevels(shoreoutline$piece)) {
        island <- paste("polySP", ii, sep="")
        island <- get(island)
        islands <- c(islands, geom_polygon(data=island, aes(x=long, y=lat), fill="white", asp = 1))
    }
}

###### Plotting commands  #####
## Generate transect plots
xsplotTEMP <- ggplot() +
  geom_tile(data = grdxsect, aes(x = x, y = y, fill = Temp)) +  ## basic plot
  scale_fill_gradientn(expression(atop("",
                                       atop(textstyle("Temperature"),
                                            atop(textstyle("(°C)"))))),  
                                              colours = rainbow(4), 
                                                breaks = c(15, 21, 27), 
                                                  limits = c(minTEMP, maxTEMP)) +
  theme(legend.text = element_text(size=14), legend.title = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  ## remove gridlines
  scale_x_continuous(expand = c(0,-0.5)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_y_reverse(expand = c(0,0), name = "Depth (m)")  ## Reverse axis on depth scale

xsplotDO <- ggplot() +
  geom_tile(data = grdxsect, aes(x = x, y = y, fill = DO)) +  ## basic plot
  scale_fill_gradientn(expression(atop("",
                                       atop(textstyle("DO"),
                                            atop(textstyle("(mg/L)"))))), 
                                              colours = rainbow(4), 
                                                breaks = c(2, 7, 12), 
                                                  limits = c(minDO, maxDO)) +
  theme(legend.text = element_text(size=14), legend.title = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  ## remove gridlines
  scale_x_continuous(expand = c(0,-0.5)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_y_reverse(expand = c(0,0), name = "Depth (m)")  ## Reverse axis on depth scale

xsplotPREY <- ggplot() +
  geom_tile(data = grdxsect, aes(x = x, y = y, fill = lnBiomass)) +  ## basic plot
  scale_fill_gradientn(expression(atop("",
                                       atop(textstyle("Prey Biomass"),
                                            atop(textstyle("[ln(g/m"^3*")]"))))), 
                                              colours = rainbow(4), 
                                                breaks = c(0, 1, 2), 
                                                  limits = c(minPREY, maxPREY)) +
  theme(legend.text = element_text(size=14), legend.title = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  ## remove gridlines
  scale_x_continuous(expand = c(0,-0.5)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_y_reverse(expand = c(0,0), name = "Depth (m)")  ## Reverse axis on depth scale

##### Generate 2D areal maps for all variables  #####
## Finds only those Easting and Northing coordinates within lake boundary
lakepts <- idwt.output
coordinates(lakepts) <- ~Easting + Northing
if(lake != "Acton") proj4string(lakepts) <- CRS("+proj=utm +zone=17") 
if(lake == "Acton") proj4string(lakepts) <- CRS("+proj=utm +zone=16")
lakepts <- gDifference(lakepts, outside)
lakepts <- as.data.frame(lakepts)

## Goes through each coordinate set in lakepts and finds the depth
for(c in 1:nrow(lakepts)){
    easting <- lakepts[c,1]
    northing <- lakepts[c,2]
    
    ## Finds bottom depth for ploting  (used idwd.oupt, changed to depthd3D)
    dnearest <- sqrt((depthd3D[, 1] - easting)^2 + (depthd3D[, 2] - northing)^2)
    lakepts[c, 3] <- round(max(depthd3D[which(dnearest == min(dnearest)), 3]))
}

## Generate 2D areal maps for temperature at each depth
TEMPlist <- list()  ## define list to hold all temperature by depth plots
DOlist <- list()    ## define list to hold all DO by depth plots
BIOMlist <- list()  ## define list to hold all biomass by depth plots
depthPLOTS <- data.frame()

if(maxz < 10) dlist = c(2, 3, 4, 5, 6, 7, 8, 9)  ## if reservoir shallow, map every 1-m
if(maxz == 10) dlist = c(2, 3, 4, 5, 6, 7, 8, 10)## if reservoir deep, skip some depths
if(maxz > 10) dlist = c(2, 3, 4, 5, 6, 7, 10, 12)

for(i in 1:8) {  ## loop through each depth plot and generate plot
  rm(depthPLOTS)   
  depthPLOTS <- data.frame()
  doi = dlist[i]
  if(doi <= maxz) {
    for(cc in 1:nrow(lakepts)){ 
      easting <- lakepts[cc, 1]
      northing <- lakepts[cc, 2]
      lakedepth <- lakepts[cc, 3]
      
      if(doi <= lakedepth) {
        poi <- subset(WQdataset, Easting == easting & Northing == northing & Depth == doi)
        depthPLOTS <- rbind(depthPLOTS, poi)
      }
    }
    
    layr <- ggplot() +   ## temperature plots
      geom_tile(data=depthPLOTS, aes(x=Easting, y=Northing, fill = Temp)) +
      scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                           limits = c(minTEMP, maxTEMP)) +
      geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
      islands +
      coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
      theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(plot.title = element_text(size=22)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
      ggtitle(paste(doi, " m", sep = ""))
    
    assign(paste0(depthsTEMP[doi+1]), layr)
    TEMPlist <- c(TEMPlist, list(layr))
    
    layr <- ggplot() +   ## DO plots
      geom_tile(data=depthPLOTS, aes(x=Easting, y=Northing, fill = DO)) +
      scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                           limits = c(minDO, maxDO)) +
      geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
      islands +
      coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
      theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(plot.title = element_text(size=22)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
      ggtitle(paste(doi, " m", sep = ""))
    
    assign(paste0(depthsDO[doi+1]), layr)
    DOlist <- c(DOlist, list(layr))
    
    layr <- ggplot() +   ## Prey plots
      geom_tile(data=depthPLOTS, aes(x=Easting, y=Northing, fill = lnBiomass)) +
      scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                           limits = c(minPREY, maxPREY)) +
      geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
      islands +
      coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
      theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(plot.title = element_text(size=22)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
      ggtitle(paste(doi, " m", sep = ""))
    
    assign(paste0(depthsBIOM[doi+1]), layr)
    BIOMlist <- c(BIOMlist, list(layr)) 
  }
}


## Code to print multiple areal temperature figures and one xsect at bottom
## Opens jpeg file and specifies size and resolution
jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/", lake,
                     dayoi, "_Temp.jpeg", sep = ""), width = 2000, 
    height = 2580, units = "px", res = 300, restoreConsole = TRUE)

## Creates grid to print multiple items on a page
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 4, heights = unit(c(0.5, 5, 5, 5, 4.5), "null"))))

## Specifies which depths to print given the depth of the reservoir
if(maxz < 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mT2, vp = vplayout(2, 1))
  if(maxz > 2) print(mT3, vp = vplayout(2, 2))
  if(maxz > 3) print(mT4, vp = vplayout(2, 3))
  if(maxz > 4) print(mT5, vp = vplayout(2, 4))
  if(maxz > 5) print(mT6, vp = vplayout(3, 1))
  if(maxz > 6) print(mT7, vp = vplayout(3, 2))
  if(maxz > 7) print(mT8, vp = vplayout(3, 3))
  if(maxz > 8) print(mT9, vp = vplayout(3, 4))
  print(xsplotTEMP, vp = vplayout(4, 1:4))
}
if(maxz == 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mT2, vp = vplayout(2, 1))
  if(maxz > 2) print(mT3, vp = vplayout(2, 2))
  if(maxz > 3) print(mT4, vp = vplayout(2, 3))
  if(maxz > 4) print(mT5, vp = vplayout(2, 4))
  if(maxz > 5) print(mT6, vp = vplayout(3, 1))
  if(maxz > 6) print(mT7, vp = vplayout(3, 2))
  if(maxz > 7) print(mT8, vp = vplayout(3, 3))
  if(maxz > 9) print(mT10, vp = vplayout(3, 4))
  print(xsplotTEMP, vp = vplayout(4, 1:4))
}
if(maxz > 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mT2, vp = vplayout(2, 1))
  if(maxz > 2) print(mT3, vp = vplayout(2, 2))
  if(maxz > 3) print(mT4, vp = vplayout(2, 3))
  if(maxz > 4) print(mT5, vp = vplayout(2, 4))
  if(maxz > 5) print(mT6, vp = vplayout(3, 1))
  if(maxz > 6) print(mT7, vp = vplayout(3, 2))
  if(maxz > 9) print(mT10, vp = vplayout(3, 3))
  if(maxz > 11) print(mT12, vp = vplayout(3, 4))
  print(xsplotTEMP, vp = vplayout(4, 1:4))
}
dev.off()


## Code to print multiple areal DO figures and one xsect at bottom
jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/", lake,
                     dayoi, "_DO.jpeg", sep = ""), width = 2000, 
    height = 2580, units = "px", res = 300, restoreConsole = TRUE)


vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 4, heights = unit(c(0.5, 5, 5, 5, 4.5), "null"))))

if(maxz < 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mD2, vp = vplayout(2, 1))
  if(maxz > 2) print(mD3, vp = vplayout(2, 2))
  if(maxz > 3) print(mD4, vp = vplayout(2, 3))
  if(maxz > 4) print(mD5, vp = vplayout(2, 4))
  if(maxz > 5) print(mD6, vp = vplayout(3, 1))
  if(maxz > 6) print(mD7, vp = vplayout(3, 2))
  if(maxz > 7) print(mD8, vp = vplayout(3, 3))
  if(maxz > 8) print(mD9, vp = vplayout(3, 4))
  print(xsplotDO, vp = vplayout(4, 1:4))
}
if(maxz == 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mD2, vp = vplayout(2, 1))
  if(maxz > 2) print(mD3, vp = vplayout(2, 2))
  if(maxz > 3) print(mD4, vp = vplayout(2, 3))
  if(maxz > 4) print(mD5, vp = vplayout(2, 4))
  if(maxz > 5) print(mD6, vp = vplayout(3, 1))
  if(maxz > 6) print(mD7, vp = vplayout(3, 2))
  if(maxz > 7) print(mD8, vp = vplayout(3, 3))
  if(maxz > 9) print(mD10, vp = vplayout(3, 4))
  print(xsplotDO, vp = vplayout(4, 1:4))
}
if(maxz > 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mD2, vp = vplayout(2, 1))
  if(maxz > 2) print(mD3, vp = vplayout(2, 2))
  if(maxz > 3) print(mD4, vp = vplayout(2, 3))
  if(maxz > 4) print(mD5, vp = vplayout(2, 4))
  if(maxz > 5) print(mD6, vp = vplayout(3, 1))
  if(maxz > 6) print(mD7, vp = vplayout(3, 2))
  if(maxz > 9) print(mD10, vp = vplayout(3, 3))
  if(maxz > 11) print(mD12, vp = vplayout(3, 4))
  print(xsplotDO, vp = vplayout(4, 1:4))
}
dev.off()


## Code to print multiple areal biomass figures and one xsect at bottom
jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/", lake,
                     dayoi, "_BIOM_.jpeg", sep = ""), width = 2000, 
    height = 2580, units = "px", res = 300, restoreConsole = TRUE)

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 4, heights = unit(c(0.5, 5, 5, 5, 4.5), "null"))))

if(maxz < 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mB2, vp = vplayout(2, 1))
  if(maxz > 2) print(mB3, vp = vplayout(2, 2))
  if(maxz > 3) print(mB4, vp = vplayout(2, 3))
  if(maxz > 4) print(mB5, vp = vplayout(2, 4))
  if(maxz > 5) print(mB6, vp = vplayout(3, 1))
  if(maxz > 6) print(mB7, vp = vplayout(3, 2))
  if(maxz > 7) print(mB8, vp = vplayout(3, 3))
  if(maxz > 8) print(mB9, vp = vplayout(3, 4))
  print(xsplotPREY, vp = vplayout(4, 1:4))
}
if(maxz == 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mB2, vp = vplayout(2, 1))
  if(maxz > 2) print(mB3, vp = vplayout(2, 2))
  if(maxz > 3) print(mB4, vp = vplayout(2, 3))
  if(maxz > 4) print(mB5, vp = vplayout(2, 4))
  if(maxz > 5) print(mB6, vp = vplayout(3, 1))
  if(maxz > 6) print(mB7, vp = vplayout(3, 2))
  if(maxz > 7) print(mB8, vp = vplayout(3, 3))
  if(maxz > 9) print(mB10, vp = vplayout(3, 4))
  print(xsplotPREY, vp = vplayout(4, 1:4))
}
if(maxz > 10) {
  grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
  if(maxz > 1) print(mB2, vp = vplayout(2, 1))
  if(maxz > 2) print(mB3, vp = vplayout(2, 2))
  if(maxz > 3) print(mB4, vp = vplayout(2, 3))
  if(maxz > 4) print(mB5, vp = vplayout(2, 4))
  if(maxz > 5) print(mB6, vp = vplayout(3, 1))
  if(maxz > 6) print(mB7, vp = vplayout(3, 2))
  if(maxz > 9) print(mB10, vp = vplayout(3, 3))
  if(maxz > 11) print(mB12, vp = vplayout(3, 4))
  print(xsplotPREY, vp = vplayout(4, 1:4))
}
dev.off()

#####  BIOENERGETIC MODEL ######
soichar <- unique(bioparam$Species)
soi <- tk_select.list(as.character(soichar), title = "Pick a Species:")
{
  ## Read in parameters from bioenergetic database (user defined params at end)
  roi <- which(bioparam$Species == soi)
  CEQ = bioparam$CEQ[roi]
  CA = bioparam$CA[roi]
  CB = bioparam$CB[roi]
  CQ = bioparam$CQ[roi]
  CTO = bioparam$CTO[roi]
  CTM = bioparam$CTM[roi]
  CTL = bioparam$CTL[roi]
  CKone = bioparam$CK1[roi]
  CKfour = bioparam$CK4[roi]
  
  REQ = bioparam$REQ[roi]
  RA = bioparam$RA[roi]
  RB = bioparam$RB[roi]
  RQ = bioparam$RQ[roi]
  RTO = bioparam$RTO[roi]
  RTM = bioparam$RTM[roi]
  RTL = bioparam$RTL[roi]
  RKone = bioparam$RK1[roi]
  RKfour = bioparam$RK4[roi]
  ACT = bioparam$ACT[roi]
  BACT = bioparam$BACT[roi]
  SDA = bioparam$SDA[roi]
  
  EGEXEQ = bioparam$EGEXEQ[roi]
  FA = bioparam$FA[roi]
  FB = bioparam$FB[roi]
  FG = bioparam$FG[roi]
  UA = bioparam$UA[roi]
  UB = bioparam$UB[roi]
  UG = bioparam$UG[roi]
  
  predW = 1000
  if(soi == "LMB3") predW = 312  #should be changed based on age of fish
  if(soi == "SAE1") predW = 593  #should be changed based on age of fish
  if(soi == "WHC2") predW = 118  #should be changed based on age of fish
  
  predED = bioparam$ED[roi]
  
  gapeonlyprey <- subset(dtrawl, Length < 100)  ## subsets prey based on gape limit of predator being modeled
  if(soi == "LMB3") gapeonlyprey <- subset(dtrawl, Length <= 110.9) # should be altered based on age (currently age-3)
  if(soi == "SAE1") gapeonlyprey <- subset(dtrawl, Length <= 122.9) # should be altered based on age (currently age-1)
  if(soi == "WHC2") gapeonlyprey <- subset(dtrawl, Length <= 73.3)# should be altered based on age (currently age-2)
  avglengthop <- mean(gapeonlyprey$Length)  ## average length of prey items from trawl
  lnshadW = 2.7875 * log(avglengthop) - 10.5461     ## General TL-W equation for Gizzard shad from Denlinger et al. 2006
  if(lake == "Acton") lnshadW = 2.25 * log(avglengthop) - 7.628 ##TL-W prey equation for Acton
  if(lake == "Burr Oak") lnshadW = 2.741 * log(avglengthop) - 10.29 ##TL-W prey equation for Burr Oak
  if(lake == "Pleasant Hill") lnshadW = 3.107 * log(avglengthop) - 11.78 ##TL-W prey equation for Pleasant Hill
  shadW=exp(lnshadW)
  
  ## Determines number of each prey species in trawl catch for weighted energy density calc  
  wGShad = sum(gapeonlyprey$Species == 20003) # number of gizzard shad
  wGShin = sum(gapeonlyprey$Species == 43003) # number of golden shiners
  wGnCat = sum(gapeonlyprey$Species == 47000) # number of general catfish
  wChCat = sum(gapeonlyprey$Species == 47002) # number of channel catfish
  wBSilv = sum(gapeonlyprey$Species == 70001) # number of brook silverside
  wWBass = sum(gapeonlyprey$Species == 74001) # number of white bass
  wHBass = sum(gapeonlyprey$Species == 74005) # number of hybrid striped bass
  wWCrap = sum(gapeonlyprey$Species == 77001) # number of white crappie
  wBCrap = sum(gapeonlyprey$Species == 77002) # number of black crappie
  wWarMt = sum(gapeonlyprey$Species == 77007) # number of warmouth
  wBlGl = sum(gapeonlyprey$Species == 77009) # number of blue gill
  wRdSn = sum(gapeonlyprey$Species == 77012) # number of redear sunfish
  wSunG = sum(gapeonlyprey$Species == 77994) # number of general sunfish
  wCraG = sum(gapeonlyprey$Species == 77996) # number of general crappie
  wLgPr = sum(gapeonlyprey$Species == 80011) # number of logperch
  wSndG = sum(gapeonlyprey$Species == 80901) # number of general sander
  wNoTa = sum(gapeonlyprey$Species == 99998) # number of non target fish
  wUnKn = sum(gapeonlyprey$Species == 99999) # number of unidentified
  
  ## Energy densities of Ohio reservoir prey items found in AC, BO, and PH trawls 2007-2016 
  ## (most values pulled from Table in Denlinger 2006)    
  wGShadED = 5103 # Miranda & Muncy 1989
  wGShinED = 4994 # Bryan et al. 1996 (value for common shiner)
  wGnCatED = 5437 # Masser et al. 1991 (value for catfishes)
  wChCatED = 5437 # Eggleton & Schramm 2002
  wBSilvED = 4392 # Pope et al. 1996
  wWBassED = 4392 # Bryan et al. 1996
  wHBassED = 5935 # average juvenile HSB from Marcek et al. 2020
  wWCrapED = 4184 # Zweifel 2000; Bajer et al. 2004
  wBCrapED = 4184 # Zweifel 2000; Bajer et al. 2004 (value for w. crappie)
  wWarMtED = 4852 # Miranda & Muncy 1989 (value for sunfishes)
  wBlGlED = 4186 # Kitchell et al. 1974
  wSunGED = 4852 # Miranda & Muncy 1989 (value for sunfishes)
  wCraGED = 4184 # Zweifel 2000; Bajer et al. 2004
  wRdSnED = 4852 # Miranda & Muncy 1989 (value for sunfishes)
  wLgPrED = 4392 # Bryan et al. 1996
  wSndGED = 4186 # zweifel et al. 2010 (value for saugeye)
  wNoTaED = 4391 # Other fish value from Denlinger 2006
  wUnKnED = 4391 # Other fish value from Denlinger 2006
  
  denom = nrow(gapeonlyprey) #calculates total number of fish in trawl
  
  #weighted energy density calc
  preyED =  (wGShad*wGShadED + wGShin*wGShinED + wChCat*wChCatED + wGnCat*wGnCatED + 
               wBSilv*wBSilvED + wWBass*wWBassED + wWCrap*wWCrapED + wBCrap*wBCrapED + 
               wWarMt*wWarMtED + wHBass*wHBassED + wBlGl*wBlGlED + wSunG*wSunGED + 
               wCraG*wCraGED + wRdSn*wRdSnED + wLgPr*wLgPrED + wSndG*wSndGED + 
               wNoTa*wNoTaED + wUnKn*wUnKnED)/denom
  
  GRPd <- data.frame()
  rownum = 0
  
  ## Generate 3D GRP data for each depth
  
  for(ccc in 1:nrow(lakepts)) {           ## loop through every cell in 3D array
    easting <- lakepts[ccc,1]
    northing <- lakepts[ccc,2]
    nearest <- sqrt((WQdataset$Easting - easting)^2 + (WQdataset$Northing - northing)^2)
    nearestdata <- WQdataset[which(nearest == min(nearest)), ]
    
    neard <- subset(lakepts, lakepts$x == easting & lakepts$y == northing)
    depthmax <- max(neard[, 3])
    
    for(z in 0:depthmax) {
      rownum = rownum +1
      
      roi <- nearestdata[which(nearestdata$Depth == z), ]
      if(z > maxz)  {
        dc = z - maxz
        roi <- nearestdata[which(nearestdata$Depth == z-dc), ]
      }
      
      wt = roi$Temp
      oxy = roi$DO
      preyd = roi$Density
      preyb = exp(roi$lnBiomass) - 1
      
      if(soi == "SMB") {         ## extra for SMB from Whitledge et al 2003
        if(wt > 22) CQ = 1.95
      }
      
      ##  Consumption Equations
      if(CEQ == 2) {     ## Consumption equation 2
        Y = log(CQ) * (CTM - CTO + 2)
        Z = log(CQ) * (CTM - CTO)
        X = (Z^2 * (1 + (1 + 40 / Y)^0.5)^2) / 400
        V = (CTM - wt) / (CTM - CTO)
        
        ftc = V^X * exp(X * (1 - V))
        if(V^X == "NaN") ftc = 0.00001  ## If water temp>CTM, use low ftc
        
        Cmax = CA * predW^CB        ## this is g prey/g pred
        
        ##  Species-specific oxygen and temperature functions
        if(soi == "LMB3") {  
          fdo = 1 / (1 + exp(-1 * (oxy - 2.281586) / 0.345078)) ## Type III Stewart et al 1967
        }
        if(soi == "MUS") {
          fdo = 1 / (1 + exp(-1 * (oxy - 2.563859) / 0.37121))  ## Type III Adelman and Smith 1970
        }
        if(soi == "WHC2") {
          fdo = 1 / (1 + exp(-1 * (oxy - 3.1485) / 0.345078))  ## Type III
        }
        if(soi == "SAE1") {          
          fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
        }
        if(soi == "SMB") {
          fdo = 1 / (1 + exp(-1 * (oxy - 3.687002) / 0.22578))  ## Type III, Bulkley 1975
        }
        if(soi == "WAE") {          
          fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
        }
        if(soi == "YEP") {     
          fdo = 1 / (1 + exp(-1 * (oxy - 2.370402) / 0.330795))  ## Type III, Roberts data
        }
        if(soi == "BCF") {     ## Torrans 2012 & 2008, Green et al 2012
          fdo = 1 / (1 + exp(-1 * (oxy - 1.19497) / 0.17284)) ## Type III
        }
        
        C = Cmax * ftc * fdo * preyb / (0.865 + preyb)  ##Constantini            
      }
      
      if(CEQ == 3) {     ## Consumption equation 3
        Gone = (1 / (CTO - CQ)) * log((0.98 * (1 - CKone)) / (CKone * 0.02)) 
        Lone = exp(Gone * (wt - CQ))
        Ka = (CKone * Lone) / (1 + CKone * (Lone - 1))
        
        Gtwo = (1 / (CTL - CTM)) * log((0.98 * (1 - CKfour)) / (CKfour * 0.02))
        Ltwo = exp(Gtwo * (CTL - wt))
        Kb = (CKfour * Ltwo) / (1 + CKfour * (Ltwo - 1))
        
        ftc = Ka * Kb
        
        Cmax = CA * predW^CB     ## this is g prey/g pred
        
        if(soi == "STB") {  ## Brandt et al 2009 and Chiba 1998
          fdo = 1 / (1 + exp(-1 * (oxy - 2.190932) / 0.290913)) ## Type III        
        }
        if(soi == "HSB") {  ## Emily's data
          fdo = 1 / (1 + exp(-1 * (oxy - 2.190932) / 0.290913)) ## Type III
        }
        
        C = Cmax * ftc * fdo * preyb / (0.865 + preyb)
      } 
      
      if(CEQ == 4) {     ## Consumption equation 4
        ftc = (CQ*wt + CK1*wt^2 + CK4*wt^3)
        
        Cmax = CA * predW^CB     ## this is g prey/g pred
        
        if(soi == "MUS") {  
          fdo = 1 / (1 + exp(-1 * (oxy - 2.563859) / 0.37121))      
        }
        
        C = Cmax * ftc * fdo * preyb / (0.865 + preyb)
      } 
      
      ## Calculate consumption for non-hypoxia simulations
      Cno = Cmax * ftc * preyb / (0.865 + preyb)
      
      ##  Respiration equations
      if(REQ == 1) {     ## Respiration equation 1
        ftr = exp(RQ * wt)
        VEL = RKone * predW^RKfour  ## no need to use equation set b/c T>RTL
        activity  = exp(RTO * VEL)
        S = SDA * (C - (FA * C))
        Sno = SDA * (Cno - (FA * Cno))  ## if no hypoxia effect
        OC = 13560 #oxycalorific coefficient
        R = (RA * predW^RB * ftr) * activity * (OC/preyED)+S #Respiration equation from Zhang et al. 2014
        Rno = (RA * predW^RB *ftr) * activity * (OC/preyED)+Sno
      }
      
      if(REQ == 2) {     ## Respiration equation 2
        V = (RTM - wt) / (RTM - RTO)
        Z = log(RQ) * (RTM - RTO)
        Y = log(RQ) * (RTM - RTO +2)
        X = (Z^2 * (1 + (1 + 40 / Y)^0.5)^2) / 400
        
        ftr = V^X * exp(X * (1 - V))
        activity  = ACT
        
        S = SDA * (C - (FA * C))
        Sno = SDA * (Cno - (FA * Cno))  ## if no hypoxia effect
        OC = 13560 #oxycalorific coefficient
        R = (RA * predW^RB * ftr) * activity * (OC/preyED)+S #Respiration equation from Zhang et al. 2014
        Rno = (RA * predW^RB *ftr) * activity * (OC/preyED)+Sno
      }
      
      ##  Egestion and excretion equations
      if(EGEXEQ == 1) {     ## Egestion and excretion equation 1
        FE = FA * C
        U = UA * (C - FE)
        
        FEno = FA * Cno  ## if no hypoxia effect
        Uno = UA * (Cno - FEno)  ## if no hypoxia effect
      }
      
      if(EGEXEQ == 2) {     ## Egestion and excretion equation 2
        p = C/Cmax
        FE = FA * wt^FB * exp(FG * p) * C
        U = UA * wt^UB * exp(UG * p) * (C - FE)
        
        pno = Cno/Cmax  ## if no hypoxia effect
        FEno = FA * wt^FB * exp(FG * pno) * Cno  ## if no hypoxia effect
        Uno = UA * wt^UB * exp(UG * pno) * (Cno - FEno)  ## if no hypoxia effect
      }
      
      ## Final Equation
      GRP = (preyED / predED)*(C-(R+FE+U))
      GRPno = (preyED / predED)*(Cno-(Rno+FEno+Uno))  ## if no hypoxia effect
      
      ## Stores data in master output dataset
      GRPd[rownum,1] <- easting
      GRPd[rownum,2] <- northing
      GRPd[rownum,3] <- z
      GRPd[rownum,4] <- wt
      GRPd[rownum,5] <- oxy
      GRPd[rownum,6] <- preyd
      GRPd[rownum,7] <- preyb
      GRPd[rownum,8] <- GRP
      GRPd[rownum,9] <- GRPno
      
    }
  } 
  
  colnames(GRPd) <- c("Easting", "Northing", "Depth", "Temp", "DO",
                      "Density", "Biomass", "GRP", "GRPnohyp")
  
  ##Subset data frame to only include depths 2 m and greater
  GRPd2 <- GRPd[which(GRPd$Depth > 1),names(GRPd) %in% c("Easting", "Northing", "Depth", "Temp", "DO",
                                                         "Density", "Biomass", "GRP", "GRPnohyp")]
    maxGRP = 0.0374729957076764
    medGRP = 0.0142377854882836
    minGRP = -0.0089974247311092

    
    ## Generate transect GRP data for each depth
    for(i in 1:nrow(xsectline)) {                ## loop through points on xsect
      easting = as.numeric(xsectline[i, 1])    ## easting point
      northing = as.numeric(xsectline[i, 2])   ## northing point
      
      ## find closest east and north in data
      distance <- sqrt((GRPd$Easting - easting)^2 + (GRPd$Northing - northing)^2)
      nearestdata <- GRPd[which(distance == min(distance)), ]
      
      ## Find bottom depth for ploting
      maxdepth3D <- max(grdxsect[which(grdxsect$x == i), 2])
      
      ## Subset commands just in case there are multiple nearest data points
      if(length(unique(nearestdata$Easting)) > 1) {
        nearestdata <- nearestdata[which(nearestdata$Easting == min(nearestdata$Easting)), ]
      }
      if(length(unique(nearestdata$Northing)) > 1) {
        nearestdata <- nearestdata[which(nearestdata$Northing == min(nearestdata$Northing)), ]
      }
      
      for(ii in 0:maxz) {           ## loop through each depth in xsect                
        xsectdata[i + nrow(xsectline) * ii, 9] <- nearestdata[ii + 1, 8]
        xsectdata[i + nrow(xsectline) * ii, 10] <- nearestdata[ii + 1, 9]
        
        if(ii >= nrow(nearestdata)) {  ## if data are lacking, find near point min for GRP (corrects for bottom anomalies)
          possibledata <- c(xsectdata[i + nrow(xsectline) * ii - nrow(xsectline), 9], 
                            xsectdata[i + nrow(xsectline) * ii - 1, 9],
                            xsectdata[i + nrow(xsectline) * ii - 2, 9],
                            xsectdata[i + nrow(xsectline) * ii - 16, 9])
          xsectdata[i + nrow(xsectline) * ii, 9] = min(na.omit(possibledata))
          possibledatano <- c(xsectdata[i + nrow(xsectline) * ii - nrow(xsectline), 10], 
                              xsectdata[i + nrow(xsectline) * ii - 1, 10],
                              xsectdata[i + nrow(xsectline) * ii - 2, 10],
                              xsectdata[i + nrow(xsectline) * ii - 16, 10])
          xsectdata[i + nrow(xsectline) * ii, 10] = min(na.omit(possibledatano))
        }
        if(ii > maxdepth3D) xsectdata[i + nrow(xsectline) * ii, 9:10] = NA
      }
    }
    colnames(xsectdata) <- c("x", "y", "Easting", "Northing", "Temp", "DO",
                             "Density", "lnBiomass", "GRP", "GRPnohyp")
    
    infval <- which(xsectdata$GRP == Inf)  ## ID infinite values because of bottom
    xsectdata[infval, 9] = NA              ## Change to NA and remove
    grdxsect <- na.omit(xsectdata)
    
    
    ## Generate GRP transect plot
    xsplotGRP <- ggplot() +
      geom_tile(data = grdxsect, aes(x = x, y = y, fill = GRP)) +  ## basic plot
      scale_fill_gradientn(expression(atop("",
                                           atop(textstyle("GRP with "*italic(f)[c]*"DO"),
                                                atop(textstyle("(g/g/d)"))))), 
                                                  colours = rainbow(4), 
                                                    breaks = c(0, 0.015, 0.03), 
                                                      limits = c(minGRP, maxGRP)) +
      theme(legend.text = element_text(size=14), legend.title = element_text(size = 18)) +
      theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  ## remove gridlines
      scale_x_continuous(expand = c(0,-0.5)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      scale_y_reverse(expand = c(0,0), name = "Depth (m)")  ## Reverse axis on depth scale
    
    xsplotGRPno <- ggplot() +
      geom_tile(data = grdxsect, aes(x = x, y = y, fill = GRPnohyp)) +  ## basic plot
      scale_fill_gradientn(expression(atop("",
                                           atop(textstyle("GRP without "*italic(f)[c]*"DO"),
                                                atop(textstyle("(g/g/d)"))))), 
                                                  colours = rainbow(4), 
                                                    breaks = c(0, 0.015, 0.03), 
                                                      limits = c(minGRP, maxGRP)) + 
      theme(legend.text = element_text(size=14), legend.title = element_text(size = 18)) +
      theme(panel.background = element_rect(fill = "black", color = "black")) +  ## black background
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  ## remove gridlines
      scale_x_continuous(expand = c(0,-0.5)) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      scale_y_reverse(expand = c(0,0), name = "Depth (m)")  ## Reverse axis on depth scale
    
    ## Generate 2D areal map for GRP
    depths = paste("m", 0:maxz, sep = "") ## make character names for 1-m plots
    GRPlist <- list()   ## define list to hold all GRP by depth plots
    
    depthsno = paste("mno", 0:maxz, sep = "") ## make character names for 1-m plots
    GRPnolist <- list()   ## define list to hold all GRP no hyp by depth plots
    
    if(maxz < 10) dlist = c(2, 3, 4, 5, 6, 7, 8, 9)
    if(maxz == 10) dlist = c(2, 3, 4, 5, 6, 7, 8, 10)
    if(maxz > 10) dlist = c(2, 3, 4, 5, 6, 7, 10, 12)
    
    for(i in 1:8) {
      doi = dlist[i]
      
      if(doi <= maxz) {
        depthGRP <- subset(GRPd, GRPd$Depth == doi, select = c(Easting, Northing,
                                                               GRP, GRPnohyp, Depth))
      }
      
      if(doi < maxz) {
        layr <- ggplot() + 
          geom_tile(data=depthGRP, aes(x=Easting, y=Northing, fill = GRP)) +
          scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                               limits = c(minGRP, maxGRP)) +
          geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
          islands +
          coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
          theme(panel.background = element_rect(fill = "black", color = "black")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(plot.title = element_text(size=22)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
          theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
          ggtitle(paste(doi, " m", sep = ""))
        
        assign(paste0(depths[doi+1]), layr)
        GRPlist <- c(GRPlist, list(layr))
        
        
        layr <- ggplot() + 
          geom_tile(data=depthGRP, aes(x=Easting, y=Northing, fill = GRPnohyp)) +
          scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                               limits = c(minGRP, maxGRP)) +
          geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
          islands +
          coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
          theme(panel.background = element_rect(fill = "black", color = "black")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(plot.title = element_text(size=22)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
          theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
          ggtitle(paste(doi, " m", sep = ""))
        
        assign(paste0(depthsno[doi+1]), layr)
        GRPnolist <- c(GRPnolist, list(layr))
      }
      
      if(doi == maxz) {
        layr <- ggplot() + 
          geom_tile(data=depthGRP, aes(x=Easting, y=Northing, fill = GRP)) +
          scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                               limits = c(minGRP, maxGRP)) +
          geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
          islands +
          coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
          theme(panel.background = element_rect(fill = "black", color = "black")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(plot.title = element_text(size=22)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
          theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
          ggtitle(paste(doi, " m", sep = ""))
        
        assign(paste0(depths[doi+1]), layr)
        GRPlist <- c(GRPlist, list(layr))
        
        layr <- ggplot() + 
          geom_tile(data=depthGRP, aes(x=Easting, y=Northing, fill = GRPnohyp)) +
          scale_fill_gradientn(colours = rainbow(4), guide = FALSE, na.value = "black", 
                               limits = c(minGRP, maxGRP)) +
          geom_polygon(data=outside, aes(x=long, y=lat), fill="white", asp = 1) +   
          islands +
          coord_fixed(ratio = 1, xlim=c(xmin,xmax), ylim=c(ymin, ymax)) +  ## aspect
          theme(panel.background = element_rect(fill = "black", color = "black")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(plot.title = element_text(size=22)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
          theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
          ggtitle(paste(doi, " m", sep = ""))
        
        assign(paste0(depthsno[doi+1]), layr)
        GRPnolist <- c(GRPnolist, list(layr))
      }
      
    }
      
    
    
    ## Code to print multiple areal GRP figures and one xsect at bottom
    jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/", lake,
                         dayoi, "_", soi, "_GRP_HYP.jpeg", sep = ""), width = 2000, 
        height = 2580, units = "px", res = 300, restoreConsole = TRUE)
    
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 4, heights = unit(c(0.5, 5, 5, 5, 4.5), "null"))))
    
    if(maxz < 10) {
      grid.text(paste(lake, dayoi, soi, "hypoxia", sep = ", "), vp = vplayout(1,2))
      if(maxz > 1) print(m2, vp = vplayout(2, 1))
      if(maxz > 2) print(m3, vp = vplayout(2, 2))
      if(maxz > 3) print(m4, vp = vplayout(2, 3))
      if(maxz > 4) print(m5, vp = vplayout(2, 4))
      if(maxz > 5) print(m6, vp = vplayout(3, 1))
      if(maxz > 6) print(m7, vp = vplayout(3, 2))
      if(maxz > 7) print(m8, vp = vplayout(3, 3))
      if(maxz > 8) print(m9, vp = vplayout(3, 4))
      print(xsplotGRP, vp = vplayout(4, 1:4))
    }
    if(maxz == 10) {
      grid.text(paste(lake, dayoi, soi, "hypoxia", sep = ", "), vp = vplayout(1,2))
      if(maxz > 1) print(m2, vp = vplayout(2, 1))
      if(maxz > 2) print(m3, vp = vplayout(2, 2))
      if(maxz > 3) print(m4, vp = vplayout(2, 3))
      if(maxz > 4) print(m5, vp = vplayout(2, 4))
      if(maxz > 5) print(m6, vp = vplayout(3, 1))
      if(maxz > 6) print(m7, vp = vplayout(3, 2))
      if(maxz > 7) print(m8, vp = vplayout(3, 3))
      if(maxz > 9) print(m10, vp = vplayout(3, 4))
      print(xsplotGRP, vp = vplayout(4, 1:4))
    }
    if(maxz > 10) {
      grid.text(paste(lake, dayoi, sep = ", "), vp = vplayout(1,2))
      if(maxz > 1) print(m2, vp = vplayout(2, 1))
      if(maxz > 2) print(m3, vp = vplayout(2, 2))
      if(maxz > 3) print(m4, vp = vplayout(2, 3))
      if(maxz > 4) print(m5, vp = vplayout(2, 4))
      if(maxz > 5) print(m6, vp = vplayout(3, 1))
      if(maxz > 6) print(m7, vp = vplayout(3, 2))
      if(maxz > 9) print(m10, vp = vplayout(3, 3))
      if(maxz > 11) print(m12, vp = vplayout(3, 4))
      print(xsplotGRP, vp = vplayout(4, 1:4))
  }
    dev.off()
    
    jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/", lake,
                         dayoi, "_", soi, "_GRP_HYPno.jpeg", sep = ""), width = 2000,
        height = 2580, units = "px", res = 300, restoreConsole = TRUE)
    
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 4, heights = unit(c(0.5, 5, 5, 5, 4.5), "null"))))
    
    if(maxz < 10) {
      grid.text(paste(lake, dayoi, soi, "no hypoxia", sep = ", "), vp = vplayout(1,2))
      if(maxz > 1) print(mno2, vp = vplayout(2, 1))
      if(maxz > 2) print(mno3, vp = vplayout(2, 2))
      if(maxz > 3) print(mno4, vp = vplayout(2, 3))
      if(maxz > 4) print(mno5, vp = vplayout(2, 4))
      if(maxz > 5) print(mno6, vp = vplayout(3, 1))
      if(maxz > 6) print(mno7, vp = vplayout(3, 2))
      if(maxz > 7) print(mno8, vp = vplayout(3, 3))
      if(maxz > 8) print(mno9, vp = vplayout(3, 4))
      print(xsplotGRPno, vp = vplayout(4, 1:4))
    }
    if(maxz == 10) {
      grid.text(paste(lake, dayoi, soi, "no hypoxia", sep = ", "), vp = vplayout(1,2))
      if(maxz > 1) print(mno2, vp = vplayout(2, 1))
      if(maxz > 2) print(mno3, vp = vplayout(2, 2))
      if(maxz > 3) print(mno4, vp = vplayout(2, 3))
      if(maxz > 4) print(mno5, vp = vplayout(2, 4))
      if(maxz > 5) print(mno6, vp = vplayout(3, 1))
      if(maxz > 6) print(mno7, vp = vplayout(3, 2))
      if(maxz > 7) print(mno8, vp = vplayout(3, 3))
      if(maxz > 9) print(mno10, vp = vplayout(3, 4))
      print(xsplotGRPno, vp = vplayout(4, 1:4))
    }
    if(maxz > 10) {
      grid.text(paste(lake, dayoi, soi, "no hypoxia", sep = ", "), vp = vplayout(1,2))
      if(maxz > 1) print(mno2, vp = vplayout(2, 1))
      if(maxz > 2) print(mno3, vp = vplayout(2, 2))
      if(maxz > 3) print(mno4, vp = vplayout(2, 3))
      if(maxz > 4) print(mno5, vp = vplayout(2, 4))
      if(maxz > 5) print(mno6, vp = vplayout(3, 1))
      if(maxz > 6) print(mno7, vp = vplayout(3, 2))
      if(maxz > 9) print(mno10, vp = vplayout(3, 3))
      if(maxz > 11) print(mno12, vp = vplayout(3, 4))
      print(xsplotGRPno, vp = vplayout(4, 1:4))
    }
    dev.off()
    


    ## Generate summary data and write to file
    Summary3D <- data.frame()
    Summary3D[1,1] <- lake
    Summary3D[1,2] <- as.character(dayoi)
    Summary3D[1,3] <- soi
    Summary3D[1,4] <- "Prop hypoxic"
    DOoi <- which(GRPd2$DO < 2.0)
    Summary3D[1,5] <- length(DOoi)/nrow(GRPd2)
    
    
    Summary3D[2,2] <- "Mean GRP"
    Summary3D[2,3] <- "SE GRP"
    Summary3D[2,4] <- "Max GRP"
    Summary3D[2,5] <- "Mean GRP pos cells"
    Summary3D[2,6] <- "Proportion pos cells"
    Summary3D[2,7] <- "Mean Temperature"
    Summary3D[2,8] <- "SE Temperature"
    Summary3D[2,9] <- "Mean Biomass"
    Summary3D[2,10] <- "SE Biomass"
    Summary3D[2,11] <- "Mean DO"
    Summary3D[2,12] <- "SE DO"
    
    
    Summary3D[3,1] <- "With hypoxia"
    summaryGRP <- na.omit(GRPd2)
    Summary3D[3,2] <- mean(GRPd2$GRP)
    Summary3D[3,3] <- sd(GRPd2$GRP) /  sqrt(length(GRPd2$GRP[!is.na(GRPd2$GRP)]))
    Summary3D[3,4] <- max(GRPd2$GRP)
    grpoi <- which(GRPd2$GRP > 0)
    Summary3D[3,5] <- mean(GRPd2$GRP[grpoi])
    Summary3D[3,6] <- length(which(GRPd2$GRP > 0)) / nrow(GRPd2)
    Summary3D[3,7] <- mean(GRPd2$Temp)
    Summary3D[3,8] <- sd(GRPd2$Temp) /  sqrt(length(GRPd2$Temp[!is.na(GRPd2$Temp)]))
    Summary3D[3,9] <- mean(GRPd2$Biomass)
    Summary3D[3,10] <- sd(GRPd2$Biomass) /  sqrt(length(GRPd2$Biomass[!is.na(GRPd2$Biomass)]))
    Summary3D[3,11] <- mean(GRPd2$DO)
    Summary3D[3,12] <- sd(GRPd2$DO) /  sqrt(length(GRPd2$DO[!is.na(GRPd2$DO)]))
    
    Summary3D[4,1] <- "Without hypoxia"
    summaryGRP <- na.omit(GRPd2)
    Summary3D[4,2] <- mean(GRPd2$GRPnohyp)
    Summary3D[4,3] <- sd(GRPd2$GRPnohyp) /  sqrt(length(GRPd2$GRPnohyp[!is.na(GRPd2$GRPnohyp)]))
    Summary3D[4,4] <- max(GRPd2$GRPnohyp)
    grpoi <- which(GRPd2$GRPnohyp > 0)
    Summary3D[4,5] <- mean(GRPd2$GRPnohyp[grpoi])
    Summary3D[4,6] <- length(which(GRPd2$GRPnohyp > 0)) / nrow(GRPd2)
    Summary3D[4,7] <- mean(GRPd2$Temp)
    Summary3D[4,8] <- sd(GRPd2$Temp) /  sqrt(length(GRPd2$Temp[!is.na(GRPd2$Temp)]))
    Summary3D[4,9] <- mean(GRPd2$Biomass)
    Summary3D[4,10] <- sd(GRPd2$Biomass) /  sqrt(length(GRPd2$Biomass[!is.na(GRPd2$Biomass)]))
    Summary3D[4,11] <- mean(GRPd2$DO)
    Summary3D[4,12] <- sd(GRPd2$DO) /  sqrt(length(GRPd2$DO[!is.na(GRPd2$DO)]))
    
    
    ## Print results
    write.table(GRPd2, paste("C:/Users/richa/OneDrive/Desktop/",
                             dayoi, "_", soi, "_", lake, ".txt", sep = ""), row.names = TRUE)
    write.table(Summary3D, paste("C:/Users/richa/OneDrive/Desktop/",
                                 dayoi, "_", soi, "_", lake, "_sum.txt", sep = ""), row.names = TRUE)
}
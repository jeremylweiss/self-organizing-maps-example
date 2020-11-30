

# EXAMPLE SELF-ORGANIZING MAP (SOM) ANALYSIS WITH NCEP-NCAR DATA


# This code downloads NCEP NCAR R2 data (see:
# https://climatedataguide.ucar.edu/climate-data/ncep-reanalysis-r2 
# and https://en.wikipedia.org/wiki/Atmospheric_reanalysis for descriptions) and 
# proceeds through a self-organizing map analysis of daily regional atmospheric 
# pressure patterns for prescribed months over the Southwestern United States 

# The following code is based on initial work by: 
# Michael Crimmins, Extension Specialist and Professor
# Department of Environmental Science
# University of Arizona

# Author:
# Jeremy Weiss, Climate and Geospatial Extension Scientist
# School of Natural Resources and the Environment
# University of Arizona
# 520-626-8063, jlweiss@email.arizona.edu

# Additional reading:
# 1. Wehrens & Buydens (2007) JStatSoft
# 2. https://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/


# SETUP --------------------


# Set months of interest
month_start <- 10
month_end <- 11

# Load needed libraries
library("maps")  # needed for "RNCEP" package
library("RNCEP")  # working with NCEP-NCAR data
library("kohonen")  # for SOM analysis

library("dplyr")
library("reshape2")
library("tidyr")
library("ggplot2")

library("PBSmapping")

# Function to help with visualization of results
colors.node.counts <- function(n, alpha = 1) {
  heat.colors(n, alpha = alpha)[n:1]
}


# DOWNLOAD AND TRANSFORM NCEP-NCAR R2 DATA --------------------


# Downloaded NCEP-NCAR R2 data is automatically brought into the workspace as a 
# 3-D array. In this case, the variable of interest is 500-millibar geopotential 
# height for the months prescribed above over the Southwestern US
gph500 <- NCEP.gather(variable = 'hgt',
                      level = 500,
                      months.minmax = c(month_start, month_end),
                      years.minmax = c(1979, 2019),
                      lat.southnorth = c(20, 50),
                      lon.westeast = c(-130, -90),
                      reanalysis2 = TRUE,
                      return.units = TRUE,
                      status.bar = TRUE)

# Save the data array to the current directory for future development
#save(gph500,file = "gph500_feb_apr_7919.Rdata")

# To load the saved array at a later time, use the following:
#load("gph500_feb_apr_7919.RData")

# Confirm area of interest by mapping a single layer (i.e., timestep) of the 
# data array
NCEP.vis.area(wx.data = gph500,
              layer = 1,
              show.pts = TRUE,
              draw.contours = TRUE,
              cols = topo.colors(n = 16, alpha = 0.5),
              transparency = 0.5,
              axis.args = NULL,
              map.args = NULL,
              grid.args= NULL,
              title.args = list(
                main = "example layer of geopotential height data",
                xlab = "longitude (degrees)",
                ylab = "latitude (degrees)"
                ),
              interp.loess.args = list(span = 0.75),
              image.plot.args = 
                list(legend.args = list(text = "m", cex = 1.0, padj = -0.5,
                                        adj = 0.0)),
               contour.args = list(labcex = 0.75),
               points.args = list(pch = 20, cex = 0.75))

# NCEP NCAR R2 data are on a six-hour timestep. Temporally aggregate the data on 
# a daily timestep, averaging values from individual days
gph500_daily <- 
  NCEP.aggregate(wx.data = gph500, YEARS = TRUE, MONTHS = TRUE, DAYS = TRUE,
                 HOURS = FALSE, fxn = "mean")
rm(gph500)

# Convert the data array to a dataframe with columns 'datetime', 'latitude', 
# 'longitude', and 'gph500' and rows as individual space-time-height value 
# combinations
gph500_daily <- NCEP.array2df(wx.data = gph500_daily, var.names = "gph500")
rm(gph500_daily)






# Reshape the dataframe from long to wide format. This is to 
# place all data for one day into one row (i.e., a sample in the
# case of SOM analysis) with following columns (i.e., gridpoints 
# in the analysis domain) as corresponding values.
gph500_df_daily_wide <- dcast(data = gph500_daily_df,
                               formula = datetime ~ latitude + longitude,
                               value.var = "gph500")
rm(gph500_daily_df)


##################################################################
## C. RUN SOM ANALYSIS, TESTING SOM GRID ALGORITHM OPTION
##################################################################


# Define the number of nodes (i.e., key patterns) to retain in 
# the SOM analysis. Start with a SOM grid size that is common
# in climate studies and increase from there.
nrows <- 5
ncols <- 7
for (i in 1:12) {
  
  # Run the self-organizing map analysis. Note that the input data
  # must be in matrix form and that the datetime column is not 
  # needed. The 'somgrid' function sets up a grid of units for the
  # analysis. This grid represents an atmospheric pattern topology
  # to which atmospheric patterns of individual days are linked, or
  # associated. 
  set.seed(7)
  gph500_som <- som(X = as.matrix(gph500_df_daily_wide[,2:ncol(gph500_df_daily_wide)]),
                     grid = somgrid(xdim = ncols,
                                     ydim = nrows,
                                     topo = "rectangular"),
                     rlen = 1000,  # number of times data presented to network
                     alpha = c(1.0,0.001),  # learning rate range
                     keep.data = TRUE)
  
  # In this case, the 'som()' function returns a kohonen-type 
  # object that is a list of 13 elements.
  summary(gph500_som)
  names(gph500_som)
  
  # Check how the SOM training progressed (progress here is defined
  # as the distance from the weights of each node to the samples
  # represented by that node becoming reduced as the number of 
  # times that the data are presented to the network increases).
  # Ideally, distances should approach minimimal values 
  # (represented on this graph as an inverse plateau).
  plot(x = gph500_som,
        type = "changes",
        main = paste0("training progress : ",
                       "nrows x ncols = ",
                       nrows,
                       " x ",
                       ncols))
  
  # Display the number of samples (i.e., daily geopotential heights)
  # per grid node. One ideal result is to have counts be relatively 
  # uniform across the grid, with 5-10 samples per node. Nodes with
  # higher (lower) sample numbers suggest that a larger (smaller)
  # map is needed.
  plot(x = gph500_som,
        type = "counts",
        main = paste0("node counts : ",
                       "nrows x ncols = ",
                       nrows,
                       " x ",
                       ncols),
        palette.name = colors.node.counts)
  
  # Increase dimensions of the SOM grid to test this algorithm
  # option.
  nrows <- nrows+2
  ncols <- ncols+2
}
rm(i,nrows,ncols)









##################################################################
## D. VISUALIZE SOM ANALYSIS RESULTS
##################################################################


# Check how the SOM training progressed (progress here is defined
# as the distance from the weights of each node to the samples
# represented by that node becoming reduced as the number of 
# times that the data are presented to the network increases).
# Ideally, distances should approach minimimal values 
# (represented on this graph as an inverse plateau).
plot(x = gph500_som,
      type = "changes",
      main = "training progress")

#
plot(x = gph500_som,
      type = "counts",
      main = "node counts",
      palette.name = colors.node.counts)





##### START EDITING HERE NEXT TIME







# From: https://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/
# 
plot(gph500_som,type="dist.neighbours",main = "SOM neighbour distances")
plot(gph500_som,type="codes")
plot(gph500_som,
      type="property",
      property=getCodes(gph500_som)[,4],
      main=colnames(getCodes(gph500_som))[4])


plot(gph500_som,
      type = "mapping",
      col = 
      main = "yeast data")



# From: 'kohonen.pdf'
# 
# Check characteristics of SOM analysis output.
# 
# Dimensions of extracted codebook vectors.
dim(getCodes(gph500_som))
# Summary 
summary(gph500_som) 




##################################################################
## D. RUN SOM ANALYSIS, TESTING ALGORITHM OPTIONS
##################################################################


# Define the number of nodes (i.e., key patterns) to retain in 
# the SOM analysis.
nrows <- 19
ncols <- 21

# Run the self-organizing map analysis. Note that the input data
# must be in matrix form and that the datetime column is not 
# needed. The 'somgrid' function sets up a grid of units for the
# analysis. This grid represents an atmospheric pattern topology
# to which atmospheric patterns of individual days are linked, or
# associated. 
set.seed(7)
gph500_som <- som(X = as.matrix(gph500_df_daily_wide[,2:ncol(gph500_df_daily_wide)]),
                   grid = somgrid(xdim = ncols,
                                   ydim = nrows,
                                   topo = "rectangular"),
                   rlen = 1000,  # number of times data presented to network
                   alpha = c(1.0,0.001),  # learning rate range
                   keep.data = TRUE)

# In this case, the 'som()' function returns a kohonen-type 
# object that is a list of 13 elements.
summary(gph500_som)
names(gph500_som)






# Convert the 'codes' of the SOM analysis to a dataframe. This
# represents 500-mb geopotential height values for each gridpoint
# for each of the SOM analysis nodes, or maps.  
gph500.codebook <- as.data.frame(gph500_som$codes)

# Extract the grid points from the SOM analysis output.
code.grid <- as.data.frame(gph500_som$grid$pts)
code.grid$mapUnit <- seq(1,nrow(code.grid))
code.grid <- code.grid %>%
  unite(y,x,col="codes",sep="_") # add mapunit back in if needed


# xCols<-as.data.frame(som.gh500$grid$pts[,1])
# colnames(xCols)<-"X_Cols"
# yRows<-as.data.frame(som.gh500$grid$pts[,2])
# colnames(yRows)<-"Y_Rows"
codebook <- cbind(code.grid,gph500.codebook)
codebook.long <- melt(codebook,id.vars=1)





## BLOCK START
# TESTING
codebook.long.mapunit <- codebook.long[1:(ncols*nrows),]

# original
#codebook.long <- separate(codebook.long,
#                          variable,
#                          convert=TRUE,
#                          into=c("lat","lon"),
#                          sep="_"
#                         )

# TESTING
codebook.long <- separate(codebook.long[(ncols*nrows+1):nrow(codebook.long),],
                           variable,
                           convert=TRUE,
                           into=c("lat","lon"),
                           sep="_"
)

## BLOCK END



codebook.long$lat <- as.numeric(gsub("X","",codebook.long$lat))

codebook.long$lon <- codebook.long$lon-360

codebook.long <- separate(codebook.long,
                           codes,
                           convert=FALSE,
                           remove=FALSE,
                           into=c("xCols","yRows"),
                           sep="_"
)

# assign days to nodes
#  NOTE THAT THE FUNCTION map() WILL REVERT TO THE 'maps' 
#  LIBRARY VERSION AND PRODUCE AN ERROR
nodes <- map(gph500_som)

somTime <- as.data.frame(cbind(gph.500.daily.df.wide$datetime,
                                 nodes$unit.classif,
                                 nodes$distances)
)
colnames(somTime) <- c("datetime","mapUnit","errorDist")
somTime <- separate(somTime,
                     datetime,
                     convert=TRUE,
                     remove=FALSE,
                     into=c("year","month","day","hr"),
                     sep="_"
)

somTime$date <- as.Date(paste(somTime$year,"-",
                                somTime$day,"-",
                                somTime$month,
                                sep=""),
                         format="%Y-%d-%m"
)

somTime$doy <- as.numeric(format(somTime$date,"%j"))

somTime$mapUnit <- as.integer(somTime$mapUnit)

somTime$errorDist <- as.numeric(as.character(somTime$errorDist))


# plot map - fix lines http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
# plot limits
#xlim = c(-125,-95)
#ylim = c(20,50)
xlim = c(-126.25,-93.75)
ylim = c(18.75,51.25)



# in ggplot
# Reintroduce this package here.
library("maps")
all_states <- map_data("state")
world<-map_data("world")



colnames(world)<-c("X","Y","PID","POS","region","subregion")
#library(PBSmapping)
world = clipPolys(world, xlim=xlim,ylim=ylim, keepExtra=TRUE)
#colnames(world)<-c("Lon","Lat","PID","POS","region","subregion") 

colnames(all_states)<-c("X","Y","PID","POS","region","subregion")
all_states = clipPolys(all_states, xlim=xlim,ylim=ylim, keepExtra=TRUE)
#colnames(all_states)<-c("Lon","Lat","PID","POS","region","subregion")


# Mike's version

#p <- ggplot() +
# geom_polygon(data=world, aes(x=X, y=Y, group = PID),colour="black", fill=NA) +
# #scale_x_continuous(breaks = c(-120,-140)) +
# facet_wrap(~codes, nrow = nrows, ncol=ncols)+theme_bw()+
# stat_contour(data=codebook.long,aes(lon,lat,z=value,
#                                       color=..level..),
#               size=1) +
# geom_raster(data=codebook.long,aes(x=lon,
#                                      y=lat,
#                                      fill=value)) +
# coord_quickmap()+
# scale_colour_distiller(palette = "Spectral", name="500mb GPH (m)")
# #scale_color_continuous(low="blue", high = "red")+ 
# #coord_map(xlim = xlim,ylim = ylim)+
# #labs(x="Lon", y="Lat")
#p





library("extrafont")
font_import()
y

p <- ggplot() +
  facet_wrap(~codes,nrow=nrows,ncol=ncols) +
  #theme_minimal() +
  #geom_raster(data=codebook.long[36:nrow(codebook.long),],
  geom_raster(data=codebook.long,
               aes(x=lon,
                    y=lat,
                    fill=value) 
 ) +
  geom_polygon(data=world,
                aes(x=X,
                     y=Y,
                     group=PID),
                color="black",
                fill=NA,
                size=0.25) +
  geom_polygon(data=all_states,
                aes(x=X,
                     y=Y,
                     group=PID),
                color="black",
                fill=NA,
                size=0.25) +
  stat_contour(data=codebook.long[36:nrow(codebook.long),],
                aes(x=lon,
                     y=lat,
                     z=value),
                color="black",
                linetype="longdash",
                size=0.25) +
  scale_fill_distiller(palette="Spectral") +
  coord_quickmap() +
  theme_minimal(base_family="Consolas") +
  theme(axis.text=element_text(color="gray30",size=8),
         axis.ticks.length=unit(0.0,"mm"),
         axis.title=element_text(color="gray30",face="bold",size=8),
         legend.key.size=unit(0.5,"cm"),
         legend.text=element_text(size=10),
         legend.title=element_text(face="bold",
                                    size=10),
         panel.grid.major.x=element_line(color="gray30",size=0.3),
         panel.grid.major.y=element_line(color="gray30",size=0.3),
         strip.background=element_rect(fill="white",color="white",size=0.1),
         strip.text=element_text(color="gray30",face="bold",size=8)
 )
p

ggsave("./som_4x6_maps_ncep_ncar_r2_southwest_jfm_7818.png",
        plot=p,
        device="png",
        path=NULL,
        scale=1,
        width=10,
        height=7.5,
        units="in",
        dpi=600)










# summary plots -- appears to plot opposite up/down from SOM plot
counts <- plot(gph500_som, type="counts", shape = "straight", labels=counts)
codes <- plot(gph500_som, type="codes", shape = "straight")
similarities <- plot(gph500_som, type="quality", palette.name = terrain.colors)
plot(gph500_som, type="dist.neighbours", main = "SOM neighbour distances")

# sammon mapping
library(MASS)
gph500.codes <- gph500_som$codes
dis <- dist(as.matrix(gph500_som$codes[[1]]))
gph500.sam <- sammon(dis)
plot(gph500.sam$points, type="n")
text(gph500.sam$points,labels=as.character(1:nrow(code.grid)))




# join codes to node table
somTime<-left_join(x=somTime,
                    y=code.grid,
                    by="mapUnit")
# plot map units
ggplot(somTime, aes(doy, year)) + 
  geom_tile(aes(fill = mapUnit), colour = "grey") + 
  geom_text(aes(label = codes), size=1.25)+
  scale_fill_gradient2(low = "red", mid = "green",
                       high = "blue", midpoint = 7, space = "Lab",
                       na.value = "grey50", guide = "colourbar")











# timeseries of map units - for facet wrap
year <- as.vector(1979:2018)
ts.facetwrap <- expand.grid(year,code.grid$codes)
colnames(ts.facetwrap) <- c("year","codes")
ts.facetwrap$ann_count <- NA

for (y in 1979:2018) {
  a <- filter(somTime,year==y)
  
  for (c in code.grid$codes) {
    b <- filter(a,codes==c)
    ts.facetwrap$ann_count[which(ts.facetwrap$year==y & ts.facetwrap$codes==c,
                                   arr.ind=FALSE)] <- nrow(b)
    rm(b)
  }
  
  rm(a)
}
rm(y)

p <- ggplot() +
  facet_wrap(~codes,nrow=nrows,ncol=ncols) +
  geom_line(data=ts.facetwrap,
             aes(x=year,
                  y=ann_count),
             color="black",
             size=0.75) +
  theme_minimal(base_family="Consolas") +
  theme(#axis.line=element_line(color="gray30",size=0.3),
    #axis.line.x.bottom=element_line(color="gray30",size=0.3),
    #axis.line.x.top=element_line(color="gray30",size=0.3),
    #axis.line.y.left=element_line(color="gray30",size=0.3),
    #axis.line.y.right=element_line(color="gray30",size=0.3),
    axis.text=element_text(color="gray30",size=8),
    #axis.ticks=element_line(color="gray30",size=0.3),
    #axis.ticks.length=unit(0.0,"mm"),
    axis.title=element_text(color="gray30",face="bold",size=8),
    legend.key.size=unit(0.5,"cm"),
    legend.text=element_text(size=10),
    legend.title=element_text(face="bold",
                               size=10),
    #panel.background=element_rect(fill="gray30"),
    #panel.border=element_rect(linetype="solid",fill=NA),
    panel.grid.major.x=element_line(color="gray50",size=0.3),
    panel.grid.major.y=element_line(color="gray50",size=0.3),
    #panel.grid.major.x=element_blank(),
    #panel.grid.major.y=element_blank(),
    strip.background=element_rect(fill="white",color="white",size=0.1),
    strip.text=element_text(color="gray30",face="bold",size=8)
 )
p

ggsave("./som_4x6_timeseries_ncep_ncar_r2_southwest_jfm_7818.png",
        plot=p,
        device="png",
        path=NULL,
        scale=1,
        width=10,
        height=7.5,
        units="in",
        dpi=600)



# Plot single time series




# OLD version
# 
# # timeseries of map units - individual nodes
year <- as.vector(1979:2018)
mapUnit <- as.vector(1:35)
ts <- expand.grid(year,mapUnit)
colnames(ts) <- c("year","mapUnit")
ts$ann_count <- NA

for (y in 1979:2018) {
  a <- filter(somTime,year==y)
  
  for (m in 1:35) {
    b <- filter(a,mapUnit==m)
    #if (nrow(b)==0) {
    # ts$ann_count[which(ts$year==y & ts$mapUnit==m,
    #                      arr.ind=FALSE)] <- 0
    #}
    #else {
    ts$ann_count[which(ts$year==y & ts$mapUnit==m,
                         arr.ind=FALSE)] <- nrow(b)
    #}
    
    rm(b)
  }
  
  rm(a)
}
rm(y)

ts2<-filter(ts,mapUnit==4)

ggplot(data=ts2,
        aes(x=year,y=ann_count,group=mapUnit)) +
  geom_line(aes(color=mapUnit))












# plot error
ggplot(somTime, aes(doy, year)) + 
  geom_tile(aes(fill = errorDist), colour = "grey") + 
  scale_fill_gradient2(low = "white", mid = "yellow",
                       high = "red", space = "Lab",
                       na.value = "grey50", guide = "colourbar")














# plot precip with local CPC precip


# plot precip with rnoaa cpc_prcp
library(rnoaa)

tmp<-cpc_prcp(date = somTime$date[i])
tmp$lon<-tmp$lon-180

for(i in 1:4){
  tmp<-cpc_prcp(date = somTime$date[i])
  if (i==1){
    tmp2 <- tmp
  }else{
    tmp2 <-bind_cols(tmp2,tmp[,3]) # brick or stack?
  }
}

tmp<-cpc_prcp(date = "1958-07-16", us=FALSE)
p <- ggplot(tmp, aes(x=lon, y=lat, fill=precip)) + theme_bw()
p + geom_tile()+
  ggtitle("cpc_prcp(date = '2016-07-16', us=FALSE)")

p <- ggplot(tmp, aes(x=lon, y=lat)) + theme_bw()
p + geom_raster(aes(fill=precip))

# plot CPC precip from netcdf files
library(ncdf4)
library(raster)

yr1<-1979
yr2<-2017

for(i in yr1:yr2){
  paste0(i)
  cpc.prcp.file <- paste0("/scratch/crimmins/cpc_global/ftp.cdc.noaa.gov/Datasets/cpc_global_precip/precip.",i,".nc")
  #"load" data into R via Raster package
  # cycle through each year, assigning each day to appropriate SOM unit
  prcp.tmp <- brick(cpc.prcp.file, varname="precip",  ncdf=TRUE)
  plot(prcp.tmp[[1]])
  
  # crop extent
  # e <- extent(-100, -66.7, 24.1, 51.2)
  # tminCrop <- crop(tmin, e)
  
  
}
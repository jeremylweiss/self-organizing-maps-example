

##################################################################
##  EXAMPLE SELF-ORGANIZING MAP (SOM) ANALYSIS WITH NCEP-NCAR DATA
##################################################################


#  This code downloads NCEP NCAR R2 data (see:
#  https://climatedataguide.ucar.edu/climate-data/ncep-reanalysis-r2 
#  and https://en.wikipedia.org/wiki/Atmospheric_reanalysis for 
#  descriptions) and proceeds through a self-organizing map 
#  analysis of daily regional atmospheric pressure patterns for 
#  the months of February, March, and April over the Southwestern
#  United States. 

#  The following code is based on that by: 
#  Michael Crimmins, Extension Specialist and Professor
#  Department of Environmental Science
#  University of Arizona

#  Author:
#  Jeremy Weiss, Climate and Geospatial Extension Scientist
#  School of Natural Resources and the Environment
#  University of Arizona
#  520-626-8063, jlweiss@email.arizona.edu


##################################################################
##  A. SETUP R ENVIRONMENT
##################################################################


#  Load the needed packages.
library( "maps" )  #  needed for "RNCEP" package
library( "RNCEP" )  #  needed for NCEP-NCAR data
library( "reshape2" )  #  for long/wide data reformatting

library( "kohonen" )  #
library( "dplyr" )  #
library( "tidyr" )  #
library( "ggplot2" )  #
library( "PBSmapping" )  #


##################################################################
##  B. DOWNLOAD AND TRANSFORM NCEP-NCAR R2 DATA
##################################################################


#  Download NCEP-NCAR R2 data. This is automatically brought into
#  the workspace as a 3-D array. In this case, the variable of 
#  interest is 500-millibar geopotential height for the months of 
#  February, March, and April over the Southwestern US.
gph_500 <- NCEP.gather( variable = 'hgt',
                        level = 500,
                        months.minmax = c( 2,4 ),
                        years.minmax = c( 1979,2019 ),
                        lat.southnorth = c( 20,50 ),
                        lon.westeast = c( -130,-90 ),
                        reanalysis2 = TRUE,
                        return.units = TRUE,
                        status.bar = TRUE )

#  Save the data array to the current directory for future code
#  development.
#save( gph_500,file = "gph_500_feb_apr_7919.Rdata" )

#  To load this array at a later time, use the following:
#load( "gph_500_feb_apr_7919.RData" )

#  Confirm area of interest by mapping a single layer (i.e.,
#  timestep) of the data array. 
NCEP.vis.area( wx.data = gph_500,
               layer = 1,
               show.pts = TRUE,
               draw.contours = TRUE,
               cols = topo.colors( 16,alpha = 0.5 ),
               transparency = 0.5,
               axis.args = NULL,
               map.args = NULL,
               grid.args= NULL,
               title.args = list( main = "example layer of geopotential height data",
                                  xlab = "longitude (degrees)",
                                  ylab = "latitude (degrees)"),
               interp.loess.args = list( span = 0.75 ),
               image.plot.args = list( legend.args = list( text = "m",
                                                           cex = 1.0,
                                                           padj = -0.5,
                                                           adj = 0.0 ) ),
               contour.args = list( labcex = 0.75 ),
               points.args = list( pch = 20,
                                   cex = 0.75 ) )


#  NCEP NCAR R2 data are on a six-hour timestep. Temporally
#  aggregate the data on a daily timestep, averaging values
#  from individual days.
gph_500_daily <- NCEP.aggregate( wx.data = gph_500,
                                 YEARS = TRUE,
                                 MONTHS = TRUE,
                                 DAYS = TRUE,
                                 HOURS = FALSE,
                                 fxn = "mean" )
rm( gph_500 )

#  Convert the data array to a dataframe. Columns will be
#  'datetime', 'latitude', 'longitude', 'gph500'. Rows will be
#  individual space-time-height value combinations.
df <- NCEP.array2df( wx.data = gph_500_daily,
                     var.names = "gph500" )
rm( gph_500_daily )


##################################################################
##  C. RUN SELF-ORGANIZING MAP (SOM) ANALYSIS
##################################################################


#  Reshape the dataframe from long to wide format. This is to 
#  place all data for one day into one row (i.e., a sample in the
#  case of SOM analysis) with following columns (i.e., gridpoints 
#  in the analysis domain) as corresponding values.
df_wide <- dcast( data = df,
                  formula = datetime ~ latitude + longitude,
                  value.var = "gph500"
)






##### START EDITING HERE NEXT TIME




#  Define the number of patterns to retain in the SOM analysis.
nrows <- 5
ncols <- 7


#  Run the Self-organizing map analysis. Note that the input data
#  must be in matrix form and that the datetime column is not 
#  needed. 'somgrid' sets up a grid of units for the analysis.
som.gph500 <- som( X=as.matrix( gph.500.daily.df.wide[ ,2:ncol( gph.500.daily.df.wide ) ] ),
                   grid=somgrid( xdim=ncols,
                                 ydim=nrows,
                                 topo="rectangular"
                   ),
                   rlen=500,
                   alpha=c( 0.05,0.01 ),  #  learning rate
                   #radius= ,  #  neighborhood radius
                   keep.data=TRUE
)



#  From: https://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/
#  
#  Check the SOM training progress As the SOM training iterations
#  progress, the distance from each node's weights to the samples
#  represented by that node is reduced. Ideally, this distance
#  should reach a minimum plateau. This plot option shows the
#  progress over time. If the curve is continually decreasing, 
#  more iterations are required.
plot( som.gph500,type="changes" )
plot( som.gph500,type="count",main="node counts" )
plot( som.gph500,type="dist.neighbours",main = "SOM neighbour distances" )
plot( som.gph500,type="codes" )
plot( som.gph500,
      type="property",
      property=getCodes( som.gph500 )[ ,4 ],
      main=colnames( getCodes( som.gph500 ) )[ 4 ] )



#  From: 'kohonen.pdf'
#  
#  Check characteristics of SOM analysis output.
#  
#  Dimensions of extracted codebook vectors.
dim( getCodes( som.gph500 ) )
#  Summary 
print( som.gph500 )
summary( som.gph500 ) 






#  Convert the 'codes' of the SOM analysis to a dataframe. This
#  represents 500-mb geopotential height values for each gridpoint
#  for each of the SOM analysis nodes, or maps.  
gph500.codebook <- as.data.frame( som.gph500$codes )

#  Extract the grid points from the SOM analysis output.
code.grid <- as.data.frame( som.gph500$grid$pts )
code.grid$mapUnit <- seq( 1,nrow( code.grid ) )
code.grid <- code.grid %>%
  unite( y,x,col="codes",sep="_" ) # add mapunit back in if needed


# xCols<-as.data.frame(som.gh500$grid$pts[,1])
# colnames(xCols)<-"X_Cols"
# yRows<-as.data.frame(som.gh500$grid$pts[,2])
# colnames(yRows)<-"Y_Rows"
codebook <- cbind( code.grid,gph500.codebook )
codebook.long <- melt( codebook,id.vars=1 )





## BLOCK START
#  TESTING
codebook.long.mapunit <- codebook.long[ 1:( ncols*nrows ), ]

#  original
#codebook.long <- separate( codebook.long,
#                           variable,
#                           convert=TRUE,
#                           into=c("lat","lon"),
#                           sep="_"
#                           )

#  TESTING
codebook.long <- separate( codebook.long[ ( ncols*nrows+1 ):nrow( codebook.long ), ],
                           variable,
                           convert=TRUE,
                           into=c("lat","lon"),
                           sep="_"
)

## BLOCK END



codebook.long$lat <- as.numeric( gsub( "X","",codebook.long$lat ) )

codebook.long$lon <- codebook.long$lon-360

codebook.long <- separate( codebook.long,
                           codes,
                           convert=FALSE,
                           remove=FALSE,
                           into=c( "xCols","yRows" ),
                           sep="_"
)

# assign days to nodes
#   NOTE THAT THE FUNCTION map() WILL REVERT TO THE 'maps' 
#   LIBRARY VERSION AND PRODUCE AN ERROR
nodes <- map( som.gph500 )

somTime <- as.data.frame( cbind( gph.500.daily.df.wide$datetime,
                                 nodes$unit.classif,
                                 nodes$distances )
)
colnames( somTime ) <- c( "datetime","mapUnit","errorDist" )
somTime <- separate( somTime,
                     datetime,
                     convert=TRUE,
                     remove=FALSE,
                     into=c( "year","month","day","hr" ),
                     sep="_"
)

somTime$date <- as.Date( paste( somTime$year,"-",
                                somTime$day,"-",
                                somTime$month,
                                sep=""),
                         format="%Y-%d-%m"
)

somTime$doy <- as.numeric( format( somTime$date,"%j" ) )

somTime$mapUnit <- as.integer( somTime$mapUnit )

somTime$errorDist <- as.numeric( as.character( somTime$errorDist ) )


# plot map - fix lines http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
# plot limits
#xlim = c(-125,-95)
#ylim = c(20,50)
xlim = c(-126.25,-93.75)
ylim = c(18.75,51.25)



#  in ggplot
#  Reintroduce this package here.
library( "maps" )
all_states <- map_data("state")
world<-map_data("world")



colnames(world)<-c("X","Y","PID","POS","region","subregion")
#library(PBSmapping)
world = clipPolys(world, xlim=xlim,ylim=ylim, keepExtra=TRUE)
#colnames(world)<-c("Lon","Lat","PID","POS","region","subregion") 

colnames(all_states)<-c("X","Y","PID","POS","region","subregion")
all_states = clipPolys(all_states, xlim=xlim,ylim=ylim, keepExtra=TRUE)
#colnames(all_states)<-c("Lon","Lat","PID","POS","region","subregion")


#  Mike's version

#p <- ggplot() +
#  geom_polygon( data=world, aes(x=X, y=Y, group = PID),colour="black", fill=NA ) +
#  #scale_x_continuous( breaks = c( -120,-140 ) ) +
#  facet_wrap(~codes, nrow = nrows, ncol=ncols)+theme_bw()+
#  stat_contour( data=codebook.long,aes( lon,lat,z=value,
#                                        color=..level..),
#                size=1 ) +
#  geom_raster( data=codebook.long,aes( x=lon,
#                                       y=lat,
#                                       fill=value ) ) +
#  coord_quickmap()+
#  scale_colour_distiller(palette = "Spectral", name="500mb GPH (m)")
#  #scale_color_continuous(low="blue", high = "red")+ 
#  #coord_map(xlim = xlim,ylim = ylim)+
#  #labs(x="Lon", y="Lat")
#p





library( "extrafont" )
font_import()
y

p <- ggplot() +
  facet_wrap( ~codes,nrow=nrows,ncol=ncols ) +
  #theme_minimal() +
  #geom_raster( data=codebook.long[ 36:nrow(codebook.long ), ],
  geom_raster( data=codebook.long,
               aes( x=lon,
                    y=lat,
                    fill=value ) 
  ) +
  geom_polygon( data=world,
                aes( x=X,
                     y=Y,
                     group=PID ),
                color="black",
                fill=NA,
                size=0.25 ) +
  geom_polygon( data=all_states,
                aes( x=X,
                     y=Y,
                     group=PID ),
                color="black",
                fill=NA,
                size=0.25 ) +
  stat_contour( data=codebook.long[ 36:nrow(codebook.long ), ],
                aes( x=lon,
                     y=lat,
                     z=value ),
                color="black",
                linetype="longdash",
                size=0.25 ) +
  scale_fill_distiller( palette="Spectral" ) +
  coord_quickmap() +
  theme_minimal( base_family="Consolas" ) +
  theme( axis.text=element_text( color="gray30",size=8 ),
         axis.ticks.length=unit( 0.0,"mm" ),
         axis.title=element_text( color="gray30",face="bold",size=8 ),
         legend.key.size=unit( 0.5,"cm" ),
         legend.text=element_text( size=10 ),
         legend.title=element_text( face="bold",
                                    size=10 ),
         panel.grid.major.x=element_line( color="gray30",size=0.3 ),
         panel.grid.major.y=element_line( color="gray30",size=0.3 ),
         strip.background=element_rect( fill="white",color="white",size=0.1 ),
         strip.text=element_text( color="gray30",face="bold",size=8 )
  )
p

ggsave( "./som_4x6_maps_ncep_ncar_r2_southwest_jfm_7818.png",
        plot=p,
        device="png",
        path=NULL,
        scale=1,
        width=10,
        height=7.5,
        units="in",
        dpi=600 )










# summary plots -- appears to plot opposite up/down from SOM plot
counts <- plot(som.gph500, type="counts", shape = "straight", labels=counts)
codes <- plot(som.gph500, type="codes", shape = "straight")
similarities <- plot(som.gph500, type="quality", palette.name = terrain.colors)
plot(som.gph500, type="dist.neighbours", main = "SOM neighbour distances")

# sammon mapping
library(MASS)
gph500.codes <- som.gph500$codes
dis <- dist(as.matrix(som.gph500$codes[[1]]))
gph500.sam <- sammon(dis)
plot(gph500.sam$points, type="n")
text(gph500.sam$points,labels=as.character(1:nrow(code.grid)))




# join codes to node table
somTime<-left_join( x=somTime,
                    y=code.grid,
                    by="mapUnit" )
# plot map units
ggplot(somTime, aes(doy, year)) + 
  geom_tile(aes(fill = mapUnit), colour = "grey") + 
  geom_text(aes(label = codes), size=1.25)+
  scale_fill_gradient2(low = "red", mid = "green",
                       high = "blue", midpoint = 7, space = "Lab",
                       na.value = "grey50", guide = "colourbar")











#  timeseries of map units - for facet wrap
year <- as.vector( 1979:2018 )
ts.facetwrap <- expand.grid( year,code.grid$codes )
colnames( ts.facetwrap ) <- c( "year","codes" )
ts.facetwrap$ann_count <- NA

for ( y in 1979:2018 ) {
  a <- filter( somTime,year==y )
  
  for ( c in code.grid$codes ) {
    b <- filter( a,codes==c )
    ts.facetwrap$ann_count[ which( ts.facetwrap$year==y & ts.facetwrap$codes==c,
                                   arr.ind=FALSE ) ] <- nrow( b )
    rm( b )
  }
  
  rm( a )
}
rm( y )

p <- ggplot() +
  facet_wrap( ~codes,nrow=nrows,ncol=ncols ) +
  geom_line( data=ts.facetwrap,
             aes( x=year,
                  y=ann_count ),
             color="black",
             size=0.75 ) +
  theme_minimal( base_family="Consolas" ) +
  theme( #axis.line=element_line( color="gray30",size=0.3 ),
    #axis.line.x.bottom=element_line( color="gray30",size=0.3 ),
    #axis.line.x.top=element_line( color="gray30",size=0.3 ),
    #axis.line.y.left=element_line( color="gray30",size=0.3 ),
    #axis.line.y.right=element_line( color="gray30",size=0.3 ),
    axis.text=element_text( color="gray30",size=8 ),
    #axis.ticks=element_line( color="gray30",size=0.3 ),
    #axis.ticks.length=unit( 0.0,"mm" ),
    axis.title=element_text( color="gray30",face="bold",size=8 ),
    legend.key.size=unit( 0.5,"cm" ),
    legend.text=element_text( size=10 ),
    legend.title=element_text( face="bold",
                               size=10 ),
    #panel.background=element_rect( fill="gray30" ),
    #panel.border=element_rect( linetype="solid",fill=NA ),
    panel.grid.major.x=element_line( color="gray50",size=0.3 ),
    panel.grid.major.y=element_line( color="gray50",size=0.3 ),
    #panel.grid.major.x=element_blank(),
    #panel.grid.major.y=element_blank(),
    strip.background=element_rect( fill="white",color="white",size=0.1 ),
    strip.text=element_text( color="gray30",face="bold",size=8 )
  )
p

ggsave( "./som_4x6_timeseries_ncep_ncar_r2_southwest_jfm_7818.png",
        plot=p,
        device="png",
        path=NULL,
        scale=1,
        width=10,
        height=7.5,
        units="in",
        dpi=600 )



#  Plot single time series




#  OLD version
#  
#  #  timeseries of map units - individual nodes
year <- as.vector( 1979:2018 )
mapUnit <- as.vector( 1:35 )
ts <- expand.grid( year,mapUnit )
colnames( ts ) <- c( "year","mapUnit" )
ts$ann_count <- NA

for ( y in 1979:2018 ) {
  a <- filter( somTime,year==y )
  
  for ( m in 1:35 ) {
    b <- filter( a,mapUnit==m )
    #if ( nrow( b )==0 ) {
    #  ts$ann_count[ which( ts$year==y & ts$mapUnit==m,
    #                       arr.ind=FALSE ) ] <- 0
    #}
    #else {
    ts$ann_count[ which( ts$year==y & ts$mapUnit==m,
                         arr.ind=FALSE ) ] <- nrow( b )
    #}
    
    rm( b )
  }
  
  rm( a )
}
rm( y )

ts2<-filter( ts,mapUnit==4 )

ggplot( data=ts2,
        aes( x=year,y=ann_count,group=mapUnit ) ) +
  geom_line( aes( color=mapUnit ) )












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
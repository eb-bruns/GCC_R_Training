### Author: Emily Beckman  ###  Date: 01/27/2021
### Funding:
	# The Morton Arboretum,
	#		Institute of Museum and Library Services (award #MA-30-18-0273-18), and
	# 	Botanic Gardens Conservation International
	# Buffer methodology based on work funded by the
	#		USDA Forest Service (Cooperative Agreement 16-CA-11132546-045); see
	# 	Beckman, E., Meyer, A., Pivorunas, D., & Westwood, M. (2021).
	#			Conservation Gap Analysis of Native U.S. Pines. Lisle, IL:
	#			The Morton Arboretum.

### DESCRIPTION: Calculate geographic and ecological coverage of ex situ
	# collections, and create a map of coverage

### INPUTS:
	## Ecoregions
	# 	Terrestrial Ecoregions of the World, via WWF (Olson et al., 2001)
	#			https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
	## Global country boundaries
	#		UIA World Countries Boundaries, UNIGIS Geospatial Education Resources, via ArcGIS Hub
	# 		Shapefile
	#			https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0
	## In situ occurrence points (latitude and longitude in decimal degrees)
	# 	Can use the output from 3-1_refine_occurrence_points.R
	# 		https://github.com/MortonArb-CollectionsValue/OccurrencePoints/tree/master/scripts
	# 		Each file has data for only one species and is named "Genus_species.csv"
	# 		You can read in data for mult. species in one file but need to edit the
	#			code to split after reading in
	## Ex situ wild localities (latitude and longitude in decimal degrees)
	# 	Can use the output from 3-1_refine_occurrence_points.R, which has a
	#			"database" column that has "Ex_situ" to distinguish the ex situ records
	#			from the rest of the in situ records

### OUTPUTS:
  ## Table with geographic and ecological coverage for each species; example:
	# | species						|	geo_20								|	eco_20					|	geo_50	|	eco_50	|	geo_100	|	eco_100	|	geo_avg	| eco_avg	|
	# |-------------------|-----------------------|-----------------|---------|---------|---------|---------|---------|---------|
	# | Taxus brevifolia	|	19,793 / 304,644 (6%)	|	53 / 161 (33%)	| ...			|					|					|					|	23%			|	47%			|
	# | Taxus canadensis	|	14,494 / 368,449 (4%)	|	...           	|     		|					|					|					|					|					|
	# ...
	## Interactive map with 50 km buffers (can change buffer size as desired)

################################################################################

#################
### LIBRARIES ###
#################

my.packages <- c("leaflet","raster","sp","rgeos","plyr","dplyr","rgdal",
	"Polychrome","cleangeo","RColorBrewer","smoothr")
#install.packages(my.packages) # turn on to install current versions
lapply(my.packages, require, character.only=TRUE)

select <- dplyr::select

#########################
### WORKING DIRECTORY ###
#########################

## set up working directories
	## for point data (in situ and ex situ)
pts_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/points"
	## for polygon data (ecoregions, states, countries)
poly_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/polygons"
	## for outputs
output_dir <- "/Volumes/GoogleDrive/My Drive/Conservation Consortia/R Training/outputs"

#################
### FUNCTIONS ###
#################

# create buffers around points, using specified projection
create.buffers <- function(df,radius,pt_proj,buff_proj,boundary){
	# select coordinate columns
	latlong <- df %>% select(decimalLongitude,decimalLatitude)
	# turn occurrence point data into a SpatialPointsDataFrame
	sp_df <- SpatialPointsDataFrame(latlong, df, proj4string = pt_proj)
	# reproject SpatialPointsDataFrame to specified projection
	proj_df <- spTransform(sp_df,buff_proj)
	# place buffer around each point
	buffers <- buffer(proj_df,width=radius,dissolve=T)
	# clip buffers by boundary (e.g., state, country)
	buffers_clip <- raster::intersect(buffers,boundary)
	# return buffer polygons
	return(buffers_clip)
}

# create buffers around points and calculate area
calc.buff.area <- function(pts,radius,pt_proj,buff_proj,boundary){
	# create buffers
	buffers <- create.buffers(pts,radius,pt_proj,buff_proj,boundary)
	# calculate buffer area
	buff_area <- buffers@polygons[[1]]@area/1000000
	#print(paste("Area covered by buffers:", round(buff_area,0),"km²"))
	return(buff_area)
}

# create data frame with ecoregion data extracted for area covered by buffers,
#		for both in situ and ex situ points, then compare count of ecoregions
count.eco <- function(pts,radius,pt_proj,buff_proj,ecoregions,boundary){
	# create buffers
	buffers <- create.buffers(pts,radius,pt_proj,buff_proj,boundary)
	# make sure ecoregions are in same projection as buffers
	eco_proj <- spTransform(ecoregions,buff_proj)
	# intersect buffers with ecoregions
	buff_join_eco <- raster::intersect(buffers,eco_proj)
	# count number of ecoregions under buffers
	count_eco <- nrow(buff_join_eco@data %>% distinct(ECO_ID))
	#print(paste0("Number of ecoregions under ex situ buffers: ",count_eco))
	return(count_eco)
}

# format text in cell for output table
format.cell <- function(ex_result,in_result,final_result){
	cell <- paste0(format(round(ex_result,0),format="d",big.mark=",")," / ",
								 format(round(in_result,0),format="d",big.mark=","),"\n",
								 "(",round(final_result,0),"%)")
	return(cell)
}

# clip points by boundary so only in target area
# (helpful if focusing on one country/region)
clip.by.boundary <- function(pts,pt_proj,boundary){
	# select coordinate columns
	latlong <- pts %>% select(decimalLongitude,decimalLatitude)
	# turn occurrence point data into a SpatialPointsDataFrame
	spatial_pts <- SpatialPointsDataFrame(latlong, pts, proj4string = pt_proj)
	# clip by boundary created earlier
	spatial_pts <- spatial_pts[boundary, ]
	# keep just the data (not the spatial info you added)
	pts_new <- spatial_pts@data
	return(pts_new)
}

################################################################################
################################################################################
## Set things up
################################################################################

### DEFINE PROJECTIONS / COORDINATE REFERENCE SYSTEM (CRS)

## define initial projection of points (usually WGS 84); also used when creating
##	 leaflet map
wgs.proj <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84
	+no_defs +towgs84=0,0,0")
## define projection for calculations (meters/km must be the unit)
aea.proj <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-110
	+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m")

### READ IN POLYGON DATA

## Countries
world_countries <- readOGR(file.path(poly_dir,"UIA_World_Countries_Boundaries-shp/World_Countries__Generalized_.shp"))
	## filter to only target countries; speeds things up and prevents errors.
	##	there is a self-intersection error when trying the aggregate function for
	##	the aea projection using all countries; tried clgeo_Clean and did not fix
sort(unique(world_countries@data$ISO))
target_iso <- c("CN","JP")
target_countries <- world_countries[world_countries@data$ISO %in% target_iso,]
	## create polygon for clipping buffers later, one in each projection
target_countries.wgs <- spTransform(target_countries,wgs.proj)
boundary.wgs <- aggregate(target_countries.wgs,dissolve = TRUE)
target_countries.aea <- spTransform(target_countries,aea.proj)
	## this is where the error occurs with certain countries.. may need to find
	##	a work-around
boundary.aea <- aggregate(target_countries.aea,dissolve = TRUE)

## Ecoregions
ecoregions <- readOGR(file.path(poly_dir,"official/wwf_terr_ecos.shp"))
#length(unique(ecoregions@data$ECO_ID)) #827
	## project the ecoregions for mapping and calculations
ecoregions.wgs <- spTransform(ecoregions,wgs.proj)
ecoregions.aea <- spTransform(ecoregions,aea.proj)
	## you can also clip the ecoregions to just be target countries (maps faster)
ecoregions_clip.wgs <- raster::intersect(ecoregions.wgs,boundary.wgs)
ecoregions_clip.aea <- raster::intersect(ecoregions.aea,boundary.aea)

### CREATE COLOR PALETTES / MAP ICONS

## Ecoregions polygon colors
eco_pal_colors <- createPalette(length(unique(ecoregions_clip.wgs@data$ECO_ID)),
	seedcolors = c("#ba3c3c","#ba7d3c","#baab3c","#3ca7ba","#3c6aba","#573cba","#943cba","#ba3ca1","#ba3c55"),
	range = c(5,42), target = "normal", M=50000)
swatch(eco_pal_colors)
eco_pal_colors <- as.vector(eco_pal_colors)
eco_pal <- colorFactor(eco_pal_colors,ecoregions_clip.wgs@data$ECO_ID)

## Ex situ point data triangle icons
triangle_sm <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/triangle-png-28.png",
 	iconWidth = 8, iconHeight = 8)
triangle_md <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/triangle-png-28.png",
 	iconWidth = 15, iconHeight = 15)
triangle_lg <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/triangle-png-28.png",
 	iconWidth = 22, iconHeight = 22)

################################################################################
## Create interactive leaflet map
################################################################################

### CREATE LIST OF TARGET SPECIES

target_sp <- c("Acer_griseum","Acer_pycnanthum")
## select species to work with now
sp <- 1

### READ IN AND PREP POINT DATA

## read in wild in situ occurrence points
insitu <- read.csv(file.path(pts_dir,paste0(target_sp[sp],
	"_distribution_points_2020.csv")),na.strings=c("","NA"),
	stringsAsFactors = F)
str(insitu)
## change column names or remove columns as needed; need at least
##	"decimalLatitude" and "decimalLongitude"
insitu <- insitu %>%
	rename(decimalLatitude = Latitude, decimalLongitude = Longitude) %>%
	filter(Category == "in_situ")
## if desired, can clip points by boundary so only in target area
## (helpful if focusing on one country/region)
insitu <- clip.by.boundary(insitu,wgs.proj,boundary.wgs)
str(insitu)

## read in ex situ wild locality points
exsitu <- read.csv(file.path(pts_dir,paste0(target_sp[sp],
	"_ALL_POSTGEO.csv")),na.strings=c("","NA"),stringsAsFactors = F)
str(exsitu)
## change column names or remove columns as needed; should have at least
##	"decimalLatitude","decimalLongitude","num_indiv"
exsitu <- exsitu %>%
	rename(decimalLatitude = lat_dd, decimalLongitude = long_dd) %>%
	select(UID,inst_short,prov_type,gps_det,decimalLatitude,decimalLongitude,
		num_indiv)
str(exsitu)
unique(exsitu$num_indiv)
## remove points that are very close together, so can see mapped better
## you can change number of digits to remove less
exsitu$lat_round <- round(as.numeric(exsitu$decimalLatitude),digits=2)
exsitu$long_round <- round(as.numeric(exsitu$decimalLongitude),digits=2)
exsitu <- exsitu %>%
	filter(!is.na(lat_round)) %>%
	group_by(lat_round,long_round) %>%
	mutate(num_indiv = sum(as.numeric(num_indiv))) %>%
	ungroup() %>%
		## add things here if you want extra columns preserved; !! could also
		##	collapse institution names so they are all listed in popup
	distinct(lat_round,long_round,.keep_all=T) %>%
	select(-lat_round,-long_round)
## if desired, can clip points by boundary so only in target area
## (helpful if focusing on one country/region)
exsitu <- clip.by.boundary(exsitu,wgs.proj,boundary.wgs)
str(exsitu)
## optionally, remove any bad ex situ points manually by ID number
#if(exsitu$species_name_acc[1] == "Juglans californica"){
#	exsitu <- exsitu %>% filter(UID != "id00092049" & UID != "id00091728")
#} else if (exsitu$species_name_acc[1] == "Juglans major"){
#	exsitu <- exsitu %>% filter(UID != "id00092280")
#}
#sum(exsitu$num_indiv)
## split by number of individuals, to use different sized symbol for each
## change as needed to get categories you want
exsitu1 <- exsitu %>% arrange(num_indiv) %>% filter(num_indiv <= 1)
exsitu2 <- exsitu %>% arrange(num_indiv) %>% filter(num_indiv > 2 & num_indiv < 4)
exsitu3 <- exsitu %>% arrange(num_indiv) %>% filter(num_indiv >= 5)

## add ex situ points to in situ points
insitu <- rbind.fill(insitu,exsitu)

### CREATE MAP

## create buffers for visualization
## change number to get different buffer size; e.g., 20000 = 20km
insitu_buff_wgs <- create.buffers(insitu,20000,wgs.proj,wgs.proj,boundary.wgs)
exsitu_buff_wgs <- create.buffers(exsitu,20000,wgs.proj,wgs.proj,boundary.wgs)

## map everthing!
	## can turn layers on or off, or switch them for other polygons, as desired
map <- leaflet(options = leafletOptions(maxZoom = 9)) %>%
  ## Base layer
  addProviderTiles(providers$CartoDB.VoyagerNoLabels) %>%
	## Species name label
	addControl(paste0("<b>",gsub("_"," ",target_sp[sp])), position = "topright") %>%
	## Ecoregions
	addPolygons(
		data = ecoregions_clip.wgs, label = ~ECO_NAME,
		fillColor = ~eco_pal(ecoregions_clip.wgs@data$ECO_ID),
		fillOpacity = 0.9, color = "#757575", weight = 1.5, opacity = 0.8) %>%
	## Country outlines
	#addPolygons(
	#	data = target_countries.wgs, label = ~COUNTRY,
	#	fillColor = "transparent", weight = 1.5, opacity = 0.5, color = "black") %>%
	## Buffers
		## In situ
	addPolygons(
		data = insitu_buff_wgs,
		fillColor = "#a3a3a3", fillOpacity = 0.45,
		weight = 1.3, opacity = 0.9, color = "#c4c4c4",
		smoothFactor = 0) %>%
		## Ex situ
	addPolygons(
		data = exsitu_buff_wgs,
		fillColor = "white", fillOpacity = 0.55,
		weight = 1.3, color = "white", opacity = 0,
		smoothFactor = 0) %>%
	## Ex situ points
  addMarkers(data = exsitu1,
		lng = ~decimalLongitude, lat = ~decimalLatitude, icon = triangle_sm,
		popup = exsitu1$inst_short) %>%
	addMarkers(data = exsitu2,
		lng = ~decimalLongitude, lat = ~decimalLatitude, icon = triangle_md,
		popup = exsitu2$inst_short) %>%
  addMarkers(data = exsitu3,
		lng = ~decimalLongitude, lat = ~decimalLatitude, icon = triangle_lg,
		popup = exsitu3$inst_short) %>%
	## (optional) In situ points
	addCircleMarkers(data = insitu,
		lng = ~decimalLongitude, lat = ~decimalLatitude,
		color = "white", radius = 3, fillOpacity = 1, stroke = F) %>%
	## Add scale bar and set view
	addScaleBar(position = "bottomright",
		options = scaleBarOptions(maxWidth = 150)) %>%
	setView(104, 32, zoom = 5)
map

## save map
		# looks like these maps are too big to save? works best to view
		#		one-by-one in browser straight from R and take screenshot
htmlwidgets::saveWidget(map, file.path(output_dir,paste0(target_sp[sp],"_leaflet_map.html")))

################################################################################
## Calculate geographic and ecological coverage of ex situ collections
################################################################################

### START SUMMARY TABLE

# we add each target species as we go along
summary_tbl <- data.frame(species = "start",
	geo_20 = "start", eco_20 = "start",
	geo_50 = "start", eco_50 = "start",
	geo_100 = "start", eco_100 = "start",
	geo_avg = "start", eco_avg = "start", stringsAsFactors=F)

### CYCLE THROUGH TARGET SPECIES TO CALCULATE EX SITU COVERAGE

for(sp in 1:length(target_sp)){

	# print progress
  cat("\tStarting ", target_sp[sp])

	### READ IN AND PREP POINT DATA

	## RIGHT NOW THIS IS JUST COPIED FROM ABOVE; WAY TO SHORTEN/NOT REPEAT?
	## read in wild in situ occurrence points
	insitu <- read.csv(file.path(pts_dir,paste0(target_sp[sp],
		"_distribution_points_2020.csv")),na.strings=c("","NA"),
		stringsAsFactors = F)
	## change column names or remove columns as needed; need at least
	##	"decimalLatitude" and "decimalLongitude"
	insitu <- insitu %>%
		rename(decimalLatitude = Latitude, decimalLongitude = Longitude) %>%
		filter(Category == "in_situ")
	## if desired, can clip points by boundary so only in target area
	## (helpful if focusing on one country/region)
	insitu <- clip.by.boundary(insitu,wgs.proj,boundary.wgs)

	## read in ex situ wild locality points
	exsitu <- read.csv(file.path(pts_dir,paste0(target_sp[sp],
		"_ALL_POSTGEO.csv")),na.strings=c("","NA"),stringsAsFactors = F)
	## change column names or remove columns as needed; should have at least
	##	"decimalLatitude","decimalLongitude","num_indiv"
	exsitu <- exsitu %>%
		rename(decimalLatitude = lat_dd, decimalLongitude = long_dd) %>%
		select(UID,inst_short,prov_type,gps_det,decimalLatitude,decimalLongitude,
			num_indiv)
	unique(exsitu$num_indiv)
	## remove points that are very close together, so can see mapped better
	## you can change number of digits to remove less
	exsitu$lat_round <- round(as.numeric(exsitu$decimalLatitude),digits=2)
	exsitu$long_round <- round(as.numeric(exsitu$decimalLongitude),digits=2)
	exsitu <- exsitu %>%
		filter(!is.na(lat_round)) %>%
		group_by(lat_round,long_round) %>%
		mutate(num_indiv = sum(as.numeric(num_indiv))) %>%
		ungroup() %>%
			## add things here if you want extra columns preserved; !! could also
			##	collapse institution names so they are all listed in popup
		distinct(lat_round,long_round,.keep_all=T) %>%
		select(-lat_round,-long_round)
	## if desired, can clip points by boundary so only in target area
	## (helpful if focusing on one country/region)
	exsitu <- clip.by.boundary(exsitu,wgs.proj,boundary.wgs)
	## optionally, remove any bad ex situ points manually by ID number
	#if(exsitu$species_name_acc[1] == "Juglans californica"){
	#	exsitu <- exsitu %>% filter(UID != "id00092049" & UID != "id00091728")
	#} else if (exsitu$species_name_acc[1] == "Juglans major"){
	#	exsitu <- exsitu %>% filter(UID != "id00092280")
	#}
	#sum(exsitu$num_indiv)

	## add ex situ points to in situ points
	insitu <- rbind.fill(insitu,exsitu)

	### CALCULATE EX SITU COVERAGE

	## Geographic coverage
		## 20 km buffers
	geo_20_insitu <- calc.buff.area(insitu,20000,wgs.proj,aea.proj,boundary.aea)
	geo_20_exsitu <- calc.buff.area(exsitu,20000,wgs.proj,aea.proj,boundary.aea)
	geo_20_percent <- (geo_20_exsitu/geo_20_insitu)*100
		## 50 km buffers
	geo_50_insitu <- calc.buff.area(insitu,50000,wgs.proj,aea.proj,boundary.aea)
	geo_50_exsitu <- calc.buff.area(exsitu,50000,wgs.proj,aea.proj,boundary.aea)
	geo_50_percent <- (geo_50_exsitu/geo_50_insitu)*100
		## 100 km buffers
	geo_100_insitu <- calc.buff.area(insitu,100000,wgs.proj,aea.proj,boundary.aea)
	geo_100_exsitu <- calc.buff.area(exsitu,100000,wgs.proj,aea.proj,boundary.aea)
	geo_100_percent <- (geo_100_exsitu/geo_100_insitu)*100

	## Ecological coverage (use count.eco.wwf function if using global ecoregions)
		## 20 km buffers
	eco_20_insitu <- count.eco(insitu,20000,wgs.proj,aea.proj,ecoregions.aea,boundary.aea)
	eco_20_exsitu <- count.eco(exsitu,20000,wgs.proj,aea.proj,ecoregions.aea,boundary.aea)
	eco_20_percent <- (eco_20_exsitu/eco_20_insitu)*100
		## 50 km buffers
	eco_50_insitu <- count.eco(insitu,50000,wgs.proj,aea.proj,ecoregions.aea,boundary.aea)
	eco_50_exsitu <- count.eco(exsitu,50000,wgs.proj,aea.proj,ecoregions.aea,boundary.aea)
	eco_50_percent <- (eco_50_exsitu/eco_50_insitu)*100
		## 100 km buffers
	eco_100_insitu <- count.eco(insitu,100000,wgs.proj,aea.proj,ecoregions.aea,boundary.aea)
	eco_100_exsitu <- count.eco(exsitu,100000,wgs.proj,aea.proj,ecoregions.aea,boundary.aea)
	eco_100_percent <- (eco_100_exsitu/eco_100_insitu)*100

	## Create text for results table
	geo_20_txt <- format.cell(geo_20_exsitu,geo_20_insitu,geo_20_percent)
	geo_50_txt <- format.cell(geo_50_exsitu,geo_50_insitu,geo_50_percent)
	geo_100_txt <- format.cell(geo_100_exsitu,geo_100_insitu,geo_100_percent)
	geo_avg_txt <- paste0(round(mean(c(geo_20_percent,geo_50_percent,geo_100_percent)),0),"%")
	eco_20_txt <- format.cell(eco_20_exsitu,eco_20_insitu,eco_20_percent)
	eco_50_txt <- format.cell(eco_50_exsitu,eco_50_insitu,eco_50_percent)
	eco_100_txt <- format.cell(eco_100_exsitu,eco_100_insitu,eco_100_percent)
	eco_avg_txt <- paste0(round(mean(c(eco_20_percent,eco_50_percent,eco_100_percent)),0),"%")

  ## Add text results to summary table
  summary_add <- data.frame(species = target_sp[sp],
			geo_20 = geo_20_txt, eco_20 = eco_20_txt,
			geo_50 = geo_50_txt, eco_50 = eco_50_txt,
			geo_100 = geo_100_txt, eco_100 = eco_100_txt,
			geo_avg = geo_avg_txt, eco_avg = eco_avg_txt, stringsAsFactors=F)
  summary_tbl[sp,] <- summary_add

}

## write summary table
summary_tbl
write.csv(summary_tbl, file.path(output_dir,"ExSituCoverage_Test_Table.csv"),
	row.names = F)

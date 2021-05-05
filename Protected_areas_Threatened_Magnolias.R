### Author: Jean Linsky  ###  Date: 04/7/2021
### Adapted From:
	# Emily Beckman et al. global_insitu_versus_exsitu_buffers script

### DESCRIPTION: Create a map of Protected areas and in situ points of threatened Magnolia

### INPUTS:
	## Protected Areas
	# 	Protected areas via ProtectedPlanet, UNEP-WCMC and IUCN (2021)
	#			https://www.protectedplanet.net/en
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

### OUTPUTS:
	## Interactive map with protected areas and insitu points visualized

################################################################################

#################
### LIBRARIES ###
#################

rm(list=ls())
my.packages <- c("leaflet","raster","sp","rgeos","plyr","dplyr","rgdal",
	"Polychrome","cleangeo","RColorBrewer","smoothr","rnaturalearth","polylabelr",
	"sf")
#install.packages(my.packages) # turn on to install current versions
lapply(my.packages, require, character.only=TRUE)

select <- dplyr::select

#########################
### WORKING DIRECTORY ###
#########################

## set up working directories
	## for point data (in situ)
pts_dir <- "C:/Users/Jean Linsky/Documents/Magnolia_Coordinator/Statistics_and_R/outputs/spp_raw_points"
	## for polygon data (ecoregions, states, countries)
poly_dir <- "C:/Users/Jean Linsky/Documents/Magnolia_Coordinator/Statistics_and_R/Magnolia/polygons"
  ## PA directory
pa_dir <- "C:/Users/Jean Linsky/Documents/Magnolia_Coordinator/GIS/Protected Areas Maps/China_PAs"
	## for outputs
output_dir <- "C:/Users/Jean Linsky/Documents/Magnolia_Coordinator/Statistics_and_R/Magnolia/Protected Area Maps for publication/Test_R_Maps"

#################
### FUNCTIONS ###
#################

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
##	leaflet map
wgs.proj <- sp::CRS(SRS_string="EPSG:4326")
	##CRS arguments: +proj=longlat +datum=WGS84 +no_defs
## define projection for calculations (meters/km must be the unit); this one
##	works best if you use projection specifically for your target region;
## 	you can search for projections and their EPSG codes here: https://epsg.org
## FOR ASIA/PACIFIC: 8859; FOR THE AMERICAS: 8858; FOR EUROPE/AFRICA: 8857;
##	FOR THE U.S. ONLY, if you want to align with USGS preference: 5070
aea.proj <- sp::CRS(SRS_string="EPSG:8859")
	##CRS arguments: +proj=eqearth +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs

### READ IN POLYGON DATA

## Countries
world_countries <- readOGR(file.path(poly_dir,"UIA_World_Countries_Boundaries-shp/World_Countries__Generalized_.shp"))
	## filter to only target countries; speeds things up and prevents errors.
	##	there is a self-intersection error when trying the aggregate function for
	##	the aea projection using all countries; tried clgeo_Clean and did not fix
sort(unique(world_countries@data$ISO))
	## Look up country codes at website below, using "Alpha 2" column:
	##	https://www.nationsonline.org/oneworld/country_code_list.htm
target_iso <- c("CN")
target_countries <- world_countries[world_countries@data$ISO %in% target_iso,]
	## create polygon for clipping buffers later, one in each projection
target_countries.wgs <- spTransform(target_countries,wgs.proj)
boundary.wgs <- aggregate(target_countries.wgs,dissolve = TRUE)
target_countries.aea <- spTransform(target_countries,aea.proj)
	## this is where the error occurs with certain countries.. may need to find
	##	a work-around
boundary.aea <- aggregate(target_countries.aea,dissolve = TRUE)

##Protected Areas
protected_areas0 <- readOGR(file.path(pa_dir,"WDPA_WDOECM_Apr2021_Public_CHN_shp-polygons0.shp"))
##clipping PAs to land only
protected_areas0_clip.wgs <- raster::intersect(protected_areas0,boundary.wgs)
protected_areas0_clip.aea <- raster::intersect(protected_areas0,boundary.aea)

protected_areas1 <- readOGR(file.path(pa_dir,"WDPA_WDOECM_Apr2021_Public_CHN_shp-polygons1.shp"))
##clipping PAs to land only
protected_areas1_clip.wgs <- raster::intersect(protected_areas1,boundary.wgs)
protected_areas1_clip.aea <- raster::intersect(protected_areas1,boundary.aea)

protected_areas2 <- readOGR(file.path(pa_dir,"WDPA_WDOECM_Apr2021_Public_CHN_shp-polygons2.shp"))
##clipping PAs to land only
protected_areas2_clip.wgs <- raster::intersect(protected_areas2,boundary.wgs)
protected_areas2_clip.aea <- raster::intersect(protected_areas2,boundary.aea)

VNMprotected_areas0 <- readOGR(file.path(pa_dir,"WDPA_WDOECM_Apr2021_Public_VNM_shp-polygons0.shp"))
##clipping PAs to land only
VNMprotected_areas0_clip.wgs <- raster::intersect(VNMprotected_areas0,boundary.wgs)
VNMprotected_areas0_clip.aea <- raster::intersect(VNMprotected_areas0,boundary.aea)

VNMprotected_areas1 <- readOGR(file.path(pa_dir,"WDPA_WDOECM_Apr2021_Public_VNM_shp-polygons1.shp"))
##clipping PAs to land only
VNMprotected_areas1_clip.wgs <- raster::intersect(VNMprotected_areas1,boundary.wgs)
VNMprotected_areas1_clip.aea <- raster::intersect(VNMprotected_areas1,boundary.aea)

VNMprotected_areas2 <- readOGR(file.path(pa_dir,"WDPA_WDOECM_Apr2021_Public_VNM_shp-polygons2.shp"))
##clipping PAs to land only
VNMprotected_areas2_clip.wgs <- raster::intersect(VNMprotected_areas2,boundary.wgs)
VNMprotected_areas2_clip.aea <- raster::intersect(VNMprotected_areas2,boundary.aea)

## States
	## read in state polygons using rnaturalearth package (or can read in other shapefile
	##		you have downloaded online)
state_bound <- ne_states(country=NULL)
	## project to WGS84
state_bound.wgs <- spTransform(state_bound,wgs.proj)
	## if desired, select states in target countries only
state_bound_clip.wgs <- state_bound.wgs[state_bound.wgs@data$iso_a2 %in% target_iso,]
	## find the visual center point of each state (still not perfect), for
	##		labeling purposes
simple_states <- st_as_sf(state_bound_clip.wgs)
state_centers <- as.data.frame(do.call(rbind, poi(simple_states, precision=0.01)))
state_centers$label <- state_bound_clip.wgs@data$name_en
state_centers$x <- as.numeric(state_centers$x)
state_centers$y <- as.numeric(state_centers$y)
	## (optional) edit state/province names to better format for map labels
	##	 remove word "Prefecture" in labels
state_centers$label <- gsub(" Prefecture","",state_centers$label)
	##	 remove words "Autonomous Region" in labels
state_centers$label <- gsub(" Autonomous Region","",state_centers$label)
	## 	 view labels
state_centers

### CREATE COLOR PALETTES / MAP ICONS

## PA polygon colors
pa_pal_colors <- createPalette(length(unique(protected_areas0@data$ECO_ID)),
	seedcolors = c("#ba3c3c","#ba7d3c","#baab3c","#3ca7ba","#3c6aba","#573cba","#943cba","#ba3ca1","#ba3c55"),
	range = c(5,42), target = "normal", M=50000)
swatch(pa_pal_colors)
pa_pal_colors <- as.vector(pa_pal_colors)
pa_pal <- colorFactor(pa_pal_colors,protected_areas0@data$ECO_ID)


################################################################################
## Create interactive leaflet map
################################################################################

### CREATE LIST OF TARGET SPECIES

target_sp <- c("Magnolia_xanthantha")
## select species to work with now
sp <- 1

### READ IN AND PREP POINT DATA

## read in wild in situ occurrence points
insitu <- read.csv(file.path(pts_dir,paste0(target_sp[sp],
	"_1.csv")),na.strings=c("","NA"),
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


### CREATE MAP


## map everything!
	## can turn layers on or off, or switch them for other polygons, as desired
map <- leaflet(options = leafletOptions(maxZoom = 9)) %>%
  ## Base layer
	##	 explore other base layer options here:
	##	 http://leaflet-extras.github.io/leaflet-providers/preview/index.html
  addProviderTiles(providers$CartoDB.VoyagerNoLabels) %>%
	## Species name label
	addControl(paste0("<b>",gsub("_"," ",target_sp[sp])),
		position = "topright") %>%
	## ProtectedAreas
	addPolygons(
		data = protected_areas0_clip.wgs, label = ~NAME,
		fillColor = ~pa_pal(protected_areas0_clip.wgs@data$ECO_ID),
		fillOpacity = 0.8, color = "#038f28", weight = 1.5, opacity = 0.8) %>%
  addPolygons(
    data = protected_areas1_clip.wgs, label = ~NAME,
    fillColor = ~pa_pal(protected_areas0_clip.wgs@data$ECO_ID),
    fillOpacity = 0.8, color = "#038f28", weight = 1.5, opacity = 0.8) %>%
  addPolygons(
    data = protected_areas2_clip.wgs, label = ~NAME,
    fillColor = ~pa_pal(protected_areas0@data$ECO_ID),
    fillOpacity = 0.8, color = "#038f28", weight = 1.5, opacity = 0.8) %>%
	## (optional) Country or state outlines
	##	when you add these outlines, ecoregion labels don't pop up anymore...
	##	not sure yet how to have both on the map
	addPolygons(
		data = state_bound_clip.wgs, label = ~name_en, fillColor = "transparent",
		weight = 1.5, opacity = 0.3, color = "black") %>%
		## (optional) Add static labels to countries/states
	addLabelOnlyMarkers(
		data = state_centers, lng = ~x, lat = ~y, label = ~label,
    labelOptions = labelOptions(noHide = TRUE, textOnly = TRUE,
			style = list("font-weight"="bold","font-size"="13px","color"="black"))) %>%
	## (optional) In situ points
	addCircleMarkers(data = insitu,
		lng = ~decimalLongitude, lat = ~decimalLatitude,
		color = "black", radius = 5, fillOpacity = 1, stroke = F) %>%
	## Add scale bar
	addScaleBar(position = "bottomright",
		options = scaleBarOptions(maxWidth = 150)) %>%

  ##Need to see if a legend with PAs, points can be added
  	## Add legend
	##	not perfect, but something! Used https://imgbb.com to host the buffer
	##	PNG images! So you could do that for any shape you'd like
	addControl(
		html = "<img src='https://i.ibb.co/1dW95pC/Insitu-buffer.png'
		style='width:40px;height:40px;'> Species' estimated native distribution<br/>
		(20 km buffer around in situ occurrence points)<br/>
		<img src='https://i.ibb.co/SR71N6k/Exsitu-buffer.png'
		style='width:40px;height:40px;'> Estimated capture of ex situ collections<br/>
		(20 km buffer around wild provenance localities)",
		position = "bottomleft") %>%
	addControl(
		html = "Source locality and number of wild provenance<br/>individuals in ex situ collections<br/>
		<img src='https://www.freeiconspng.com/uploads/triangle-png-28.png'
		style='width:8px;height:8px;'> 1-4
		<img src='https://www.freeiconspng.com/uploads/triangle-png-28.png'
		style='width:15px;height:15px;'> 5-14
		<img src='https://www.freeiconspng.com/uploads/triangle-png-28.png'
		style='width:22px;height:22px;'> 15+",
		position = "bottomleft") %>%
	## Set view (long and lat) and zoom level, for when map initially opens
	setView(104, 32, zoom = 5)
map

## save map as html file, so you can embed on a webpage or share with others
		# looks like some of these maps are too big to save? works best to view
		#		one-by-one in browser straight from R and take screenshot
htmlwidgets::saveWidget(map, file.path(output_dir,paste0(target_sp[sp],"_leaflet_map.html")))


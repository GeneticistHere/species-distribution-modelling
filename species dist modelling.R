setwd('/Users/macbook/Downloads/modeling')

library(dismo)

# library(rgbif)

sp <- gbif("Pistacia","vera",download = T, geo=T,sp=F)
sp
# Save the GBIF data to a CSV file
write.csv(sp, file = "Pistacia_vera_gbif_data.csv", row.names = FALSE)

# Verify the file
print("File saved as 'Pistacia_vera_gbif_data.csv'")

class(sp)
str(sp)

sp <- sp[,c('lon','lat')]
head(sp)
w <- which(is.na(sp$lon))
sp <- sp[-w,]
w <- which(is.na(sp$lat))

sp$species <- 1
head(sp)

coordinates(sp) <- ~ lon + lat
class(sp)

plot(sp)
#### Download Bioclim data for the current time

install.packages("geodata")
library(geodata)

bio <- raster::getData('worldclim', var='bio', res=10)
# Specify the path to save the data (e.g., "data" folder in your working directory)
bio <- geodata::worldclim_global(var = "bio", res = 10, path = "/Users/macbook/Downloads/modeling")

plot(bio[[1]])
names(bio)
#############extra
v1<-vifstep(bio)
v2<-vifstep(bio)
biom<-exclude(bio,v2)
plot(biom[[1]])
############
points(sp)
# Download CMIP6 bioclimatic data
# Get CMIP6 data
# Note: CMIP6 is the newer version of CMIP5
biof <- cmip6_world(
  model = "ACCESS-ESM1-5",  # This is equivalent to the AC model
  ssp = "585",             # This is equivalent to RCP8.5
  time = "2061-2080",      # This is equivalent to year 70
  var = "bio",
  res = 10,
  path = "/Users/macbook/Downloads/modeling"
)
# Check the downloaded data
print(biof)

plot(biof[[1]],col=bpy.colors(100))

names(biof) <- names(bio)
#--------
head(sp)
#---- remove collinear variables

# Convert SpatialPointsDataFrame to terra vect object
library(usdm)
sp_vect <- vect(sp)
spx <- extract(bio, sp_vect)

head(spx)
class(spx)
spx <- data.frame(spx)

########
library(sdm)
library(terra)
library(raster)

# Convert SpatRaster to RasterStack
bio_stack <- stack(bio)
d <- sdmData(species~., train=sp, predictors=bio_stack, bg=list(n=1000))
d
m <- sdm(species ~ . , d, methods=c('glm','brt','rf','svm','bioclim.dismo'),
         replication='sub', test.p=30)

m

gui(m)

en <- ensemble(m, bio, filename = 'ens_current.img', setting=list(method='weighted',stat='AUC'))

enf <- ensemble(m, biof, filename = 'ens_future.img', setting=list(method='weighted',stat='AUC'))

library(mapview)

plot(en)
plot(enf)

# Convert SpatRaster to raster stack
en_raster <- raster(en)
enf_raster <- raster(enf)

# Now create the stack and view it
combined_stack <- stack(en_raster, enf_raster)
mapview(combined_stack)
mapview(stack(en, enf))

# Convert the SpatRaster layers to raster format
bio13_raster <- raster(bio[[13]])
biof12_raster <- raster(biof[[12]])

# Create the stack and view it
combined_stack <- stack(bio13_raster, biof12_raster)
mapview(combined_stack)

#------

ch <- enf - en
cl <- colorRampPalette(c('red','white','green','darkblue'))
plot(ch,col=cl(100))

#--------------------

th <- getEvaluation(m,stat='threshold')
mean(th[,2])

df <- data.frame(as.data.frame(d),coordinates(d))
head(df)

pr <- extract(en, df[,c('lon','lat')])

head(pr)
ev <- evaluates(df$species, pr)
ev@threshold_based

th <- 0.51247

pa <- raster(en)

pa[] <- ifelse(en[] >= 0.51247, 1, 0)
plot(pa)

paf <- raster(enf)

paf[] <- ifelse(enf[] >= 0.51247, 1, 0)

plot(paf)

pa.ch <- paf - pa
plot(pa.ch)

cl <- colorRampPalette(c('red','gray80','darkgreen'))
plot(pa.ch,col=cl(3))

##########
# Create a directory for plots if it doesn't exist
dir.create("pistacia_plots", showWarnings = FALSE)

# 1. Current Species Occurrences with Bioclimatic Data
png("pistacia_plots/1_current_occurrences.png", width=1000, height=800)
plot(bio[[1]])
points(sp)
title("Current Distribution of Pistacia vera with Bioclimatic Variables")
dev.off()

# 2. Current Habitat Suitability
png("pistacia_plots/2_current_suitability.png", width=1000, height=800)
plot(en, main="Current Habitat Suitability for Pistacia vera")
dev.off()

# 3. Future Habitat Suitability
png("pistacia_plots/3_future_suitability.png", width=1000, height=800)
plot(enf, main="Projected Future Habitat Suitability (2061-2080)")
dev.off()

# 4. Change in Habitat Suitability
png("pistacia_plots/4_suitability_change.png", width=1000, height=800)
plot(ch, col=cl(100), main="Change in Habitat Suitability (Future - Current)")
dev.off()

# 5. Current Presence/Absence
png("pistacia_plots/5_current_presence_absence.png", width=1000, height=800)
plot(pa, main="Current Binary Presence/Absence (Threshold = 0.51247)")
dev.off()

# 6. Future Presence/Absence
png("pistacia_plots/6_future_presence_absence.png", width=1000, height=800)
plot(paf, main="Future Binary Presence/Absence (2061-2080)")
dev.off()

# 7. Change in Presence/Absence
png("pistacia_plots/7_presence_absence_change.png", width=1000, height=800)
plot(pa.ch, col=cl(3), main="Changes in Species Distribution (Current - Future)")
dev.off()














# Set up required libraries
library(gridExtra)
library(ggplot2)
library(grid)

# Function to save individual plots in high quality
save_high_quality_plots <- function() {
  # 1. Current Species Occurrences
  png("pistacia_plots/1_current_occurrences_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(bio[[1]])
  points(sp)
  title("Current Distribution of Pistacia vera with Bioclimatic Variables")
  dev.off()
  
  # 2. Current Habitat Suitability
  png("pistacia_plots/2_current_suitability_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(en, main="Current Habitat Suitability for Pistacia vera")
  dev.off()
  
  # 3. Future Habitat Suitability
  png("pistacia_plots/3_future_suitability_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(enf, main="Projected Future Habitat Suitability (2061-2080)")
  dev.off()
  
  # 4. Change in Habitat Suitability
  png("pistacia_plots/4_suitability_change_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(ch, col=cl(100), main="Change in Habitat Suitability (Current - Future)")
  dev.off()
  
  # 5. Current Presence/Absence
  png("pistacia_plots/5_current_presence_absence_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(pa, main="Current Binary Presence/Absence (Threshold = 0.51247)")
  dev.off()
  
  # 6. Future Presence/Absence
  png("pistacia_plots/6_future_presence_absence_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(paf, main="Future Binary Presence/Absence (2061-2080)")
  dev.off()
  
  # 7. Change in Presence/Absence
  png("pistacia_plots/7_presence_absence_change_hq.png", 
      width = 2000, height = 1600, res = 300)
  plot(pa.ch, col=cl(3), main="Changes in Species Distribution (Current - Future)")
  dev.off()
}

# Function to compile plots into a single figure
compile_plots <- function() {
  # Read all high-quality plots
  plots <- list.files("pistacia_plots", pattern = "_hq.png", full.names = TRUE)
  plot_list <- lapply(plots, png::readPNG)
  
  # Create a new device for the compiled figure
  png("pistacia_plots/compiled_analysis.png", 
      width = 3000, height = 3500, res = 300)
  
  # Set up the layout matrix for 4x2 grid (7 plots)
  layout_matrix <- matrix(c(1:8), nrow = 4, ncol = 2, byrow = TRUE)
  layout(layout_matrix)
  
  # Plot each image
  for(i in 1:7) {
    plot.new()
    rasterImage(plot_list[[i]], 0, 0, 1, 1)
  }
  
  dev.off()
}

# Execute the functions
save_high_quality_plots()
compile_plots()


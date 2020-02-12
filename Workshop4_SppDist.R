install.packages("dismo")
install.packages("maptools")
install.packages("rgdal")
install.packages("raster")
install.packages("sp")
#Installation for utilizing BIOCLIM 

#create folders in workspace
dir.create(path = "data")
dir.create(path = "output")

#data was called using such obs.data <- getData since it would not import properly
obs.data <- read.csv('C:/Users/katca/Desktop/Reproducible_Science/data/SDM_Data.csv')

library("sp")
library("raster")
library("maptools")

library("rgdal")

library("dismo")




#downloading bioclimatic variable data with the get data function

bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5,
                        path = "C:/Users/katca/Desktop/Reproducible_Science")

# got path using getwd() 

#Read in saguaro observations
obs.data <- read.csv(file = "C:/Users/katca/Desktop/Reproducible_Science/data/SDM_Data.csv")

#check the data to make sure it loaded correctly
summary(obs.data)

#removing NAs in dataset
obs.data <- obs.data[!is.na(obs.data$latitude), ]
summary(obs.data)

#to determine geographic extent of data
max.lat <- ceiling(max(obs.data$latitude))
min.lat <- floor(min(obs.data$latitude))
max.lon <- ceiling(max(obs.data$longitude))
min.lon <- floor(min(obs.data$longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

#before doing any modeling, we are going to plot the points on a map.
#load the data to use for our base map
data("wrld_simpl")

#plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

#add points for individual observation
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

#draw a little box around the graph
box()

###Building a model and visualizing results

#all coding above is occurrence data. Now it is time to use bioclimatic variables to create a model
#we are restricting the bioclimatic variable data to the geographic extent of our occurrence data, the SW Region

#crop bioclim data to geographic extent of saguaro
bioclim.data <- crop(x = bioclim.data, y = geographic.extent)

#Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)

#error message given after running the above. Data is not in the right format. 
#In order to solve this issues we are going to
# drop unused column.. drops everything else but lat and long
obs.data <- obs.data[, c("latitude", "longitude")]

#Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)

#Error message recieved, insufficient records.
#To check dataframe use obs.data, going to see how it looks like now
head(obs.data)
 #First column Lat, Second Long
 #Think about how R deals with coordinates. 
    #When plotting something, syntax used is plot(x,y)
 #the bioclim data, lat and long is not set up in the manner that R can plot
 #Latitude is corresponding to the y-axis, Longitude to the x-axis, 
 #R takes the lat long as it is set up, in this case the 1st column x, 2nd y
 #we need to swap the lat and long so R can correspong to the correct data for plotting before passing obs.data to bioclim

#Reverse order of columns
obs.data <- obs.data[, c("longitude", "latitude")]

#Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)

#Before proceeding to plot model on map. Generate object that has the models probability of occurrence for saguros
#Use the predict model from the dismo package

#Predict presence from model
predict.presence <- dismo::predict(object = bc.model, x = bioclim.data, ext = geographic.extent)
     #dismo::predict is used rather than predict. NEEDs to be clarified as there are 3 predict functions loaded into memory

#plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

#Add model probabilities
plot(predict.presence, add = TRUE)

#Redraw those county boarders
plot(wrld_simpl, add = TRUE, border = "grey5")

#Add original Observations
points(obs.data$longitude, obs.data$latitude, col = "olivedrab", pch = 20, cex = 0.75)
box()

##__the plot shows probability of occurence of saguros, though max probability anywhere on the map is 0.78, less than 1
  #Saguaros are found a pretty broad area of the Sonoran Desert, we do have observations to prove it.
  #For the map to better reflect this, analyses needs to be re-run but this time include absence points where they are known not to occur. Though, we only have presence data for saguaros


###The Pseudo-Absence Point
  #Psudeo-absence means one randomly samples points from a given geographic area, treat them like locations where the species of interest is absent. 
  #We are going to create a set of background (PSudoabs) points at random with as many points as we have observations.
  #using bioclim data files for determining spataial resolution of the points, restricting the sampling area to the general region of the observ of saguaros.

#Use the bioclim data files for sampling resolution (remember to have path correct)
bil.files <- list.files(path = "C:/Users/katca/Desktop/Reproducible_Science/wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)
#We only need one file, so use the first one in the list of .bil files
mask <- raster(bil.files[1])

#Randomly sample points (same number as our observed points)
background <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(obs.data),      # Number of random points
                           ext = geographic.extent, # Spatially restricts sampling
                           extf = 1.25)             # Expands sampling a little bit

#to see selected random points:
head(background)

#to visualize on a map, as we had done so for observed points
#plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

#add background points
points(background, col = "grey30", pch = 1, cex = 0.75)

#Add the observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

box()

#now that we have psudo-absence points, we require a more traditional-range-map-lookingfigure requires post hov evaluation of the model
#to do this post hoc evaluation we are to build the model using only part of our data (the training data), reserving a portion of the data for ecaluation of the model after it is build (the testing data).
#We are reserving 20% of the data for testing, so we could use kfold function in the dismo package to evenly assign each observation to a random group

#Arbitrarily assign group 1 as the testing data group
testing.group <- 1

#create vector of group memberships
group.presence <- kfold(x = obs.data, k = 5) # kfold is in dismo package

  #take a look at group presence vector just created utilizing formula:
head(group.presence)

#should see even representation in each group
table(group.presence)
    #the output of table shows how many points have been assigned to each of the five groups, in this case, the points have been evenly distributed, with 20% of points in group 1, out testing group

#Using group.presence vector with the observed data to seperate our observations into a training data set and a testing data set:
#Seperate observations into training and testing groups
presence.train <- obs.data[group.presence != testing.group, ]
presence.test <- obs.data[group.presence == testing.group, ]

#Repeat the process for psuedo-ansence points
group.background <- kfold(x = background, k = 5)
background.train <- background[group.background != testing.group, ]
background.test <- background[group.background == testing.group, ]


###Training and Testing the Model
#we have 1 psuedo-absence points, and 2 seperate training and testing data, we can rebuild model and evaluate its performance, draw a more aesthetically pleasing map
#we build model with bioclim function as before, but instead of using all observations in obs.data we only use the training data stored in presence.train:

# Build a model using training data
bc.model <- bioclim(x = bioclim.data, p = presence.train)

# Predict presence from model (same as previously, but with the update model)
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data, 
                                   ext = geographic.extent)

#we take model, and evaluate it using the observation data and the psuedo-absence points we reserved for model testing. We then use this test to establish a cutoff of occrrence probability to determine the boundaries of the saguaro range
# Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test,   # The presence testing data
                    a = background.test, # The absence testing data
                    model = bc.model,    # The model we are evaluating
                    x = bioclim.data)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "spec_sens")

#threshold function offers a number of means of determining the threshols cutoff through the stat parameter. In this we chose "spac_sens", which sets "the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest."

#We can use that threshold to paint a map with the predicted range of sagauro
# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")
# Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = "olivedrab")
# And add those observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)
# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()
#"doesnt look right, plotted a large portion of the map green. Lets see what we actually asked R to plot, that is, we plot the value of predict.presence > bc.threshold. SO what is that?
predict.presence > bc.threshold

#the comparison of these two rasters produces another raster with values of only 0 or 1:0 where the comparison evaluates as FASLE (ie when the value in a grid cell of predict.presence is less than or equal to the value in the corressponding grid cell of bc.threshold)
#and 1 where the comparison evaluates as TRUE. 
#Since there are 2 values in this comparison (the 0 and 1 in the values field), we need to update what we pass to the col parameter in our plot call. Instead of just passing a single value, we provide a color for 0 (NA) and a color for 1 ("olivedrab"):

# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "olivedrab"))


# And add those observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

#Conlusion: A final note on our approach: the map we have drawn presents a categorical classification of whether a particular point on the landscape will be suitable or not for the species of interest. 
  #This classification relies quite heavily on the value of the threshold (see bc.threshold and the documentation for threshold) and the pseudo-absence points. 
  #Given that we used random sampling to generate those pseudo-absence points, there is potential for variation in the predicted range if you run this code more than once (try it! if you re-run the code from the point of creating the pseudo-absence points, you are almost guaranteed a different map.). 
  #There are a number of approaches to dealing with this variation, and the paper by Barbet-Massin et al. (2012) is a great resource. I’ll leave it as homework for you to determine which approach is most appropriate here!


* Setup/configuration
  * imports
  * paths
  * constants
  * target species
* Extract presence data (done several times for different sources)
  * Query API by target species
  * Download occurrences
  * Clean/QAQC
  * Standardize columns
  * Intersect with buffered lake polygons
* Extract pseudoabsence data
  * Get observer list from presence data
  * Query API for other observations reported by those observers
  * Remove observations of target species
  * Spatially filter to lakes
* NHD addressing
  * Snap presence observations to the nearest NHD flowline
  * Iterate through all observations chronologically so you're only calculating nearest-neighbor distances to presences that were observed before the focal observation
  * Trace the target observation upstream
  * Count how many species are present upstream
  * Get the alongstream and Euclidean distances to the closest previously reported presence
* Generate XGBoost models
  * Clean environmental variables
  * Sample from the pseudoabsence observations so that for each presence observation there is a corresponding pseudoabsence observation
  * Split off most recent records as an external test set
  * Use TimeSeriesSplit to split the observations for a time-aware cross-validation
  * Use a grid search to identify optimal hyperparameter values
  * Fit a boosted regression model to the presence, pseudoabsence and environmental data
  * Produce some metrics and charts for model evaluation
  * Predict presence probability for all lakes in region

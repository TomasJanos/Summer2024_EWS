################################################################################
#
# Heat-related mortality in Europe during 2024 and health emergency forecasting
# to reduce preventable deaths
#
# R Code illustrating the analyses - part 1
#
# This is a simplified code with sample data.
# Results are expected to differ from those published in the article.
# 
# Contact: Tomáš Janoš (tomas.janos@isglobal.org)
# Last updated: 09 July 2025
#
################################################################################

# Clear the environment and the console
rm( list = ls() );
cat("\014");

# Required libraries
suppressMessages( library(dlnm) ); # Required Functions: logknots, crossbasis
suppressMessages( library(splines) ); # Required Functions: ns, bs
suppressMessages( library(lubridate) ); # Required Functions: wday, month
suppressMessages( library(mixmeta) ); # Required Functions: mixmeta
suppressMessages( library(ISOweek) ); # ISOweek2date


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
#### Definition of the parameters ####
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#

# Period for the predictions: summer 2024
DATE1_PRED = as.Date( "2024-06-01" );
DATE2_PRED = as.Date( "2024-09-30" );

# Formula for the Best Linear Unbiased Predictions and pooled coefficients
FORMULA_META = COEF_MODEL ~ MP_TEMP_ANN + MP_TEMP_IQR;

# Formula of the model
FORMULA = mort ~ ns( date, df = round( DF_SEAS * length(date) / 365.25 ) ) + dow + CROSSBASIS_TEMP;

# Exposure-Response: type
FUNC_TEMP = "ns"

# Exposure-Response: percentiles of the temperature knots
PERC_TEMP = c(10,75,90) / 100;

# Lag-Response: minimum and maximum lag
MIN_LAG =  0;
MAX_LAG = 21;

# Lag-Response: lag knots
LAG_KNOTS = logknots( c( MIN_LAG, MAX_LAG ), 3 );

# Degrees of freedom per year for the seasonal and long-term trends
DF_SEAS = 8;

# Percentiles for the predictions of the cumulative Exposure-Response
vPERC_CROSSPRED = sort( unique( c( seq(  0.0,   1.0, 0.1 ),
                                   seq(  1.5,   5.0, 0.5 ),
                                   seq(  6.0,  94.0, 1.0 ),
                                   seq( 95.0,  98.5, 0.5 ),
                                   seq( 99.0, 100.0, 0.1 ) ) / 100 ) );

# Number of simulations
nSIM = 1000;

# Vector of time periods for all regions with mortality data: whole period
vPER = c( "Whole Period" );
nPER = length(vPER);

# Vector of temperature ranges
vRNG = c( "Total", "Total Cold", "Total Heat", "Moderate Cold", "Moderate Heat", "Extreme Cold", "Extreme Heat" );
nRNG = length(vRNG);

# 95% Confidence Interval
vCI95 = c(0.025,0.975);

# Value for seeding random variate generators
SET_SEED_VALUE = 76284739;


#xxxxxxxxxxxxxxxxxxxxxxxx#
#### Data Preparation ####
#xxxxxxxxxxxxxxxxxxxxxxxx#

# Reading the daily mortality data table
DATA_MORT_DAILY = read.csv( "./data_mort_daily.csv" );
# Ensuring date is in a date format
DATA_MORT_DAILY$date = as.Date(DATA_MORT_DAILY$date)

# Reading the daily temperature data table
DATA_TEMP_DAILY = read.csv( "./data_temp_daily.csv" );
# Ensuring date is in a date format
DATA_TEMP_DAILY$date = as.Date(DATA_TEMP_DAILY$date)

# Adding temperature data to daily mortality data table
DATA_MORT_DAILY = plyr::join( DATA_MORT_DAILY, DATA_TEMP_DAILY, by = c( "location", "date" ), type = "left" );

# Adding day of the week indicator to daily mortality data table
DATA_MORT_DAILY$dow = lubridate::wday( DATA_MORT_DAILY$date, week_start = 1 );

# Reading the weekly mortality data table
DATA_MORT_WEEKLY = read.csv( "./data_mort_weekly.csv" );

# Transformation of the weekly mortality time series into daily mortality by dividing by seven
# creating time series of days within the period and vector of regions
date = seq( as.Date( DATE1_PRED, "%Y%m%d"), as.Date( DATE2_PRED + MAX_LAG, "%Y%m%d"), by = "days");
location = unique(DATA_MORT_WEEKLY$location);

# Creating daily data table
DATA_MORT_DW = expand.grid(location, date); colnames(DATA_MORT_DW) <- c("location","date");
rm(date, location);

# adding week-year variable
DATA_MORT_DW$week = ISOweek::ISOweek(DATA_MORT_DW$date);

# Joining tables and dividing weekly mortality by seven to get daily mortality
DATA_MORT_DW = plyr::join( DATA_MORT_DW, DATA_MORT_WEEKLY, by = c( "location", "week" ), type = "left" );
DATA_MORT_DW$mort = DATA_MORT_DW$mort / 7;
rm(DATA_MORT_WEEKLY)

# Adding temperature data to daily-weekly mortality data table
DATA_MORT_DW = plyr::join( DATA_MORT_DW, DATA_TEMP_DAILY, by = c( "location", "date" ), type = "left" );
rm(DATA_TEMP_DAILY)

# Creating the vector of codes of the regions in the temperature and mortality Table
nREG = as.numeric(length(unique(DATA_MORT_DAILY$location)));
vREG = unique(DATA_MORT_DAILY$location);

# Creating the temperature and mortality lists with all the regions for the calibration and prediction periods
DATALIST_DATA_CALI = lapply( vREG, function(x) DATA_MORT_DAILY[ DATA_MORT_DAILY$location == x, ] );
DATALIST_DATA_PRED = lapply( vREG, function(x) DATA_MORT_DW[ DATA_MORT_DW$location == x, ] );
names(DATALIST_DATA_CALI) = names(DATALIST_DATA_PRED) = vREG;
rm(DATA_MORT_DAILY, DATA_MORT_DW);


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
#### Location-Specific Associations ####
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#

# Reduced coefficients
COEF_MODEL = matrix( data = NA, nREG, length(PERC_TEMP) +    1   , dimnames = list(vREG) );

# Reduced covariance matrices
VCOV_MODEL = vector( "list", nREG );
names(VCOV_MODEL) = vREG;

for( iREG in 1:nREG ){
  
  # Cross-Basis of temperature
  CROSSBASIS_TEMP = crossbasis( DATALIST_DATA_CALI[[iREG]]$temp,
                                c( MIN_LAG, MAX_LAG ),
                                arglag = list( knots = LAG_KNOTS ),
                                argvar = list( knots = quantile( DATALIST_DATA_CALI[[iREG]]$temp, PERC_TEMP, na.rm = TRUE ),
                                               fun = FUNC_TEMP,
                                               Bound = range( DATALIST_DATA_CALI[[iREG]]$temp, na.rm = TRUE ) ) );

  # Fitting the Cross-Basis model
  GLM_MODEL = glm( formula = FORMULA,
                   DATALIST_DATA_CALI[[iREG]],
                   family = quasipoisson,
                   na.action = "na.exclude" );

  # For now, centering at median temperature
  CEN = median(DATALIST_DATA_CALI[[iREG]]$temp)
  
  # Regional reduced coefficients and covariance matrix: cumulative Exposure-Response
  REDUCED = crossreduce( CROSSBASIS_TEMP,
                         GLM_MODEL,
                         model.link = "log",
                         cen = CEN );
  COEF_MODEL[iREG,] = coef( REDUCED );
  VCOV_MODEL[[iREG]] = vcov( REDUCED );
  rm(REDUCED);
  
  rm(CROSSBASIS_TEMP, GLM_MODEL);
  
}
rm(iREG);


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
#### Multivariate Meta-Analysis and Best Linear Unbiased Predictions ####
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#

# Creating objects
META = BLUP = vector( "list", 1 );

# Regional cumulative Exposure-Response after the meta-analysis
CROSSPRED = vector( "list", nREG );
names(CROSSPRED) = vREG;

# Regional Minimum Mortality Temperature after the meta-analysis
MMT = array( NA, dim = nREG, dimnames = list(vREG) );

# Regional meta-predictors (see FORMULA_META):
MP_TEMP_ANN = sapply( DATALIST_DATA_CALI, function(x) mean( x$temp, na.rm = TRUE ) );
MP_TEMP_IQR = sapply( DATALIST_DATA_CALI, function(x) IQR( x$temp, na.rm = TRUE ) );

# Multivariate meta-analysis of the cumulative Exposure-Response #
META = mixmeta( FORMULA_META,
                VCOV_MODEL,
                data = data.frame( vREG = vREG ),
                control = list( showiter = TRUE, igls.inititer = 10, maxiter = 1000 ),
                method = "reml" );

# Regional Best Linear Unbiased Predictions
BLUP = blup( META, vcov = TRUE );
rm(COEF_MODEL,VCOV_MODEL);

for( iREG in 1:nREG ){
  
  # Vector of temperatures for the cross-predictions
  TEMP_CROSSPRED = quantile( DATALIST_DATA_CALI[[iREG]]$temp, vPERC_CROSSPRED, na.rm = TRUE );
  
  # One-Basis of temperature
  ONEBASIS_TEMP = onebasis( TEMP_CROSSPRED,
                            knots = TEMP_CROSSPRED[ paste0( 100 * PERC_TEMP, ".0%" ) ],
                            fun = FUNC_TEMP,
                            Bound = range( TEMP_CROSSPRED, na.rm = TRUE ) );
  
  # Very exhaustive ordered vector with all possible temperatures to calculate the Minimum Mortality Temperature
  # This is done to minimize the occurrence of negative predicted Attributable Fractions due to lack of precision of the Minimum Mortality Temperature
  TIMESERIES_CROSSPRED = sort( unique( c( DATALIST_DATA_CALI[[iREG]]$temp, DATALIST_DATA_PRED[[iREG]]$temp ) ) );
  
  # Regional cumulative Exposure-Response without centering
  MORT_CROSSPRED = crosspred( ONEBASIS_TEMP,
                              coef = BLUP[[iREG]]$blup,
                              vcov = BLUP[[iREG]]$vcov,
                              model.link = "log",
                              at = TIMESERIES_CROSSPRED,
                              cen = mean( TEMP_CROSSPRED, na.rm = TRUE ) );

  # Regional Minimum Mortality Temperature after the meta-analysis
  # The Minimum Mortality Temperature is the lowest Relative Risk
  MMT[iREG] = MORT_CROSSPRED$predvar[ which.min( MORT_CROSSPRED$allRRfit ) ];
  
  rm(TIMESERIES_CROSSPRED, MORT_CROSSPRED);
  
  # Regional cumulative Exposure-Response with centering
  CROSSPRED[[iREG]] = crosspred( ONEBASIS_TEMP,
                                 coef = BLUP[[iREG]]$blup,
                                 vcov = BLUP[[iREG]]$vcov,
                                 model.link = "log",
                                 at = TEMP_CROSSPRED,
                                 cen = MMT[iREG] );

  }
rm(iREG);


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
#### Attributable Mortality for the Period of the Predictions ####
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#

# Values and simulations of the regional attributable number (in total deaths)
ATTNUM_VAL = array( NA, dim = c( nREG, nPER, nRNG       ), dimnames = list( vREG, vPER, vRNG         ) );
ATTNUM_SIM = array( NA, dim = c( nREG, nPER, nRNG, nSIM ), dimnames = list( vREG, vPER, vRNG, 1:nSIM ) );
ATTNUM     = array( NA, dim = c( nREG, nPER, nRNG, 3    ), dimnames = list( vREG, vPER, vRNG, c( "att_val", "att_low", "att_upp" ) ) );

# Thresholds of the regional temperature ranges
P025_REG = sapply( DATALIST_DATA_PRED, function(x) quantile( x$temp, 0.025, na.rm = TRUE ) ); names(P025_REG) = vREG;
P975_REG = sapply( DATALIST_DATA_PRED, function(x) quantile( x$temp, 0.975, na.rm = TRUE ) ); names(P975_REG) = vREG;

for( iREG in 1:nREG ){
  
  # One-Basis of temperature
    ONEBASIS_TEMP = onebasis( DATALIST_DATA_PRED[[iREG]]$temp,
                              knots = quantile( DATALIST_DATA_CALI[[iREG]]$temp, PERC_TEMP, na.rm = TRUE ),
                              fun = FUNC_TEMP,
                              Bound = range( DATALIST_DATA_CALI[[iREG]]$temp, na.rm = TRUE ) );

  # One-Basis of the Minimum Mortality Temperature to center the One-Basis of temperature
    ONEBASIS_MMT = onebasis( MMT[iREG],
                             knots = quantile( DATALIST_DATA_CALI[[iREG]]$temp, PERC_TEMP, na.rm = TRUE ),
                             fun = FUNC_TEMP,
                             Bound = range( DATALIST_DATA_CALI[[iREG]]$temp, na.rm = TRUE ) );
  
  # One-Basis of temperature centered at the Minimum Mortality Temperature
  ONEBASIS_CEN = scale( ONEBASIS_TEMP, center = ONEBASIS_MMT, scale = FALSE );
  rm(ONEBASIS_TEMP,ONEBASIS_MMT);
  
  # Matrix of deaths with lags (dimensions: time x lags)
  LAGGED_MORT_MATRIX = tsModel::Lag( DATALIST_DATA_PRED[[iREG]]$mort, -seq( MIN_LAG, MAX_LAG ) );
  
  # Average of deaths across lags (dimensions: time x 1)
  LAGGED_MORT_VECTOR = rowMeans( LAGGED_MORT_MATRIX, na.rm = FALSE );

  # Daily time series of the regional Attributable Number (dimensions: time x 1)
  ATTNUM_TS_REF = ( 1 - exp( -ONEBASIS_CEN %*% BLUP[[iREG]]$blup ) ) * LAGGED_MORT_VECTOR;

  # Perturbed Best Linear Unbiased Predictions (dimensions: coefficients x simulations)
  set.seed(SET_SEED_VALUE);
  COEF_SIM = t( MASS::mvrnorm( nSIM, BLUP[[iREG]]$blup, BLUP[[iREG]]$vcov ) );
  
  # Matrix of perturbed daily time series of the regional Attributable Number (dimensions: time x simulations)
  ATTNUM_TS_SIM = ( 1 - exp( -ONEBASIS_CEN %*% COEF_SIM ) ) * LAGGED_MORT_VECTOR;
  rm(COEF_SIM);
  
  for( iPER in 1:length(vPER) ){
    
    # Data selection of the time period
    if( vPER[iPER] == "Whole Period" ){ vTIM = 1:length( DATALIST_DATA_PRED[[iREG]]$date ); }

      for( iRNG in 1:nRNG ){
        
        # Data selection of the temperature range      
        if     ( vRNG[iRNG] == "Total"         ){ vRNG_THRES = 1:length(vTIM); }
        else if( vRNG[iRNG] == "Total Cold"    ){ vRNG_THRES = which( DATALIST_DATA_PRED[[iREG]]$temp[vTIM] < MMT[iREG] ); }
        else if( vRNG[iRNG] == "Total Heat"    ){ vRNG_THRES = which( MMT[iREG] < DATALIST_DATA_PRED[[iREG]]$temp[vTIM] ); }
        else if( vRNG[iRNG] == "Moderate Cold" ){ vRNG_THRES = which( P025_REG[iREG] < DATALIST_DATA_PRED[[iREG]]$temp[vTIM] & DATALIST_DATA_PRED[[iREG]]$temp[vTIM] < MMT[iREG] ); }
        else if( vRNG[iRNG] == "Moderate Heat" ){ vRNG_THRES = which( MMT[iREG] < DATALIST_DATA_PRED[[iREG]]$temp[vTIM] & DATALIST_DATA_PRED[[iREG]]$temp[vTIM] < P975_REG[iREG] ); }
        else if( vRNG[iRNG] == "Extreme Cold"  ){ vRNG_THRES = which( DATALIST_DATA_PRED[[iREG]]$temp[vTIM] < P025_REG[iREG] ); }
        else if( vRNG[iRNG] == "Extreme Heat"  ){ vRNG_THRES = which( P975_REG[iREG] < DATALIST_DATA_PRED[[iREG]]$temp[vTIM]                                                                       ); }
        else if( vRNG[iRNG] == "Summer Heat"   ){ vRNG_THRES = which( MMT[iREG] < DATALIST_DATA_PRED[[iREG]]$temp[vTIM] & DATALIST_DATA_PRED[[iREG]]$month[vTIM] %in% c(6,7,8,9) ); }
        else                                    { stop("ERROR: Invalid Temperature Range !!!"); }
        if( length(vRNG_THRES) > 0 ){
          
        # Attributable Number: Value
          ATTNUM_VAL[iREG,iPER,iRNG ] = sum( ATTNUM_TS_REF[ vTIM[vRNG_THRES] ], na.rm = TRUE ) *
            sum( rowMeans( LAGGED_MORT_MATRIX[ vTIM, , drop=FALSE ], na.rm = TRUE ), na.rm = TRUE ) /
            sum(           LAGGED_MORT_VECTOR[ vTIM               ]                , na.rm = TRUE );

        # Attributable Number: Simulations
          ATTNUM_SIM[iREG,iPER,iRNG,] = colSums( ATTNUM_TS_SIM[ vTIM[vRNG_THRES], , drop=FALSE ], na.rm = TRUE ) *
            sum( rowMeans( LAGGED_MORT_MATRIX[ vTIM, , drop=FALSE ], na.rm = TRUE ), na.rm = TRUE ) /
            sum(           LAGGED_MORT_VECTOR[ vTIM               ]                , na.rm = TRUE );

        # Attributable Number: Value and Confidence Interval
          ATTNUM[iREG,iPER,iRNG, 1 ] =           ATTNUM_VAL[iREG,iPER,iRNG ]                  ;
          ATTNUM[iREG,iPER,iRNG,2:3] = quantile( ATTNUM_SIM[iREG,iPER,iRNG,], vCI95 );
          
        }else{
          ATTNUM[iREG,iPER,iRNG,] = 0;
          ATTNUM_VAL[iREG,iPER,iRNG ] = 0;
          ATTNUM_SIM[iREG,iPER,iRNG,] = 0;
        }
        rm(vRNG_THRES);
      }
      rm(iRNG);
    }
  rm(ONEBASIS_CEN, LAGGED_MORT_MATRIX, LAGGED_MORT_VECTOR, ATTNUM_TS_REF, ATTNUM_TS_SIM, iPER);
}
rm(P025_REG,P975_REG, iREG);







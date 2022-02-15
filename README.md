# AnalogVolcanoForecast

The computational code and data are for the paper: Wang, T., Bebbington, M., Cronin, S., and Carman, J. (2022) Forecasting eruptions at poorly known volcanoes using analogues and multivariate renewal processes, under review.

All .R and .jags files can be downloaded and saved in the same folder. Instructions are given if they depend on code in other sections, in which case one simply needs to run the R code in that section first.

### GVPmod0WeibsepZuse.jags 
This file contains the JAGS "model" for the Weibull model in equation (1) in the paper Wang et al. (2022).

### GVPmod1WeibsepZuse.jags
This file contains the JAGS "model" for model M1 in equation (2) in the paper Wang et al. (2022).

### GVPmod2WeibsepZuse.jags
This file contains the JAGS "model" for model M2 in equation (3) in the paper Wang et al. (2022).

### GVPmod3WeibsepZuse.jags
This file contains the JAGS "model" for model M3 in equation (4) in the paper Wang et al. (2022).

### ReadGVPdata.R
This file reads all the data files for the empirical analogues saved in the folder "HoloceneRecords". These data files are obtained from the Smithsonian Institution's GVP catalogue.

## Fitting models M1, M2, and M3 to the empirical analogues 
#### Run the R code in "ReadGVPdata.R" to read all the GVP data first


## Using eruption records with minimum VEI 3

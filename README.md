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
### GVPEmpAnalog.R
##### Make sure to save all .R and .jags files and the two data folders "HoloceneRecords" and "QuaternaryRecords" in the same main folder.
Run the file GVPEmpAnalog.R, and the MCMC samples of the posterior distributions will be saved in .image files.

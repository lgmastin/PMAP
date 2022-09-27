# PMAP software:  installation and usage

The PMAP tool is a one-dimensional plume model that predicts ash aggregation within an erupting volcanic plume.  The tool was developed by D. W. Hoffman at Stanford, extending from the Plumeria model written by L. G. Mastin at USGS CVO.  The code operates using Matlab.

The purpose of this repository is to present the code used in the following publication:

Hoffman, D.W., Mastin, L.G., Van Eaton, A.R., Solovitz, S.A., Cal, R.B., and Eaton, J.K., in review, A 1-D volcanic plume model for predicting ash aggregation

## OPERATION

The code can be run by opening the file "example_usage.m", and hitting the "run" button in Matlab.

* Input properties (such as plume radius, velocity, mass fraction gas, etc.) are specified in the "Model Setup" section of example_usage.m. 
* Ambient conditions (such as from a meteorological sounding) are specified in the "input:ambient conditions" subsection of example_usage.m.
* Grain size distribution (such as provided from a field measurement or estimate) and aggregation parameters are specified in the "Input: aggregation parameters" subsection of the "Model Setup" section of example_usage.m.

## TESTING

Running the script “example_usage.m” will run the model.  If no parameters are changed before execution, it will generate outputs for an eruption of Redoubt Volcano in 2009 identical to those presented in Hoffman et al. (in review). It will generate the following outputs:

* Figure 1: plume mass balance (%) versus height asl
* Figure 2: plume velocity versus height asl
* Figure 3: plume radius (km) versus height asl
* Figure 4: mass fractions of water (all phases), air, and magma (particles) versus height asl
* Figure 5: air pressure versus height asl
* Figure 6: temperature of the plume and ambient air versus height asl
* Figure 7: density of the plume and ambient air versus height asl
* Figure 8: a plot of the mass fractions of water vapor, liquid water, and ice versus height above sea level in the plume
* Figure 9: original and aggregated size distributions.  This plot should be a mirror image of Figure 7 of Hoffman et al. (in review), but contain only two of the four curves: the original grain-size distribution (labeled "GSD" in H. et al.), and the aggregated size distribution (labeled "AGSD - PMAP/baseline" in Hoffman et al. (in reveiew) To reproduce the curve labeled "AGSD - PMAP/Costa", change the parameter "input.StickModel" in example_usage.m (line 44) from 'Hoffman' to 'Costa', and rerun the model.
* Figure 10:  a plot that should be identical to Figure 6 of Hoffman et al. (in review)

Figures 1-7 in the model output have no analogous figures in Hoffman et al. (in review)

**--note:** the code gives file paths in a Microsoft Windows format (e.g. 'c:\filepath\filename').  If you are using a Macintosh, Linux, or some other system, you will have to modify the paths to conform with your operating system.

## THEORY

The primary code, PMAP.m, uses a series of fluid, thermodynamic, and aggregation models to predict the resulting plume properties and aggregated grain size distribution.  Conservation equations and constitutive properties are presented in Mastin (2007, 2014).  Methodology for aggregation calculations are provided in Hoffman et al. (in review).  Individual M-files apply the various models.  Notably, the user may select two different aggregation efficiency models, based on earlier work by Costa et al. (2010) and more recent experimental studies by Hoffman (2022).  Aggregation settings are specified in the "Input: aggregation parameters" subsection of the "Model setup" section of example_usage.m.

## References:
* Costa A, Folch A, Macedonio G (2010) A model for wet aggregation of ash particles in volcanic plumes and clouds: 1. Theoretical formulation. Journal of Geophysical Research 115:B09201. https://doi.org/10.1029/2009JB007175
* Hoffman DW (2021) The Role of Turbulence in the Aggregation of Ash Particles in Volcanic Plumes.  Department of Mechanical Engineering. Stanford University, Stanford, CA, p 177, https://purl.stanford.edu/dz177dw3065
* Hoffman, D.W., Mastin, L.G., Van Eaton, A.R., Solovitz, S.A., Cal, R.B., and Eaton, J.K., in review, A 1-D volcanic plume model for predicting ash aggregation
* Mastin LG (2007) A user-friendly one-dimensional model for wet volcanic plumes. Geochemistry, Geophysics, Geosystems 8(Q03014):doi:10.1029/2006GC001455. https://doi.org/10.1029/2006GC001455
* Mastin LG (2014) Testing the accuracy of a 1-D volcanic plume model in estimating mass eruption rate. Journal of Geophysical Research: Atmospheres 119(5). http://dx.doi.org/10.1002/2013JD020604









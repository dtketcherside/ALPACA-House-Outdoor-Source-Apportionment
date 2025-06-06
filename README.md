This repository contains the an R project dedicated to reproducing the analysis and plots for 
"Wintertime Abundance and Sources of Key Trace Gas and Particle Species in Fairbanks, Alaska": 
a paper submitted to the Journal for Geophysical Research: Atmospheres in February of 2025.

ALPACA 2022 House Site Indoor and Outdoor Trace Gas and Particulate Data

Prepared by: Damien Ketcherside (damien.ketcherside@umontana.edu)

Additional Researchers Involved in Data Collection Efforts: 
            
    Dr. Vanessa Selimovic, Dr. Lu Hu, Dr. Bob Yokelson (The University of Montana)
    Andrew Holen, Judy Wu, Dr. Kerri Pratt (The University of Michigan)
    Dr. Ellis Robinson and Dr. Peter DeCarlo (Johns Hopkins University)
    Dr. Jonas Khun and Dr. Jochen Stutz (University of California - Los Angeles)
    Dr. Meeta Cesler-Maloney, Dr. Jingqui Mao, and Dr. William Simpson (The University of Alaska - Fairbanks)
			
R Project Description:

    The R project in this repository, titled "ALPACA House Outdoor Source Apportionment.Rproj",
    contains three distinct R scripts:
    
    	1) OSTRTA.R - Control Script: This single script will rerun the entire analysis, reading
     	   data from the data folder and sourcing the next two scripts.
     
     	2) Custom Functions.R - Function Repository: This script contains all of the functions
      	   used to perform analyses and generate figures used in the manuscript.
	  
      	3) Make and Save Figures.R - Data Visualization Script: Generates and saves figures to the
       	   output folder as .tiff format images.

    The R project employs the renv package to ensure stability of the project regardless of
    updates made by developers to their respective packages used to perform these analyses.

Data Description:

    The included data was collected at an unoccupied Airbnb in Fairbanks, Alaska as 
    part of the ALPACA 2022 field intensive. Other data include measurements taken at the
    CTC and Birch Hill sites. Data have been background corrected.

    VOC, NOx, and O3 measurements shared an inlet. An inlet switch was used to transition 
    between two heated inlets, one indoor and another outdoor, every 10 minutes. This switch 
    consisted of two electronic 1/4" PTFE three-way valves. Additionally, a bypass was 
    included in the design to facilitate quick instrument response during switching. Data on 
    either side of an inlet transition (+/- 1 minute) was removed to limit measurement of 
    mixed indoor and outdoor samples. Data points that were found to be abnormally high for a 
    short period of time (1 minute) were removed, as these are likely not "real" or important 
    (e.g., artifacts, transient source, etc.), nor are they useful for source attribution.
  
    All other measurements were taken via a separate, non-heated inlet switch system. While 
    researchers attempted to sync indoor/outdoor switching of the two inlet systems, this was 
    not always possible. Therefore, there are time periods in which VOCs, NOx, and O3 were 
    being sampled indoors while formaldehyde, CO, CO2, CH4, and particulates were being 
    sampled outdoors. All viable outdoor samples are included in this file, regardless of inlet 
    overlap.

Special Notes: 

    1) The 'Experiment_Type' flags begin 1 hour prior to experiment start and three hours
       after experiment end.
         
    3) A butanol CPC was installed at the house site on February 11th, 2022. This
       instrument uses butanol and diethylene glycol as solvents during operation. The
       instrument was located inside the house, with a liquid exhaust line routed to the
       outside. Due to fragmentation of butanol and diethylene glycol, there is a
       substantial interference with two VOC species included in this dataset:
       propyne and butene.
         
    5) Several VOC species were calibrated in the field using four certified gas standards.
       Formic and acetic acid were quantified using a humidity-dependent sensitivity
       function derived from liquid calibrations performed prior to deployment.
       D5-siloxane mixing ratios were quantified using a sensitivity derived from a gas
       standard calibration post-deployment. Mixing ratios for uncalibrated species
       were estimated using calculated sensitivities.

Data File Descriptions:

    The file 'data/ALPACA_Outdoor_1min_20240607.txt' contains 1-minute resolution outdoor data 
    collected at the house site.
    
    The file 'data/PACT1D_ALPACA_wSnow3m.txt' contains model output used for NO3 radical plots.
    
    The file 'data/PAFA.csv' contains meteorological data from the Fairbanks International Airport.
    
    The file 'data/PMF_species_perc_dist.csv' contains the PMF factor contributions to individual
    species' fit by the model solution.
    
    The file 'data/PMF_TS_contributions.csv' contains the PMF factor time series.

Variable Definitions (Variable, Unit, Instrument/Description):
		
    Time_AKST,AKST,Local Measurement Time 
    Methanol, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 33.033)
    Propyne, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 41.039)
    Acetonitrile, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 42.034)
    Acetaldehyde, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 45.033)
    Formic_Acid, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 47.013)
    Ethanol, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 47.049)
    Butene, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 57.070
    Acetone, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 59.049)
    Acetic_Acid, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 61.028)
    Furan, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 69.033)
    Isoprene, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 69.070)
    MVK_MACR, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 71.049)
    Methyl_Glyoxal, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 73.028)
    MEK, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 73.065)
    Hydroxyacetone, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 75.044)
    Benzene, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 79.054)
    Methylfuran, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 83.049)
    Pentanone, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 87.079)
    Toluene, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 93.070)
    Furfural, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 97.028)
    Maleic_Anhydride, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 99.011)
    Naphthalene, ppbv,  Ionicon Analytik PTR-ToF-MS 4000 (m/z 129.070)
    Hexanone, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 101.096)
    C8_Aromatics, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 107.086)
    Methylfurfural, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 111.044)
    C9_Aromatics, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 121.101)
    C10_Aromatics, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 135.117)
    Monoterpenes, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 137.132)
    D5_Siloxane, ppbv, Ionicon Analytik PTR-ToF-MS 4000 (m/z 371.101)
    Ozone, ppbv, 2B Tech Model 211 Scrubberless Ozone Monitor
    NO, ppbv, 2B Tech Model 405 NO/NO2/NOx Monitor
    NO2, ppbv, 2B Tech Model 405 NO/NO2/NOx Monitor
    NOx, ppbv, 2B Tech Model 405 NO/NO2/NOx Monitor
    Formaldehyde, ppbv, Picarro G2401
    CO, ppmv, Picarro G2401
    CO2, ppmv, Picarro G2401
    CH4, ppmv, Picarro G2401
    NH3, ppbv, Picarro G2103
    Org, ug/m3, AMS (Organic Fraction)
    SO4, ug/m3, AMS (Sulfate Fraction)
    NO3, ug/m3, AMS (Nitrate Fraction)
    NH4, ug/m3, AMS (Ammonium Fraction)
    Chl, ug/m3, AMS (Chelated Species Fraction)
    TNRM, ug/m3, AMS (Total non-refractory mass)
    PM1, ug/m3, TRNM+BC
    BC, ug/m3, Aerosol Magee Scientific AE-33 Aethalometer
    Experiment_Type, NA, Indicates what kind of house experiment, if any, was being conducted.
    House_Occupied, NA, Indicates whether the study house was occupied.
    CTC_temp_3m,deg C,Temperature measured 3 meters a.g.l. at CTC site.
    CTC_temp_6m,deg C,Temperature measured 6 meters a.g.l. at CTC site.
    CTC_temp_11m,deg C,Temperature measured 11 meters a.g.l. at CTC site.
    CTC_temp_23m,deg C,Temperature measured 23 meters a.g.l. at CTC site.
    CTC_CO2_3m,ppm,CO2 concentration measured 3 meters a.g.l. at CTC site.
    CTC_CO2_23m,ppm,CO2 concentration measured 23 meters a.g.l. at CTC site.
    CTC_WS,m/s,Wind speed measured 23 meters a.g.l. at CTC site.
    CTC_WD,degrees,Wind direction measured 23 meters a.g.l. at CTC site.
    CTC_SO2,ppb,SO2 concentration measured at CTC site.
    CTC_CO,ppm,CO concentration measured at CTC site.
    CTC_NO2,ppb,NO2 concentration measured at CTC site.
    CTC_NO,ppb,NO2 concentration measured at CTC site.
    CTC_O3,ppb,Ozone concentration measured at CTC site.
    Birch_Hill_O3,ppb,Ozone concentration measured 155 meters a.g.l. at Birch Hill site.

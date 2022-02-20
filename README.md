# PulseWaveform
An automated pipeline for the extraction of PPG waveform features, including the Hybrid Excess and Decay (HED) model for extraction of component wave features.


```{r setup, include=FALSE}
```

# Description

This is a ReadMe for the PulseWaveform package and complementary analysis scripts as published in []. PulseWaveform contains a collection of functions for the processing and analysis of photoplethysmography (PPG) data. The General purpose script, located within the 'Scripts' folder in this repository, provides a foundational pipeline structure for implementing the package's functions to generate both morphological and HED model outputs that characterise the PPG pulse waveform. What follows is an expanded overview of the general analysis pipeline, in greater depth than was possible in the paper. Further details on assumptions made in constraining model behaviour are also included. Further clarification on code specifics can be found by examining the PulseWaveform.R file, which includes all functions with detailed comments. For any further queries beyond this, contact stw32@cam.ac.uk. 

What is the Pulse Waveform? The fundamental signal inherent to PPG, reflecting the peripheral changes in blood volume and pressure corresponding to each heart beat. Waveforms vary in their shape depending on a number of psychological and physiological factors; a detailed characterisation of these changes in morphology is the primary purpose of this repository. 

What is the HED Model? The Hybrid Excess and Decay Model is an attempt to parameterize the waveform such that it may bestow more useful information than simple descriptive measures of waveform morphology alone. It assumes each waveform is composed of component waves that suprimpose to produce the overall morphology (excess element), and that each component wave also has an exponential decay (decay element). 

## Prerequisites
Install PulseWaveform: devtools::install_github(repo = 'stw32/PulseWaveform')

The PulseWaveform package and associated pipelines make use of the following packages:  
- library(tidyverse)
- library(splines2)  
- library(pracma)  
- library(SplinesUtils)    
- library(spectral)
- library(zoo)
- library(readr)
- library(PulseWaveform)

# Analysis Pipeline:
The general analysis pipeline's structure is illustrated below. Pulse decomposition modelling (and subsequent outputs) are optional and can be bypassed if only descriptive morphological features are desired. 

<img width="909" alt="Screenshot 2022-02-20 at 16 28 36" src="https://user-images.githubusercontent.com/63592847/154853033-c604131b-e1b1-47b9-bde7-d6398217d298.png">

## Starting Parameters:
Once all relevant packages are installed, starting parameters must be specified before running the general purpose script. These include: 

1. Information about the PPG time series data to be inputted:
  - Sample rate
  - Peak finding threshold (follow instructions within GeneralPurposeScript to set this correctly)
  
2. Preferences for pulse decomposition modelling:
  - To run the HED model or not
  - The degree of model accuracy desired (number of simplex iterations)
  - To model each waveform individually, or to fix certain parameters across waveforms (by specifying batch numbers of >1)
  - To model only a subsection of the time series, or the entire time series (`all_beats <- FALSE` or `TRUE`)

## Read in Data:
Specify the file path of the PPG time series data to be analysed.

## Preprocessing:
Two different preprocessing routines are possible with PulseWaveform:
  - Bioradio preprocessing: downsampling, correction of detrending, and removal of baseline wander (`preproc`)
  - GE MR750 MRI scanner (ISO study) preprocessing: adjustment of successive decays in amplitude between datapoints, correction of entire time series gradient, and removal of baseline wander* (`FindUndetrendingParams` and `UnDetrend`)
  
An example of an individual waveform from the ISO dataset, before and after preprocessing, is given below. New time series data from other hardware sources will require unique preprocessing solutions, though the above routines may be useful as starting points. 

<img width="889" alt="Screenshot 2022-02-20 at 16 56 23" src="https://user-images.githubusercontent.com/63592847/154854317-8460e4f9-6c16-4b0b-96c3-0f47c14a582c.png">

*baseline wander, whilst considered preprocessing, actually occurs later in the general purpose script, once trough points are known, using `baseline`. The below figure illustrates removal of baseline wander:

<img width="940" alt="Screenshot 2022-02-20 at 17 20 49" src="https://user-images.githubusercontent.com/63592847/154855404-d6c94a1a-ef02-478c-975e-dab8bc0da43d.png">

## Beat Segmentation
We interpolate a cubic spline of the provided data points 
`sfunction <- splinefun(1:length(undetrended), undetrended, method = "natural")`
and take its first derivative
`deriv1 <- sfunction(seq(1, length(undetrended)), deriv = 1)`. 
Then the maximum of the first derivative, called the W point(s) are found:  
`w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly, sr = samplingRate)`
![](Wpoint.png)  
Beat segmentation is based on the correct identification of the W points.  
Next U, V and O are detected.
Based on the O points we remove a baseline in order to filter out low frequency modulation of the signal without altering the morphology of individual waveforms. 
![](baseline.png)  
 
The `sep_beats` function then takes care of some rudimentary cleaning procedures - plotting options of the outliers can be added by setting `q = TRUE`:   
- long waves (where a distance of more than two waves (Os) was counted as one) are removed  provided the waves to either side of them are not similarly long (if so there is a plotting option to make a manual decision)
- abnormally high amplitude
- waves that have a second systolic upstroke are removed  
- waves that fall below the O-O threshold (currently hardcoded to 4?) are excluded  
- waves that deviate more than 2 SDs from the average in their residuals
- waves that deviate more than 5 standard deviatiosn from the average wave (calculated after all of the above are already excluded) are removed  
- when `subset == T` is only necessary for the identification of a subset of data with an autonomic arousal manipulation (based on inter.beat intervals)  

In the next step, in the `FindStartParams` function each beat is modelled as a composition of three peaks (S, D and N) with parameters determining the amplitude, timing and width of each within the dataframe `beat`: `STime`, `SAmplitude` `SWidth` and so forth. Initial parameters are estimated using the excess of each beat which is what is left once the assumed decay towards a variable baseline is factored out:     
`excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))`

The width of each peak is initialised at 0.25 i.e. `par[4] <- 0.25`. 


## Improving fit using downhill simplex

To estimate the timing of the first and second reflectance peaks, we take the median of the time values generated by the initial parameter estimation.

`  renal_param <- median(beat$NTime)`

`  dias_param <- median(beat$DTime)`

Then we begin to refine the parameters originally generated by `FindStartParams`. We extract both within and across beat parameters. 
The across beat parameters extracted are those that are expected to be held constant across beats. These are:

- The width of S, R1 and R2
- The timing for R1 and R2
- The config rate

These parameters are fed into the `simplex.MakeSimplex2` function. This function slightly improves each of the across beat parameters, and outputs a suitable matrix to feed into the downhill simplex algorithm. For further explanation about the simplex creation and algorithm, please refer to Numerical Recipes (Press, Teukolsky, Vetterling & Flannery, 2007).

The within beat parameters produced are those which are able to vary from beat to beat. These are:

- The timing of S
- The amplitude of S, R1 and R2
- The two baseline values 

These parameters are fed into the `FindWithinParams` function, which subsequently feeds into a simplex. This improves each of the within beat parameters and outputs a suitable matrix for the downhill simplex algorithm. The distinction between this and the across beat parameters is that the within beat parameter values are calculated for each individual beat, whereas the across beat parameters are calculated for each batch.

The within and across beat parameters are combined into a matrix, and then fed into the downhill simplex algorithm, which iteratively improves the parameter estimation by minimizing the chi-sq value. The improved parameter values are extracted using `extractOutput`. Finally, `fixOutput` is utilised to ensure that there are no excessive deviations of the new parameters from our original estimations.

This process is then iterated upon four times, to ensure a the simplex does not get stuck in local minima. 


## Penalties applied to the chi-sq value
Penalties are applied to the chi-sq value if any of the below conditions are met: 

- If any peak amplitudes are estimated as negative 
- If the peak width is estimated as <0.05s or >0.5s (units?) 
- If the width of the R1 peak is <0.1 or >0.25
- If the diastolic width is >0.45
- If the renal peak timing strays too far from initial parameter estimation
- If the S-peak deviates too far from the maxima of the beat
- If the timing of the diastolic is < 0.2
- If there is less than 0,1s between the peaks
- If there is a large shift between baselines (Currently not applied!)
- If the configuration rate is above 0.95
- The renal peak is gradually penalised as the amplitude increase

# Contributing 
We would love to hear from people who would like to contribute or have ideas for developing out model further. 

# PulseWaveform
An automated pipeline for the extraction of PPG waveform features, including the Hybrid Excess and Decay (HED) model for extraction of component wave features.


```{r setup, include=FALSE}
```

# Description

This is a ReadMe for the PulseWaveform package and complementary analysis scripts as published in []. PulseWaveform contains a collection of functions for the processing and analysis of photoplethysmography (PPG) data. The General purpose script, located within the 'Scripts' folder in this repository, provides a foundational pipeline structure for implementing the package's functions to generate both morphological and HED model outputs that characterise the PPG pulse waveform. What follows is an expanded overview of the general analysis pipeline, in greater depth than was possible in the paper. Further clarification on code specifics can be found by examining the PulseWaveform.R file, which includes all functions with detailed comments. For any further queries, please contact stw32@cam.ac.uk. 

What is the Pulse Waveform? The fundamental signal inherent to PPG, reflecting the peripheral changes in blood volume and pressure corresponding to each heart beat. Waveforms vary in their shape depending on a number of psychological and physiological factors; a detailed characterisation of these changes in morphology is the primary purpose of this repository. 

What is the HED Model? The Hybrid Excess and Decay Model is an attempt to parameterize the waveform such that it may bestow more useful information than simple descriptive measures of waveform morphology alone (Williamson et al). It assumes each waveform is composed of component waves that suprimpose to produce the overall morphology (excess element), and that each component wave also has an exponential decay (decay element). 

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

To extract fiducial point data the pre-processed time series must first be segmented into individual waveforms. This entails peak identification followed by identification of troughs, which mark the start and end of each beat. To ensure accuracy, inflection points crucial to peak and trough detection are identified from cubic splines of the data. The steps are as follows:

  - Cubic Spline Interpolation: 
    - A cubic spline is generated from existing data points `sfunction <- splinefun(1:length(undetrended), undetrended, method = "natural")`
    - The first derivative of the spline is then calculated `deriv1 <- sfunction(seq(1, length(undetrended)), deriv = 1)`. 

  - Peak Detection: The `find_w` algorithm identifies peaks in the first derivative of the time series, denoted w. An adaptive rolling window is employed to identify successive peaks in the time series. Within the window, inflection points are confirmed as peaks if they exceed thresholds relative to values within the window and across the entire time series. The size of the window is determined from previous inter-beat intervals and is therefore adaptive to changes in heart rate. It is also adaptive in its ability to skip forward (over artefactual regions) or widen in the forward direction (if inadequate beats detected). Artefacts are also detected by the window, again according to local and global thresholds. Key assumptions made by the algorithm are as follows: 
    - (i) The duration of an inter-beat interval can be reasonably approximated from the previous inter-beat interval. 
    - (ii) All peaks are inflection points in the data above the 95th percentile for amplitude within the window. 
    - (iii) Peaks with extreme amplitude values relative to adjacent peaks and the entire time series should be labelled as artefactual. 
    - (iv) Extended windows with no identifiable peaks are likely low-amplitude artefactual regions. 
    - (v) Peaks immediately preceding / following artefacts are likely to also be compromised and should therefore be removed.
    
Identified peaks in the first derivative of a PPG time series are shown below. Peaks detected by the algorithm are shown as black circles. Prominent secondary peaks, precluding the use of a simple objective threshold, are indicated by red arrows.

<img width="818" alt="Screenshot 2022-02-21 at 12 45 09" src="https://user-images.githubusercontent.com/63592847/154957759-042ca7b0-3470-4804-a5fe-f46318b915c6.png">

  - Trough Detection: The minima preceding w, denoted O, are calculated as the inflection points preceding w, using `find_o`. In the absence of an inflection point, O points are derived from inflection points preceding w in the first derivative. The removal of baseline wander as mentioned above is achieved through cubic spline interpolation of O points. The constructed spline is then subtracted from the signal.
  
On completion of peak and trough detection, inter-beat intervals between adjacent w peaks are used to calculate heart rate and heart rate variability. O points defining the start and end of each beat are then used to derive individual beat segments. Once segmented, waveforms are scaled for amplitude (an optional feature) and aligned in time by w. The figure below shows a sample of aligned beats with the average (mean) waveform overlayed in red. 

<img width="710" alt="Screenshot 2022-02-21 at 18 26 31" src="https://user-images.githubusercontent.com/63592847/155009496-55803e90-c3bd-4114-b5c9-f10a4244acd7.png">

## Cleaning

Waveforms heavily influenced by noise are more likely to yield inaccurate features and should be removed. The peak detection algorithm (`find_w`) outlined above excludes major artefacts but is unlikely to detect subtler deformations in morphology due to less pronounced movements (e.g speaking). Therefore, a cleaning step to assess each waveform’s quality is included after beat segmentation and before feature extraction.

The approach to assessing waveform quality follows closely that of Orphanidou et al. The main principle is to utilize the average waveform of a sample, to which all individual waves can be compared*. The `sep_beats` function is applied to exclude beats with the following features: 
  - (i) Wave length: wave segments that are excessively long or short relative to adjacent waveforms and the average
  - (ii) Two peaks: segments with two maxima (or a rising tail portion) suggestive of two systolic peaks (e.g ectopic beats)
  - (iii) Values below baseline: segments with values significantly below the corrected baseline
  - (iv) Standard Deviations from the Average: Segments with values exceeding a defined number of standard deviations from the average
  - (v) Standard Deviation of Residuals: A set of residual values is first calculated by subtracting a given waveform from the average. The standard deviation of all residual values is then calculated. Beat segments with high residual values due to physiological change tend to have consistent residual values, whereas beat segments with high residual values due to noise tend to have inconsistent residual values. Therefore, waves with noise can be removed more selectively.
  
The cutoff thresholds for the above criteria were determined empirically and may need adjusting for new datasets. `PlotRejects` plots the number of beats rejected for each criterion and can be used at the end of the pipeline to assess if adjustments to thresholds are needed. Below is an example of an aligned set of segmented beats before and after cleaning.

<img width="881" alt="Screenshot 2022-02-21 at 18 27 07" src="https://user-images.githubusercontent.com/63592847/155009585-2ce4b92e-e576-491b-a8d5-d421c0f9c651.png">

* setting `subset == T` for `sep_beats` allows only a subset of beats to be taken forward in the analysis pipeline. It is only necessary for autonomic arousal manipulations where inter-beat intervals are altered (such as the ISO study). 

## Morphological Feature Extraction

  - Fiducial Point Identification: To obtain morphological features, key fiducial points OSND must first be identified on each waveform. The first derivative is utilized for this purpose. To increase robustness, OSND values for the average wave of a sample are identified first (`osnd_of_average`) and inform identification of points on individual waves. For each wave, values are identified in the following order: 
     - (i) N: The notch is assumed to occur between S and D as identified on the average wave. Within this section, inflection points in the first derivative are identified. The point with the greatest amplitude is selected, and the immediately preceding inflection point on the original trace is labelled as N. In the absence of a notch, the inflection point closest to zero in the first derivative is instead identified, and the corresponding point on the original waveform is taken as N.
     - (ii) D: The peak following N is taken as D. If no peak follows N, D is assumed to have the same value as N (i.e. the waveform is presumed to be class 3 or above).
     - (iii) O: The first value of the waveform, as defined during beat segmentation. In case of error during beat segmentation, a check for inflection points preceding w of lower amplitude is performed. O is reassigned if one is present.
     - (iv) S: The maximum of the waveform is taken as S. In the case of a maximum that does not correspond to S (i.e. due to reflectance wave augmentation in late systole), the timing from w to the maximum is calculated. If the timing exceeds an empirical threshold, S is instead defined as the point occurring a set time after w.
  
  - Feature Extraction: Once OSND points have been identified, morphological features derived from them can be extracted using `feature_extract`. Table 1 outlines these, as well as additional features that are not solely calculated from fiducial points. The list is not extensive but does offer a comprehensive view of waveform morphology sufficient to quantify morphological differences between groups. Since fiducial points are already identified at this stage, additional measures can be incorporated with relative ease.
  
<img width="542" alt="Screenshot 2022-02-21 at 18 31 14" src="https://user-images.githubusercontent.com/63592847/155010028-ffaa21f3-cf31-4277-9063-5f7a10ec27b8.png">
  
## Pulse Decomposition Modelling: The HED Model

### Theoretical Foundations and Model Components:

A single prevailing physiological view on the nature of arterial waves and reflections has not yet been reached, largely due to the inherent complexity of the arterial tree. Nonetheless it is generally accepted that backward travelling reflectance waves exist and are likely to be measurable as composites from multiple reflection sites, irrespective of origin. Hitherto, three component waves of the pulse waveform have been described:
 - (i) A systolic wave caused by the increase in pressure arising from left ventricular contraction and ejection of stroke volume. 
 - (ii) A first reflectance wave, also described as the ‘renal’ wave’.
 - (iii) A second reflectance wave, often referred to as the ‘diastolic’ wave.
  
PDA models have taken a largely data-driven approach to decomposing the waveform, using combinations of component waves of varying number and shape. Comparisons of these suggest three to be the optimal number for capturing the waveform’s morphology (Tigges et al, 2017). The current model starts from this point, modelling the waveform as a composite of three component waves: 
 - (i) systolic; 
 - (ii) first reflectance; 
 - (iii) second reflectance. 
    
Nine parameters are used to model the timing, width and amplitude of each component wave, as shown below. The current model, however, differs from existing PDA models in its approach to modelling the diastolic decay. A tendency for waves of prolonged duration to decay below baseline, an established phenomenon in aortic studies (Parker, 2017), was noted in the PPG signal. This is not considered by current PDA models. To account for this and better elucidate the underlying factors driving waveform morphology, a further three parameters are incorporated:
 - (i) Decay rate: The rate at which the signal decays exponentially to baseline in the absence of component wave influence.
 - (ii) Baseline 1: The baseline towards which the signal decays during the systolic portion of the waveform. 
 - (iii) Baseline 2: The baseline towards which the signal decays during the diastolic portion of the waveform.
    
Overall, therefore, the waveform is modelled as a composite of signal due to the initial systolic pressure wave, an exponential decay, and two reflectance waves. For simplicity, the systolic and reflectance waves are henceforth referred to as the excess element, and the exponential decay as the decay element. The output of the model is a 12-parameter vector for each waveform, which can be used to construct a modelled wave fitting the original PPG data.
    
<img width="898" alt="Screenshot 2022-02-21 at 18 40 30" src="https://user-images.githubusercontent.com/63592847/155011100-33f58d5e-bf06-4997-84b8-59cabf55e647.png">

A. The first 9 parameters of the model compose the excess element. TS = Timing of systolic wave; TR1 = Timing of 1st reflectance wave; TR2 = Timing of 2nd reflectance wave; AS = Amplitude of systolic wave; AR1 = Amplitude of 1st reflectance wave; AR2 = Amplitude of 2nd reflectance wave; WS = Width of systolic wave; WR1 = Width of 1st reflectance wave; WR2 = Width of 2nd reflectance wave. B. The final 3 parameters compose the decay element, which when incorporated with the excess yield the final modelled waveform. D R= Decay rate; B1 = 1st baseline; B2 = 2nd baseline.

Below are three Bioradio-acquired waves fitted with the HED model (using `GGplotFits`). Each wave is recorded from a different individual and is of a different class. 

<img width="903" alt="Screenshot 2022-02-21 at 18 41 47" src="https://user-images.githubusercontent.com/63592847/155011236-06b140ba-2a33-4ece-b206-197e3602ea6a.png">

A. class 1, NRMSE 0.93; B. class 2, NRMSE 0.94; C. class 3, NRMSE 0.85 (values to two decimal places).

### Model refinement: 
In refining the model's behaviour certain constraints were imposed. These were done according to the following assumptions:
 - (i) component waves cannot have the same timing and must occur in order
 - (ii) component waves cannot have negative amplitudes
 - (iii) the first reflectance wave plays a relatively minor role in shaping the waveform and must therefore have the smallest amplitude and width
 - (iv) the number of component waves need not be three if a more parsimonious fit can be achieved with fewer.

In the code, this manifests as penalties applied to goodness of fit (ChiSq) values when any of the following conditions are met within the `model2.FIX_PAR3` function:
 - If any peak amplitudes are estimated as negative 
 - If the peak width is estimated as <0.05s or >0.5s 
 - If the width of the R1 peak is <0.1 or >0.25
 - If the diastolic width is >0.45
 - If the renal peak timing strays too far from initial parameter estimation
 - If the S-peak deviates too far from the maxima of the beat
 - If the timing of the diastolic is < 0.2
 - If there is less than 0.1s between the peaks
 - If there is a large shift between baselines
 - If the configuration rate is above 0.95
 - The renal peak is gradually penalised as the amplitude increase
 
### Parameter Optimisation:
In modeling each waveform, a set of starting parameters are required to begin the fitting process. These are first estimated using `FindStartParams`. An estimated decay element is subtracted from each waveform to yield an 'excess' (see below). Once the excess is defined, peaks in the excess are searched for within empirically defined ranges. Values for width, timing and amplitude for each peak (the first 9 parameters) can then be estimated. The final 3 parameters are inputted as constants. 

<img width="1091" alt="Screenshot 2022-02-21 at 19 06 34" src="https://user-images.githubusercontent.com/63592847/155013849-219d353f-6f5e-48b5-9fca-79642f623008.png">

Starting parameters are then iteratively improved upon using `simplex.MakeSimplex2` for parameters fixed across beats:
- The width of S, R1 and R2
- The timing for R1 and R2
- The config rate

and `simplex.MakeSimplex3` and `FindWithinParams` for parameters with unique values for each beat:
- The timing of S
- The amplitude of S, R1 and R2
- The two baseline values 

These improved parameters are combined into a matrix (`make_matrix`) for the next stage. The measure of goodness of fit by which parameters are improved is the reduced ChiSq function:

<img width="788" alt="Screenshot 2022-02-21 at 19 20 00" src="https://user-images.githubusercontent.com/63592847/155015221-c031a5f6-ed2b-47c7-85f8-a9cbd2e10241.png">

Residuals are calculated by subtracting a modelled version of the waveform, constructed with the parameter set, from the original. Central morphology (from w to D) residuals are weighted more heavily, such that the model tends to fit key features like the systolic decline and notch well without overfitting noise in the tail of diastolic decay.

Finally, the parameters are optimised using a downhill simplex routine (`simplex.Run2`). A simplex is a geometric structure comprised of n + 1 vertices when generated in n-dimensional space. Thus it takes the form of a line in one dimension, a triangle in two dimensions, and a 13-vertex structure in 12 dimensions. In the case of the current model, each vertex of the simplex is defined by a different set of the 12 model parameters. 

For each iteration of the simplex, ChiSq values are first determined for each vertex. The vertex with the highest value (poorest fit) is identified, and the structure of the simplex altered so as to reflect (or, alternatively, shrink) away from it. By reflecting away from the highest vertex (worst set of parameter values) over several iterations, the simplex moves ‘downhill’ until a minimum point is reached, where further reflections do not result in significant improvements in fit. The end result is a set of parameters optimized for goodness of fit. The parameter values are then extracted using `extractOutput`. Finally, `fixOutput` is utilised to ensure that there are no excessive deviations of the new parameters from initial estimations and established constraints. The simplex is restarted four times to mitigate the risk of convergence on local minima in multi-dimensional space. 

### Assessment of Model Performance:
Following the completion of the downhill simplex routine, goodness of fit is calculated according to additional measures using `model2.ChiSq4`:

Normalized Root Mean Square Error (NRMSE) (Sorelli et al):

<img width="887" alt="Screenshot 2022-02-21 at 19 24 20" src="https://user-images.githubusercontent.com/63592847/155015675-cf987aea-136a-4e19-a5aa-8f5c4131cb4d.png">

An NRMSE value of 1 indicates a perfect fit to the data, whilst a value of 0 indicates no improved performance over the null model.

An alternative calculation of NRMSE (aNRMSE) (Wang et al):

<img width="927" alt="Screenshot 2022-02-21 at 19 25 03" src="https://user-images.githubusercontent.com/63592847/155015760-0f711ec7-78de-46f2-9de2-44e1fcc7e9d6.png">

The calculated aNRMSE value is an expression of the residual error as a percentage of the data. Values less than 2% are generally considered acceptable.

Furthermore, `osnd_fit` can be used to assess the HED model's ability to recapitulate relevant morphology as fiducial points O, S, N and D.

## Outputs
Outputs from the analysis pipeline can be broadly classified as:
  - Morphological feature outputs
  - HED model parameter outputs
  - HED Model Performance 
  - Inter-beat intervals 
  - Metadata (number of beats in sample, number of beats rejected through cleaning)
  - Processed PPG data (individual waveforms, OSND values) 
  
They are stored within the `AllOutputs` list, which acts as the input to the ISO_study_analysis_script within the Scripts folder, though custom analysis scripts will likely be required for novel data. 


# Contributing 
We would love to hear from anyone who would like to contribute to PulseWaveform, or who has ideas for developing the HED model further. 

# Version Updates
This repository is maintained by stw32. Further methods for extracting meaningful features from the PPG pulse waveform may be added in due course. 

# References
Orphanidou, C. et al. Signal-Quality Indices for the Electrocardiogram and Photoplethysmogram: Derivation and Applications to Wireless Monitoring. IEEE J. Biomed. Health Inform. 19, 832–838 (2015).

Parker, K. The reservoir-wave model | Atlantis Press. https://www.atlantis-press.com/journals/artres/125924989 (2017).

Sorelli, M., Perrella, A. & Bocchi, L. Detecting Vascular Age Using the Analysis of Peripheral Pulse. IEEE Trans. Biomed. Eng. 65, 2742–2750 (2018).

Tigges, T. et al. Model selection for the Pulse Decomposition Analysis of fingertip photoplethysmograms. in 2017 39th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC) 4014–4017 (2017). doi:10.1109/EMBC.2017.8037736.

Wang, L., Xu, L., Feng, S., Meng, M. Q.-H. & Wang, K. Multi-Gaussian fitting for pulse waveform using Weighted Least Squares and multi-criteria decision making method. Comput. Biol. Med. 43, 1661–1672 (2013).

Williamson ST, Daniel-Watanabe L, Finnemann J, Powell C, Teed A, Paulus M, Khalsa SS, Fletcher PC. 2021 Characterising the Photoplethysmography Pulse Waveform for Use in Human Neuroscience: The Hybrid Excess and Decay (HED)
Model. bioRxiv: http://doi.org/10.1101/2021.08.19.456935.


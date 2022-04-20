# 4DFlowPWVTool
a Matlab-based app for 4D flow MRI measurement of pulse wave velocity

This tool is based off a previously published post-processing software. Use of this tool should be associated and cited with the following reference: 
E. Schrauben, A. Wåhlin, K. Ambarki, E. Spaak, J. Malm, O. Wieben, A. Eklund. “Fast 4D flow MRI intracranial segmentation and quantification in tortuous arteries.” JMRI 2015, doi: 10.1002/jmri.24900.

**Installation:**
This app is built using Matlab's appdesigner functionality. It was built using Matlab version 2019a; it is therefore recommended to use this version or **newer** for running this app. Download and unzip the files. All functionality of the app is built into PulseWaveVelocityTool.mlapp and the corresponding utilities/ subfolder.

Additionally, the use of this software requires the following Matlab toolboxes:
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

**Data needed:**
This version of the app works directly on reconstructed dicom files from a GE 4D flow acquisition. It assumes these are all located within a top folder. The subfolders that contain the dicoms are assumed to have the following names: 
ANATOMY (containing the time-resolved **magnitude** data)
AP FLOW (containing the time-resolved **ap velocity** data)
LR FLOW (containing the time-resolved **lr velocity** data)
SI FLOW (containing the time-resolved **si velocity** data)

**Getting Started:**
Within Matlab, navigate to the folder containing PulseWaveVelocityTool.mlapp. Type this into your command window:
PulseWaveVelocityTool

Further details for each step within the app can be found in the corresponding user manual: 4DFlowPWVTool User Manual.pdf

# 4DFlowPWVTool
a Matlab-based app for 4D flow MRI measurement of pulse wave velocity in the aorta

This tool is based off a previously published post-processing software. Use of this tool should be associated and cited with the following reference: 
E. Schrauben, A. Wåhlin, K. Ambarki, E. Spaak, J. Malm, O. Wieben, A. Eklund. “Fast 4D flow MRI intracranial segmentation and quantification in tortuous arteries.” JMRI 2015, doi: 10.1002/jmri.24900.

Adding tool to Matlab path All files for the Centerline GUI are located in within a single folder that need to on the local machine. To add to Matlab’s current path, simply locate this folder, and type this into you command window:

Data needed

The app works directly on reconstructed .par / .rec files from a Philips 4D flow acquisition. It assumes these are all located within a single folder. The files needed are as follows: SUBJECTID or other identifier_1.par SUBJECTID or other identifier_1.rec SUBJECTID or other identifier_2.par SUBJECTID or other identifier_2.rec SUBJECTID or other identifier_3.par SUBJECTID or other identifier_3.rec

Getting Started
Type into Matlab’s command window:

Click ‘Load Data’. You will be prompted to find the folder your data are located in. After loading completion, the tool will automatically display axial, coronal, and sagittal calculated maximum intensity correction. It also automatically thresholds this image to produce a 3D isosurface of vessels within the field of view.

If predefined manual segmentation of the aorta has already been performed, click 'Load Segmentation'. Navigate to the folder containing segmentation dicoms

Data can be cropped and images can be updated with appropriate buttons.

Additional automatic velocity unwrapping is performed with 'Unwrap Velocity'

After pressing ‘Calculate Centerline and PWV tool’, the vasculature is skeletonized and individual vessel centerlines are extracted and labelled with unique identifiers. A new window will pop-up, which corresponds to all remaining steps for pulse wave velocity calculation.

** TO DO **

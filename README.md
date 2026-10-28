# Aggregated Data and Code for Fraser Fir and Machine Learning Article 
Repository curated by William Parrott, data and code access.

## Project Overview
*Proof of Concept*  
Following methods from Fearer et al. 2022, supervised classification pipelines employing support vector machine (SVM) and random forest (RF) models were used to identify Fraser fir that were resistant to infestation by the balsam woolly adelgid using NIR spectra. Three measurements of reflectance at 257 bands in the NIR range were taken from a total of 99 trees with known resistance levels, then were averaged by tree ID to ensure independence in training-testing splits<sup>1</sup>. 33 trees showed no resistance, 32 showed medium resistance, and 34 showed high resistance. Resulting  spectra were transformed using a second derivative Savitzky-Golay transformation, then were split into 70% training and 30% testing sets. SVM and RF were conducted on the full dataset and important variables for these classifications were identified; reflectance from important wavelengths was then used for a new training/testing split and subsequent RF and SVM models. A selection of every fifth wavelength was used to establish another training/testing split for an SVM. Finally, sparse partial least squares discriminant analysis was used to identify the loadings and stability of important wavelengths<sup>1</sup>.

<sup>1</sup> This was not originally done and was added/corrected in later runs. 

*Protein Analysis*  
Wavelengths related to proteins were found to be important in model performance. Phloem samples were also collected from trees with known resistance during sample collection in June and August of 2024. These samples were used for protein analysis with the Bradford assay. Resulting data were used to compare protein concentrations between resistance levels across the growing season.

*Field Model Performance*  
Partial least squares discriminant analysis (PLS-DA) models were trained with collected NIR spectra and used to attempt classification of both training and naive samples *in situ* in October 2026. Resulting calibration, cross validated accuracies, and evaluation accuracies were calculated from confusion matrices. Raw performance data and preprocessing hyperparameters are documented.


## Repository Structure
**1.** Archived: Data, code, and statistical analyses from the original runs that featured all spectra, unaveraged.  
**2.** Data: .csv files containing the phloem and needle spectral datasets.  
**3.** Scripts: .Rmd files containing the supervised classification pipelines. Filepaths are relative (./filepath) for ease of access.  
**4.** Protein Analysis: files containing raw and processed protein data, as well as the R script for statistical analysis on proteins.  
**5.** Field Validation Analysis: files containing data about field deployed models, including preprocessing values and performance.  
**6.** References: file containing all references for this project.  
**7.** README.md: The file containing this text.  
**8.** wavelengthref.txt: Reference file containing the 257 wavelengths that were measured to develop the spectra.  


## References (for the README):
Conrad, A. O., Li, W., Lee, D.-Y., Wang, G.-L., Rodriguez-Saona, L., and Bonello, P. (2020). Machine Learning-Based Presymptomatic Detection of Rice Sheath Blight Using Spectral Profiles. Plant Phenomics 2020, 8954085. doi: 10.34133/2020/8954085

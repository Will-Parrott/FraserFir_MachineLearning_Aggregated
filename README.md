# Aggregated Data and Code for Fraser Fir and Machine Learning Article 
Repository curated by William Parrott, data and code access for Frontiers in Plant Science.

## Project Overview
Following methods from Conrad et al. 2020, supervised classification pipelines employing support vector machine (SVM) and random forest (RF) models were used to identify Fraser fir that were resistant to infestation by the balsam woolly adelgid using NIR spectra. Three measurements of reflectance at 257 bands in the NIR range were taken from a total of 99 trees with known resistance levels, then were averaged by tree ID to ensure independence in training-testing splits<sup>1</sup>. 33 trees showed no resistance, 32 showed medium resistance, and 34 showed high resistance. Resulting  spectra were transformed using a second derivative transformation, then were split into 70% training and 30% testing sets.

<sup>1</sup> This was not originally done and was corrected in later runs. 

## Repository Structure
**1.** Archived: Data, code, and statistical analyses from the original runs that featured all spectra, unaveraged.
**2.** Data: .csv files containing the phloem and needle spectral datasets.
**3.** Scripts: .Rmd files containing the supervised classification pipelines. Filepaths are relative (./filepath) for ease of access.
**4.** README.md: The file containing this text.
**5.** wavelengthref.txt: Reference file containing the 257 wavelengths that were measured to develop the spectra.


References:
Conrad, A. O., Li, W., Lee, D.-Y., Wang, G.-L., Rodriguez-Saona, L., and Bonello, P. (2020). Machine Learning-Based Presymptomatic Detection of Rice Sheath Blight Using Spectral Profiles. Plant Phenomics 2020, 8954085. doi: 10.34133/2020/8954085

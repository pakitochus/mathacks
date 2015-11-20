MatHacks
========

MatHacks is a collection of miscellaneous files generated to help in 
Matlab (R) and OCTAVE computing. It currently comprises:

* **FROMSTRUCTTOLATEX** Creates an accuracy-sensitivity-specificity latex table from ```perfStruct``` (an struccture constructed with ```structConfMatErr```), which can be a single strucutre or an array of structures, along with their standard errors. 
* **MATRIXTOLATEX** converts a 1D, 2D or 3D matrix into a latex table that can be imported into whichever latex document. 
* **PCA_FS** performs a feature selection using PCA
* **RANKFEATURESINTERFACE** provides an interface to integrate filtering methods other than the five included in ```RANKFEATURES```, such as FDR or the p-value estimating method using ```ANOVA1```, and expanding the t-test using the two-tailed, two-sample t-test. 
* **SAVEFIGUREASGIF** Saves an animated 3D figure as a GIF 
* **SELECTC** Selects optimum C value in a linear kernel (using ```LIBLINEAR```) or another SVM kernel (using ```LIBSVM```). 
* **SELECTOPTIMUMC** Selection of optimum C by x-validation  Requires: ```LIBSVM```
* **SLICESDISPLAY** takes a 4-D image (with the first singleton dimension, such as the ones needed for the use within ```MONTAGE```) and returns an image showing the selected slices side by side. 
* **STRUCTCONFMAT** returns a performance structure, with the values of Correct Rate, Error Rate, Sensitivity, Specificity, npv, ppv, Positive Likelihood and Negative Likelihood. It uses a 2x2 confusion matrix, so no higher dimension or non-binary classifiers are allowed. 
* **STRUCTCONFMATERR** same as the previous one, but returning also the standard error of each field. 
* **SUMMARY** summarizes a matrix telling how many elements of each value there are. 
* **SUPERIMPOSEIMAGE** superimposes a color map ```MAP``` to a grayscale image ```IM``` using ```IM``` as an alpha map (transparency index), so that the final result is fake coloured. 
* **WRITE2FILE** Writes and appends lines to file. 

MatHacks
========

MatHacks is a collection of miscellaneous files generated to help in 
Matlab (R) and OCTAVE computing. It currently comprises:

* **MAPBRAIN** performs a proyection of a 3D brain (from a struct I that contains an field img) into a two dimensional plane, using the spherical coordinates and different approaches, as in <http://www.ncbi.nlm.nih.gov/pubmed/25488228>
* **MATRIXTOLATEX** converts a 1D, 2D or 3D matrix into a latex table that can be imported into whichever latex document. 
* **PCA_FS** performs a feature selection using PCA
* **RANKFEATURESINTERFACE** provides an interface to integrate filtering methods other than the five included in ```RANKFEATURES```, such as FDR or the p-value estimating method using ```ANOVA1```, and expanding the t-test using the two-tailed, two-sample t-test. 
* **SAVEFIGUREASGIF** Saves an animated 3D figure as a GIF 
* **SELECTC** Selects optimum C value in a linear kernel (using ```LIBLINEAR```) or another SVM kernel (using ```LIBSVM```). 
* **SELECTOPTIMUMC** Selection of optimum C by x-validation  Requires: ```LIBSVM```
* **SLICESDISPLAY** takes a 4-D image (with the first singleton dimension, such as the ones needed for the use within ```MONTAGE```) and returns an image showing the selected slices side by side. 
* **SUMMARY** summarizes a matrix telling how many elements of each value there are. 
* **SUPERIMPOSEIMAGE** superimposes a color map ```MAP``` to a grayscale image ```IM``` using ```IM``` as an alpha map (transparency index), so that the final result is fake coloured. 
* **WRITE2FILE** Writes and appends lines to file. 

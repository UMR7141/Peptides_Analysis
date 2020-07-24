# Peptides_Analysis

Code for paper 'Evidence supporting an antimicrobial origin of targeting peptides to endosymbiotic organelles', Clotilde Garrido, Oliver D Caspari, Yves Choquet, Francis-André Wollman, Ingrid Lafontaine

Analysing short peptides using z-scale factor (Hellberg 1987) that represent the physicochimical properties. 
Process auto-cross covariance, dimensional reduction and classification for vizulising and compare differents groups of peptides.

@autor: Ingrid Lafontaine, Oliver Caspari and Clotilde Garrido, IBPC UMR7141

#DEPENDENCIES
Python 3.5.2
The package requirements can be found in requirements.txt
You can use this file with virtalenv

First install virtualenv
> pip install virtualenv

Creat a new environnement
> virtualenv [name of new environement] --python=python3

Open environnement
> source [name of new environement]/bin/activate

Install all package need
> pip install -r requirements.txt

#INSTALLATION

Downlod repertory Peptides_Alalysis
Run Main.py
You can replicate our paper results or change data used in data.py
Data of our paper are in In_File folder

#FILE OVERVIEW
* Blots: original uncropped image files for Western Blots presented in the paper
* correlation quant: individual cell images and background fluorescence values used in the imagestat_PCC-bg.R script
* Epifluorescence originals: original epifluorescence microscopy image files encompassing the entire fields of view from which cells for presentation in the paper were cropped out
* ACC: empty folder, will contain backup csv files of auto-cross covariance calculated by the program
* In_File: Folder with data of the paper
	- HA-RAMP: Folder with csv files of HA-RAMP peptides by family and to class establish in paper
	- Other: Folder with csv files of globular AMP and random peptides.
	- SP: Folder with csv files of signal peptides 
	- TP: Folder with csv files of targeting peptides
	- ResultHeliquest : Foldes with csv files of amphiohilic helix prediction and script for generate them (Look start of script to view dependencies) 
* New3Helix: Folder with csv file of helix and peptide features and script for make distribution of theses features
* New_csv: empty folder, will contain backup csv files of peptides sequences used by the program
* Out_File: Folder were result file of program are registred
	- Out_Freq: result file of frequency amino acid analys
	- Out_Kmean: result file of k-mean analys
	- Out_Length: result file of peptid length analys
	- Out_PCA: result file of principal component analys
	- Out_Tree: result file of ?? analys
	- Out_Zscale: result file of mean z-scale analys
* Z_scale_Wold: empty folder, will contain backup csv files of mean z-scale calculated by the program
* confocal originals: original confocal microscopy image files encompassing the entire fields of view from which cells for presentation in the paper were cropped out
* acc.R: R program that calculate auto-cross-covariance
Can be use independly Command line : Rscript --vanilla acc.R -f <csv-file> -a <column-name> -l <lag> -o <out-file>
* data.py: Python script that indicate data use in the analyse (you can edit this file for use your data)
* functions.py : Python script contain function used in program 
* imagestat_PCC-bg.R: R sript using to calculate Pearson Correllation Coefficients for Venus compared to chlorophyll and mitotracker channels in batch
* Main.py: Python script for helping to use the needed function for analyse
* MultChi2.r: R script that make 
* requierments.txt: file contain all package need for the program
* settings.py: Python script containing some backup setting. 
* Z_scale_Wold.R: R script that calculate mean of z-scale
* 


# Methylation data analysis
R script for basic methylation data analysis

I arrange the methylation data analysis into 4 steps:

1. target file preparation

The target file is used to import raw data (idat data: ***Grn.idat and ***Red.idat) with R package "minfi".
The target file is a data frame including the following information:
Sample_Name, Sample_Well, Sample_Plate, Sample_Group, Pool_ID, Sentrix_ID, Sentrix_Position, and Basename. 
These information are usually included in sample sheet file from the lab taking the experiment. We can also include other phenotype information such as sex, age.

2. load methylation from raw data

The key function in this step is "loadData". The function import methylation data from raw file and take a basic QC check. 
The output value of the function is a list including a RGChannelSet(methylation information), predicted sex from methylation data, and detection p-value (small p-value indicates high quality).

3. quality control

Quality control and quantile normalization of the data. Low quality probes and samples are removed from the analysis. 
Check whether the predicted sex and age are the same or similar to the information in phenotype data. 

4. methylation analysis

Analyze the methylation data in probe, region and pathway level. 


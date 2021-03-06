The our training and testing datasets, and some datasets from other approaches are included in this folder.

training set:

1. genuine_premirna.fa:  
Pre-miRNA sequences of Viridiplantae were downloaded from the miRBase database (Release 21). 
After removing sequences containing non-AGCU characters, 3,044 known pre-miRNAs were collected as positive training dataset, 
which derived from 9 species, including Arabidopsis thaliana, Glycine max, Oryza sativa, Physcomitrella patens, Medicago truncatula, 
Sorghum bicolor, Arabidopsis lyrata, Zea mays and Solanum lycopersicum. 

2. negData.fa:  
Pseudo pre-miRNAs we selected from cDNAs.  
The CDS sequences of Arabidopsis thaliana, Medicago truncatula, Oryza sativa, Glycine max, and Populus trichocarpa were obtained from
 EnsemblPlants database (Release 18), and then fragmented into non-overlapped segments under a constraint condition 
that the length distribution of extracted segments was identical with that of known plant pre-miRNAs. 
Further, they should satisfy some criteria. The criteria are determined by observing real  plant pre-miRNAs. 
The criteria for selecting the pseud pre-miRNA are: minimum of 14 base pairings in the hairpins 
and maximum of −7.4 kcal/mol free energy of secondary structures (including GU wobble pairs)


testing set:

1. positive_testing_data.fa: 
Pre-miRNA sequences of Viridiplantae were downloaded from the miRBase database (Release 21). 
The remaining pre-miRNA sequences of Viridiplantae (excluding pre-miRNAs in "genuine_premirna.fa") were collected. 

2. negative_testing_data.fa:
The details about this dataset construction are the same with "negData.fa". But The sequences in two files are different. 


datasets_from_miPlantPreMat:
The training and testing datasets of miPlantPreMat, which included in the software of miPlantPreMat, were downloaded from website: https://github.com/kobe-liudong/miPlantPreMat.


datasets_from_PlantMiRNAPred:
The training and testing datasets of PlantMiRNAPred were downloaded from web site of PlantMiRNAPred: http://nclab.hit.edu.cn/PlantMiRNAPred/.

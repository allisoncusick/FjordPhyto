# FjordPhyto
FjordPhyto (https://fjordphyto.ucsd.edu/) started in 2015 to investigate how melting glaciers are affecting the phytoplankton along the western Antarctic Peninsula.

## Data Description

Data includes information on seawater temperature, conductivity, salinity, Secchi depth readings, euphotic depth, ocean color (RGB), glacial meltwater (oxygen isotope), phytoplankton abundance, carbon biomass, and species diversity (microscopy, metabarcoding). We have detected hundreds of species of phytoplankton in samples collected by travelers and we are in the process of creating a phytoplankton identification guide for use in the field.

### Data Processing

#### Metadata
If you would like to utilize the clean and processed metadata for use in other scripts see the file beginning with "fjord_phyto_processed" follwed by the most recent date in "mmddyyyy" format in the data folder.

If you would like to view the data in it's raw format and see the raw processing code then see "metadata_process-2024.r". 
The full metadata file (.xlsx) may be updated and redownloaded periodically so look for the file beginning with "FjordPhyto-MASTER-accessed_" follwed by the most recent date in "mmddyyyy" format. 

#### Metagenomic Data
##### DNA Extraction and Sequencing
Seawater collected using a phytoplankton net (20-micron) surface tow, from the western Antarctica Peninsula, filtered on 0.2-micron 47mm filters then preserved in RNAlater on ship. Samples are processed at J Craig Venter Institute in San Diego, CA, USA where genomic-DNA is extracted using the Machery-Nigel plant extraction kit and prepped and sequenced using Allen Lab 18SV9 and 16SV4-5 primer sets then submitted to the UCSD IGM Genomics center for Illumina MiSeq PE150 sequencing. 

##### Demultiplexing and Processing
Files are demultiplexed by IGM for i7, then further demultiplexed for i5 following this protocol: https://github.com/allenlab/i5_index_demultiplex
Fastq files are processed through the Qiime2 pipeline using these protocols: 
16S https://github.com/allenlab/QIIME2_16S_ASV_protocol 
18SV9 https://github.com/allenlab/QIIME2_18Sv9_ASV_protocol 






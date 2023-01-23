# FjordPhyto
Seawater collected using a phytoplankton net (20-micron) surface tow, from the western Antarctica Peninsula, filtered on 0.2-micron 47mm filters then preserved in RNAlater on ship. Samples are processed at J Craig Venter Institute in San Diego, CA, USA where genomic-DNA is extracted using the Machery-Nigel plant extraction kit and prepped and sequenced using Allen Lab 18SV9 and 16SV4-5 primer sets then submitted to the UCSD IGM Genomics center for Illumina MiSeq PE150 sequencing. 

Files are demultiplexed by IGM for i7, then further demultiplexed for i5 following this protocol: https://github.com/allenlab/i5_index_demultiplex
Fastq files are processed through the Qiime2 pipeline using these protocols: 
16S https://github.com/allenlab/QIIME2_16S_ASV_protocol 
18SV9 https://github.com/allenlab/QIIME2_18Sv9_ASV_protocol 

# Metagenome-analysis-pipeline


#ASSEMBLY_PIPELINE.sh
##############################################################################################################################################################
INPUT: sbatch ASSEMBLY_PIPELINE.sh assembly_output sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq ...
This script is meant to be a comprehensive solution for producing the best metagenomic assembly given paired end reads from one or more samples.
Ideally it should take in fully QC'd reads. First, the reads are assembled with metaSPAdes3.10, then all the reads that did not map back to the
contigs are re-assembled with MEGAHIT (which works better on lower coverage contigs. The resulting assemblies are combined, sorted, and short 
contigs are removed. The finall assembly is then QCed by QUAST. (Please indicate path to quast at the sart of script)

Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses.
For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.

##############################################################################################################################################################




#BINNING_PIPELINE.sh
##############################################################################################################################################################
INPUT: sbatch BINNING_PIPELINE.sh assembly.fa sampleA_1.fastq sampleA_2.fastq [ sampleX_1.fastq sampleX_2.fastq ]
This script is meant to be run on the outputs of the QC_PIPELINE and ASSSEMBLY_PIPELINE to split the assembly contigs into metagenomic bins.
Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.
The more samples, the better the binning. 

The script uses metaBAT to bin the contigs, then uses bwa to recruit reads back to the bins and ressembles them with SPAdes. It then uses KRAKEN 
to assign taxonomy to each bin and producing a joint kronagram of all the bins. Finally, it uses CheckM to test the contamination and completion
of the bins, and sepperates bins with \>20% completion and \<10% contamintation. 

NOTE: Samtools, CheckM, and KRAKEN require instalation, and be sure the configure the right path to the KRAKEN folder and the KRAKEN database.
Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses.
For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
##############################################################################################################################################################


#BLOBOLOGY_PIPELINE.sh
##############################################################################################################################################################

INPUT: sbatch BLOBOLOGY_PIPELINE.sh assembly.fa reads_1.fastq reads_2.fastq
This script is a modified pipeline from the program 'BLOBOLOGY'. The original script has been modified to be better suited for running on clusters.

Author of original pipeline: Sujai Kumar (https://github.com/blaxterlab). Author modifications: German Uritskiy. I do not take any credit for the original pipeline.
For questions, bugs, and suggestions, contact German Uritskiy at guritsk1@jhu.edu.

##############################################################################################################################################################



#KRAKEN_PIPELINE.sh
##############################################################################################################################################################
READS INPUT: sbatch KRAKEN_PIPELINE.sh output_folder sampleA_1.fastq sampleA_2.fastq [ sampleX_1.fastq sampleX_2.fastq ]
ASSEMBLY INPUT: sbatch KRAKEN_PIPELINE.sh output_flder sampleA.fasta sampleA.fa
This script is meant to be run on paired end reads (with extensions *_1.fastq and *_2.fastq) or assembled contigs (*.fa or *.fasta).
The script runs KRAKEN on the sequences, then translates them to taxonomy form with kraken-translate. Then in-house scripts are used to 
parse out the taxonomy into the format for KRONA-TOOLS, colapse the file to save memory, and finally produce a prety kronagram with all the files.

NOTE: KRAKEN KronaTools requires instalation, and be sure the configure the right path to the KRAKEN folder and the KRAKEN database.

Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses.
For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.

##############################################################################################################################################################


#QC_PIPELINE.sh
##############################################################################################################################################################

INPUT: sbatch QC_PIPELINE.sh reads_1.fastq reads_2.fastq
This script is meant to be a comprehensive solution to QC new HiSeq reads in preparation for assembly, and other operations.
The main things this pipeline accomplishes are read trimming based on quality scores, and removal of human sequences.
The script also produces a FASTQC report before and after the procedures.
Requires bmtagger to be installed in the PATH! (by conda for example).

Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses."
For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.

##############################################################################################################################################################


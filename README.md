# Metagenome-analysis-pipeline


  This script is meant to be a comprehensive solution for producing the best metagenomic assembly, taxonomic profiling, and binning from paired end reads from one or more samples. Of course this is not the best pipeline out there. But this is the pipeline taht worked for me so far.
  
  There are 5 major parts to the pipeline:
  
    1) Read QC
    
    2) Assembly
    
    3) Binning
    
    4) Taxonomic profiling
    
    5) Visualization with Blobplots
    
    
  ![Detailed pipeline walkthrough](http://i.imgur.com/NTmjezn.png)

  
  At this moment there are seperate scripts to accomplish these steps, as it is best to make sure each step(s) is correctly completed before going on to the next. In future iterations I wil also include a master script to run everything at once (at your own risk). I highly advise you to run these pipelines with a tiny subsections of your reads first to make sure it works to the end. There are many checkpoints throught the pipelines that will exit with an error message if the previous step failed.
  
  The biggest challenge will likely be to configure all the softwares correctly. The folder "SOFTWARE" has numerous scripts required for running this pipeline, so make sure to download it and set the right path to it in each script. The pipeline scripts are currently configured for MARCC cluster (as seen by the headers) but they are written in bash so they should run on any linux-based environment.
  
  System requirements: The resourse requirements for this pipeline wil vary greatly based on the number of your reads, but I would advise against attempting to run it on anything less than 10 cores and 100GB RAM.

  Author of pipeline: German Uritskiy. I do not clain to have any authorship of the many programs this pipeline uses. For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.

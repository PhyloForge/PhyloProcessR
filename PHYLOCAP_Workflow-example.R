#Libraries
library(rdrop2)


#################################################
## Step 0: Data download. Optional step.
##################

### Example usage
sample.spread = "/home/c111h652/scratch/MitoGenomes/Mitogenome_study.csv"
out.dir = "/home/c111h652/scratch/MitoGenomes/raw-reads-frogs"
dropbox.dir = "/Research/3_Sequence-Database/Raw-Reads"
dropbox.tok = "/home/c111h652/dropbox-token.RDS"

#Run download function
dropboxDownload(sample.spreadsheet = sample.spread,
                dropbox.directory = dropbox.dir,
                dropbox.token = dropbox.tok,
                out.directory = out.dir,
                overwrite = TRUE)


#################################################
## Step 1: Preprocess reads with read processing function
##################






#################################################
## Step 2: Clean out contamination
##################






#################################################
## Step 3: iteratively assemble mt-genomes
##################







#################################################
## Step 4: do all the alignment stuff
##################











# END SCRIPT


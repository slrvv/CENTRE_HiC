########################################################################
# Download HiC script
########################################################################

##Download all the HiC data in subsetted_HiC.tsv using the \
## encode_downloader.py script from Kundaje Lab
import os
import sys

name=sys.argv[1]
file = open(name, 'r')
dir =sys.argv[2]

for line in file.readlines()[1:]:
    accession_id = line.split(',')[0]
    print(accession_id)
    bio_name = line.split(',')[4]
    print(bio_name)
    print("Making the directories for: "+bio_name) 
    #Making the directory for each sample 
    #os.system("mkdir "+dir+bio_name)
    
    dir0 = dir+bio_name

    cmd1 = "python encode_downloader.py "+accession_id+" --file-types ""hic"" --dir "+dir0

    print("Downloading HiC maps ...")
    os.system(cmd1)
    print(" ")
    print("Finished downloading all files for : "+bio_name)
    print("Moving on to next entry")

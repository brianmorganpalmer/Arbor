#####################################################################################
#
# Author: Brian Palmer
# Date: Dec 10, 2018
#
# Abstract: Simple pipeline to evaluate Cas13d crRNA processing
# Provided by David Scott
#
#
####################################################################################

file317_1 = "../data/SRR6800317_1.fastq"
file317_2 = "../data/SRR6800317_2.fastq"
file318_1 = "../data/SRR6800318_1.fastq"
file318_2 = "../data/SRR6800318_2.fastq"

# read in the files
from Bio import SeqIO
records1 = list(SeqIO.parse(file317_1, "fastq"))
records2 = list(SeqIO.parse(file317_2, "fastq"))

# get the size of the sequence files
size1 = len(records1)
size2 = len(records2)

# make sure they're the same size
if size1 != size2:
    print("Error: File sizes are different...")

# go through and make sure there ids are the same
for i in range(size1):

    # need to normalize the names, just first 2 numeric digits
    sId1 = records1[i].id                   # get the id string
    idArray = sId1.split(".")               # split into an array
    sId1 = idArray[0] + "." + idArray[1]    # assign first 2 vals to id / name

    sId2 = records2[i].id                   # same for file 2
    idArray = sId2.split(".")
    sId2 = idArray[0] + "." + idArray[1]

    if sId1 != sId2:
        print("Error: Identifiers don't match...")

print("Execution complete for SRR6800317")
print("Read Ids are the same")
print(size1, " - Records processed")

###############################
#
# TBD: Duplicate code needs to be combined
#
# read in the files SRR6800318
#
from Bio import SeqIO
records1 = list(SeqIO.parse(file318_1, "fastq"))
records2 = list(SeqIO.parse(file318_2, "fastq"))

# get the size of the sequence files
size1 = len(records1)
size2 = len(records2)

# make sure they're the same size
if size1 != size2:
    print("Error: File sizes are different...")

# go through and make sure there ids are the same
for i in range(size1):

    # need to normalize the names, just first 2 numeric digits
    sId1 = records1[i].id                   # get the id string
    idArray = sId1.split(".")               # split into an array
    sId1 = idArray[0] + "." + idArray[1]    # assign first 2 vals to id / name

    sId2 = records2[i].id                   # same for file 2
    idArray = sId2.split(".")
    sId2 = idArray[0] + "." + idArray[1]

    if sId1 != sId2:
        print("Error: Identifiers don't match...")

print("Execution complete for SRR6800318")
print("Record Ids are the same")
print(size1, " - Records processed")

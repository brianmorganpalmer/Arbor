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

"""
A SeqRecord object holds a sequence and information about it.

Main attributes:
    id          - Identifier such as a locus tag (string)
    seq         - The sequence itself (Seq object or similar)

Additional attributes:
    name        - Sequence name, e.g. gene name (string)
    description - Additional text (string)
    dbxrefs     - List of database cross references (list of strings)
    features    - Any (sub)features defined (list of SeqFeature objects)
    annotations - Further information about the whole sequence (dictionary).
    letter_annotations - Per letter/symbol annotation (restricted dictionary).
      This holds Python sequences (lists, strings or tuples) whose length matches that of the sequence.
      A typical use would be to hold a list of integers representing sequencing quality scores, or a string
      representing the secondary structure. """




# read in the files
from Bio import SeqIO
records1 = list(SeqIO.parse("../data/SRR6800317_1.fastq", "fastq"))
records2 = list(SeqIO.parse("../data/SRR6800317_2.fastq", "fastq"))

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

print("Execution complete")
print("Record Ids are the same")
print(size1, " - Records processed")

#####################################################################################
#
# Author: Brian Palmer
# Date: Dec 10, 2018
#
# Abstract: Simple pipeline to evaluate Cas13d crRNA processing
# Provided by David Scott
#
"""
SequenceRecord Object
        seq - Sequence, required (Seq, MutableSeq or UnknownSeq)
        id - Sequence identifier, recommended (string)
        name - Sequence name, optional (string)
        description - Sequence description, optional (string)
        dbxrefs - Database cross references, optional (list of strings)
        features - Any (sub)features, optional (list of SeqFeature objects)
        annotations - Dictionary of annotations for the whole sequence
        letter_annotations - Dictionary of per-letter-annotations,
         values should be strings, list or tuples of the same length as the full sequence.
"""
####################################################################################
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def modifyReadId(seq_record):
    seq_id = seq_record.id
    seq_id_array = seq_id.split(".")               # split into an array
    new_seq_name = seq_id_array[0] + "." + seq_id_array[1]
    return new_seq_name

def reduceSeqString(seq_record):
    seq = seq_record.seq[0:30]
    if len(seq) < 30:           # discard small reads
        return None
    return seq

def createNewSquenceRecord(seq_record):
    # get new sequence values
    new_id = modifyReadId(seq_record)
    new_seq = reduceSeqString(seq_record)

    if new_seq is None:
        print("Warning: Sequence less than 30 characters")
        return None

    phred_quality = list(seq_record.letter_annotations["phred_quality"])
    phred_quality = phred_quality[0:30]

    # Phred score less than 30
    for score in phred_quality:
        if score < 30:
            print("Warning: PHRED scroe lower than 30 - discarding - ", score)
            return None

    # construct a new sequence record
    new_seq_record = SeqRecord(Seq(str(new_seq)), id=new_id, name=seq_record.name,
                               description=seq_record.description, dbxrefs=seq_record.dbxrefs,
                               features=seq_record.features, annotations=seq_record.annotations)

    print("Dude")
    return new_seq_record


# read in the files
records1 = list(SeqIO.parse("../data/SRR6800317_1.fastq", "fastq"))
records2 = list(SeqIO.parse("../data/SRR6800317_2.fastq", "fastq"))

# new sequence record arrays
new_seq_records1 = []
new_seq_records2 = []

# get the size of the sequence files
size1 = len(records1)
size2 = len(records2)

# make sure they're the same size
if size1 != size2:
    print("Error: File sizes are different...")

# go through and make sure there ids are the same
for i in range(size1):

    new_seq_id = modifyReadId(records1[i])              # get new id value
    new_seq_str = reduceSeqString(records1[i])          # get new sequence string
    # create new sequence record

    seq_record = createNewSquenceRecord(records1[i])
    # add to new record array
    if seq_record is not None:
        new_seq_records1.append(seq_record)

with open("317Update.fastq", "w") as output_handle:
    SeqIO.write(new_seq_records1, output_handle, "fastq")


print("Execution complete")
print("Record Ids are the same")
print(size1, " - Records processed")
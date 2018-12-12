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
from Bio.Alphabet import SingleLetterAlphabet

def mod_read_id(seq_record):
    seq_id = seq_record.id
    seq_id_array = seq_id.split(".")               # split into an array
    return seq_id_array[0] + "." + seq_id_array[1]

def reduce_seq_length(seq_record):
    seq = seq_record.seq[0:30]
    if len(seq) < 30:           # discard small reads
        return None
    return seq

def create_new_sequence_record(seq_record):
    # get new sequence values
    new_id = mod_read_id(seq_record)
    new_seq = reduce_seq_length(seq_record)

    if new_seq is None:
        print("Warning: Sequence less than 30 characters")
        return None

    phred_quality = list(seq_record.letter_annotations["phred_quality"])
    phred_quality = phred_quality[0:30]

    # Phred score less than 30
    for score in phred_quality:
        if score < 30:
            # print("Warning: PHRED scroe lower than 30 - discarding - ", score)
            return None

    # construct a new sequence record
    new_seq_record = SeqRecord(Seq(str(new_seq), SingleLetterAlphabet()), id=new_id, name=seq_record.name,
                               description=seq_record.description, dbxrefs=seq_record.dbxrefs,
                               features=seq_record.features, annotations=seq_record.annotations)

    new_seq_record.letter_annotations["phred_quality"] = phred_quality
    return new_seq_record


def create_new_record_list(record_list):

    new_seq_records = []
    for record in record_list:

        # create new sequence record
        seq_record = create_new_sequence_record(record)

        # add to new record array
        if seq_record is not None:
            new_seq_records.append(seq_record)

    return new_seq_records

# End of functions


# file name defines
file317_1 = "../data/SRR6800317_1.fastq"
file317_2 = "../data/SRR6800317_2.fastq"
file318_1 = "../data/SRR6800318_1.fastq"
file318_2 = "../data/SRR6800318_2.fastq"

new_file317_1 = "./New_SRR6800317_1.fastq"
new_file317_2 = "./New_SRR6800317_2.fastq"
new_file318_1 = "./New_SRR6800318_1.fastq"
new_file318_2 = "./New_SRR6800318_2.fastq"

# read in the files
records1 = list(SeqIO.parse(file317_1, "fastq"))
records2 = list(SeqIO.parse(file317_2, "fastq"))
records3 = list(SeqIO.parse(file318_1, "fastq"))
records4 = list(SeqIO.parse(file318_2, "fastq"))

# get the size of the sequence files
size1 = len(records1)
size2 = len(records2)
size3 = len(records3)
size4 = len(records4)

# make sure they're the same size
if size1 != size2:
    print("Error: File sizes are different for: SRR6800317 ")

# make sure they're the same size
if size3 != size4:
    print("Error: File sizes are different for: SRR6800318 ")

# new sequence record arrays
# with modified read id, seq strings and phred filter
new_records1 = create_new_record_list(records1)
new_records2 = create_new_record_list(records2)
new_records3 = create_new_record_list(records3)
new_records4 = create_new_record_list(records4)

# Write out 317_1
with open(new_file317_1, "w") as output_handle:
    SeqIO.write(new_records1, output_handle, "fastq")
rec_count = len(new_records1)

print("Execution complete - SRR6800317_1")
print(rec_count, " - Records processed")

# Write out 317_2
with open(new_file317_2, "w") as output_handle:
    SeqIO.write(new_records2, output_handle, "fastq")
rec_count = len(new_records2)

print("Execution complete - SRR6800317_2")
print(rec_count, " - Records processed")

# Write out 318_1
with open(new_file318_1, "w") as output_handle:
    SeqIO.write(new_records3, output_handle, "fastq")
rec_count = len(new_records3)

print("Execution complete - SRR6800318_1")
print(rec_count, " - Records processed")

# Write out 318_2
with open(new_file318_2, "w") as output_handle:
    SeqIO.write(new_records4, output_handle, "fastq")
rec_count = len(new_records4)

print("Execution complete - SRR6800318_2")
print(rec_count, " - Records processed")
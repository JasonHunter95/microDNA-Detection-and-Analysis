import pysam

bam = pysam.AlignmentFile("data/SRR413984.querysorted.bam", "rb")

for read in bam:
    if read.cigartuples and (read.cigartuples[0][0] == 4 or read.cigartuples[-1][0] == 4):
        print("Read name:", read.query_name)
        print("Position:", read.reference_start)
        print("CIGAR:", read.cigarstring)
        print("MAPQ:", read.mapping_quality)
        print("Read sequence:", read.query_sequence)
        print("Flags:", read.flag)
        break
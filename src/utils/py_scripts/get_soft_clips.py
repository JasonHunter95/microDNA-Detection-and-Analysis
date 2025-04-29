import pysam

bam_path = "data/SRR413984.querysorted.bam"
bam = pysam.AlignmentFile(bam_path, "rb")

start_softclip = 0
end_softclip = 0
total_reads = 0

for read in bam:
    total_reads += 1
    if read.cigartuples is None:
        continue
    # soft-clip at start
    if read.cigartuples[0][0] == 4:
        start_softclip += 1
    # soft-clip at end
    if read.cigartuples[-1][0] == 4:
        end_softclip += 1

bam.close()

print(f"Total reads: {total_reads}")
print(f"Soft-clipped at start: {start_softclip}")
print(f"Soft-clipped at end:   {end_softclip}")
print(f"Total soft-clipped:    {start_softclip + end_softclip}")

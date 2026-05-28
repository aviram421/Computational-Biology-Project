from pathlib import Path
from pyfaidx import Fasta

# Before starting this step, we downloaded hg18 genome files from UCSC as explained in the readme file.
# Step 1: Combining all the hg18 genome files to one fasta file of the whole hg18 genome.
chromosome_path = Path("Separate_chromosome_files")
combined_file = Path("Input_for_fixed_windows") / "hg18.fa"
fasta_files = sorted(chromosome_path.glob("*.fa"))
with open(combined_file, "w") as combined:
    for fasta_file in fasta_files:
        with open(fasta_file, "r") as single_fasta:
            combined.write(single_fasta.read())

# Step 2: Create fixed 50 nucleotide windows around the summit of each binding site in BED format and convert to Fasta format.
input_bed = "Output_for_3UTR_bed/par_clip_in_3utr.bed" # The original PAR-CLIP entries from the experiment, but only those that are within the 3'UTR.
genome_fasta = "Input_for_fixed_windows/hg18.fa" # Fasta file of the whole hg18 genome without random, chromosome M and the 2 variants of chromosome 6.

windows_bed = "Output_for_fixed_windows/par_clip_windows_bed.bed" # BED file of the 50nt windows.
windows_fasta = "Output_for_fixed_windows/par_clip_windows_fasta.fa" # Fasta file of the 50nt windows.

half_window = 25  # 25 nt on each side of the summit.

# Create fixed windows around all the summits as a BED6 file.
with open(input_bed, "r") as full, open(windows_bed, "w") as final:
    for line in full:
        cols = line.strip().split("\t")

        chrom = cols[0]
        name = cols[3]
        score = cols[4]
        strand = cols[5]

        summit_start = int(cols[6])
        summit_end = int(cols[7])

        summit_center = (summit_start + summit_end) // 2

        if summit_center - half_window < 0: # Protecting against negative values.
            start = 0
        else:
            start = summit_center - half_window

        end = summit_center + half_window

        final.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

# Convert BED to FASTA while preserving the exact order of windows_bed.
genome = Fasta(genome_fasta)

complement = str.maketrans("ACGTacgt", "TGCAtgca")

with open(windows_bed, "r") as bed_file, open(windows_fasta, "w") as final_fasta:
    for line in bed_file:
        cols = line.strip().split("\t")

        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        name = cols[3]
        strand = cols[5]

        seq = genome[chrom][start:end].seq

        if strand == "-":
            seq = seq.translate(complement)[::-1]

        final_fasta.write(f">{name}\n{seq}\n")

# Sanity check - all sequences should be 50nt long.
all_ok = True

with open(windows_fasta) as f:
    for line in f:
        if line.startswith(">"):
            continue
        if len(line.strip()) != 50:
            all_ok = False
            break

if all_ok:
    print('All sequences are 50nt long. All ok!')
else:
    print('Not all sequences are 50nt long!')
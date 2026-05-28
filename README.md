# Project structure:

## Directories:
- Input_for_3UTR_bed/
- Output_for_3UTR_bed/
- Input_for_fixed_windows/
- Output_for_fixed_windows/
- Separate_chromosome_files/

## Files:
- create_3utr_bed.py
- create_fixed_windows_and_fasta.py
- main.py
- general_prediction.py
- requirements.txt
- README.md


### DATA FILES AND SOURCES


1. PAR-CLIP binding sites in a BED file
File:
Input_for_3UTR_bed/GSM1144507_RC3H1_sites.bed

Description:
This file contains the PAR-CLIP binding sites of Roquin-1 in BED8 format (hg18 assembly).

Source:
GEO accession GSM1144507:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1144507

Associated article:
Murakawa, Y., Hinz, M., Mothes, J. et al. RC3H1 post-transcriptionally regulates A20 mRNA and modulates the activity of the IKK/NF-κB pathway. Nat Commun 6, 7367 (2015). 

Article link: https://doi.org/10.1038/ncomms8367

Note:
The compressed BED file was downloaded from GEO and then extracted using Archive Utility on a Mac (or by extracting a zip file on Windows). Then, we erased the first lines which are not coordinates but only describe the file itself.

--------------------------------------------------

2. 3'UTR annotations in a BED file
File:
Input_for_3UTR_bed/hg18_3UTR annotations.bed

Description:
This file contains hg18 3'UTR exon annotations in BED format.

The annotations were obtained from UCSC Table Browser by selecting:
- genome: Human
- assembly: hg18 (Wasn't chosen in the table browser itself, but in the following link: https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu.).
- track: NCBI RefSeq
- Output format: BED - browser extensible data
- Output filename: hg18_3UTR annotations

Then we pressed on "Get output", in the next window (Output refGene as BED) we chose 3' UTR Exons under the label of "Create one BED record per:" and then we pressed on "get BED".

Note: For the remaining parameters, we chose the default values and didn't change them.

Source link:
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3787459551_wJpAUP3jHgf5HSf4p1GAIMA81srd&db=hg18&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C799%2C529-155%2C812%2C871&hgta_outputType=primaryTable&hgta_outFileName=hg18_3UTR annotations

--------------------------------------------------

3. hg18 genome FASTA file
File (Expected after running the create_fixed_windows_and_fasta.py file):
Input_for_fixed_windows/hg18.fa

-- Note: You won't see this file on Moodle because it weighs more than 500 MB. To get this file, please download the 4 zipped folders named "chr_1to6", "chr_7to13", "chr_14to20", and "chr_21to22_xy" from our Moodle submission and extract them. Then move all files (24 chr files in total) from the 4 extracted folders to the "Separate_chromosome_files" folder within the main project directory. Then, after running the create_fixed_windows_and_fasta.py file, you'll get the "hg18.fa" file automatically.

Description:
This file contains the combined FASTA sequences of the hg18 human genome used to extract 50nt windows.

Source for all the separate chromosome files that we include in our submission:
https://hgdownload.soe.ucsc.edu/goldenPath/hg18/chromosomes/

Notes:
We downloaded all files named chr*.fa.gz (chr1-22, X, Y), excluding:
- random chromosomes
- chromosome M
- the two chromosome 6 variant files

After extraction, the chromosome FASTA files are merged into a single file named hg18.fa by running the script create_fixed_windows_and_fasta.py.


### CODE FILES


1. create_3utr_bed.py

Relevant project part:
Preprocessing step for restricting the PAR-CLIP analysis only to binding sites located in 3'UTR.

Input:
- Input_for_3UTR_bed/GSM1144507_RC3H1_sites.bed
- Input_for_3UTR_bed/hg18_3UTR annotations.bed

Output:
- Output_for_3UTR_bed/par_clip_in_3utr.bed

What this code does:
This script keeps only PAR-CLIP binding sites that overlap 3'UTR.

How it works:
The script loads the PAR-CLIP and 3'UTR BED files using pandas and adds an index column (OriginalOrderIndex) to preserve the original order of the experimental data.
Entries on random chromosomes are removed from both datasets.

The datasets are merged by chromosome, and genomic overlap is determined using coordinate conditions (Start < UTR3_End and End > UTR3_Start).
Duplicate entries are removed to ensure each binding site is counted once, and the original order is restored before saving the filtered BED file.

--------------------------------------------------

2. create_fixed_windows_and_fasta.py

Relevant project part:
Preprocessing step for generating 50nt sequence windows around binding-site summits and converting them to FASTA format.

Input:
- Output_for_3UTR_bed/par_clip_in_3utr.bed
- Input_for_fixed_windows/hg18.fa

Output:
- Output_for_fixed_windows/par_clip_windows_bed.bed
- Output_for_fixed_windows/par_clip_windows_fasta.fa

What this code does:
This script combines all the chromosome files into a combined hg18 genome file, creates fixed 50nt windows around the summit of each PAR-CLIP binding site located in the 3'UTR, and then converts these windows from BED format to FASTA format.

How it works:
The script first merges all chromosome FASTA files into a single genome FASTA file.
It then reads the filtered PAR-CLIP BED file and creates fixed 50nt windows around each binding-site summit.

Next, the genome FASTA file is accessed using the pyfaidx libraryץ
For each window, the corresponding nucleotide sequence is extracted based on chromosome, start, and end positions.

Strand orientation is taken into account: for negative-strand entries, the reverse complement is generated using a translation table and sequence reversal.

Finally, all sequences are written to a FASTA file in the same order as the BED file, and a sanity check verifies that all sequences are exactly 50 nucleotides long.

--------------------------------------------------

3. main.py

Relevant project part:
Main data analysis of enriched motifs, threshold evaluation, background comparison, significance filtering, Bayes-based probability estimation and normalizing the results to 0.2 threshold scale.

Input:
- Output_for_fixed_windows/par_clip_windows_fasta.fa

Output:
This script does not create output files. Its outputs are plots, relevant motifs and printed dictionaries.

Graphical output:
- Plots of the top 5 enriched motifs compared to background for thresholds 0.1, 0.2, and 0.3.
- Plot of motif mean enrichment versus threshold (for thresholds 0.1-0.5).
- Plots of the top 5 enriched motifs in strong versus weak subsets for thresholds 0.1 and 0.2.
- Plots of binding probability for exact motifs, one-mismatch motifs, and their background for thresholds 0.1 and 0.2.
- Plots of binding probability for the new degenerate motif pattern for thresholds 0.1 and 0.2.
- Plot of binding probability for all relevant motifs for the threshold = 0.2 subset.

Printed output:
- Significant exact and degenerative motifs.
- Dictionary of significant motifs for threshold = 0.1.
- Dictionary of significant motifs for threshold = 0.2.
- Final dictionary of all relevant motifs for the threshold = 0.2 scale.

What this code does:
This script performs the main motif analysis on the 50nt binding-site windows.

How it works:
First, the script converts DNA sequences from the FASTA file into RNA sequences by replacing T with U. It then evaluates motif mean enrichment for several thresholds (0.1–0.5) and compares the top enriched motifs to a shuffled background in order to determine the relevant thresholds for further analysis. Enrichment vs. background plots are generated for thresholds 0.1, 0.2 and 0.3.

Based on this initial evaluation, the analysis proceeds with thresholds 0.1 and 0.2. For each threshold, the script identifies the top 5 most enriched 7-mer motifs in the strong-binding subset compared to the weak subset.

Next, one-mismatch degenerate motifs and additional degenerate motif patterns are generated and analyzed. A Bayes probability model is applied to estimate the probability that a sequence belongs to the strong-binding subset given that it contains a specific motif, and motifs that exceed the background are considered significant.

Finally, the relevant motifs from thresholds 0.1 and 0.2 are combined and normalized to the threshold = 0.2 subset, and the final set of motifs is presented in a summary plot and dictionary.

--------------------------------------------------

4. general_prediction.py

Relevant project part:
Final prediction tool for testing a new 50nt sequence against the relevant motifs found in the main analysis.

Input:
There are no input files.
The script receives an input sequence from the user interactively.

Output:
There are no output files.
The script prints all the exact and degenerate motifs found in the user sequence, together with their strong-binding probabilities for the threshold = 0.2 subset.

What this code does:
This script serves as a simple prediction tool for new sequences.

How it works:
The script contains the final probability dictionary obtained from the main analysis for the threshold = 0.2 subset. The user is asked to enter a sequence of exactly 50 nucleotides. The script first validates that the sequence length is correct and that it contains only valid DNA/RNA letters (A, U, G, C, or T). It then converts T to U and scans the sequence for all exact motifs and degenerate motifs stored in the final dictionary. Finally, it prints all matching motifs together with their corresponding strong-binding probabilities.

** Important note:
In the main section, we put 4 different examples of sequences that you can run in the interactive platform (as the docstring for the main). Of course, you can put any other sequence as long as it is 50nt long. We also put them here for your convenience:

1. AUGCAUUUAUAGCUACGUGACACAUCUUACGGAUUACGGAUUACGGAAUC
2. ACAAUACGGAUCUAAUACAGCUUACUCGAUUACGGAUACAACUUAGGAUC
3. UAAUACAGGAACUAUUUUACUCGGAACAAUACGGUACAUUACGGAACGUA
4. GGAUACAAUACCUUUACUCAUACUAAGGAACAACUUGGAUUACGAUACAA


### PYTHON ENVIRONMENT AND DEPENDENCIES


This project was developed inside a Python virtual environment.

All required Python packages and their exact versions are listed in the file:

requirements.txt

To install all dependencies, run the following command from the project
directory (the directory containing the requirements.txt file):

pip install -r requirements.txt

This command will install all the packages required to run the project.


### HOW TO RUN THE PROJECT


Run the scripts in the following order:

1. create_3utr_bed.py
2. create_fixed_windows_and_fasta.py
3. main.py
4. general_prediction.py


#### ADDITIONAL NOTES


- Binding sites located on unplaced or random chromosomes were excluded from the analysis.
- The sequence order in the FASTA file reflects binding strength ranking from strongest to weakest, as provided by the input data.

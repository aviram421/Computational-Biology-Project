import pandas as pd

# Before starting this step, we downloaded 3'UTR annotations from UCSC
# and the PAR-CLIP BED file from GEO as explained in the readme file.

# Step 1: Keeping only PAR-CLIP binding sites that are associated with 3'UTR regions.
par_clip_bed = "Input_for_3UTR_bed/GSM1144507_RC3H1_sites.bed"
utr3_bed = "Input_for_3UTR_bed/hg18_3UTR annotations.bed"
output_bed = "Output_for_3UTR_bed/par_clip_in_3utr.bed"

# Read the original PAR-CLIP BED file with explicit column names.
par_clip_df = pd.read_csv(par_clip_bed, sep="\t", header=None, names=["Chromosome", "Start", "End", "Name", "Score", "Strand", "SummitStart", "SummitEnd"])

# Preserve the original order from the experiment file.
par_clip_df["OriginalOrderIndex"] = range(len(par_clip_df))

# Remove random chromosomes before further analysis.
par_clip_df = par_clip_df[~par_clip_df["Chromosome"].str.contains("random", na=False)].copy()


# Read the 3'UTR BED file.
# We only need the first 3 BED columns for overlap checking.
utr3_df = pd.read_csv(utr3_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["Chromosome", "UTR3_Start", "UTR3_End"])

# Remove random chromosomes from the annotation file as well.
utr3_df = utr3_df[~utr3_df["Chromosome"].str.contains("random", na=False)].copy()

# Merge by chromosome so that each PAR-CLIP site is compared only to 3'UTR regions on the same chromosome.
merged_df = par_clip_df.merge(utr3_df, on="Chromosome", how="inner")

# Keep only rows where the intervals overlap.
overlap_df = merged_df[(merged_df["Start"] < merged_df["UTR3_End"]) & (merged_df["End"] > merged_df["UTR3_Start"])].copy()

# Keep each PAR-CLIP site only once, even if it overlaps multiple 3'UTRs.
df_utr3 = overlap_df.drop_duplicates(subset=["OriginalOrderIndex"])

# Restore the original order from the experiment file.
df_utr3 = df_utr3.sort_values("OriginalOrderIndex")

# Save in BED format without the temporary order column.
df_utr3[["Chromosome", "Start", "End", "Name", "Score", "Strand", "SummitStart", "SummitEnd"]].to_csv(output_bed, sep="\t", header=False, index=False)

# Sanity check - most binding sites should be in the 3'UTR region.
num_total = len(par_clip_df)
num_3utr = len(df_utr3)

percentage_in_3utr = 100 * (num_3utr / num_total)
print(f"{percentage_in_3utr:.0f}% of the total binding sites are in the 3'UTR.") # The output is 83%.
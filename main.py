import matplotlib.pyplot as plt
import random
import numpy as np
import re

def main():
    """
    Run the full motif analysis pipeline on the 50nt FASTA windows.

    The function converts the input DNA sequences to RNA, evaluates suitable
    thresholds for strong/weak subsets, identifies enriched exact and degenerate motifs,
    compares them to shuffled background, applies a Bayes-based probability model,
    and summarizes the relevant motifs for the threshold=0.2 subset.

    :return: None
    """

    # Step 1: converting T to U and generating a list of RNA sequences.
    fasta_windows = "Output_for_fixed_windows/par_clip_windows_fasta.fa"  # Input file: Fasta file of the 50nt windows.
    rna_windows = rna_conversion(fasta_windows)
    # num_sequences = len(rna_windows)  # 13399

    # print(rna_windows[:3])  # Sanity check: we compared the first elements in the list to the first sequences in the Fasta file.

    # Step 2: Finding the ideal threshold for the strong and weak subsets.
    ideal_threshold(rna_windows)

    # Step 3: finding and presenting the top 5 enriched motifs in the strongest 10% and 20% binding sequences.
    dict_prob_1 = {}
    dict_prob_2 = {}
    label_dict_1 = {}
    label_dict_2 = {}

    for threshold in [0.1, 0.2]:
        strong_sequences, top_motifs_list, enrichment_list, strong_top_counts = top_motifs_generator(threshold, rna_windows)
        shuffled_strong_all = generate_shuffled_sets(strong_sequences, 200)

        plt.figure()
        plt.bar(top_motifs_list, enrichment_list)
        plt.xlabel('Motif', fontsize=14, labelpad=14)
        plt.ylabel('Enrichment (Strong / Weak)', fontsize=14, labelpad=14)
        plt.title(f'Top 5 Enriched Motifs for threshold={threshold:.1f}', fontsize=14, fontweight='bold', pad=14)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()

        one_mismatch_list = one_mismatch_generator(top_motifs_list)
        label_dict_current = motif_label_generator(top_motifs_list, one_mismatch_list)

        # Step 4: generation of random sequences multiple times and comparison of the strong vs shuffled sequences.
        # Generating background allows us to verify that the motifs are not random,
        # and that their enrichment in comparison to the weak binding sequences is not because of a nucleotide composition bias of the genes.
        bg_enrichment_list, bg_std_dev_list, P_ab_background_mean_exact, P_ab_background_std_exact, P_ab_background_mean_regex, P_ab_background_std_regex = \
            background_generator(top_motifs_list, strong_top_counts, rna_windows, shuffled_strong_all, one_mismatch_list)

        # Step 5: Applying Bayes probability model.
        P_ab_exact = bayes_model(strong_sequences, rna_windows, list(top_motifs_list), motif_type='exact')
        P_ab_1mm = bayes_model(strong_sequences, rna_windows, one_mismatch_list, motif_type='regex')

        plt.figure()
        x = np.arange(5)
        width = 0.2
        plt.bar(x - 1.5 * width, P_ab_exact, width=width, label='Exact', color='tab:blue')
        plt.bar(x - 0.5 * width, P_ab_background_mean_exact, yerr=P_ab_background_std_exact, capsize=8, width=width, label='Background (Exact)', color='tab:green')
        plt.bar(x + 0.5 * width, P_ab_1mm, width=width, label='1 mismatch', color='tab:orange')
        plt.bar(x + 1.5 * width, P_ab_background_mean_regex, yerr=P_ab_background_std_regex, capsize=8, width=width, label='Background (Regex)', color='tab:red')
        plt.xticks(x, top_motifs_list, rotation=90)
        plt.xlabel('Motif')
        plt.ylabel(f'P ({int(threshold * 100)}% most strong | motif) [%]')
        plt.title(f'Binding Probability For Exact, Degenerate and Background Motifs for threshold={threshold:.1f}', fontsize=10)
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Step 6: Applying Bayes probability model for the new degenerate patterns.
        if threshold == 0.2:
            new_pattern = "(AC[AU]A[CU][UA][CU])"
        else:
            new_pattern = "([AU][AU]A[UC][AU][AC]A)"

        label_dict_current[new_pattern] = new_pattern

        P_ab_new_regex = bayes_model(strong_sequences, rna_windows, [new_pattern], motif_type='regex')
        P_ab_new_regex_background_mean, P_ab_new_regex_background_std = shuffled_background_bayes_model(rna_windows, shuffled_strong_all, [new_pattern], motif_type='regex')

        plt.figure()
        x = [0, 0.07]
        width = 0.05
        plt.bar(x[0], P_ab_new_regex[0], width=width, color='tab:orange', label='Strong')
        plt.bar(x[1], P_ab_new_regex_background_mean[0], yerr=P_ab_new_regex_background_std[0], capsize=8, width=width, color='tab:red', label='Background')
        plt.xticks(x, ['Strong', 'Background'])
        plt.xlabel('Subset')
        plt.ylabel(f'P ({int(threshold * 100)}% most strong | regex pattern) [%]')
        plt.title(f'Binding Probability For The New {new_pattern} Pattern for threshold={threshold:.1f}', fontsize=10)
        plt.tight_layout()
        plt.show()

        # Step 7: creating a probability dictionary for each threshold.
        one_mismatch_list.append(new_pattern)
        P_ab_1mm.append(P_ab_new_regex[0])
        P_ab_background_mean_regex.append(P_ab_new_regex_background_mean[0])
        P_ab_background_std_regex.append(P_ab_new_regex_background_std[0])

        dict_prob_exact, dict_prob_regex = relevant_motifs(
            top_motifs_list,
            one_mismatch_list,
            P_ab_exact,
            P_ab_1mm,
            P_ab_background_mean_exact,
            P_ab_background_std_exact,
            P_ab_background_mean_regex,
            P_ab_background_std_regex,
        )

        if threshold == 0.1:
            dict_prob_1 = dict_prob_exact | dict_prob_regex
            label_dict_1 = label_dict_current
        else:
            dict_prob_2 = dict_prob_exact | dict_prob_regex
            label_dict_2 = label_dict_current

        print()

    print("Dictionary for threshold=0.1:")
    print(dict_prob_1)
    print("\nDictionary for threshold=0.2:")
    print(dict_prob_2)

    # Step 8: taking all the relevant motifs from both thresholds,
    # recalculating only the threshold=0.1 motifs for threshold=0.2,
    # and keeping the threshold=0.2 motifs as they are.
    dict_prob_all_for_02 = recalculate_relevant_motifs_to_threshold_02(
        dict_prob_1, dict_prob_2, label_dict_1, label_dict_2, rna_windows
    )

    print("\nDictionary for all relevant motifs for the threshold=0.2 subset:")
    print(dict_prob_all_for_02)

def rna_conversion(fasta_file):
    """
     Convert T to U and generate a list of the RNA sequences.
     The function reads the FASTA file line by line, skips header lines starting with '>',
     replaces all T nucleotides with U in each sequence line, and returns a list of the converted RNA sequences.
     (First sequences in the list correspond to the highest binding score, while the last ones correspond to the lowest binding score.)

    :param fasta_file: File in a fasta format.
    :return: List with all the sequences from the fasta file in RNA representation.
    """

    rna_seqs = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                continue
            rna_seqs.append(line.replace("T", "U"))

    return rna_seqs


def find_k_mers(sequences_list, k):
    """
     Find all possible k-mers in all the sequences.
     The function iterates over each sequence and extracts every substring of length k
     using a sliding window, adding each k-mer to a set so that duplicate k-mers are stored only once.

    :param sequences_list: A list of sequences.
    :param k: Defines the size of the k-mers.
    :return: A set of all possible k-mers.
    """

    k_mers = set()
    for seq in sequences_list:
        for i in range(0, len(seq) - k + 1):
            k_mers.add(seq[i:i + k])
    return k_mers


def k_mer_counting(sequences_list, motifs, k):
    """
     Count the fraction of sequences in which each motif appears.
     The function scans each sequence using a sliding window of size k and checks
     whether each k-mer is in the given motif set. Each motif is counted only once
     per sequence, even if it appears multiple times. Finally, counts are normalized by the total number of
     sequences and returned as a dictionary.

    :param sequences_list: A list of sequences.
    :param motifs: A set of potential motifs to search for.
    :param k: Length of the motif.
    :return: A dictionary of the motifs found in the sequences_list as keys and their fractions as values.
    """

    fractions = {}
    n = len(sequences_list)
    for seq in sequences_list:
        k_mers = set()
        for i in range(len(seq) - k + 1):
            k_mer = seq[i:i + k]
            if k_mer in motifs:
                k_mers.add(k_mer)

        for k_mer in k_mers:
            if k_mer in fractions:
                fractions[k_mer] += 1
            else:
                fractions[k_mer] = 1

    for k_mer in fractions:
        fractions[k_mer] = fractions[k_mer] / n

    return fractions


def top_enriched_motifs(strong_counts, weak_counts, n):
    """
     Find the top n enriched motifs based on enrichment ratio.
     The function calculates enrichment for motifs present in both strong and weak subsets
     by dividing their proportional frequency in the strong set by their
     proportional frequency in the weak set. Motifs are then sorted in
     descending order of enrichment, and the top n motifs are returned.

    :param strong_counts: A dictionary of the strong binding motifs and their proportional counts.
    :param weak_counts: A dictionary of the weak binding motifs and their proportional counts.
    :param n: The number of top enriched motifs to return.
    :return: A list of the top n motifs and their enrichment value as tuples (motif, enrichment), sorted by enrichment (highest first).
    """

    enrichment = {}
    for motif in strong_counts:
        if motif not in weak_counts:
            continue

        s = strong_counts[motif]
        w = weak_counts[motif]
        enrichment[motif] = s / w

    motifs_list = list(enrichment.items())
    motifs_list.sort(key=lambda x: x[1], reverse=True)

    top_motifs = motifs_list[:n]
    return top_motifs


def shuffle_sequences(sequences_list):
    """
     Randomly shuffle each sequence in the given list.
     The function converts each sequence into a list of characters, applies
     random shuffling directly to the character list, and then joins the characters back into a string.
     A new list containing the shuffled sequences is returned.

    :param sequences_list: A list of sequences to shuffle.
    :return: A list of shuffled sequences.
    """

    shuffled = []
    for seq in sequences_list:
        chars_seq = list(seq)
        random.shuffle(chars_seq)
        shuffled.append("".join(chars_seq))

    return shuffled


def generate_shuffled_sets(sequences_list, repeats):
    """
     Generate multiple shuffled versions of the same sequence list.
     The function shuffles the given sequence list several times and stores
     all the shuffled sequence lists in one outer list.

    :param sequences_list: A list of sequences to shuffle.
    :param repeats: Number of shuffled repeats to generate.
    :return: A list containing multiple shuffled sequence lists.
    """

    shuffled_all = []

    for _ in range(repeats):
        shuffled_all.append(shuffle_sequences(sequences_list))

    return shuffled_all


def regex_fraction(re_pattern, sequences_list):
    """
     Calculate the fraction of sequences containing a given regular expression pattern.
     The function compiles the provided regex pattern and scans each sequence.
     If the pattern is found at least once in a sequence, it is counted.
     The final count is normalized by the total number of sequences.

    :param re_pattern: A regular expression pattern.
    :param sequences_list: A list of sequences.
    :return: The fraction of sequences in which the pattern appears.
    """

    compiled_re = re.compile(re_pattern)
    counts = 0
    for seq in sequences_list:
        if compiled_re.search(seq):
            counts += 1

    return counts / len(sequences_list)


def one_mismatch_generator(top_motifs_list):
    """
     Generate regex patterns that represent one mismatch from each exact motif.
     The function creates all possible one-mismatch options for every position in each motif,
     joins them into one regex pattern, and returns the full list of regex motifs.

    :param top_motifs_list: A list of exact motifs.
    :return: A list of regex patterns with one mismatch for each exact motif.
    """

    nucleotides = "AUGC"
    one_mismatch_list = []

    for motif in top_motifs_list:
        one_mismatch_patterns = []
        for i in range(len(motif)):
            mismatch = nucleotides.replace(motif[i], "")
            one_mismatch_patterns.append(motif[:i] + f"[{mismatch}]" + motif[i + 1:])

        one_mismatch_regex = "(" + "|".join(one_mismatch_patterns) + ")"
        one_mismatch_list.append(one_mismatch_regex)

    return one_mismatch_list


def motif_label_generator(top_motifs_list, one_mismatch_list):
    """
     Create a dictionary of labels for the exact motifs and their matching regex motifs.

    :param top_motifs_list: A list of exact motifs.
    :param one_mismatch_list: A list of one-mismatch regex motifs.
    :return: A dictionary of labels for plotting and presentation.
    """

    label_dict = {}

    for i in range(len(top_motifs_list)):
        label_dict[top_motifs_list[i]] = top_motifs_list[i]
        label_dict[one_mismatch_list[i]] = f"regex of {top_motifs_list[i]} motif"

    return label_dict


def motif_type_detector(motif):
    """
     Detect whether a motif is an exact motif or a regex motif.

    :param motif: A motif string.
    :return: 'regex' if the motif is a regex pattern, otherwise 'exact'.
    """

    if "[" in motif or "(" in motif:
        return 'regex'
    return 'exact'


def single_bayes_probability(motif, sequences_list, rna_windows, motif_type='exact', P_b=None):
    """
     Calculate the Bayes probability for one motif.
     The function calculates motif fractions in the subset and in the entire dataset,
     and applies the Bayes equation P(A|B) = (P(B|A) * P(A)) / P(B).
     The motif can be an exact 7-mer or a degenerate motif represented as a regex pattern.
     P_b can be pre-computed and passed in to avoid redundant calculations.

    :param motif: An exact motif or a regex pattern.
    :param sequences_list: A subset of sequences (e.g., strong binding sequences).
    :param rna_windows: A list of all sequences in the dataset.
    :param motif_type: Defines whether the motif is 'exact' or 'regex'.
    :param P_b: Pre-computed fraction of the motif in rna_windows (optional).
    :return: The probability (in %) that a sequence belongs to the subset given that it contains the motif.
    """

    P_a = len(sequences_list) / len(rna_windows)

    if motif_type == 'exact':
        fractions_subset = k_mer_counting(sequences_list, {motif}, 7)
        P_ba = fractions_subset.get(motif, 0)
        if P_b is None:
            fractions_all = k_mer_counting(rna_windows, {motif}, 7)
            P_b = fractions_all.get(motif, 0)

    else:  # motif_type == 'regex'
        P_ba = regex_fraction(motif, sequences_list)
        if P_b is None:
            P_b = regex_fraction(motif, rna_windows)

    if P_b == 0:
        return 0

    return ((P_ba * P_a) / P_b) * 100


def shuffled_background_bayes_model(rna_windows, shuffled_strong_all, motifs_list, motif_type='exact'):
    """
     Calculate the background Bayes probability for a list of motifs of one type.
     The function pre-computes P_b for all motifs once before the loop,
     then goes over all the pre-generated shuffled sequence lists,
     calculates the Bayes probability for all the motifs in each shuffle,
     and finally returns the mean and the standard deviation for each motif.

    :param rna_windows: A list of all sequences in the dataset.
    :param shuffled_strong_all: A list containing multiple shuffled sequence lists.
    :param motifs_list: A list of motifs of one type (exact or regex).
    :param motif_type: Defines whether the motifs are 'exact' or 'regex'.
    :return: Two lists of the background mean and the background std for each motif.
    """

    # Calculate P_b once for all motifs before the loop.
    if motif_type == 'exact':
        fractions_all = k_mer_counting(rna_windows, set(motifs_list), 7)
        P_b_list = [fractions_all.get(motif, 0) for motif in motifs_list]
    else:
        P_b_list = [regex_fraction(motif, rna_windows) for motif in motifs_list]

    P_ab_background_all = []

    for shuffled_strong in shuffled_strong_all:
        current_list = [
            single_bayes_probability(motif, shuffled_strong, rna_windows, motif_type, P_b=P_b_list[i])
            for i, motif in enumerate(motifs_list)
        ]
        P_ab_background_all.append(current_list)

    P_ab_background_mean = list(np.mean(P_ab_background_all, axis=0))
    P_ab_background_std = list(np.std(P_ab_background_all, axis=0))

    return P_ab_background_mean, P_ab_background_std


def mixed_background_bayes_model(rna_windows, shuffled_strong_all, motifs_list):
    """
     Calculate the background probability for a mixed list of exact and regex motifs.
     The function pre-computes P_b for all motifs once before the loop,
     then goes over all the pre-generated shuffled sequence lists,
     calculates the Bayes probability for each motif according to its type,
     and finally returns the mean and the standard deviation for each motif.

    :param rna_windows: A list of all sequences in the dataset.
    :param shuffled_strong_all: A list containing multiple shuffled sequence lists.
    :param motifs_list: A mixed list of exact motifs and regex motifs.
    :return: Two lists of the background mean and the background std for each motif.
    """

    # Calculate P_b once for all motifs before the loop.
    P_b_list = []
    for motif in motifs_list:
        if motif_type_detector(motif) == 'exact':
            fractions_all = k_mer_counting(rna_windows, {motif}, 7)
            P_b_list.append(fractions_all.get(motif, 0))
        else:
            P_b_list.append(regex_fraction(motif, rna_windows))

    P_ab_background_all = []

    for shuffled_strong in shuffled_strong_all:
        current_list = [
            single_bayes_probability(motif, shuffled_strong, rna_windows, motif_type_detector(motif), P_b=P_b_list[i])
            for i, motif in enumerate(motifs_list)
        ]
        P_ab_background_all.append(current_list)

    P_ab_background_mean = list(np.mean(P_ab_background_all, axis=0))
    P_ab_background_std = list(np.std(P_ab_background_all, axis=0))

    return P_ab_background_mean, P_ab_background_std


def background_generator(top_motifs_list, strong_top_counts, rna_windows, shuffled_strong_all, regex_motifs_list=None):
    """
     Generate background motif fractions by using pre-generated shuffled strong sequences.
     The function counts exact motif fractions in each shuffled repeat, and calculates the mean and
     standard deviation of the motif fractions across all the shuffled repeats.
     Enrichment vs. background is computed as strong_fraction / mean_background_fraction.
     If regex motifs are also given, the function additionally returns background Bayes mean/std
     for the exact and regex motifs.

    :param top_motifs_list: A list of exact motifs to evaluate.
    :param strong_top_counts: A dictionary of motif fractions in the strong subset.
    :param rna_windows: A list of all sequences in the dataset.
    :param shuffled_strong_all: A list containing multiple shuffled sequence lists.
    :param regex_motifs_list: A list of regex motifs to evaluate.
    :return: Two lists containing motif enrichment vs. background and the corresponding standard deviations,
             and Bayes background mean/std for exact and regex motifs.
    """

    bg_values = {motif: [] for motif in top_motifs_list}

    for shuffled_strong in shuffled_strong_all:
        fractions_subset_exact = k_mer_counting(shuffled_strong, set(top_motifs_list), 7)

        for motif in top_motifs_list:
            bg_values[motif].append(fractions_subset_exact.get(motif, 0))

    bg_mean = {motif: np.mean(bg_values[motif]) for motif in top_motifs_list}  # Mean value of all the shuffled repeats for each top motif.
    bg_std_dev = {motif: np.std(bg_values[motif]) for motif in top_motifs_list}  # Standard deviation value of all the shuffled repeats for each top motif.

    bg_enrichment = {}
    for motif in top_motifs_list:
        if bg_mean[motif] > 0:
            bg_enrichment[motif] = strong_top_counts[motif] / bg_mean[motif]  # Enrichment vs. Background for each top motif.

    bg_enrichment_list = list(bg_enrichment.values())
    bg_std_dev_list = [bg_std_dev[motif] for motif in bg_enrichment.keys()]

    P_ab_background_mean_exact = []
    P_ab_background_std_exact = []
    P_ab_background_mean_regex = []
    P_ab_background_std_regex = []

    if regex_motifs_list is not None:
        P_ab_background_mean_exact, P_ab_background_std_exact = shuffled_background_bayes_model(
            rna_windows, shuffled_strong_all, top_motifs_list, motif_type='exact'
        )

        P_ab_background_mean_regex, P_ab_background_std_regex = shuffled_background_bayes_model(
            rna_windows, shuffled_strong_all, regex_motifs_list, motif_type='regex'
        )

    return bg_enrichment_list, bg_std_dev_list, P_ab_background_mean_exact, P_ab_background_std_exact, P_ab_background_mean_regex, P_ab_background_std_regex


def top_motifs_generator(threshold, rna_windows):
    """
     Generate the top enriched motifs for a given threshold.
     The function splits the sequences into strong and weak subsets according
     to the threshold, extracts all possible k-mers from the strong subset,
     counts motif fractions in both subsets, and calculates enrichment
     (strong_fraction / weak_fraction) to identify the top motifs.

    :param threshold: Fraction defining the size of the strong and weak subsets.
    :param rna_windows: A list of RNA sequences ordered by binding score.
    :return: The strong subset, the top motifs,
             their enrichment values, and the motif fractions in the strong subset.
    """

    num_sequences = len(rna_windows)
    cut_threshold = int(threshold * num_sequences)
    strong_sequences = rna_windows[:cut_threshold]
    weak_sequences = rna_windows[-cut_threshold:]
    potential_k_mers = find_k_mers(strong_sequences, 7)
    strong_counts = k_mer_counting(strong_sequences, potential_k_mers, 7)
    weak_counts = k_mer_counting(weak_sequences, potential_k_mers, 7)

    top_motifs = top_enriched_motifs(strong_counts, weak_counts, 5)
    top_motifs_list, enrichment_list = zip(*top_motifs)
    return strong_sequences, top_motifs_list, enrichment_list, strong_counts


def ideal_threshold(rna_windows):
    """
     Evaluate motif enrichment for different thresholds.
     The function tests several thresholds, calculates motif enrichment
     vs. shuffled background for each threshold, and plots the enrichment
     of the top motifs as well as the mean enrichment across thresholds.

    :param rna_windows: A list of RNA sequences ordered by binding score.
    :return: None (plots are generated).
    """

    threshold_list = []
    mean_enrichment = []
    mean_std_dev_list = []

    for i in range(1, 6):
        threshold = 0.1 * i
        strong_sequences, top_motifs_list, _, strong_top_counts = top_motifs_generator(threshold, rna_windows)
        shuffled_strong_all = generate_shuffled_sets(strong_sequences, 200)

        bg_enrichment_list, bg_std_dev_list, _, _, _, _ = background_generator(
            top_motifs_list, strong_top_counts, rna_windows, shuffled_strong_all
        )

        threshold_list.append(threshold)
        mean_enrichment.append(np.mean(bg_enrichment_list))
        mean_std_dev_list.append(np.mean(bg_std_dev_list))

        if mean_enrichment[i - 1] >= 1:
            plt.figure()
            plt.bar(top_motifs_list, bg_enrichment_list, yerr=bg_std_dev_list, capsize=8)
            plt.axhline(y=1, linestyle='--', color='gray')
            plt.xlabel('Motif', fontsize=14, labelpad=14)
            plt.ylabel('Enrichment vs. Background', fontsize=14, labelpad=14)
            plt.title(f'Motif Enrichment Compared to Background for threshold={threshold:.1f}', fontsize=10, fontweight='bold', pad=14)
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.show()

    plt.figure()
    plt.bar(threshold_list, mean_enrichment, yerr=mean_std_dev_list, capsize=8, width=0.07)
    plt.axhline(y=1, linestyle='--', color='gray')
    plt.xlabel('Threshold', fontsize=14, labelpad=14)
    plt.ylabel('Mean Enrichment vs. Background', fontsize=14, labelpad=14)
    plt.title('Motif Mean Enrichment vs. Threshold', fontsize=14, fontweight='bold', pad=14)
    plt.xticks(threshold_list, [f"{t:.1f}" for t in threshold_list])
    plt.tight_layout()
    plt.show()


def bayes_model(sequences_list, rna_windows, motifs_list, motif_type='exact'):
    """
     Calculate the probability that a sequence belongs to a given subset if it contains specific motifs using Bayes theorem.
     The function goes over each motif in the given list, calculates its Bayes probability separately,
     and returns all the probabilities in the same order as the input motifs list.
     Motifs can be exact 7-mers or degenerate motifs represented as regex patterns.

    :param sequences_list: A subset of sequences (e.g., strong binding sequences).
    :param rna_windows: A list of all sequences in the dataset.
    :param motifs_list: A list of motifs (exact motifs or regex patterns).
    :param motif_type: Defines whether the motifs are 'exact' or 'regex'.
    :return: A list of probabilities (in %) that a sequence belongs to the subset given that it contains each motif.
    """

    P_ab_list = []

    for motif in motifs_list:
        P_ab_list.append(single_bayes_probability(motif, sequences_list, rna_windows, motif_type))

    return P_ab_list


def relevant_motifs(motif_list_exact, motif_list_regex, P_ab_exact, P_ab_regex, P_ab_bg_mean_exact, P_ab_bg_std_exact, P_ab_bg_mean_regex, P_ab_bg_std_regex):
    """
     Evaluate whether exact and regex motifs are significant relative to their matching background.

    :param motif_list_exact: A list of the exact motifs.
    :param motif_list_regex: A list of the regex motifs.
    :param P_ab_exact: A list of the Bayes probabilities for the exact motifs.
    :param P_ab_regex: A list of the Bayes probabilities for the regex motifs.
    :param P_ab_bg_mean_exact: A list of the background mean values for the exact motifs.
    :param P_ab_bg_std_exact: A list of the background std values for the exact motifs.
    :param P_ab_bg_mean_regex: A list of the background mean values for the regex motifs.
    :param P_ab_bg_std_regex: A list of the background std values for the regex motifs.
    :return: Two dictionaries of the significant exact motifs and regex motifs.
    """

    dict_prob_exact = {}
    dict_prob_regex = {}

    # Exact motifs
    for i in range(len(motif_list_exact)):
        bg_height_exact = P_ab_bg_mean_exact[i] + P_ab_bg_std_exact[i]

        if P_ab_exact[i] > bg_height_exact:
            print(f"The exact motif {motif_list_exact[i]} is significant.")
            dict_prob_exact[motif_list_exact[i]] = P_ab_exact[i]
        else:
            print(f"The exact motif {motif_list_exact[i]} is not significant.")

    # Regex motifs, including the new pattern
    for i in range(len(motif_list_regex)):
        bg_height_regex = P_ab_bg_mean_regex[i] + P_ab_bg_std_regex[i]

        if P_ab_regex[i] > bg_height_regex:
            print(f"The regex pattern {motif_list_regex[i]} is significant.")
            dict_prob_regex[motif_list_regex[i]] = P_ab_regex[i]
        else:
            print(f"The regex pattern {motif_list_regex[i]} is not significant.")

    return dict_prob_exact, dict_prob_regex


def recalculate_relevant_motifs_to_threshold_02(dict_prob_1, dict_prob_2, label_dict_1, label_dict_2, rna_windows):
    """
     Recalculate only the relevant motifs from threshold=0.1 according to threshold=0.2.
     The function takes all the motifs that were found significant for threshold=0.1,
     calculates their Bayes probability again for the strongest 20% subset,
     and then adds the motifs from threshold=0.2 with their original probabilities.
     Finally, it compares all of them to shuffled background and presents the result in a plot.

    :param dict_prob_1: Dictionary of the significant motifs found for threshold=0.1.
    :param dict_prob_2: Dictionary of the significant motifs found for threshold=0.2.
    :param label_dict_1: Dictionary of labels for threshold=0.1 motifs.
    :param label_dict_2: Dictionary of labels for threshold=0.2 motifs.
    :param rna_windows: A list of all sequences in the dataset.
    :return: A dictionary of all the relevant motifs for the threshold=0.2 subset.
    """

    strong_sequences_02, _, _, _ = top_motifs_generator(0.2, rna_windows)
    shuffled_strong_all_02 = generate_shuffled_sets(strong_sequences_02, 200)

    motifs_list = []
    for motif in dict_prob_1:
        motifs_list.append(motif)
    for motif in dict_prob_2:
        if motif not in motifs_list:
            motifs_list.append(motif)

    dict_prob_all_for_02 = {}

    for motif in dict_prob_1:
        current_type = motif_type_detector(motif)
        dict_prob_all_for_02[motif] = round(single_bayes_probability(motif, strong_sequences_02, rna_windows, current_type), 2)

    for motif in dict_prob_2:
        dict_prob_all_for_02[motif] = round(dict_prob_2[motif], 2)

    P_ab_list = []
    for motif in motifs_list:
        P_ab_list.append(dict_prob_all_for_02[motif])

    P_ab_background_mean, P_ab_background_std = mixed_background_bayes_model(rna_windows, shuffled_strong_all_02, motifs_list)

    # Creating short labels for the x-axis.
    short_labels = []
    for motif in motifs_list:
        if motif in label_dict_2:
            short_labels.append(label_dict_2[motif])
        elif motif in label_dict_1:
            short_labels.append(label_dict_1[motif])
        else:
            short_labels.append(motif)

    plt.figure(figsize=(max(10, len(motifs_list) * 1.1), 8))
    x = np.arange(len(motifs_list))
    width = 0.35
    plt.bar(x - width / 2, P_ab_list, width=width, label='Strong', color='tab:blue')
    plt.bar(x + width / 2, P_ab_background_mean, yerr=P_ab_background_std, capsize=8, width=width, label='Background', color='tab:red')
    plt.xticks(x, short_labels, rotation=45, ha='right')
    plt.xlabel('Relevant motifs')
    plt.ylabel('P (20% most strong | motif) [%]')
    plt.title('Binding Probability For All Relevant Motifs for the threshold=0.2 subset', fontsize=10)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return dict_prob_all_for_02


if __name__ == '__main__':
    main()
import re

# Final probability dictionary for the threshold=0.2 subset.
dict_prob_all_for_02 = {
    'UAAUACA': 40.85,
    'AUACUAA': 39.19,
    '([AUG]UAUUUA|C[AGC]AUUUA|CU[UGC]UUUA|CUA[AGC]UUA|CUAU[AGC]UA|CUAUU[AGC]A|CUAUUU[UGC])': 28.06,
    '([AUG]UUUAUA|C[AGC]UUAUA|CU[AGC]UAUA|CUU[AGC]AUA|CUUU[UGC]UA|CUUUA[AGC]A|CUUUAU[UGC])': 25.34,
    '([AU][AU]A[UC][AU][AC]A)': 28.31,
    'ACAACUU': 40.35,
    'ACUAUUC': 47.22,
    'ACAAUAC': 64.0,
    'UUUACUC': 32.88,
    '([UGC]CAACUU|A[AUG]AACUU|AC[UGC]ACUU|ACA[UGC]CUU|ACAA[AUG]UU|ACAAC[AGC]U|ACAACU[AGC])': 26.23,
    '([UGC]CAAUAC|A[AUG]AAUAC|AC[UGC]AUAC|ACA[UGC]UAC|ACAA[AGC]AC|ACAAU[UGC]C|ACAAUA[AUG])': 29.57,
    '([AGC]UUACUC|U[AGC]UACUC|UU[AGC]ACUC|UUU[UGC]CUC|UUUA[AUG]UC|UUUAC[AGC]C|UUUACU[AUG])': 24.51,
    '(AC[AU]A[CU][UA][CU])': 34.85
}

def find_all_matching_motifs(sequence, dict_probabilities):
    """
    Find all relevant motifs that appear in a given sequence.
    The function scans the sequence for all the exact and degenerate motifs
    in the final probability dictionary, and returns a dictionary containing
    all the motifs that were found and their strong-binding probabilities.

    :param sequence: A given sequence of 50nt.
    :param dict_probabilities: A dictionary of the probabilities of each exact or degenerate motif on threshold=0.2 subset.
    :return: A dictionary of all the motifs found in the sequence and their binding probabilities.
    """

    sequence = sequence.upper().replace("T", "U")
    found = {}

    for motif in dict_probabilities:
        if "[" in motif or "(" in motif:
            if re.search(motif, sequence):
                found[motif] = dict_probabilities[motif]
        else:
            if motif in sequence:
                found[motif] = dict_probabilities[motif]

    return found


def sequence_validation(sequence):
    """
    Validate that the given sequence is exactly 50 nucleotides long
    and contains only valid RNA/DNA characters.

    :param sequence: Input sequence from the user.
    :return: True if the sequence is valid, otherwise False.
    """

    sequence = sequence.upper()

    if len(sequence) != 50:
        print("Invalid input: the sequence must contain exactly 50 nucleotides.")
        return False

    for nucleotide in sequence:
        if nucleotide not in "AUGCT":
            print("Invalid input: the sequence must contain only A, U, G, C or T.")
            return False

    return True


def print_prediction_results(found_motifs):
    """
    Print all the motifs found in the sequence and their strong-binding probabilities.

    :param found_motifs: A dictionary of motifs found in the sequence and their probabilities.
    :return: None
    """

    if len(found_motifs) == 0:
        print("No significant strong-binding motif was found in the sequence.")
        return

    print("\nMotifs found in the sequence and their strong-binding probabilities")
    print("(for threshold = 0.2 subset, i.e. the top 20% most strong-binding sequences):\n")

    for motif, probability in found_motifs.items():
        print(f"{motif} -> {probability:.2f}%")


def main():
    """
    Interactive prediction system for a new sequence.

    Example sequences of 50 nt that can be copied as input:
    1. AUGCAUUUAUAGCUACGUGACACAUCUUACGGAUUACGGAUUACGGAAUC
    2. ACAAUACGGAUCUAAUACAGCUUACUCGAUUACGGAUACAACUUAGGAUC
    3. UAAUACAGGAACUAUUUUACUCGGAACAAUACGGUACAUUACGGAACGUA
    4. GGAUACAAUACCUUUACUCAUACUAAGGAACAACUUGGAUUACGAUACAA
    """

    user_sequence = input("Please enter a sequence of exactly 50 nucleotides:\n").strip()

    while not sequence_validation(user_sequence):
        user_sequence = input("Please enter a valid sequence of exactly 50 nucleotides:\n").strip() # strip() removes extra spaces, tabs, or newline characters from the user input.

    found_motifs = find_all_matching_motifs(user_sequence, dict_prob_all_for_02)
    print_prediction_results(found_motifs)


if __name__ == '__main__':
    main()
from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs, print_alignment


if __name__ == '__main__':
    # Read sequences from file
    list_of_pairs= read_fasta_file_as_list_of_pairs(r'OxbenchClustalW\1.fa')
    sequences_id = [pair[0] for pair in list_of_pairs]
    aligned_sequences = [pair[1] for pair in list_of_pairs]
    

    msa = MSA(aligned_sequences, sequences_id)
    print_alignment(msa)

    # Percentage of non-gaps and totally conserved columns
    non_gaps = PercentageOfNonGaps(msa)
    totally_conserved_columns = PercentageOfTotallyConservedColumns(msa)

    percentage = non_gaps.compute()
    print("Percentage of non-gaps: {0} %".format(percentage))

    conserved = totally_conserved_columns.compute()
    print("Percentage of totally conserved columns: {0}".format(conserved))

    # Entropy
    value = Entropy(msa).compute()
    print("Entropy score: {0}".format(value))

    # Sum of pairs
    value = SumOfPairs(msa, Blosum62()).compute()
    print("Sum of Pairs score (Blosum62): {0}".format(value))

    value = SumOfPairs(msa, PAM250()).compute()
    print("Sum of Pairs score (PAM250): {0}".format(value))

    # Star
    value = Star(msa, Blosum62()).compute()
    print("Star score (Blosum62): {0}".format(value))

    value = Star(msa, PAM250()).compute()
    print("Star score (PAM250): {0}".format(value))
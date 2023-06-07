from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs
from pprint import pprint


if __name__ == '__main__':
    # Read sequences from file
    list_of_pairs= read_fasta_file_as_list_of_pairs(r'OxbenchClustalW\1.fa')
    sequences_id , aligned_sequences =  zip(*list_of_pairs)
    
    msa = MSA(aligned_sequences, sequences_id)
    
    stats = {
        'Entropy': Entropy(msa).compute(),
        'PercentageOfNonGaps': PercentageOfNonGaps(msa).compute(),
        'PercentageOfTotallyConservedColumns': PercentageOfTotallyConservedColumns(msa).compute(),
        'SumOfPairsBlosum62': SumOfPairs(msa, Blosum62()).compute(),
        'SumOfPairsPAM250': SumOfPairs(msa, PAM250()).compute(),
        'StarBlosum62': Star(msa, Blosum62()).compute(),
        'StarPAM250': Star(msa, PAM250()).compute()
    }
    pprint(stats)
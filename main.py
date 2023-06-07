from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs
from pathlib import Path
import pandas as pd
from tqdm import tqdm

def eval_msa(fasta_file_path:str) -> dict:
    list_of_pairs= read_fasta_file_as_list_of_pairs(fasta_file_path)
    sequences_id , aligned_sequences =  zip(*list_of_pairs)    

    msa = MSA(aligned_sequences, sequences_id)

    return {
        'Entropy': Entropy(msa).compute(),
        'PercentageOfNonGaps': PercentageOfNonGaps(msa).compute(),
        'PercentageOfTotallyConservedColumns': PercentageOfTotallyConservedColumns(msa).compute(),
        'SumOfPairsBlosum62': SumOfPairs(msa, Blosum62()).compute(),
        'SumOfPairsPAM250': SumOfPairs(msa, PAM250()).compute(),
        'StarBlosum62': Star(msa, Blosum62()).compute(),
        'StarPAM250': Star(msa, PAM250()).compute(),
    }
    
def bulk_eval_msa(fasta_files_path:Path, output_path:Path) -> pd.DataFrame:
    # check if input and output paths exist
    if not fasta_files_path.exists():
        FileNotFoundError(f'Input path {fasta_files_path} does not exist')
    if not output_path.exists():
        FileNotFoundError(f'Output path {output_path} does not exist')
    # check if input path contains fasta files
    fasta_files = list(fasta_files_path.glob('*.fa'))
    if not fasta_files:
        FileNotFoundError(f'No fasta files found in {fasta_files_path}')
    frame_data = [
        {'file_id': fasta_file.stem} | eval_msa(str(fasta_file))
        for fasta_file in tqdm(fasta_files)
    ]
    return pd.DataFrame(frame_data)
    
       


if __name__ == '__main__':
    oxbench_dataset_path = Path('OXBench')
    aligned_clustalw_path = Path('OXBenchClustalW')
    aligned_clustalo_path = Path('OXBenchClustalO')
    output_path = Path('OutputStats')
    #statstics_unaligned = bulk_eval_msa(oxbench_dataset_path, output_path)
    statstics_clustalw = bulk_eval_msa(aligned_clustalw_path, output_path)
    statstics_clustalo = bulk_eval_msa(aligned_clustalo_path, output_path)
    #statstics_unaligned.to_csv(output_path / 'unaligned.csv')
    statstics_clustalw.to_csv(output_path / 'clustalw.csv')
    statstics_clustalo.to_csv(output_path / 'clustalo.csv')
    
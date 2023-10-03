import argparse
import os
import subprocess
from Bio import SeqIO
import logging
import pandas as pd
import tempfile
import random

logging.basicConfig(level=logging.INFO)

def process_fasta(fasta_path, max_subseqs=20, class_="", kmer_size=7, temp_dir="", random_mode=False):
    output_files = []
    total_seq = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq)
        subseq_start_indices = range(0, len(seq_str), 2000)
        
        if random_mode:
            subseq_start_indices = random.sample(subseq_start_indices, min(max_subseqs, len(subseq_start_indices)))
        
        subseqs = [seq_str[i:i+2000] for i in subseq_start_indices][:max_subseqs]
        total_seq += len(subseqs)
        for idx, subseq in enumerate(subseqs):
            subseq_file = os.path.join(temp_dir, f"subseq_{idx}.fasta")
            with open(subseq_file, "w") as f:
                f.write(f">subseq_{idx}\n{subseq}")
            run_jellyfish(subseq_file, kmer_size, temp_dir)
            output_files.append(os.path.join(temp_dir, f"{os.path.basename(subseq_file).replace('.fasta', '')}_jf_formatted.csv"))
    logging.info(f"Processed {total_seq} sequences from {fasta_path}")
    return output_files

def run_jellyfish(fasta_path, kmer_size, temp_dir):
    output_prefix = os.path.join(temp_dir, os.path.basename(fasta_path).replace('.fasta', '') + '_jf')
    command = f"jellyfish count -m {kmer_size} -o {output_prefix} -s 10000000 -t 20 {fasta_path}"
    try:
        subprocess.run(command, shell=True, check=True)
        logging.info(f'Jellyfish completed for {fasta_path}')
    except subprocess.CalledProcessError:
        logging.error(f'Jellyfish failed for {fasta_path}')
    dump_command = f"jellyfish dump {output_prefix} > {output_prefix}_dump.fa"
    subprocess.run(dump_command, shell=True)
    reformat_jellyfish_dump(f"{output_prefix}_dump.fa")

def reformat_jellyfish_dump(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    output_file = file_path.replace('_dump.fa', '_formatted.csv')
    with open(output_file, 'w') as file:
        for i in range(0, len(lines), 2):
            count, sequence = lines[i].strip('>\n'), lines[i+1].strip()
            file.write(f'{sequence},{count}\n')

def aggregate_jellyfish_output(fasta_path, output_files, class_):
    aggregated_data = []
    for idx, file in enumerate(output_files):
        df = pd.read_csv(file, header=None, names=['kmer', 'count'])
        unique_id = f"{os.path.splitext(os.path.basename(fasta_path))[0]}_chunk_{idx+1}"
        data = {'unique_id': unique_id, 'file_name': os.path.basename(fasta_path), 'sample_id': idx+1, 'class': class_}
        data.update(df.set_index('kmer')['count'].to_dict())
        aggregated_data.append(data)
    result_df = pd.DataFrame(aggregated_data).fillna(0)
    result_df.to_csv(f"{os.path.splitext(os.path.basename(fasta_path))[0]}_aggregated_jellyfish_output.csv", index=False)

def main(args):
    all_aggregated_data = []
    if args.folder_input:
        fasta_files = [os.path.join(args.folder_input, f) for f in os.listdir(args.folder_input) if f.endswith('.fasta')]
    else:
        fasta_files = [args.fasta_path]
    
    for fasta_file in fasta_files:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_files = process_fasta(fasta_file, args.max, args.class_, args.kmer_size, temp_dir)
            aggregate_jellyfish_output(fasta_file, output_files, args.class_)
            result_file = f"{os.path.splitext(os.path.basename(fasta_file))[0]}_aggregated_jellyfish_output.csv"
            all_aggregated_data.append(pd.read_csv(result_file))

    if args.folder_input and all_aggregated_data:
        combined_df = pd.concat(all_aggregated_data, ignore_index=True)
        combined_df.to_csv(f"{os.path.basename(args.folder_input)}_all_aggregated_jellyfish_output.csv", index=False)


def main(args):
    all_aggregated_data = []
    folder_inputs = args.folder_input.split(',')
    class_labels = args.class_.split(',')

    if len(folder_inputs) != len(class_labels):
        logging.error("Number of folder inputs must match number of class labels")
        return

    for folder_input, class_label in zip(folder_inputs, class_labels):
        fasta_files = [os.path.join(folder_input, f) for f in os.listdir(folder_input) if f.endswith('.fasta')]
        for fasta_file in fasta_files:
            with tempfile.TemporaryDirectory() as temp_dir:
                output_files = process_fasta(fasta_file, args.max, class_label, args.kmer_size, temp_dir, args.random_mode)
                aggregate_jellyfish_output(fasta_file, output_files, class_label)
                result_file = f"{os.path.splitext(os.path.basename(fasta_file))[0]}_aggregated_jellyfish_output.csv"
                all_aggregated_data.append(pd.read_csv(result_file))

    if all_aggregated_data:
        combined_df = pd.concat(all_aggregated_data, ignore_index=True)
        output_file = "all_aggregated_jellyfish_output.csv"
        combined_df.to_csv(output_file, index=False)

        # Calling the R script
        r_script_path = "pcoa_plot.R"  # Adjust the path if necessary
        output_plot_path = "output_plot.pdf"  # Adjust the path if necessary
        subprocess.run(f"Rscript {r_script_path} -f {output_file} -o {output_plot_path} -c {args.color_by}", shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a fasta file with Jellyfish.')
    parser.add_argument('fasta_path', type=str, nargs='?', default='', help='Path to the fasta file')
    parser.add_argument('--max', type=int, default=20, help='Maximum number of 2000 nt subsequences to process (default: 20)')
    parser.add_argument('--class_', type=str, required=True, help='Comma-separated class labels to be added to each row in the output data')
    parser.add_argument('--kmer_size', type=int, default=7, help='Size of k-mers to be used by Jellyfish (default: 7)')
    parser.add_argument('--folder_input', type=str, required=True, help='Comma-separated folders containing fasta files to be processed')
    parser.add_argument('--random_mode', action='store_true', help='Enable random mode to randomly select subsequences')
    parser.add_argument('--color_by', type=str, default='class', help='Column name to color points by in the PCoA plot (default: class)')

    args = parser.parse_args()
    
    main(args)

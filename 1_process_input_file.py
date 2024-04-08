import subprocess
from Bio import SeqIO
import os
import random

def create_file(file_dir):
    try:
        subprocess.call(["touch", file_dir])
        with open(file_dir, 'w', encoding='utf-8') as f:
            pass
    except:
        print("An exception occurred, solving")
        problematic_file = os.path.basename(file_dir)
        print(file_dir)
        subprocess.call(["touch", file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv")])
        print(file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv"))
        if os.path.isfile(file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv")):
            with open(file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv"), 'w', encoding='utf-8') as f:
                pass
            print("Exception file created:"+file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv"))
        subprocess.call(["mv", file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv"), file_dir])
        if os.path.isfile(file_dir):
            print("solved")
        else:
            print("not solved")

def GC(seq):
    seqlen = len(seq)
    if seqlen == 0:
        return 0
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    gc_percent = (gc_count / seqlen) * 100
    return gc_percent

def parse_fasta_to_dict(fasta_file):
    sequences_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description.split('|')
        seq_id = None
        variant_H = 0
        variant_N = 0
        gene = None
        for field in header:
            if field.startswith('EPI_ISL'):
                seq_id = field
            elif len(field) >= 4:
                if 'H' in field and any(c.isdigit() for c in field) and 'N' in field:
                    variant_H = field[field.find('H')+1]
                    variant_N = field[field.find('N')+1]
                    
            elif field in ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS', 'HE', 'P3']:
                gene = field
        key = seq_id
        if key not in sequences_dict:
            sequences_dict[key] = []
        sequences_dict[key].append([gene, variant_H, variant_N, str(record.seq), len(record.seq), GC(record.seq)])
    return sequences_dict

def dict_to_tsv(dictionary, output_tsv):
    if not os.path.exists(output_tsv):
        create_file(output_tsv)

    with open(output_tsv, 'w') as tsv_file:
        for seq_id, gene_sequences in dictionary.items():
            for values in gene_sequences:
                tsv_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq_id, values[0], values[1], values[2], values[3], values[4], values[5]))

    print("Archivo {output_tsv} guardado exitosamente.")

def split_tsv2(input_file, output_dir, lines_per_file):
    print("Splitting: " + input_file)
    if os.path.getsize(input_file) == 0:
        print("Input file is empty")
        return
    tsv_file_name = os.path.basename(input_file).replace(".tsv", "")
    with open(input_file, 'r', encoding='utf-8') as f:
        file_count = 1
        line_count = 0
        output_file = str(output_dir) + str(tsv_file_name) + "_" + str(file_count) + ".tsv"
        create_file(output_file)
        out = open(output_file, 'w', encoding='utf-8')
        
        for line in f:
            out.write(line)
            line_count += 1
            if line_count >= lines_per_file:
                out.close()
                file_count += 1
                line_count = 0
                output_file = str(output_dir) + str(tsv_file_name) + "_" + str(file_count) + ".tsv"
                create_file(output_file)
                out = open(output_file, 'w', encoding='utf-8')
    if not out.closed:
        out.close()

def split_tsv(input_file, output_dir, output_prefix, train_ratio):
    print("Splitting:", input_file)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    random.shuffle(lines)
    
    total_lines = len(lines)
    train_lines = int(total_lines * train_ratio)
    
    train_data = lines[:train_lines]
    test_data = lines[train_lines:]
    
    train_output_file = os.path.join(output_dir, output_prefix + "_train.tsv")
    with open(train_output_file, 'w', encoding='utf-8') as train_file:
        train_file.writelines(train_data)
    print("Train file saved:", train_output_file)
    
    subprocess.call(['rm', input_file])

    test_output_file = os.path.join(output_dir, output_prefix + "_test.tsv")
    with open(test_output_file, 'w', encoding='utf-8') as test_file:
        test_file.writelines(test_data)
    print("Test file saved:", test_output_file)

def update_input(input_dir, tsv_database_dir):
    print("Updating database")
    fasta_files = [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith('.fasta')]
    
    for fasta_file in fasta_files:
        if True:
            sequences_dict = parse_fasta_to_dict(fasta_file)
            fasta_name = os.path.basename(fasta_file)
            output_tsv = os.path.join(tsv_database_dir, fasta_name.replace(".fasta", ".tsv"))

            create_file(output_tsv)
            with open(output_tsv, 'w') as tsv_file:
                for seq_id, gene_sequences in sequences_dict.items():
                    for values in gene_sequences:
                        tsv_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq_id, values[0], values[1], values[2], values[3], values[4], values[5]))
            print("File "+output_tsv+" was created successfully.")
        else:
            print('File '+fasta_file+' was already processed.')
        print("Database updated successfully.")

def update_datasets(tsv_database_dir, dataset_dir, lines_per_file, train_ratio):
    print("Updating datasets")
    tsv_files = [os.path.join(tsv_database_dir, file) for file in os.listdir(tsv_database_dir) if file.endswith('.tsv')]
    for tsv_file in tsv_files:
        tsv_file_base = os.path.splitext(os.path.basename(tsv_file))[0]
        split_tsv(tsv_file, tsv_database_dir, tsv_file_base, train_ratio)

        output_prefix = os.path.join(dataset_dir, tsv_file_base)
    tsv_files = [os.path.join(tsv_database_dir, file) for file in os.listdir(tsv_database_dir) if file.endswith('.tsv')]
    for tsv_file in tsv_files:
        split_tsv2(tsv_file, dataset_dir, lines_per_file)

input_dir = 'Path to the input files directory'
tsv_database_dir = 'Path to the output TSV directory'
dataset_dir = 'Path to the output divided data TSV directory'
lines_per_file = 50
train_ratio = 0.8
update_input(input_dir, tsv_database_dir)
update_datasets(tsv_database_dir, dataset_dir, lines_per_file, train_ratio)
tsv_database_dir = 'Path to the output TSV directory'
dataset_dir = 'Path to the output divided data TSV directory'
lines_per_file = 50
train_ratio = 0.8
update_input(input_dir, tsv_database_dir)
update_datasets(tsv_database_dir, dataset_dir, lines_per_file, train_ratio)

import os

datasets_dir = "input TSV files"
encoded_dir = "output encoded TSV files"
encoding_type = "Encoding type (label, OH or OHE)"
labels = [0.25, 0.5, 0.75, 1, 0] #Just for label encoding (encoding for [A, C, G, T, N])
fill_rows = False
window_size = 100
step_size = 10

def create_file(file_dir):
    try:
        subprocess.call(["touch", file_dir])
        with open(file_dir, 'w', encoding='utf-8') as f:
            pass
    except Exception as e:
        print("An exception occurred:", e)
        problematic_file = os.path.basename(file_dir)
        new_file_dir = file_dir.replace(problematic_file, "EMPTY_FILE_EXCEPTION.tsv")
        subprocess.call(["touch", new_file_dir])
        print("Exception file created:", new_file_dir)
        if os.path.isfile(new_file_dir):
            print("Exception file created successfully.")
        else:
            print("Failed to create exception file.")

def new_gc(seq, encoding_type, labels):
    total_bases = 0
    gc_bases = 0
    if encoding_type == "label":
        for base in seq:
            if base != '0':  
                if base == str(labels[1]) or base == str(labels[2]):
                    gc_bases += 1
                total_bases += 1
    
    elif encoding_type == "OH":
        for base in seq:
            if base != '[0,0,0,0]':  
                if base == '[0,1,0,0]' or base == '[0,0,1,0]':
                    gc_bases += 1
                total_bases += 1
    elif encoding_type == "OHE":
        for base in seq:
            if base != '[0,0,0,0]':  
                if base == '[0,1,0,0]' or base == '[0,0,1,0]':
                    gc_bases += 1
                total_bases += 1
    gc_percent = (gc_bases / total_bases) * 100
    return gc_percent


def encode_sequence(sequence, encoding_type, labels):
    encoded_sequence = ''
    if encoding_type == "label":
        encoding = {'A': str(labels[0]), 'C': str(labels[1]), 'G': str(labels[2]), 'T': str(labels[3])}
    elif encoding_type == "OH":
        encoding = {'A': '[1,0,0,0]', 'C': '[0,1,0,0]', 'G': '[0,0,1,0]', 'T': '[0,0,0,1]'}
    elif encoding_type == "OHE":
        encoding = {'A': '1000', 'C': '0100', 'G': '0010', 'T': '0001'}

    for base in sequence.upper():
        if encoding_type == "label":
            encoded_sequence += encoding.get(base, str(labels[4]))
        if encoding_type == "OH":
            encoded_sequence += encoding.get(base, '[0,0,0,0]')
        if encoding_type == "OHE":
            encoded_sequence += encoding.get(base, '0000')    
    return encoded_sequence


def sliding_window(sequence, window_size, step_size):
    windows = []
    for i in range(0, len(sequence) - window_size + 1, step_size):
        windows.append(sequence[i:i+window_size])
    return windows

def encode_file(dataset_tsv, encoded_dataset_tsv, window_size, fill_rows, step_size, encoding_type, labels):
    data = []
    longest_seq = 0
    with open(dataset_tsv, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            seq = parts[4]
            if len(seq) > window_size:
                subseqs = [seq[i:i+window_size] for i in range(0, len(seq), step_size)]
            else:
                subseqs = [seq.ljust(window_size, 'N')]
            for i in range(len(subseqs)):
                subseqs[i] = encode_sequence(subseqs[i], encoding_type, labels)
            if len(subseqs) > longest_seq:
                longest_seq = len(subseqs)
            if not fill_rows:
                rows = sliding_window(subseqs, window_size, step_size)
                for row in rows:
                    encoded_gc = new_gc(row, encoding_type, labels)
                    data.append([parts[0], parts[2], parts[3], parts[5], encoded_gc] + row)
            else:
                data.append([parts[0], parts[2], parts[3], parts[5], int(parts[6].split('.')[0])] + subseqs)

    if fill_rows:
        n_subseqs = longest_seq - window_size
        n_columns = n_subseqs + 5
        empty_subseq = '0' * window_size
        for i in range(len(data)):
            if len(data[i]) < n_columns:
                data[i] += [empty_subseq] * (n_columns - len(data[i]))

    label_map = {'HA': 0, 'NA': 1, 'other': 2}
    with open(encoded_dataset_tsv, 'w') as output_file:
        header = "ID\tgene\tH\tN\tlen\tgc_percent\t" + "\t".join(f"seq_{i+1}" for i in range(window_size))
        output_file.write(header + "\n")
        for row in data:
            label = row[1] if row[1] in label_map else 'other'
            row[1] = label_map[label]
            output_file.write('\t'.join(map(str, row)) + '\n')

def update_encoded_datasets(datasets_dir, encoded_dir, encoding_type, labels, fill_rows, window_size, step_size):
    tsv_files = [os.path.join(datasets_dir, file) for file in os.listdir(datasets_dir) if file.endswith('.tsv')]
    if not os.path.exists(encoded_dir):
        os.makedirs(encoded_dir)
    for tsv_file in tsv_files:
        dataset_filename = os.path.basename(tsv_file)
        encoded_dataset_tsv_path = os.path.join(encoded_dir, dataset_filename.replace('.tsv', '_encoded.tsv'))
        encode_file(tsv_file, encoded_dataset_tsv_path, window_size, fill_rows, step_size, encoding_type, labels)

update_encoded_datasets(datasets_dir, encoded_dir, encoding_type, labels, fill_rows, window_size, step_size)


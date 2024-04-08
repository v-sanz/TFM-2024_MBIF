import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import load
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout

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

def sliding_window(sequence, window_size, step_size):
    subseqs = []
    for i in range(0, len(sequence) - window_size + 1, step_size):
        subseqs.append(sequence[i:i+window_size])
    return subseqs

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

def preprocess_seq(sequence, window_size, step_size, encoding_type, labels, gc_inclusion):
    encoded_subseqs = []
    preprocessed_seq = []
    subseqs = sliding_window(sequence, window_size, step_size)
    for subseq in subseqs: 
        encoded_subseq = encode_sequence(subseq, encoding_type, labels)
        encoded_subseqs.append(encoded_subseq)
    for row in encoded_subseqs:
        if gc_inclusion:
            gc_percent = new_gc(row, encoding_type, labels)
            row.insert(0, gc_percent)
        if encoding_type == "OHE": 
            new_row = []
            for i in range(len(row)):
                if i == 0:
                    new_row.append(row[i])
                else:
                    for number in row[i]:
                        new_row.append(number)
            row = new_row
        preprocessed_seq.append(row)
    return preprocessed_seq

def predict_subseq(subseq, model_path, model_type, gc_inclusion, encoding_type):
    if os.path.exists(model_path):
        model = load(model_path)
    else:
        print("Model not found.")
        exit()
    if gc_inclusion:
        cols = ['gc_percent']
    else:
        cols = []
    for i in range(1, len(subseq)+1):
        cols.append(str('seq_'+str(i)))

    X = pd.DataFrame([subseq], columns=cols)
    if encoding_type == "OH":
        seq_columns = [col for col in X.columns if col.startswith('seq_')]
        X = X[seq_columns].apply(lambda x: np.array(x)).values
        X = np.array([np.array(x.tolist()) for x in X])
    
    prediction = model.predict(X)
    return prediction

def predict_seq(sequence, window_size, step_size, encoding_type, labels, gc_inclusion, model_path, model_type):
    subseqs = preprocess_seq(sequence, window_size, step_size, encoding_type, labels, gc_inclusion)
    prediction_dict = {}
    for i in range(len(sequence)):
        prediction_dict[i] = [0, 0, 0]  # Initialize with zeros for each class

    for n, subseq in enumerate(subseqs):
        prediction = predict_subseq(subseq, model_path, model_type, gc_inclusion, encoding_type)
        for i in range(n * step_size, n * step_size + window_size):
            prediction_dict[i][prediction] += 1
    
    for pos in prediction_dict:
        total_freq = sum(prediction_dict[pos])
        if total_freq > 0:
            prediction_dict[pos] = [freq / total_freq for freq in prediction_dict[pos]]

    return prediction_dict

def save_probability_distribution_plot(predictions_dict, png_path):
    positions = list(predictions_dict.keys())
    probs = np.array(list(predictions_dict.values()))

    plt.figure(figsize=(10, 6))
    plt.plot(positions, probs[:, 0], label='HA', color='blue')
    plt.plot(positions, probs[:, 1], label='NA', color='red')
    #plt.plot(positions, probs[:, 2], label='Other', color='green')

    plt.xlabel('Sequence position')
    plt.ylabel('Probability')
    plt.title('Probability Distribution of Gene Prediction')
    plt.legend()
    plt.savefig(png_path)
    plt.close()


sequence = "ATGC"  # Input DNA sequence
model_path = "Path to the trained model"
model_type = "CNN"
encoding_type = "OH"
gc_inclusion = True
window_size = 100
step_size = 10
png_path = "Path to save the output PNG"

predictions_dict = predict_seq(sequence, window_size, step_size, encoding_type, gc_inclusion, model_path, model_type)
save_probability_distribution_plot(predictions_dict, png_path)

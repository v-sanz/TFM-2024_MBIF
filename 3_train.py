import os
import pandas as pd
import numpy as np
from joblib import dump, load

train_dir = "Path of the TSV input files dir"
encoding_type = "Encoding type (label, OH or OHE)"
gc_inclusion = True  # True or False
model_path = "Path to the training model"
model_type = "LR"  # "LR" or "CNN"

if model_type == "LR":
    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score
    from joblib import dump, load
    from sklearn.preprocessing import StandardScaler
if model_type == "CNN":
    from sklearn.preprocessing import LabelEncoder
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout

def load_model(model_path, model_type, encoding_type, gc_inclusion):
    if os.path.exists(model_path):
        model = load(model_path)
    elif model_type == "LR":
        model = LogisticRegression(multi_class='ovr')
    elif model_type == "CNN":
        if encoding_type == "label":
            dimension1 = 100
            dimension2 = 1
        elif encoding_type == "OH":
            dimension1 = 100
            dimension2 = 4
        elif encoding_type == "OHE":
            dimension1 = 400
            dimension2 = 1
        if gc_inclusion:
            dimension1 += 1
        model = Sequential([
                Conv1D(filters=32, kernel_size=5, activation='relu', input_shape=(dimension1, dimension2)),
                MaxPooling1D(pool_size=4, strides=1),
                Conv1D(filters=64, kernel_size=5, activation='relu'),
                MaxPooling1D(pool_size=4, strides=1),
                Conv1D(filters=128, kernel_size=5, activation='relu'),
                MaxPooling1D(pool_size=4, strides=1),
                Flatten(),
                Dense(1025, activation='relu'),
                Dropout(0.2),
                Dense(512, activation='relu'),
                Dropout(0.2),
                Dense(128, activation='relu'),
                Dropout(0.2),
                Dense(3, activation='softmax')
            ])
    return model

def load_data(file_path, encoding_type, gc_inclusion, model_type):
    print("Training input:", file_path)
    df = pd.read_csv(file_path, delimiter="\t")
    parts = df.columns.tolist()
    if gc_inclusion:
        X = df[parts[5:]]
    else:
        X = df[parts[6:]]
    if model_type == "LR":
        Y = df['gene']
    if model_type == "CNN":
        encoder = LabelEncoder()
        Y = encoder.fit_transform(df['gene'])
        if encoding_type != "OH":
            X = X.values.reshape(X.shape[0], X.shape[1], 1)
        else:
            seq_columns = parts[6:]
            X = df[seq_columns].map(lambda x: np.array(eval(x))).values
            X = np.array([np.array(x.tolist()) for x in X])  # Converto to NumPy matrix
    return X, Y

def train(train_dir, model_path, model_type, encoding_type, gc_inclusion):
    for tsv_file in os.listdir(train_dir):
        if tsv_file.endswith(".tsv"):
            file_path = os.path.join(train_dir, tsv_file)
            print("Training input:", file_path)
            model = load_model(model_path, model_type, encoding_type, gc_inclusion)
            X, Y = load_data(file_path, encoding_type, gc_inclusion)
            if model_type == "LR":
                model.fit(X, Y)
            elif model_type == "CNN":
                model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
                model.fit(X, Y, epochs=10, batch_size=64)
            print("Training successful")
            dump(model, model_path)

train(train_dir, model_path, model_type, encoding_type, gc_inclusion)

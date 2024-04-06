import os
import pandas as pd
from joblib import load
import numpy as np
from sklearn.preprocessing import LabelEncoder

test_dir = "Path of the TSV input files dir"
encoding_type = "Encoding type (label, OH or OHE)"
gc_inclusion = True  # True or False
model_path = "Path to the training model"
model_type = "LR"  # "LR" or "CNN"
predictions_path = "Path of the output predictions file"
report_path = "Path of the output predictions report"

if model_type == "CNN":
    import keras

def load_data(file_path, encoding_type, gc_inclusion, model_type):
    print("Loading data from:", file_path)
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
            X = np.array([np.array(x.tolist()) for x in X])  # Convertir a matriz NumPy tridimensional
    return X, Y

def predict(model_path, model_type, X):
    model = load(model_path)
    predictions = model.predict(X)
    if model_type == "CNN":
        predictions = np.argmax(predictions, axis=1)
    return predictions

def write_predictions(predictions, predictions_path):
    with open(predictions_path, 'w') as file:
        file.write("Gene\tPrediction\n")
        for gene, pred in predictions:
            file.write(f"{gene}\t{pred}\n")
    print("Predictions saved in:", predictions_path)

def write_report(predictions, report_path):
    total_predictions = len(predictions)
    correct_predictions = sum(1 for gene, pred in predictions if gene == pred)
    precision = (correct_predictions / total_predictions) * 100

    confusion_matrix = {0: [0, 0, 0], 1: [0, 0, 0], 2: [0, 0, 0]}
    for gene, pred in predictions:
        confusion_matrix[gene][pred] += 1

    with open(report_path, 'w') as file:
        file.write("\t0\t1\t2\tTotal\n")
        for i in range(3):
            file.write(f"{i}\t")
            for j in range(3):
                file.write(f"{confusion_matrix[i][j]}\t")
            file.write(f"{sum(confusion_matrix[i])}\n")
        file.write("Total\t")
        for i in range(3):
            confusion_matrix_total = sum(confusion_matrix[j][i] for j in range(3))
            file.write(f"{confusion_matrix_total}\t")
        file.write(f"{total_predictions}\n")
        file.write("Accuracy\t\t\t\t")
        file.write(f"{precision}%\n")
    print("Report saved in:", report_path)

def evaluate(test_dir, encoding_type, gc_inclusion, model_path, model_type, predictions_path, report_path):
    predictions = []
    if not os.path.exists(model_path):
        print("The model is not in the specified path:", model_path)
        return
    
    for file_name in os.listdir(test_dir):
        if file_name.endswith(".tsv"):
            file_path = os.path.join(test_dir, file_name)
            print("Predicting classes for:", file_path)
            data = load_data(file_path, encoding_type, gc_inclusion, model_type)
            file_predictions = predict(model_path, model_type, data[0])
            predictions.extend(zip(data[1], file_predictions))
    
    write_predictions(predictions, predictions_path)
    write_report(predictions, report_path)

evaluate(test_dir, encoding_type, gc_inclusion, model_path, model_type, predictions_path, report_path)

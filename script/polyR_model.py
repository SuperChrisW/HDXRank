## polynomial regression model ##
## HHblits+HDMD+SASA average along sequence as input ##

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
import numpy as np
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os
import pandas as pd
from tensorflow.keras.preprocessing.sequence import pad_sequences

correction = 0
def seq_embedding(HDX_file, vector_file):
  # read HDX data file
    datafile = pd.read_csv(HDX_file)
    #datafile = HDX_file
    # filter peptides > 30 AA
    max_len = 30
    df = datafile[datafile['Sequence'].str.len() < max_len] # sequence column

    start_pos = [x + correction for x in df['Start'].tolist()]
    end_pos = [x + correction for x in df['End'].tolist()]


    seq = df['Sequence'].tolist()
    #print(f"Number of peptides: {len(seq)}")
    # read embedding file
    embedding = pd.read_table(vector_file, header=None)
    embed_array = embedding.loc[:,2:].to_numpy()
    embed_seq = embedding.loc[:,1].to_numpy()

    row_size = len(start_pos)
    nfactor = 36

    input_array = np.empty((row_size,nfactor,max_len)) # create an empty array

    for i, (start, end, sequence) in enumerate(zip(start_pos, end_pos, seq)):
        embed_sequence = ''.join(embed_seq[start - 1:end])
        if sequence != embed_sequence:
            #print(f"Sequence mismatch at index {i}: {sequence} != {embed_sequence}")
            continue
        else:
            #print(f"Sequence match at index {i}: {sequence} == {embed_sequence}")
            seq_array = embed_array[start - 1:end]
            #seq_array_float = [[round(float(item), 3) for item in row] for row in seq_array]
            #if np.isnan(seq_array).any():
            #    continue 
            seq_arrayT = np.transpose(seq_array)
            padded_seq = pad_sequences(seq_arrayT, maxlen=max_len, padding="post",dtype="float64")
            input_array[i,:,:] = padded_seq

    return input_array, df

root_dir = "/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/dataset/filtered_results"
manager_file = f"{root_dir}/manager_finished.txt"
manager = open(manager_file, "r")
input_array = []
truth_array = []
error_id = []


file_num = 28
count = 0
for line in manager:
    line = line.split(':')
    if len(line) > 1:
        csv_file = line[0]
        identifier = line[1].strip()
        prediction_time = line[2].strip()
        Uni_ID = line[3].strip()
        correction_value = line[4].strip()
        #print(csv_file, identifier, prediction_time, Uni_ID, correction_value)
        if correction_value == 'skip':
            continue
    else:
        continue
    example_embedding = f"/Users/liyao/Desktop/Tsuda_Lab/Source_code/AI-HDX-main/HDX_MS_dataset/embedding_files/{Uni_ID}.embedding.txt"
    if os.path.isfile(example_embedding) == False:
        continue
    if 'CUL' in csv_file or 'ELOC' in csv_file:
        continue
    if 'state_19810' in csv_file:
        continue
    elif 'PRAS40' in csv_file:
        continue

    print(csv_file, identifier, prediction_time, Uni_ID, correction_value)
    count += 1
    # import pre-processed sequence embedding vector and the peptide fragment tables. Save them into two variables, 'prot' is the peptide embedding array, 'df1' is the table to store outputs.
    example_input = f"{root_dir}/{csv_file}_filtered.csv"
    prot1, df1 = seq_embedding(example_input, example_embedding)
    input_array.append(prot1)
    truth = df1["%D"].to_list()
    truth_array.extend(truth)
    if count >= file_num:
        break

input_array = np.concatenate(input_array, axis=0)
truth_array = np.array(truth_array)

# Find indices where NaNs are present in input_array and truth_array
nan_indices_input = np.where(np.isnan(input_array).any(axis=(1, 2)))[0]
nan_indices_truth = np.where(np.isnan(truth_array))[0]

# Combine the indices
nan_indices = np.union1d(nan_indices_input, nan_indices_truth)

# Filter out rows with NaNs in both arrays
input_array_filtered = np.delete(input_array, nan_indices, axis=0)
truth_array_filtered = np.delete(truth_array, nan_indices)

print(f"Original input_array shape: {input_array.shape}")
print(f"Filtered input_array shape: {input_array_filtered.shape}")
print(f"Original truth_array shape: {truth_array.shape}")
print(f"Filtered truth_array shape: {truth_array_filtered.shape}")

# featrue format[1-36]: 1: SASA, 2-6: HDMD, 7-36: HHblits
input_array_filtered = input_array_filtered[:,1:,:] 
x = input_array_filtered.reshape(input_array_filtered.shape[0], -1, input_array_filtered.shape[1])
y = truth_array_filtered

x = np.mean(x, axis=2)

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=42)

print('training x:', x_train.shape)
print('training y:', y_train.shape)
print('test x:', x_test.shape)
print('test y:',y_test.shape)

# Split the dataset into training and testing sets
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=42)
try:
# Initialize the Random Forest Regressor
    random_forest_model = RandomForestRegressor(n_estimators=100, random_state=42, max_features='sqrt')
    random_forest_model.fit(x_train, y_train)
    y_pred = random_forest_model.predict(x_test)

    print(y_pred)
    print('R2 score:', r2_score(y_test, y_pred))
    print('RMSE:', np.sqrt(mean_squared_error(y_test, y_pred)))
    print('Pearson:', pearsonr(y_test, y_pred))

    plt.scatter(y_test, y_pred)
    plt.xlabel('True Values [y_test]')
    plt.ylabel('Predictions [y_pred]')
    plt.title('True vs Predicted Values')

    # Optional: add a line representing perfect predictions
    max_value = max(y_test.max(), y_pred.max())
    min_value = min(y_test.min(), y_pred.min())
    plt.plot([min_value, max_value], [min_value, max_value], color='red', linestyle='--', lw=2)

    plt.show()      
except Exception as e:
    print(e)
    error_id.append(csv_file)  

print(error_id)
import pandas as pd
import json

data_tgs = []

with open('1_tagesschau.jsonl', 'r') as file:
    for line in file:
        data_tgs.append(json.loads(line))

df_tgs = pd.DataFrame(data_tgs)

# Define the output TSV file path
tsv_file_path = 'tgs_txts.tsv'

# Create a new DataFrame with the desired structure
new_data = {
    'ID': range(1, len(df_tgs) + 1),
    'SecondColumn': 'unknown',
    'ExtractedColumn': df_tgs['shorttext']
}

new_df = pd.DataFrame(new_data)

# Save the new DataFrame to a single TSV file without header
new_df.to_csv(tsv_file_path, sep='\t', header=False, index=False)

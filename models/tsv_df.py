import pandas as pd
import os

# Specify the parent directory where TSV files are located
parent_directory = 'gerVader/GerVADER/results/mode_All/sentiments_tgs'

# Get a list of all TSV files in the directory
tsv_files = [file for file in os.listdir(parent_directory) if file.endswith('.tsv')]

# Initialize an empty list to store DataFrames
dfs = []

# Process each TSV file in the directory
for tsv_file in tsv_files:
    tsv_file_path = os.path.join(parent_directory, tsv_file)

    is_negative = 'negative' in tsv_file.lower()
    is_positive = 'positive' in tsv_file.lower()
    is_neutral = 'neutral' in tsv_file.lower()

    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_file_path, sep='\t', header=None)  # Set header=None if your TSV file has no header

    if is_negative:
        df['Sentiment'] = 'negative'
    elif is_positive:
        df['Sentiment'] = 'positive'
    elif is_neutral:
        df['Sentiment'] = 'neutral'

    # Append the DataFrame to the list
    dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
final_df = pd.concat(dfs, ignore_index=True)

csv_file_path = 'sentiments_tgs.csv'
final_df.to_csv(csv_file_path, index=False)

# Display the final DataFrame
print("sentiments_tgs:\n", final_df)

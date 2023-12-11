import pandas as pd
import json

data_tgs = []

with open('1_tagesschau.jsonl', 'r') as file:
    for line in file:
        data_tgs.append(json.loads(line))

df_tgs = pd.DataFrame(data_tgs)
print(df_tgs)

df_sent = pd.read_csv("sentiments_tgs.csv")
df_sent.columns = ["shorttext", "unknown", "?", "?", "sentiment"]
print(df_sent)

count_lines_in_sent = df_tgs['shorttext'].isin(df_sent['shorttext']).sum()

print(f"Number of lines in df_tgs present in df_sent: {count_lines_in_sent}")

# Merge df_tgs with df_sent based on matching shorttext
merged_df = pd.merge(df_tgs, df_sent, how='left', on='shorttext')

# Update df_tgs with sentiment values from df_sent
df_tgs['sentiment'] = merged_df['sentiment']

# Display the updated DataFrame
print("Updated df_tgs:\n", df_tgs)


# Convert DataFrame to JSON Lines (JSONL) format
jsonl_string = df_tgs.to_json(orient='records', lines=True)

# Save the JSONL string to a file
with open('1_tagesschau_sentiment.jsonl', 'w') as file:
    file.write(jsonl_string)

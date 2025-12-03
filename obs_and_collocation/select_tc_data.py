import glob
import pandas as pd

# Get all Data_*.txt files
files = glob.glob('Data_*.txt')

# Initialize empty list to store filtered data
filtered_data = []

# Process each file
for file in files:
    print(f"Processing {file}...")
    df = pd.read_csv(file, sep='\t')
    
    # Filter for cmap > 0 and csec > 0
    filtered = df[(df['cmap'] > 0) & (df['csec'] > 0)]
    
    if not filtered.empty:
        filtered_data.append(filtered)
        print(f"  Found {len(filtered)} matching rows")

# Combine all filtered data
if filtered_data:
    combined_df = pd.concat(filtered_data, ignore_index=True)
    
    # Save to Data_TC.txt
    combined_df.to_csv('Data_TC.txt', sep='\t', index=False)
    print(f"\nCombined dataset saved as Data_TC.txt with {len(combined_df)} total rows")
else:
    print("No matching data found")


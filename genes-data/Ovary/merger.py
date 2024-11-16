import os
import pandas as pd

root_dir = os.getcwd()

for subdir, dirs, files in os.walk(root_dir):
    if subdir == root_dir:
        continue
    tsv_files = sorted((f for f in files if f.endswith('.tsv')), key=lambda x: int(os.path.splitext(x)[0]))
    if not tsv_files:
        continue

    merged_df = pd.read_csv(os.path.join(subdir, tsv_files[0]), sep='\t')

    for tsv in tsv_files[1:]:
        df = pd.read_csv(os.path.join(subdir, tsv), sep='\t')
        merged_df = pd.concat([merged_df, df])

    merged_df.to_csv(os.path.join(subdir, f"{os.path.basename(subdir)}.tsv"), sep='\t', index=False)
    print(f"Merged TSV files have been saved to {os.path.basename(subdir)}.tsv")

print("Finished merging all TSV files.")
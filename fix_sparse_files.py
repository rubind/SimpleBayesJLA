import os
import subprocess
from tqdm import tqdm
from numpy import random

def is_sparse(filepath):
    """Check if a file is sparse."""
    stat = os.stat(filepath)
    logical_size = stat.st_size            # Bytes
    blocks = stat.st_blocks                # Number of 512-byte blocks
    real_size = blocks * 512                # Actual storage used

    return real_size < logical_size

def fix_sparse(filepath):
    """Copy the file non-sparsely."""
    if random.random() < 0.01:
        print(f"Fixing sparse file: {filepath}")
    temp_file = filepath + ".nonsparse"
    # Use 'cp --sparse=never'
    subprocess.run(["cp", "--sparse=never", filepath, temp_file], check=True)
    os.replace(temp_file, filepath)  # Safe overwrite

def walk_and_fix(directory="."):
    """Walk through all files and fix sparse ones."""
    all_files = []
    for root, dirs, files in os.walk(directory):
        for name in files:
            filepath = os.path.join(root, name)
            all_files.append(filepath)

    for filepath in tqdm(all_files, desc="Checking files", unit="file"):
        try:
            if is_sparse(filepath):
                fix_sparse(filepath)
        except Exception as e:
            #print(f"Skipping {filepath}: {e}")
            pass

if __name__ == "__main__":
    walk_and_fix(".")

import os
import subprocess
from tqdm import tqdm

def is_sparse(filepath):
    """Check if a file is sparse."""
    stat = os.stat(filepath)
    logical_size = stat.st_size
    blocks = stat.st_blocks
    real_size = blocks * 512
    return real_size < logical_size

def fix_sparse(filepath):
    """Copy the file non-sparsely."""
    temp_file = filepath + ".nonsparse"
    subprocess.run(["cp", "--sparse=never", filepath, temp_file], check=True)
    os.replace(temp_file, filepath)  # Safe overwrite

def walk_and_fix(directory=".", min_file_size=1024):
    """Walk through all files and fix sparse ones."""

    # Generator version — don't store all files in RAM!
    files = (os.path.join(root, name)
             for root, dirs, filenames in os.walk(directory)
             for name in filenames)

    files = (f for f in files if os.path.isfile(f))  # Only regular files

    # Count files lazily
    files = list(files)
    total_files = len(files)

    with tqdm(total=total_files, desc="Checking files", unit="file") as pbar:
        for filepath in files:
            try:
                stat = os.stat(filepath)
                if stat.st_size < min_file_size:
                    # Skip tiny files that can't meaningfully be sparse
                    pbar.update(1)
                    continue

                if stat.st_blocks * 512 < stat.st_size:
                    try:
                        fix_sparse(filepath)
                        tqdm.write(f"Fixed sparse file: {filepath}")
                    except Exception as e:
                        tqdm.write(f"Failed to fix {filepath}: {e}")
            except Exception as e:
                tqdm.write(f"Skipping {filepath}: {e}")
            pbar.update(1)

if __name__ == "__main__":
    walk_and_fix(".", min_file_size=1024)  # Ignore files <1 KB

#### cleanup

#!/usr/bin/env python3.6
# %%
import os
import sys

# List of files to keep (everything else will be deleted)
files_to_keep = {
    "KPOINTS","POTCAR","QPOINTS",
    "cleanup.py"
}

def cleanup(files_keep):
    # Resolve current script path (robust for most run modes)
    script_path = os.path.abspath(sys.argv[0]) if sys.argv and sys.argv[0] else None
    if script_path and os.path.isfile(script_path):
        # Ensure the script itself is kept
        files_keep = set(files_keep)
        files_keep.add(os.path.basename(script_path))
    else:
        files_keep = set(files_keep)

    root = os.getcwd()

    # Traverse current directory and all subdirectories
    for dirpath, dirnames, filenames in os.walk(root):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)

            # Safety: keep hidden files (e.g., .gitignore, .DS_Store)
            if filename.startswith("."):
                continue

            # Safety: never delete this script itself
            if script_path and os.path.abspath(file_path) == script_path:
                continue

            # Delete everything not in the keep list
            if filename not in files_keep:
                try:
                    os.remove(file_path)
                    print(f"Deleted {file_path}")
                except (FileNotFoundError, PermissionError, IsADirectoryError) as e:
                    print(f"Error deleting {file_path}: {e}")

    print("cleanup complete")

#%%
# Execute cleanup function with the list of files to delete

cleanup(files_to_keep)

# %%
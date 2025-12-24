import os
import shutil

# CONFIGURATION
source_dir = "af3_inputs_ancient"
local_dir = "batch_local_domains"
cloud_dir = "batch_cloud_fulllength"

# Create directories
os.makedirs(local_dir, exist_ok=True)
os.makedirs(cloud_dir, exist_ok=True)

# Counters
count_local = 0
count_cloud = 0

files = [f for f in os.listdir(source_dir) if f.endswith(".json")]

print(f"Sorting {len(files)} files...")

for filename in files:
    src_path = os.path.join(source_dir, filename)
    
    # Logic: If filename contains "full", it goes to Cloud. Else (cat/edit) goes Local.
    if "full" in filename:
        dst_path = os.path.join(cloud_dir, filename)
        shutil.move(src_path, dst_path)
        count_cloud += 1
    else:
        dst_path = os.path.join(local_dir, filename)
        shutil.move(src_path, dst_path)
        count_local += 1

print("-" * 40)
print(f"✅ Moved {count_local} files to '{local_dir}/' (Run on Boron)")
print(f"☁️ Moved {count_cloud} files to '{cloud_dir}/' (Upload to RunPod)")
print("-" * 40)

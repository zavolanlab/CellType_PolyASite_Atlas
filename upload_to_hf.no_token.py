from huggingface_hub import HfApi, login

# 1. Authenticate using your write token
login(token="hf_YOUR_COPIED_TOKEN_HERE")

# 2. Initialize the API
api = HfApi()
repo_id = "magmir71/ucsc_tracks"

print("Starting massive upload to Hugging Face...")

# 3. Upload your entire directory natively
api.upload_folder(
    folder_path=".",               # Uploads the current directory, change to subdirs["UCSCtracks_dir"]
    path_in_repo="",               # Puts it at the root of the repo
    repo_id=repo_id,
    repo_type="dataset",
    ignore_patterns=[
        ".ipynb_checkpoints", 
        "__pycache__", 
        "*.pyc", 
        ".git", 
        ".env",
        "upload_to_hf.py"          # Prevents the script from uploading itself
    ], 
    commit_message="Initial clean upload of trackhubs and large binaries via API"
)

print("Upload complete!")
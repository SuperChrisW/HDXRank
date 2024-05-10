import os
import requests
import urllib.request
import urllib.error

def get_pdbModel(pdb_ids, save_dir):
    # URL for downloading PDB files
    url_template = "https://files.rcsb.org/download/{}.pdb"
    
    # Store failed downloads
    failed_downloads = []

    for pdb_id in pdb_ids:
        # Create the URL for this PDB ID
        url = url_template.format(pdb_id)
        print('downloading {}'.format(pdb_id))
        # Send a GET request to the URL
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            # Write the response content to a file
            fpath = os.path.join(save_dir, f"{pdb_id}.pdb")
            with open(fpath, 'wb') as f:
                f.write(response.content)
        else:
            print(f"Failed to download {pdb_id}. Status code: {response.status_code}")
            failed_downloads.append(pdb_id)

    # Write failed downloads to a file
    fpath = os.path.join(save_dir, "fail_download.txt")
    with open(fpath, 'w') as f:
        for pdb_id in failed_downloads:
            f.write(f"{pdb_id}\n")

def get_AFmodel(uniprot_list, download_dir):
    #download the corresponding model from alphafold database
    os.chdir(download_dir)
    download_list = []
    count = 0
    for key in uniprot_list:
        file = key.split('.')
        key = file[0]
        try:
            url = f"https://alphafold.ebi.ac.uk/files/{key}.pdb" ## need uniprot accession
            filename = f"{key}.pdb"
            check_path = os.path.join(download_dir, filename)
            if not os.path.isfile(check_path):
                urllib.request.urlretrieve(url, filename)
            '''
            url = f"https://alphafold.ebi.ac.uk/files/AF-{key}-F1-predicted_aligned_error_v4.json" ## error matrix
            filename = f"{key}_AF_error.json"
            check_path = os.path.join(download_dir, filename)
            if not os.path.isfile(check_path):
                urllib.request.urlretrieve(url, filename)
            '''
            download_list.append(key)
            count += 1

        except urllib.error.URLError as e:
            print(f"Error occurred while downloading {key} file: {e}")
    print("obtain", count, "PDB files from alphafold database")
    return download_list
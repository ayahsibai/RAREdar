import os
os.environ['SSL_CERT_FILE'] = '/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/certifi/cacert.pem'
#import required packages
from Bio import Entrez
import pandas as pd
import time
import ssl

def get_gene_name_from_refseq(refseq_id):
    Entrez.email = "xususan@reed.edu" # Always tell NCBI who you are. this is required for accessing their server. 
    ssl_context = ssl._create_unverified_context()
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, retmode="xml", context=ssl_context)
    records = Entrez.read(handle)
    return records[0]['GBSeq_definition']

# Load your data
data = pd.read_csv('output/RAREdar_Results.txt', sep='\t')  # adjust the file name and separator to match your file

# Extract the column with RefSeq IDs. Replace 'column_name' with the name of your column
raw_data = data['column_name']

N = 21

# Split the strings by underscore and keep only the last two parts (the identifier)
refseq_ids = ["_".join(s.split('_')[-2:]) for s in raw_data]

# Initialize a dictionary to hold RefSeq ID - gene name pairs
gene_names = {}

for i, refseq_id in enumerate(refseq_ids):
    try:
        gene_name = get_gene_name_from_refseq(refseq_id)
        data.loc[i, 'Gene Name'] = gene_name
        print(f"RefSeq ID: {refseq_id}, Gene Name: {gene_name}")
        time.sleep(0.5)  # a pause of half a second between requests to avoid overloading NCBI servers
    except Exception as e:
        print(f"Failed to get gene name for RefSeq ID {refseq_id}: {e}")
        continue

# Write the results back to the original file
data.to_csv('new_file.txt', sep='\t', index=False)
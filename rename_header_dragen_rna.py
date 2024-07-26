import os
import sys
import csv
from Bio import SeqIO

date_dict = dict()
with open(os.path.join(sys.argv[1], "sample_date.tsv")) as tab:
    reader = csv.DictReader(tab, delimiter='\t')
    for row in reader:
        sample = row['sample'].split(".")[0]
        date = row['date'].split(".")[0]
        date_dict[sample] = date

with open(os.path.join(sys.argv[1], "all_sequences_rna_pathogen.fasta")) as f:
    with open(os.path.join(sys.argv[1], "all_sequences_header_gisaid_ok.fasta"), "w") as out1, open(os.path.join(sys.argv[1], "all_sequences_header_icogen_ok.fasta"), "w") as out2:
        for record in SeqIO.parse(f, "fasta"):
            seq_id = record.id.split("|")[0]
            print(seq_id)
            if seq_id in date_dict:
                id_gisaid = "hCoV-19/Italy/LAZ-IFO-" + seq_id + "/" + date_dict[seq_id].split("-")[0]
                # hCoV-19/Italy/LAZ-IFO-30202/2021
                id_icogen = "LAZ-IFO-" + seq_id
                # LAZ-IFO-30202
                print(id_gisaid)
                print(id_icogen)
                record.description = ""
                record.id = id_gisaid
                SeqIO.write(record, out1, "fasta")
                record.id = id_icogen
                SeqIO.write(record, out2, "fasta")

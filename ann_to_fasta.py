from bioservices import UniProt
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from tqdm import tqdm
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd


def parse_organism_name(str):
    tax_items = str.split('\n')[1]
    name = ''
    if len(tax_items) > 1:
        name = tax_items[-1].split('(')[0]
    return name

def parse_taxonomy(taxonomy):
    needed_levels = ['superkingdom', 'kingdom', 'phylum', 'class', 'genus']
    taxonomy_ = taxonomy.split('\n')
    result = []
    if len(taxonomy_) == 2:
        for level in taxonomy_[1].split(','):
            for name in needed_levels:
                if name == level.split('(')[-1][:-1]:
                    result.append(' '.join(level.split()[:-1]) + '_' + name)
    else:
        return '' 
    return '/'.join(result)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--output', type=str)
    parser.add_argument('--input', type=str)
    parser.add_argument('--save-taxonomy', action='store_true')
    parser.add_argument('--unique', action='store_true')
    args = parser.parse_args()

    stk_file = args.input
    fasta_file = args.output

    records = AlignIO.read(open(stk_file, "r"), "stockholm")
    print(f"{len(records)} records extracted from {stk_file}")

    u = UniProt()
    if args.unique:
        updated_records = defaultdict(str)
    else:
        updated_records = []
    for i, record in enumerate(tqdm(records)):
        organism = ''
        try:
            find = record.id.split("/")[0] + " fragment:false"
            uniprot_taxonomy = u.search(find, columns="lineage") 

            if args.save_taxonomy:
                organism = parse_taxonomy(uniprot_taxonomy)
            else:
                organism = parse_organism_name(uniprot_taxonomy)
                
            if organism == '':
                continue                    
            
            if args.unique and organism not in updated_records:
                updated_records[organism] = SeqRecord(
                    record.seq,
                    id=organism + '/' + record.id,
                    name='',
                    description="",
                )
            elif not args.unique:
                updated_records.append(
                    SeqRecord(
                    record.seq,
                    id=organism + '/' + record.id,
                    name=record.id,
                    description="",
                ))
        except Exception as e:
            print(f"Error processing record {record.id} : {e.with_traceback()}")

    if args.unique:
        updated_records = list(updated_records.values())

    sorted_records = MultipleSeqAlignment(updated_records)
    sorted_records.sort()

    AlignIO.write(sorted_records, open(fasta_file, "w"), "fasta")
    print(f"{len(sorted_records)} records saved to {fasta_file}")


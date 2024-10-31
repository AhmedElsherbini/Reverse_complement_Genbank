#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 15:26:40 2024

@author: ahmed
"""
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import warnings

def reverse_complement_genbank(input_file):
    # Read the GenBank file
    record = SeqIO.read(input_file, "genbank")
    
    # Reverse complement the sequence
    reversed_seq = record.seq.reverse_complement()
    
    # Create a new SeqRecord with the reversed sequence
    reversed_record = SeqRecord(
        reversed_seq,
        id=record.id,
        name=record.name,
        description=record.description,
        annotations=record.annotations,
        dbxrefs=record.dbxrefs,
        features=[]
    )
    
    # Update the features' locations to match the reversed sequence
    for feature in record.features:
        start = len(record) - feature.location.end
        end = len(record) - feature.location.start
        strand = -feature.strand if feature.strand is not None else None
        
        new_location = FeatureLocation(start, end, strand)
        new_feature = SeqFeature(
            location=new_location,
            type=feature.type,
            id=feature.id,
            qualifiers=feature.qualifiers
        )
        reversed_record.features.append(new_feature)
    
    # Write the reversed record to a new GenBank file
    SeqIO.write(reversed_record, "rev_%s.gbk"%(str(input_file).split(".gbk")[0]), "genbank")

def main():
    try:
        parser = argparse.ArgumentParser(description="Reverse the orientation of a GenBank file.")
        parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to your genebank file')
        args = parser.parse_args()
        warnings.filterwarnings('ignore')
        reverse_complement_genbank(args.input_file)
        print("Job is done, pal!")
    except:
        print("Something is wrong with your file")
        
if __name__ == "__main__":
    main()

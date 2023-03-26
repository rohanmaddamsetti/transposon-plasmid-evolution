#!/usr/bin/env python

## check-barcodes-for-duplicates.py

import pprint

input_csv = "../data/tetA-GFP-barcode-clones.csv"

barcodes_list = []

barcode_to_strain_dict = {}

with open(input_csv) as input_fh:
    for i,line in enumerate(input_fh):
        if i == 0: continue ## skip the header
        line = line.strip()
        strain, plasmid, barcode, barcode_verification = line.split(',')
        barcodes_list.append(barcode_verification)
        if barcode_verification not in barcode_to_strain_dict:
            barcode_to_strain_dict[barcode_verification] = [strain]
        else:
            barcode_to_strain_dict[barcode_verification].append(strain)

print("barcode list has " + str(len(barcodes_list)) + " elements")
print("barcode set has " + str(len(set(barcodes_list))) + " elements")

pprint.pprint(barcode_to_strain_dict)

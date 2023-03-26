#!/usr/bin/env python

"""
verify-checksums.py by Rohan Maddamsetti

This script runs md5 on each file in a directory containing a checksum file, based on contexts of the checksum file.

Usage: python verify-checksums.py ../data/genome-sequencing/Maddamsetti_8336_230317A8/Maddamsetti_8336_230317A8.checksum
Usage: python verify-checksums.py ../data/genome-sequencing/Maddamsetti_8337_230322A8/Maddamsetti_8337_230322A8.checksum
"""

import argparse
import os.path
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Run checksums for each file in directory.")
    parser.add_argument('checksum_file', type=str, help='path to checksum file.')
    args = parser.parse_args()

    datadir = os.path.dirname(args.checksum_file)
    
    with open(args.checksum_file, "r") as checksum_fh:
        for line in checksum_fh:
            line = line.strip()
            my_base_filename, my_checksum = line.split()
            my_file_path = os.path.join(datadir, my_base_filename)
            ## run md5 on my_file_path, using the directory for the checksum file,
            ## and get the output.
            md5_call = subprocess.run(["md5", my_file_path], capture_output=True, text=True)
            my_md5_checksum = md5_call.stdout.split("=")[-1].strip()
            ## verify that the checksums match.
            if (my_md5_checksum == my_checksum):
                print(my_md5_checksum, "matches", my_checksum)
            else:
                error_message = "ERROR: " + my_md5_checksum + " does not match " + my_checksum + " for file " + my_file_path
                raise AssertionError(error_message)


main()


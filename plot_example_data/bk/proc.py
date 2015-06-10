#! /usr/local/bin/python

import sys

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')

key2info = {}
for line in hIN:
    F = line.rstrip('\n').split('\t')
    key = F[0] + '\t' + F[1]

    if key not in key2info:
        key2info[key] = ["---", "---"]

    if F[2] in ["point mutation,indel", "structural variation", "HBV integration"]:
        key2info[key][0] = F[2]

    if F[2] in ["HBV fusion", "gene fusion", "over expression", "splicing aberration"]:
        key2info[key][1] = F[2]

hIN.close()


for key in sorted(key2info):
    print key + '\t' + '\t'.join(key2info[key])


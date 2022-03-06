import os
import re
from math import sqrt
import argparse

parser = argparse.ArgumentParser(description="Calculate the scalar moment of a given event")
parser.add_argument("event_code", type=str, help="GCMT event code")
parser.add_argument("outdir", type=str, help="Folder in which to save output file")
args = parser.parse_args()

filename = os.path.join("..", "database", args.event_code, args.event_code)

# all of this is just pulling out the moments infomration from the CMTsolution file ( I THINK THIS IS A STUPID WAY TO DO THIS
with open(filename, "r") as f:
    lines = f.readlines()
full_list = []
for i in range(7, 13):
    s = lines[i]
    match_number = re.compile("-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?")
    line_list = [float(x) for x in re.findall(match_number, s)]
    full_list.append(line_list)

Mrr = full_list[0]
Mtt = full_list[1]
Mpp = full_list[2]
Mrt = full_list[3]
Mrp = full_list[4]
Mtp = full_list[5]

# this is just putting the exponetials back together becuase it did not like the way the exponetials were written becuase they had a plus in front so I removed it and then put the bits back together
Mrr = Mrr[0] * (10 ** Mrr[1]) * 1e-7  # converted to Nm
Mtt = Mtt[0] * (10 ** Mtt[1]) * 1e-7
Mpp = Mpp[0] * (10 ** Mpp[1]) * 1e-7
Mrt = Mrt[0] * (10 ** Mrt[1]) * 1e-7
Mrp = Mrp[0] * (10 ** Mrp[1]) * 1e-7
Mtp = Mtp[0] * (10 ** Mtp[1]) * 1e-7

# now I am going to attempt to work out the moment magnitude using the formula from the shearer textbook

Mo = (1 / sqrt(2)) * (
    sqrt((Mrr ** 2) + (Mtt ** 2) + (Mpp ** 2) + (Mrt ** 2) + (Mrp ** 2) + (Mtp ** 2))
)

with open(os.path.join(args.outdir, "scalar_moment.txt"), "w") as f:
    f.write(f"{Mo}")

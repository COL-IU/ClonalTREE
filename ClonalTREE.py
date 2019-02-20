from predict import *
import sys

if len(sys.argv) <= 4:
    sys.exit("Usage: python3 ClonalTree.py <F matrix> <algo> <fail_threshold> <out_prefix> <optional: clones>")

algo = int(sys.argv[2])
fail_threshold = Decimal(sys.argv[3])
out_prefix = sys.argv[4]

f1 = open(sys.argv[1])
F_in, var_in = read_F(f1)
f1.close()

clones = []
if len(sys.argv) == 6:
    f2 = open(sys.argv[5])
    lines = f2.readlines()
    f2.close()
    for line in lines:
        clone = list(map(int, line.split(None)))
        clones.append(clone)

var1, parents1, p1, c1 = non_square_predict(F_in, var_in, algo, fail_threshold, clones=clones)

o1 = open(out_prefix + ".tree", "w")
print_parents(var1, parents1, o1)
o1.close()

c_norm = []
for row in c1:
    for i in range(0, len(row)):
        if row[i] < 0:
            row[i] = 0
    s = sum(row)
    for i in range(0, len(row)):
        row[i] = row[i] / s
    c_norm.append(row)

o2 = open(out_prefix + ".C", "w")
o2.write(" ".join(var1) + "\n")
write_dm(c_norm, o2)
o2.close()


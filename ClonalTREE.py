from predict import *
import sys

if len(sys.argv) <= 4:
    sys.exit("Usage: python3 ClonalTREE.py <VAF file> <algorithm> <fail threshold> <out prefix> <clones (optional)>\n\n"
             "VAF file:\t[String] Input file containing the variant allele frequencies matrix (F).\n"
             "algorithm:\t[Int] 0 (RP-RT); 1 (EP-ET); 2 (EP-GT); 3 (GP-GT).\n"
             "fail threshold:\t[Float] Minimum value allowed in the C matrix to define failure.\n"
             "out prefix:\t[String] File path to prefix all output files.\n"
             "clones:\t\t[String] A file containing the composition of known clones. (Optional argument)\n")

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

if len(parents1) == 0:
    print("\nUnable to find a valid solution under the given parameters. Please try lowering the fail threshold.\n")
    sys.exit()

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

print("\nSolution output to files: \n" + out_prefix + ".tree\n" + out_prefix + ".C\n")

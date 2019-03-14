import sys
from decimal import *
from copy import *


zero = Decimal('0.0')


def reorder(l, order):
    return [l[i] for i in order]


def rearrange(F, variants, remove_redundancy=False):
    F_out = deepcopy(F)
    m = len(F_out)
    if m == 0:
        return F, variants
    if len(F[0]) == 0:
        return F, variants
    # Make F_out a diagonal (step) matrix
    F_out = list(map(list, zip(*F_out)))
    order = []
    i = 0
    for row in F_out:
        for i in range(0, len(row)):
            if row[i] > zero:
                break
        order.append(i)
    indices = list(range(0, len(F_out)))
    [order, indices, F_out] = list(zip(*sorted(zip(order, indices, F_out), key=lambda x: x[0])))
    variants_out = reorder(variants, indices)

    non_founder_F = []
    non_founder_var = []

    for i in range(0, len(variants_out)):
        row = F_out[i]
        is_founder = True
        for j in range(0, len(row)):
            if row[j] < Decimal('0.9'):
                is_founder = False
                break
        if not is_founder:
            non_founder_F.append(row)
            non_founder_var.append(variants_out[i])

    non_founder_F = list(map(list, zip(*non_founder_F)))

    if remove_redundancy:
        # Remove unnecessary time points
        z = 0
        z_prev = -1
        new_F = []
        for j in range(0, len(non_founder_F)):
            row = non_founder_F[j]
            for z in range(len(row)-1, -2, -1):
                if row[z] > zero:
                    break
            if z - z_prev > 0:
                new_F.append(non_founder_F[j])
            z_prev = z
        non_founder_F = new_F

    return non_founder_F, non_founder_var


def is_valid(t):
    if len(t) < 2:
        return False
    ts = sorted(t)
    # if max(ts) != 6:
    #     return False
    for i in range(1, len(ts)):
        if ts[i] - ts[i-1] > 1:
            return False
    return True


def write_dm(mat, out):
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write("{0:.5f}".format(mat[i][j]) + "\t")
        out.write("\n")


if len(sys.argv) < 3:
    sys.exit("Usage: filterFreqs.py <freqs> <outprefix>")

f = open(sys.argv[1])
lines = f.readlines()
f.close()

times = {}

for line in lines:
    words = line.split(None)
    locus = words[1]
    time = int(words[0][1])
    if locus in times:
        times[locus].append(time)
    else:
        times[locus] = [time]

filtered = []

for locus in sorted(times.keys()):
    if is_valid(times[locus]):
        filtered.append(locus)

f1 = open(sys.argv[2]+".filtered.freqs", "w")
F = [[zero for _ in range(len(filtered))] for _ in range(6)]
for line in lines:
    words = line.split(None)
    locus = words[1]
    time = int(words[0][1])-1
    if locus in filtered:
        F[time][filtered.index(locus)] = Decimal(words[5])
        f1.write(line.strip() + "\n")
f1.close()

F_new, var_new = rearrange(F, filtered, False)

f2 = open(sys.argv[2]+".vaf", "w")
if len(F_new) > 0:
    n = len(F_new[0])
    f2.write("\t".join(list(map(str, range(1,n+1)))))
    f2.write("\n")
    write_dm(F_new, f2)
f2.close()

f3 = open(sys.argv[2]+".var", "w")
for var in var_new:
    f3.write(var + "\n")
f3.close()


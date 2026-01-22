

using DelimitedFiles
using Printf

input = ARGS[1] #"results/nosnoc-ipopt-scholtes.csv"

data = readdlm(input)

m, n = size(data)
@printf(
    "%-28s %5s %5s %5s %3s %10s %5s %10s\n",
    "instance",
    "n",
    "m",
    "ncc",
    "st.",
    "obj.",
    "it",
    "time",
)
@printf(
    "%-28s %5s %5s %5s %3s %10s %5s %10s\n",
    "--------",
    "----",
    "----",
    "----",
    "---",
    "---------",
    "----",
    "---------",
)

for k in 2:m
    @printf("%-28s %5i %5i %5i ", data[k, 1], data[k, 2], data[k, 3], data[k, 4])
    status = (data[k, 5] == 1 || data[k, 5] == 2) ? 1 : 0
    @printf("%3i %10.3f %5i %10.2f\n", status, data[k, 6], data[k, 7], data[k, 8])
end






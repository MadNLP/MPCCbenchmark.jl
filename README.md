# MPCCBenchmark

A benchmark for large-scale mathematical programs with complementarity constraints (MPCCs).

The instances are implemented in JuMP and reformulated with [ComplementOpt](https://github.com/blegat/ComplementOpt.jl).

The package implements the following benchmarks:

- Nonsmooth optimal control problems formulated with complementarity constraints (the instances are a subset of [NOSNOC](https://arxiv.org/html/2312.11022v2));
- Power-flow problems with PV/PQ switches, following [this article](https://link.springer.com/article/10.1007/s10589-015-9745-5);
- SCOPF problems with nonsmooth recourse, following [this article](https://arxiv.org/abs/2510.13333).

## Basic usage

We recommend using the script `main.jl` to launch the benchmark in parallel on `p` processors using Distributed.jl. It proceeds as follows:

```shell
julia -p 6 --project=. main.jl --solver={ipopt,madnlpc} --benchmark={pscc-pf,pscc-scopf,nosnoc}

```
The results are dumped in a CSV file in the directory `./results/`.


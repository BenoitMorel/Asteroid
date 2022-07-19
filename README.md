# Asteroid (Accurate Species Tree Estimation RObust to Incomplete Data sampling)

Asteroid is a method that infers an unrooted species tree from a set of unrooted gene trees. It is robust to both incomplete lineage sorting and high proportions of missing genes. 

Asteroid computes for each input gene tree a distance matrix based on the gene internode distance. Then, it computes a species tree from this set of distance matrices under the minimum balanced evolution principle.

Asteroid is parallelized, and can take as input multi-furcating and multi-copy gene trees. It can computes bootstrap support values, either from the input gene trees, or from an alternative set of bootstrap gene trees.

For more information about the method, please check our manuscript (TODO!).

## Installation

To download Asteroid, please use git,  and clone with --recursive!!!
```
git clone --recursive https://github.com/BenoitMorel/Asteroid
```

To install Asteroid, run:
```
./install.sh
```
The resulting executable is: `build/bin/asteroid`


## Running Asteroid

The minimum syntax is:
```
./asteroid -i GENE_TREES 
```

To run with parallelization:
```
mpiexec -np CORES asteroid -i GENE_TREES
```

Here is a list of the available options:


|    Command                |  Comment  |
|---------------------------|-----------|
| `-h`, `--help`  | Print the help message |
| `-i`, `--input-gene-trees <STRING>`     | 	 Path to file containing one gene tree per line, in newick format |
|  `-r`,  `--random-starting-trees <INTEGER>` 	   |   Number of starting random trees. `0` to start from an ASTRID tree. Default value is `1`. |
|`-b`, `--bs-trees <INTEGER>`             |	 Number of bootstrap trees to compute. Default is 0. |
|`-p`, `--prefix <STRING>`               |  Prefix for the output files.|
|`-m`, `--gene-species-mapping <STRING>`  	  | Path to the mapping file. If unset, Asteroid assumes that the gene tree leaf labels correspond to the species names. See our wiki for the format.|
|`--seed <INTEGER>`                      |	 Random seed. Default is 1.| 
|`-n`, `--no-correction`                  	| When set, Asteroid disables its missing data correction, and the algorithm is then equivalent to ASTRID.|
|`--min-bl <REAL>`                       	| Branches with a length under this value will be contracted into polytomies before computing the distance matrices. Default value is `-1.0` (no contraction) | 
| `--use-gene-bl`               |           	 Use the gene tree branch lengths instead of their internode distance to build the distance matrices |


## Report bugs, asking questions




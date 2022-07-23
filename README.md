# Asteroid (Accurate Species Tree Estimation RObust to Incomplete Data sampling)

Asteroid is a method that infers an unrooted species tree from a set of unrooted gene trees. It is robust to both incomplete lineage sorting and high proportions of missing genes. 

Asteroid computes for each input gene tree a distance matrix based on the gene internode distance. Then, it computes a species tree from this set of distance matrices under the minimum balanced evolution principle. Asteroid is robust to high proportions of missing genes because, unlike ASTRID, it does not merge all distance matrices into one single distance matrix, but rather keeps all the distance matrices distinct and applies a bias correction to each of them. For a more detailed explanation of the method, please read our manuscript (Todo: add link).

Asteroid is parallelized, and can take as input multi-furcating and multi-copy gene trees. It can computes bootstrap support values, either from the input gene trees, or from an alternative set of bootstrap gene trees.


## Installation

If you do not need the parallelized version, you can get the linux and mac (asteroid_linux and asteroid_mac) binaries here: https://github.com/BenoitMorel/Asteroid/releases/tag/1.0
If it doesn't work or if you want to parallelize with MPI, please read the following instructions:


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

We provide several example command lines in the examples folder. To run the examples, you have to be at the root of the repository (and run, for instance, `./examples/run_asteroid_simple.sh`).

Here is a list of the available options:


|    Command                |  Comment  |
|---------------------------|-----------|
| `-h`, `--help`  | Print the help message |
| `-i`, `--input-gene-trees <STRING>`     | 	 Path to file containing one gene tree per line, in newick format |
|  `-r`,  `--random-starting-trees <INTEGER>` 	   |   Number of starting random trees. `0` to start from an ASTRID tree. Default value is `1`. |
|`-b`, `--bs-replicates <INTEGER>`             |	 Number of bootstrap trees to compute. Default is 0. |
|`--input-bs-gene-trees <STRING>`   | The set of gene trees used to perform boostrapping. If not set, the intput gene trees (`-i`) are used instead.|
|`-p`, `--prefix <STRING>`               |  Prefix for the output files.|
|`-m`, `--gene-species-mapping <STRING>`  	  | Path to the mapping file. If unset, Asteroid assumes that the gene tree leaf labels correspond to the species names. The syntax is described below.|
|`--seed <INTEGER>`                      |	 Random seed. Default is 1.| 
|`-n`, `--no-correction`                  	| When set, Asteroid disables its missing data correction, and the algorithm is then equivalent to ASTRID.|
|`--min-bl <REAL>`                       	| Branches with a length under this value will be contracted into polytomies before computing the distance matrices. Default value is `-1.0` (no contraction) | 
| `--use-gene-bl`               |           	 Use the gene tree branch lengths instead of their internode distance to build the distance matrices |

## Mapping file

A gene-species mapping file is necessary if you have multi-copy gene trees (but also works for single-copy gene trees). There are two supported syntaxes:


### Species-genes mapping file (PHYLDOG format)
```
species1:gene1;gene2;gene3
species2:gene4;gene6
species3:gene5
```


### Genes-species mapping file (TreeRecs format)
```
gene1 species1
gene2 species1
gene3 species1
gene4 species2
gene6 species2
gene5 species3
```


## Report bugs, ask questions

If you encounter a bug or have any question regarding Asteroid, please open an issue at https://github.com/BenoitMorel/Asteroid/issues




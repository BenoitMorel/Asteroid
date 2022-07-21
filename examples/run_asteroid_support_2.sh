
# Run asteroid + boostrapping
# The bootstrapping is done by sampling the gene trees
# from an alternative set of gene trees

data="examples/data"
genetrees="$data/gene_trees_1.newick"
bsgenetrees="$data/bs_gene_trees_1.newick"
outputdir="examples/output"
mkdir -p $outputdir
prefix="$outputdir/support_2"
rm $prefix*

./build/bin/asteroid -i $genetrees -p $prefix --bs-replicates 100 --input-bs-gene-trees $bsgenetrees


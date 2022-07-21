
# Run asteroid + boostrapping
# The bootstrapping is done by sampling the gene trees
# from the input gene trees ($genetrees)

data="examples/data"
genetrees="$data/gene_trees_1.newick"
outputdir="examples/output"
mkdir -p $outputdir
prefix="$outputdir/support"
rm $prefix*

./build/bin/asteroid -i $genetrees -p $prefix --bs-replicates 100



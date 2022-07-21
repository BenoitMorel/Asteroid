data="examples/data"
genetrees="$data/gene_trees_1.newick"
outputdir="examples/output"
mkdir -p $outputdir
prefix="$outputdir/starting_trees"
rm $prefix*

./build/bin/asteroid -i $genetrees -p $prefix -r 20



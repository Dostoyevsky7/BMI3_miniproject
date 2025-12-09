nucmer -p out synthetic_genome.fasta synthetic_genome.fasta
mummerplot --fat --png -p out out.delta
show-coords -rcl out.delta > out.coords
show-coords -rcl -T out.delta > syn_out.tsv
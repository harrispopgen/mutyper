#!/usr/bin/env bash
set -eo pipefail

if [ $1 == "1" ]; then
	
	snakemake \
		--configfile .test/aln_config.yaml \
		--cores $@

elif [ $1 == "2" ]; then

	snakemake \
		--configfile .test/chain_config.yaml \
		--cores $@

elif [ $1 == "3" ]; then

	snakemake \
		--configfile .test/config.yaml \
		--cores $@

fi

echo 
# check lint 
snakemake \
	--configfile .test/config.yaml \
	--cores 4 --lint 

snakefmt . 
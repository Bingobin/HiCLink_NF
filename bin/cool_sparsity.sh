#!/usr/bin/bash

COOL=$1

nbins=$(cooler info $COOL | jq .nbins)
nnz=$(cooler info $COOL | jq .nnz)
total=$(( nbins * (nbins + 1) / 2 ))
density=$(echo "scale=4; $nnz / $total" | bc)
sparsity=$(echo "scale=4; 1 - $density" | bc)
echo "$COOL: density=$density, sparsity=$sparsity"

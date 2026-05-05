#!/bin/bash

date
echo "Starting simulation screens..."

for i in {1..10}
do
    screen -dmS sims$i bash -c "Rscript sims$i.R"
    echo "Started sims$i"
done

echo "All simulations launched."
echo "Check running screens with: screen -ls"

echo 'job done'
date

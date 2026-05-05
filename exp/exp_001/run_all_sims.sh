#!/bin/bash

echo "$(date) Starting simulation screens..."

DIR="$(pwd)"

for i in {2..10}
do
    screen -dmS sims$i bash -c "cd $DIR && Rscript sims$i.R > sims$i.out 2>&1"
    echo "Started sims$i"
done

echo "All simulations launched."
echo "Check running screens with: screen -ls"
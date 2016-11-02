#!/bin/bash
grep new_matches ../tmp.log | awk '{print $2}' > new_matches
grep better_matches ../tmp.log | awk '{print $2}' > better_matches
join <(cat -n better_matches) <(cat -n new_matches) | awk '{print $2 + $3}' > sum
echo 'set term png; plot "sum" using ($0):(d($1)) smooth csplines' | gnuplot > derivative.png

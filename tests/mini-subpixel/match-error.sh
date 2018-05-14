#!/bin/bash

num=$(ls matches-*ppm | tr -d '[:alpha:].-' | sort --numeric-sort | tail -n 1)
../../validate matches-$num.ppm | awk '{ print $3 - $2 }' | sort -n | uniq -c

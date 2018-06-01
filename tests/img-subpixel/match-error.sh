#!/bin/bash

echo "X: "
../../validate $1 | awk '/^X/ { print $3 - $2 }' | sort -n | uniq -c

echo "Y: "
../../validate $1 | awk '/^Y/ { print $3 - $2 }' | sort -n | uniq -c

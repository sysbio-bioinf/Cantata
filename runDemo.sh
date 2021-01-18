#!/bin/bash

# Run demo after installation of CANTATA using make

./bin/cantata --optimize -n ./fissionyeast_trunc_1.txt -r ./fissionyeast-rules.txt -o ./fissionyeast_results.txt -on ./fissionyeast_result_%d.txt -ps 100 -no 200 -ni 100 -ns 1
#!/bin/bash

for i in $(seq 500 1001); do python3.4 main.py $i; done


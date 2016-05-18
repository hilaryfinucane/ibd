#!/bin/bash

data_dir=/groups/price/hilary/ibd/data

# doubleton frequency = 1/N = 1/489 = 0.002045
grep "0\.002045" $data_dir/1000G.EUR.QC.22.frq | 
    sed 's/\s\+/\t/g' | cut -f3 > $data_dir/1000G.EUR.QC.22.doubletons

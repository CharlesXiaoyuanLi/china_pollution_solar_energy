#! /bin/bash

for i in $(ls |grep .gz)
do
    gunzip $i
done

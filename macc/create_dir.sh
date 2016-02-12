#!/bin/bash 

while read line
do
    name=$line
    mkdir $name
done < ../../../../cn_province_nmlist.txt

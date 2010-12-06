#!/bin/bash

list=""
File=files.txt

{
while read line; do
      list="$list $line"
done
} < $File

echo $list
hadd /tmp/hwoehri/merge.root $list

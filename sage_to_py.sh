#!/bin/sh

rm *.py
sage --preparse *.sage
for file in *.sage.py; do
    mv -- "$file" "${file%.sage.py}.py"
done

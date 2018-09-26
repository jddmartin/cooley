#!/bin/sh
coverage3 run --source cooley.py -m pytest -s
coverage3 report
coverage3 html --directory htmlcov_generated 
pycodestyle --verbose --ignore=E402,W504  *.py examples_and_tests/*.py > pycodestyle_output_generated.txt

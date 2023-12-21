#!/bin/bash
set -o errexit

find . -name '*.py' \
  | xargs pylint \
  --errors-only 

echo 'Passes pylint check'

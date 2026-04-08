#!/bin/bash
set -o errexit

ruff check gtfparse/ tests/

echo 'Passes ruff check'

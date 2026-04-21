#!/usr/bin/env bash

set -e

SOURCES="gtfparse tests"

echo "Running ruff format..."
ruff format $SOURCES

echo "Formatting complete!"

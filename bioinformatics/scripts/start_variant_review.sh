#!/bin/bash
set -e
cd /home/ravanox/codex/bioinformatics
exec tools/wes-lab-env/bin/python -m wes_lab.igv_review "$@"

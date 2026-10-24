#!/usr/bin/env bash
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
exec tools/wes-lab-env/bin/python -m wes_lab "$@"

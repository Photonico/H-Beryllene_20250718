#!/bin/bash
# Lightweight login-node downloads and isolated Python environment setup.
set -euo pipefail
tools_dir=$(cd "${1:?Expected the private tools directory}" && pwd -P)
test ! -e "$tools_dir/DOWNLOADS_READY"
cd "$tools_dir"
for target in venv wannier90-3.1.0.tar.gz python-packages.txt irrep-help.txt; do
  if [ -e "$target" ]; then
    printf 'Refusing to overwrite existing private tool data: %s\n' "$target" >&2
    exit 2
  fi
done
python3 -c 'import sys; assert (3, 9) <= sys.version_info[:2] < (3, 13), "Pinned packages require Python 3.9-3.12"'
curl --fail --location --retry 3 https://github.com/wannier-developers/wannier90/archive/refs/tags/v3.1.0.tar.gz -o wannier90-3.1.0.tar.gz
tar tzf wannier90-3.1.0.tar.gz > /dev/null
python3 -m venv venv
venv/bin/python -m pip install 'pip<26' 'setuptools<81' wheel
venv/bin/python -m pip install 'numpy==1.26.4' 'scipy==1.13.1' 'spglib==2.5.0' 'irreptables==2.0.0' 'irrep==2.1.3' 'z2pack==2.2.1'
venv/bin/python -m pip check
venv/bin/python -m pip freeze > python-packages.txt
venv/bin/irrep --help > irrep-help.txt
venv/bin/python -c 'import z2pack, irrep, spglib; print("Python packages import successfully")'
sha256sum wannier90-3.1.0.tar.gz > DOWNLOADS_READY

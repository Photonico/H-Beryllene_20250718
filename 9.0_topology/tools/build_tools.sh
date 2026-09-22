#!/bin/bash
# Build only inside a fresh private prefix; never alter the site VASP installation.
set -euo pipefail
tools_dir=$(cd "${1:?Expected the private tools directory}" && pwd -P)
cd "$tools_dir"
test -s DOWNLOADS_READY
sha256sum --check DOWNLOADS_READY
test ! -e READY
test ! -e vasp
test ! -e wannier90-3.1.0
command -v ifx
command -v mpiifort
mkdir vasp
cp /cmt2/ocon2505/VASP/vasp.6.5.0/makefile vasp/
cp /cmt2/ocon2505/VASP/vasp.6.5.0/makefile.include vasp/
cp -r /cmt2/ocon2505/VASP/vasp.6.5.0/src vasp/
cp -r /cmt2/ocon2505/VASP/vasp.6.5.0/arch vasp/
mkdir vasp/bin vasp/build
# Use the same Intel Fortran family as the site's VASP 6.5 makefile.
tar xzf wannier90-3.1.0.tar.gz
cd wannier90-3.1.0
cat > make.inc <<'EOF'
F90 = ifx
FCOPTS = -O2
LDOPTS = -O2
COMMS = serial
LIBS = -qmkl=sequential
EOF
# The two targets share Fortran objects/modules in src/obj. Build serially
# and separately so recursive make processes cannot race on those files.
make -j 1 lib
make -j 1 wannier
cd "$tools_dir/vasp"
cat >> makefile.include <<EOF

# Private W90 3.1 serial interface; no changes to the site installation.
CPP_OPTIONS += -DVASP2WANNIER90
LLIBS += $tools_dir/wannier90-3.1.0/libwannier.a
EOF
make DEPS=1 -j 6 ncl
ldd bin/vasp_ncl > "$tools_dir/vasp-linked-libraries.txt"
if grep -q 'not found' "$tools_dir/vasp-linked-libraries.txt"; then exit 8; fi
test -x bin/vasp_ncl
test -s "$tools_dir/wannier90-3.1.0/libwannier.a"
# Verify the static Wannier interface was actually linked, not only compiled.
nm -g bin/vasp_ncl | grep -E '(^|[[:space:]])wannier_setup_$' > "$tools_dir/wannier-link-symbols.txt"
"$tools_dir/venv/bin/python" -m pip check
"$tools_dir/venv/bin/irrep" --help > "$tools_dir/irrep-help.txt"
"$tools_dir/venv/bin/python" -c 'import z2pack, irrep, spglib; print("Python packages import successfully")'
sha256sum bin/vasp_ncl "$tools_dir/wannier90-3.1.0/libwannier.a" > "$tools_dir/READY"
printf 'Tool build completed; WCC job performs the runtime overlap pilot.\n'

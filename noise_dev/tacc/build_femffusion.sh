#!/bin/bash
# Build FEMFFUSION (https://github.com/Zonni/FEMFFUSION, GPL-3.0) outside
# the repo, for the deterministic cross-checks (Stage 9a C5G7, and the
# informational 4b/6c/8c comparisons). Run on a login node:
#
#   bash tacc/build_femffusion.sh            # from noise_dev/
#
# FEMFFUSION needs deal.II (>= 8.0) configured with PETSc and SLEPc, real
# scalars: its noise solver works on 2 x 2 real blocks. It is tested with
# deal.II 9.5.1 / PETSc 3.19 (Ubuntu 24.04 packages, the local build). LS6
# has no 9.5; 9.4.0-real is the closest older release (9.6 removed many
# functions 9.5 deprecated). TACC's PETSc modules include SLEPc.
#
# Environment:
#   FEMFFUSION_DIR       clone and build here (default: <parent of
#                        $OPENMC_DEV_SRC>/FEMFFUSION, else $SCRATCH/FEMFFUSION)
#   FEMFFUSION_PREREQ    modules to load before deal.II (compiler, MPI), as
#                        listed by `module spider <deal.II module>`; for
#                        dealii/9.4.0-real on LS6: "intel/19.1.1 impi/19.0.9"
#   FEMFFUSION_DEALII    deal.II module (default dealii/9.4.0-real)
#
# Pinned to commit 1611f90 (design doc) plus one local fix: std::filesystem
# replaced by POSIX stat/mkdir for a pre-C++17 deal.II build (see below).
# The link line is also adjusted where deal.II names a system library that
# LS6 lacks (/usr/lib64/libtbb.so; see below), and PyTrilinos, whose
# constructor crashes femffusion.exe at load, is dropped from the linked
# libraries (patchelf; see the end). None of this changes behaviour.
# TACC's deal.II has no LAPACK: FEMFFUSION is compiled with
# DEAL_II_WITH_LAPACK (its own LAPACK calls go to MKL) and the three ROM
# unit tests, which need LAPACK inside libdeal_II, are skipped (see below).
#
# The same modules must be loaded to run femffusion.exe; the script writes
# them to $FEMFFUSION_DIR/modules.sh.
set -eo pipefail

DEALII=${FEMFFUSION_DEALII:-dealii/9.4.0-real}
PREREQ=${FEMFFUSION_PREREQ:-}
if [ -n "${OPENMC_DEV_SRC:-}" ]; then
  DEST=${FEMFFUSION_DIR:-$(dirname "$OPENMC_DEV_SRC")/FEMFFUSION}
else
  DEST=${FEMFFUSION_DIR:-$SCRATCH/FEMFFUSION}
fi

# `module` is a shell function; make sure it exists in this non-interactive
# shell (Lmod exports LMOD_PKG)
if ! type module > /dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "${LMOD_PKG:-/opt/apps/lmod/lmod}/init/bash"
fi

module reset > /dev/null 2>&1 || true
for m in $PREREQ; do
  module load "$m"
done
if ! module load "$DEALII"; then
  echo
  echo "Could not load $DEALII. Its prerequisites are:"
  module spider "$DEALII" 2>&1 | sed -n '/You will need/,/^$/p'
  echo "Set FEMFFUSION_PREREQ to one of those lines (space-separated) and rerun."
  exit 1
fi
module list 2>&1 | sed -n '1,20p'

if [ ! -d "$DEST/.git" ]; then
  git clone https://github.com/Zonni/FEMFFUSION.git "$DEST"
fi
cd "$DEST"
git fetch -q origin || true
git checkout -q -f 1611f90       # -f: start from clean sources (patch below)
rm -f CMakeCache.txt.tmp*

# Local fix: TACC's deal.II 9.4 (Intel 19.1) compiles in a pre-C++17 mode
# and puts its -std flag after ours, so <filesystem> has no
# std::filesystem. FEMFFUSION uses it only to create an output folder;
# replace those calls with POSIX stat/mkdir. Behaviour is unchanged.
python3 - <<'PYEOF'
from pathlib import Path
helpers = '''#include <sys/stat.h>  // local fix (noise_dev/tacc/build_femffusion.sh)
static bool femff_dir_exists(const std::string &p)
{ struct stat s; return stat(p.c_str(), &s) == 0; }
static void femff_mkdir(const std::string &p) { mkdir(p.c_str(), 0755); }
'''
for name in ("static_diffusion.cc", "static_full_spn.cc", "static_spn.cc",
             "static_sdpn.cc"):
    f = Path("src") / name
    s = f.read_text()
    if "femff_mkdir" in s:
        continue
    assert s.count("#include <filesystem>") == 1, name
    s = s.replace("#include <filesystem>", "#include <string>\n" + helpers)
    s = s.replace("std::filesystem::exists(folderPath)",
                  "femff_dir_exists(folderPath)")
    s = s.replace("std::filesystem::create_directory(folderPath)",
                  "femff_mkdir(folderPath)")
    assert "std::filesystem" not in s, name
    f.write_text(s)

# The three ROM (POD/LUPOD) unit tests need LAPACK inside libdeal_II
# (LAPACKFullMatrix::compute_svd, invert); a deal.II built without LAPACK
# runs them as silent no-ops. Guard them (no-op unless
# FEMFF_SKIP_ROM_TESTS is defined; see below).
f = Path("src/test.cc")
s = f.read_text()
rom = ("  test_POD_groupwise();\n  test_LUPOD_extended();\n"
       "  test_LUPOD_extended_group_wise();\n")
if "FEMFF_SKIP_ROM_TESTS" not in s:
    assert s.count(rom) == 1
    s = s.replace(rom, "#ifndef FEMFF_SKIP_ROM_TESTS  // local fix\n" + rom
                  + "#endif\n")
    f.write_text(s)
PYEOF

# Use deal.II's own compiler. FEMFFUSION's CMakeLists calls project()
# before deal.II's autopilot, so CMake would otherwise pick c++ (g++) and
# then get deal.II's Intel flags (-w2, -wd..., -qopenmp-simd, -ip).
DEALII_ROOT=${TACC_DEALII_DIR:-${DEAL_II_DIR:-}}
CFG=$(find "$DEALII_ROOT" -name deal.IIConfig.cmake 2> /dev/null | head -1)
if [ -z "$CFG" ]; then
  echo "deal.IIConfig.cmake not found under '$DEALII_ROOT'"
  exit 1
fi
DEALII_CXX=$(sed -n 's/^[Ss][Ee][Tt](DEAL_II_CXX_COMPILER "\(.*\)")/\1/p' \
  "$CFG" | head -1)
echo "deal.II config: $CFG"
echo "deal.II C++ compiler: $DEALII_CXX"

# TACC's dealii/9.4.0-real is configured without LAPACK (config.h:
# "#undef DEAL_II_WITH_LAPACK"), although MKL is on its link line through
# PETSc/Trilinos. deal.II's inline LAPACK wrappers (lapack_templates.h)
# then do nothing in Release mode, so FEMFFUSION's own geev calls
# (compute_eigs, used by the 2-group eigensolver) return zeros. Defining
# DEAL_II_WITH_LAPACK for FEMFFUSION makes those calls real (MKL); the only
# other header that reads it, full_matrix.templates.h, then uses BLAS
# instead of its own loops, which is equivalent. LAPACK inside libdeal_II
# stays off, which only the ROM features use: skip their three tests.
EXTRA_CXX=""
CFGH=$(find "$DEALII_ROOT" -path '*deal.II/base/config.h' 2> /dev/null |
  head -1)
if [ -n "$CFGH" ] && ! grep -q '^#define DEAL_II_WITH_LAPACK' "$CFGH"; then
  EXTRA_CXX="-DDEAL_II_WITH_LAPACK -DFEMFF_SKIP_ROM_TESTS"
  echo "deal.II has no LAPACK: compiling with $EXTRA_CXX"
fi

# Fresh configure: a cache from an earlier attempt keeps its compiler
rm -rf CMakeCache.txt CMakeFiles

# deal.II's autopilot build, as in the FEMFFUSION README.
# CMAKE_POLICY_VERSION_MINIMUM: CMake >= 4 rejects FEMFFUSION's
# cmake_minimum_required(VERSION 2.8.8) (as for Abeille's docopt.cpp)
cmake . -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DDEAL_II_DIR="$(dirname "$CFG")" \
  ${DEALII_CXX:+-DCMAKE_CXX_COMPILER="$DEALII_CXX"} \
  -DCMAKE_CXX_FLAGS="$EXTRA_CXX"
make release
if [ -n "$EXTRA_CXX" ] && ! grep -qiE 'mkl|lapack' \
    CMakeFiles/femffusion.exe.dir/link.txt; then
  echo "DEAL_II_WITH_LAPACK forced, but no LAPACK on the link line"
  exit 1
fi

# deal.II's exported link line may name unversioned system libraries that
# existed on TACC's build host but not on LS6 (e.g. /usr/lib64/libtbb.so;
# only the runtime libtbb.so.2 is installed). Point each missing one at the
# file libdeal_II itself loads (ldd), else at a versioned file next to it.
DEALII_LIBS=$(find "$DEALII_ROOT" -name 'libdeal_II*.so' 2> /dev/null | head -2)
LDD_OUT=$(for l in $DEALII_LIBS; do ldd "$l" 2> /dev/null; done)
LDD_OUT="$LDD_OUT" python3 - <<'PYEOF'
import glob, os, re
from pathlib import Path
ldd = {}
for line in os.environ.get("LDD_OUT", "").splitlines():
    m = re.match(r"\s*(\S+)\s+=>\s+(/\S+)", line)
    if m:
        ldd[m.group(1)] = m.group(2)
def resolve(path):
    base = os.path.basename(path)
    for soname, real in sorted(ldd.items()):
        if soname.startswith(base + ".") and os.path.exists(real):
            return real
    cands = sorted(glob.glob(path + ".*"))
    return cands[0] if cands else None
d = Path("CMakeFiles/femffusion.exe.dir")
pat = re.compile(r"(/[^\s;\"']+\.so)(?=[\s;\"']|$)")
fixes = {}
for name in ("link.txt", "build.make"):
    f = d / name
    s = f.read_text()
    for p in set(pat.findall(s)):
        if os.path.exists(p) or p in fixes:
            continue
        r = resolve(p)
        if r is None:
            raise SystemExit(f"missing library {p} and no versioned copy found")
        fixes[p] = r
    for p, r in fixes.items():
        s = re.sub(re.escape(p) + r"(?=[\s;\"']|$)", r, s)
    f.write_text(s)
for p, r in fixes.items():
    print(f"link fix: {p} -> {r}")
PYEOF

make -j "${BUILD_JOBS:-16}"

# TACC's deal.II links all of Trilinos 13, including PyTrilinos, whose
# library constructor segfaults in a C++ program (femffusion.exe died in
# libpytrilinos.so.13's init, before main; diag_femffusion.sh). Nothing uses
# PyTrilinos: drop it from femffusion.exe and from a private copy of
# libdeal_II (TACC's install is untouched), found first through
# femffusion.exe's RPATH.
DROP=libpytrilinos.so.13
if ldd ./femffusion.exe | grep -q "$DROP"; then
  if ! command -v patchelf > /dev/null; then
    if [ ! -x "$DEST/.tools/bin/patchelf" ]; then
      python3 -m venv "$DEST/.tools"
      "$DEST/.tools/bin/pip" install -q patchelf
    fi
    PATH="$DEST/.tools/bin:$PATH"
  fi
  DLIB=$(ldd ./femffusion.exe | awk '/libdeal_II[^ ]*\.so/ {print $3; exit}')
  mkdir -p "$DEST/lib"
  cp -f "$DLIB" "$DEST/lib/"
  COPY="$DEST/lib/$(basename "$DLIB")"
  # the copy keeps finding its dependencies where the original does
  OLD=$(patchelf --print-rpath "$DLIB" | sed "s|\$ORIGIN|$(dirname "$DLIB")|g")
  patchelf --set-rpath "$(dirname "$DLIB")${OLD:+:$OLD}" "$COPY"
  patchelf --remove-needed "$DROP" "$COPY"
  patchelf --remove-needed "$DROP" ./femffusion.exe
  OLD=$(patchelf --print-rpath ./femffusion.exe)
  patchelf --force-rpath --set-rpath "$DEST/lib${OLD:+:$OLD}" ./femffusion.exe
  echo "dropped $DROP; libdeal_II from $COPY"
  if ldd ./femffusion.exe | grep -q "$DROP"; then
    echo "$DROP is still loaded, through:"
    for l in $(ldd ./femffusion.exe | awk '$3 ~ /^\// {print $3}'); do
      readelf -d "$l" 2> /dev/null | grep -q "NEEDED.*$DROP" && echo "  $l"
    done
    exit 1
  fi
fi

{
  echo "# modules for $DEST/femffusion.exe (written by build_femffusion.sh)"
  echo "module reset > /dev/null 2>&1 || true"
  for m in $PREREQ; do echo "module load $m"; done
  echo "module load $DEALII"
  echo "export OMP_NUM_THREADS=1"
} > "$DEST/modules.sh"

echo "FEMFFUSION built: $DEST/femffusion.exe"
echo "Running its test suite (about a minute) ..."
# Line-buffered and logged, so a crash does not swallow the output
set +e
OMP_NUM_THREADS=1 stdbuf -oL -eL ./femffusion.exe -t > "$DEST/test.log" 2>&1
rc=$?
set -e
tail -4 "$DEST/test.log"
if [ $rc -ne 0 ]; then
  echo "femffusion.exe -t exited with status $rc; full log: $DEST/test.log"
  exit $rc
fi

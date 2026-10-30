#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./debug_slump.sh [plain|gdb|blocking] [simulation arguments...]
  plain     Run normally; attach GDB only if it stalls (default).
  gdb       Launch under GDB with normal address randomization.
  blocking  Run with CUDA_LAUNCH_BLOCKING=1, without GDB.
Set DEME_PYTHON to override /home/ruochunz/miniconda3/envs/ML/bin/python.
Each invocation creates a fresh directory under debug-runs/.
The launcher owns --out-dir; other arguments go to pyDEME_DFCSlump.py.
EOF
}

mode=${1:-plain}
if [[ $# -gt 0 ]]; then shift; fi
case "$mode" in
    -h|--help) usage; exit 0 ;;
    plain|gdb|blocking) ;;
    *) usage >&2; exit 2 ;;
esac
for arg in "$@"; do
    case "$arg" in
        --out-dir|--out-dir=*|-d|-d?*)
            echo 'Output directories are assigned by this launcher.' >&2
            exit 2 ;;
    esac
done

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$repo_dir"
python_bin=${DEME_PYTHON:-/home/ruochunz/miniconda3/envs/ML/bin/python}
[[ -x "$python_bin" ]] || { echo "Python not executable: $python_bin" >&2; exit 2; }
if [[ "$mode" == gdb ]]; then command -v gdb >/dev/null; fi
if [[ "$mode" == blocking ]]; then export CUDA_LAUNCH_BLOCKING=1; fi
mkdir -p debug-runs
run_dir=$(mktemp -d "$repo_dir/debug-runs/${mode}-XXXXXXXX")
export DEME_DEBUG_RUN_DIR="$run_dir"

# Record provenance in a separate short-lived process; don't alter solver settings.
"$python_bin" - <<'PY' | tee "$run_dir/environment.log"
import importlib.metadata
import os
import sys
import DEME
print('Python:', sys.executable)
print('DEME:', DEME.__file__)
print('DEME version:', getattr(DEME, '__version__', 'unknown'))
for name in ('deme', 'pyDEME'):
    try:
        print('Distribution:', name, importlib.metadata.version(name))
    except importlib.metadata.PackageNotFoundError:
        pass
for name, module in list(sys.modules.items()):
    path = getattr(module, '__file__', '') or ''
    if path.endswith('.so') and 'deme' in (name + path).lower():
        print('Extension:', path)
for key in ('CONDA_PREFIX', 'PYTHONPATH', 'LD_LIBRARY_PATH',
            'CUDA_VISIBLE_DEVICES', 'CUDA_LAUNCH_BLOCKING'):
    print(f'{key}={os.environ.get(key)!r}')
print('Run directory:', os.environ['DEME_DEBUG_RUN_DIR'])
PY
git rev-parse HEAD >> "$run_dir/environment.log"
git status --short >> "$run_dir/environment.log"
sha256sum pyDEME_DFCSlump.py DFCModel.cu >> "$run_dir/environment.log"
cmd=("$python_bin" -u "$repo_dir/pyDEME_DFCSlump.py" "$@" --out-dir "$run_dir/output")
printf '%q ' "${cmd[@]}" >> "$run_dir/environment.log"
printf '\n' >> "$run_dir/environment.log"
echo "Launcher PID: $$ (becomes Python in plain/blocking mode)"
if [[ "$mode" == gdb ]]; then
    exec gdb -ex 'set pagination off' -ex 'set breakpoint pending on' \
        -ex 'set disable-randomization off' --args "${cmd[@]}"
fi
exec "${cmd[@]}"

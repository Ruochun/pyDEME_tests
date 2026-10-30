# Debugging the pyDEME stall

Use the same interpreter, input files, and arguments as the failing run. The
launcher defaults to the ML interpreter and prints the loaded extension and
environment. It leaves solver settings unchanged. Each invocation saves provenance
and simulation output in a new ignored `debug-runs/` directory.

## Reproduce first, then attach

```bash
cd /home/ruochunz/pyDEME_tests
conda activate ML
./debug_slump.sh plain
# Optional interpreter override and simulation arguments:
DEME_PYTHON="$CONDA_PREFIX/bin/python" ./debug_slump.sh plain --settle
```

`--settle` actually **disables** settling in this script; omit it for the default
run. Do not run multiple experiments concurrently when comparing stall behavior.
The launcher uses `python -u` for immediate Python output. It does not save console
output automatically; redirecting or piping it can change timing.

If it stalls, use the printed launcher PID (which becomes Python's PID):

```bash
gdb -p PID
```

In GDB, capture all threads before changing anything:

```gdb
set pagination off
set logging file /tmp/deme-stall.txt
set logging overwrite on
set logging enabled on
info threads
thread apply all bt full
info sharedlibrary
set logging enabled off
continue
```

Pause with Ctrl+C after a few seconds and compare another stack snapshot. A main
thread waiting on a condition variable alone is not proof of deadlock. `detach`
releases an attached process. If ptrace permissions deny attachment, launch under
GDB instead:

```bash
./debug_slump.sh gdb
```

Then enter `run`; press Ctrl+C at the stall and use the commands above. This
launcher preserves normal address randomization (`set disable-randomization off`).
Debugger pauses and startup can still change thread timing and hide races.

## VS Code

Open this repository and install Microsoft's C/C++ extension (`ms-vscode.cpptools`)
in the same Linux/WSL/SSH environment. In Run and Debug choose (Ctrl+Shift+D):

- **Slump: attach GDB to running Python** to capture an existing stall.
- **Slump: launch Python under GDB** to start a new run with F5.

Supply the reproducing interpreter when prompted (default: ML). Launch mode uses
`OUT_debug_vscode`; the simulation clears that directory on each launch. Do not
launch two sessions with that same output directory. For unique directories use
the shell launcher and attach.

Click Pause, then enter these in the Debug Console:

```text
-exec info threads
-exec thread apply all bt full
-exec info sharedlibrary
```

Prefix the other GDB commands above with `-exec` as well. Select threads and frames
in CALL STACK. A useful native function breakpoint is
`deme::DEMSolver::DoDynamicsThenSync`. Native source stepping requires matching
C++ source and debug symbols; a release wheel may only provide partial stacks.
Python breakpoints and CUDA kernel stepping require their respective debuggers.
If using a separately staged symbols-enabled build, set `PYTHONPATH` to its install
directory before using the launcher, or add it to the VS Code configuration's
`environment` array. Verify the extension path; changing builds can change the bug.

## Compare blocking CUDA launches

```bash
./debug_slump.sh blocking
```

This sets `CUDA_LAUNCH_BLOCKING=1` before Python starts; no rebuild is required.
For VS Code launch mode, add:

```json
"environment": [{ "name": "CUDA_LAUNCH_BLOCKING", "value": "1" }]
```

Plain/GDB modes inherit the environment, so use `unset CUDA_LAUNCH_BLOCKING` before
a nonblocking comparison if you previously exported it. Blocking launches can be
much slower and can hide races; success is diagnostic evidence, not a fix.
`DoDynamicsThenSync()` synchronizes simulation workers at a call boundary, not
after every kernel. Some DEME versions also expose `SetCudaDebugSync(True)`;
availability and defaults depend on the installed build. This setup does not call it.

References: [VS Code GDB configuration](https://code.visualstudio.com/docs/cpp/launch-json-reference),
[GDB address randomization](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Starting.html),
[CUDA launch blocking](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/environment-variables.html).

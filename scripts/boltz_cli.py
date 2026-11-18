"""Local launcher for the `boltz` CLI, used instead of the installed shim.

Why this exists: uv writes console-script shims as .exe trampolines with the
absolute interpreter path baked in, so moving the working folder breaks every
invocation with `error: uv trampoline failed to canonicalize script path`.
`python -m boltz` is not an escape hatch either -- the wheel ships no
__main__.py, so it answers "boltz is a package and cannot be directly executed".

The installed entry point (boltz-2.0.3.dist-info/entry_points.txt) is
`boltz = boltz.main:cli`, which is exactly what that shim wraps, so calling it
through any python that can import boltz reproduces the shim without a
hardcoded path.
"""

import sys

from boltz.main import cli

if __name__ == "__main__":
    cli()

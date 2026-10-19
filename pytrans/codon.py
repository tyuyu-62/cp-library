from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .translator import check_source


def emit_codon(source: str, filename: str = "main.py") -> str:
    check_source(source, filename, target="codon")
    return source


def build_codon(source_path: Path, output: Path, codon: str = "codon") -> None:
    executable = shutil.which(codon)
    if executable is None:
        raise RuntimeError(f"Codon executable was not found: {codon}")
    command = [executable, "build", "-exe", "-release", "-o", str(output), str(source_path)]
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        raise RuntimeError(f"Codon build failed with exit code {completed.returncode}")

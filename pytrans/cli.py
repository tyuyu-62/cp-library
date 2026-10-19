from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from . import __version__
from .codon import build_codon, emit_codon
from .native import transpile_native_cpp
from .translator import TranspileError, check_source, parse_source, transpile_cpp


def read_source(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        raise RuntimeError(f"cannot read {path}: {error}") from error


def write_output(path: Path, text: str) -> None:
    try:
        path.write_text(text, encoding="utf-8")
    except OSError as error:
        raise RuntimeError(f"cannot write {path}: {error}") from error


def default_executable(source: Path) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return source.with_suffix(suffix)


def compile_cpp(source: str, source_path: Path, output: Path, compiler: str, opt: str, target: str = "cpp", int_type: str = "Int") -> None:
    executable = shutil.which(compiler) or compiler
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pytrans-") as temporary:
        cpp_path = Path(temporary) / (source_path.stem + ".cpp")
        translator = transpile_native_cpp if target == "native" else transpile_cpp
        generated = translator(source, str(source_path), int_type=int_type) if target == "native" else translator(source, str(source_path))
        cpp_path.write_text(generated, encoding="utf-8")
        command = [executable, "-std=c++20", opt, "-pipe", str(cpp_path), "-o", str(output)]
        if '#include <boost/' in generated and Path('/opt/homebrew/include/boost').is_dir():
            command += ['-isystem','/opt/homebrew/include']
        completed = subprocess.run(command, check=False)
        if completed.returncode:
            raise RuntimeError(f"C++ compilation failed with exit code {completed.returncode}")


def run_program(command: list[str], input_data: bytes, timeout: float) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(command, input=input_data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, check=False)


def build_for_diff(source: str, source_path: Path, target: str, compiler: str, codon: str, temporary: Path, int_type: str = "Int") -> list[str]:
    if target in {"cpp", "native"}:
        executable = temporary / "a.out"
        compile_cpp(source, source_path, executable, compiler, "-O2", target, int_type)
        return [str(executable)]
    generated = temporary / source_path.name
    write_output(generated, emit_codon(source, str(source_path)))
    executable = temporary / "a.out"
    build_codon(generated, executable, codon)
    return [str(executable)]


def command_check(args: argparse.Namespace) -> int:
    source_path = Path(args.source)
    source = read_source(source_path)
    if args.target == "native":
        transpile_native_cpp(source, str(source_path), int_type=args.native_int)
    else:
        check_source(source, str(source_path), target=args.target)
    print(f"OK: {source_path} ({args.target}, no annotations required)")
    return 0


def command_emit(args: argparse.Namespace) -> int:
    source_path = Path(args.source)
    source = read_source(source_path)
    if args.target == "native":
        output = transpile_native_cpp(source, str(source_path), int_type=args.native_int)
    elif args.target == "cpp":
        output = transpile_cpp(source, str(source_path))
    else:
        output = emit_codon(source, str(source_path))
    if args.output:
        write_output(Path(args.output), output)
    else:
        sys.stdout.write(output)
    return 0


def command_build(args: argparse.Namespace) -> int:
    source_path = Path(args.source)
    source = read_source(source_path)
    output = Path(args.output) if args.output else default_executable(source_path)
    if args.target in {"cpp", "native"}:
        compile_cpp(source, source_path, output, args.compiler, args.opt, args.target, args.native_int)
    else:
        generated = source_path
        if args.emit_source:
            generated = output.with_suffix(".py")
            write_output(generated, emit_codon(source, str(source_path)))
        else:
            check_source(source, str(source_path), target="codon")
        build_codon(generated, output, args.codon)
    print(output)
    return 0


def command_diff(args: argparse.Namespace) -> int:
    source_path = Path(args.source)
    input_path = Path(args.input)
    source = read_source(source_path)
    input_data = input_path.read_bytes()
    with tempfile.TemporaryDirectory(prefix="pytrans-diff-") as directory:
        target_command = build_for_diff(source, source_path, args.target, args.compiler, args.codon, Path(directory), args.native_int)
        python_command = [args.python, str(source_path)]
        python_result = run_program(python_command, input_data, args.timeout)
        target_result = run_program(target_command, input_data, args.timeout)
    if python_result.returncode != 0:
        sys.stderr.buffer.write(python_result.stderr)
        raise RuntimeError(f"source program failed with exit code {python_result.returncode}")
    if target_result.returncode != 0:
        sys.stderr.buffer.write(target_result.stderr)
        raise RuntimeError(f"translated program failed with exit code {target_result.returncode}")
    left = python_result.stdout.split() if args.tokens else python_result.stdout
    right = target_result.stdout.split() if args.tokens else target_result.stdout
    if left != right:
        print("mismatch", file=sys.stderr)
        print("--- source ---", file=sys.stderr)
        sys.stderr.buffer.write(python_result.stdout)
        print("\n--- target ---", file=sys.stderr)
        sys.stderr.buffer.write(target_result.stdout)
        return 1
    print("OK: source and translated output match")
    return 0


def command_explain(args: argparse.Namespace) -> int:
    source_path = Path(args.source)
    source = read_source(source_path)
    tree, analyzer = parse_source(source, str(source_path))
    counts: dict[str, int] = {}
    for node in tree.body:
        key = type(node).__name__
        counts[key] = counts.get(key, 0) + 1
    print(f"file: {source_path}")
    print(f"functions: {', '.join(analyzer.functions) or '(none)'}")
    print(f"imports: {', '.join(sorted(analyzer.imports)) or '(none)'}")
    print("module statements:")
    for key, count in sorted(counts.items()):
        print(f"  {key}: {count}")
    print("inferred signatures:")
    for name, node in analyzer.functions.items():
        params, return_type = analyzer.signature(node)
        values = ", ".join(f"{key}: {value.kind}" for key, value in params.items())
        print(f"  {name}({values}) -> {return_type.kind}")
    return 0


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pytrans", description="Deterministic competitive-programming Python transpiler")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check", help="check a source file")
    check.add_argument("source")
    check.add_argument("--target", choices=("cpp", "native", "codon"), default="cpp")
    check.set_defaults(function=command_check)

    emit = subparsers.add_parser("emit", help="emit translated source")
    emit.add_argument("source")
    emit.add_argument("--target", choices=("cpp", "native", "codon"), default="cpp")
    emit.add_argument("-o", "--output")
    emit.set_defaults(function=command_emit)

    build = subparsers.add_parser("build", help="translate and compile")
    build.add_argument("source")
    build.add_argument("--target", choices=("cpp", "native", "codon"), default="cpp")
    build.add_argument("-o", "--output")
    build.add_argument("--compiler", default="c++")
    build.add_argument("--opt", default="-O2")
    build.add_argument("--codon", default="codon")
    build.add_argument("--emit-source", action="store_true")
    build.set_defaults(function=command_build)

    diff = subparsers.add_parser("diff", help="compare source and translated output")
    diff.add_argument("source")
    diff.add_argument("--input", required=True)
    diff.add_argument("--target", choices=("cpp", "native", "codon"), default="cpp")
    diff.add_argument("--python", default=sys.executable)
    diff.add_argument("--compiler", default="c++")
    diff.add_argument("--codon", default="codon")
    diff.add_argument("--timeout", type=float, default=5.0)
    diff.add_argument("--tokens", action="store_true")
    diff.set_defaults(function=command_diff)

    for command in (check, emit, build, diff):
        command.add_argument("--native-int", choices=("Int", "int64", "int", "long long"), default="Int", help="native integer mode: adaptive Int, full fixed int64, or legacy scalar-only int/long long")

    explain = subparsers.add_parser("explain", help="show inferred types and source summary")
    explain.add_argument("source")
    explain.set_defaults(function=command_explain)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        return args.function(args)
    except (OSError, RuntimeError, TranspileError, ValueError) as error:
        print(f"pytrans: {error}", file=sys.stderr)
        return 2

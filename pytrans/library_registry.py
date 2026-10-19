from __future__ import annotations

import ast
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from .translator import TranspileError


def definition_fingerprint(source):
    tree = ast.parse(source) if isinstance(source, str) else source
    node = tree.body[0] if isinstance(tree, ast.Module) and len(tree.body) == 1 else tree
    return hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()


def bundle_cpp(path, seen=None):
    path = Path(path).resolve()
    seen = set() if seen is None else seen
    if path in seen:
        return ''
    seen.add(path)
    lines = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped == '#pragma once':
            continue
        if stripped.startswith('#include "'):
            relative = stripped.split('"')[1]
            child = path.parent / relative
            if not child.is_file():
                raise ValueError(f'unresolved C++ library include: {child}')
            lines.append(bundle_cpp(child, seen))
        else:
            lines.append(line)
    return '\n'.join(lines) + '\n'


@dataclass(frozen=True)
class NativeLibrary:
    kind: str
    name: str
    python_sha256: str
    cpp_source: str


class LibraryRegistry:
    def __init__(self, libraries=()):
        self.libraries = tuple(libraries)
        keys = [(entry.kind, entry.name, entry.python_sha256) for entry in self.libraries]
        if len(set(keys)) != len(keys):
            raise ValueError('duplicate native library registration')
        for entry in self.libraries:
            if entry.kind not in {'class','function'} or not entry.name.isidentifier() or len(entry.python_sha256) != 64:
                raise ValueError('invalid native library registration')

    @classmethod
    def load(cls, path=None):
        if path is None:
            path = Path(__file__).with_name('native_libraries.json')
            if not path.is_file():
                return cls()
        path = Path(path)
        try:
            data = json.loads(path.read_text())
            if data.get('version') != 1:
                raise ValueError('unsupported native library registry version')
            entries = [NativeLibrary(item['kind'],item['name'],item['python_sha256'],bundle_cpp(path.parent / item['cpp_path'])) for item in data['libraries']]
            return cls(entries)
        except (OSError,ValueError,KeyError,TypeError) as error:
            raise TranspileError(f'cannot load native libraries: {error}',filename=str(path)) from error

    def match(self, tree):
        definitions = {}
        for node in tree.body:
            if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
                kind = 'class' if isinstance(node,ast.ClassDef) else 'function'
                definitions[(kind,node.name)] = definition_fingerprint(node)
        return {(entry.kind,entry.name):entry for entry in self.libraries if definitions.get((entry.kind,entry.name)) == entry.python_sha256}

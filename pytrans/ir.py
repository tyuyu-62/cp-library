from __future__ import annotations

import ast
from dataclasses import dataclass
from enum import Enum
from typing import Any, Mapping


class Ownership(Enum):
    VALUE = "value"
    BORROWED = "borrowed"
    MUTABLE_BORROW = "mutable_borrow"
    COPIED = "copied"


class Storage(Enum):
    DIRECT = "direct"
    LOCAL_LIST = "local_list"
    CALLABLE_ALIAS = "callable_alias"


@dataclass(frozen=True)
class LocalIR:
    name: str
    value_type: Any
    ownership: Ownership
    storage: Storage
    declared: bool


@dataclass(frozen=True)
class FunctionIR:
    name: str
    node: ast.FunctionDef
    parameter_types: Mapping[str, Any]
    parameter_ownership: Mapping[str, Ownership]
    mutable_parameters: frozenset[str]
    return_type: Any
    environment: Mapping[str, Any]
    value_types: Mapping[str, Any]
    locals: Mapping[str, LocalIR]
    copied_parameters: frozenset[str]
    rebound_types: Mapping[str, Any]
    rebound_names: frozenset[str]
    rebound_targets: Mapping[str, tuple[Any, ...]]
    reference_aliases: frozenset[str]
    native_callable_aliases: frozenset[str]
    private_lists: frozenset[str]
    local_lists: frozenset[str]
    generator: bool


@dataclass(frozen=True)
class ProgramIR:
    source: str
    filename: str
    tree: ast.Module
    functions: Mapping[str, FunctionIR]
    methods: Mapping[tuple[str, str], FunctionIR]

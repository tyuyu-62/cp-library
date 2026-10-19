#!/bin/zsh

SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
/usr/bin/env python3 "$SCRIPT_DIR/run.py"
translation_status=$?

echo
read -k 1 "?終了するにはキーを押してください: "
exit "$translation_status"

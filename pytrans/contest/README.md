# 変換サンプル

`input.py` は、要素数 N と N 個の整数を読み、その合計を出力します。

入力例:

```text
5
1 -2 3 4 5
```

出力例:

```text
11
```

このフォルダで次のコマンドを実行すると、`output.cpp` を再生成します。

```sh
python3 -B run.py
```

macOSでは `run.command` からも実行できます。変換にはPython 3.10以降、生成C++のコンパイルにはC++20対応コンパイラとBoostヘッダが必要です。

```sh
pypy3 -B make_testcase.py
```

このコマンドはseed固定の小さな入力100件でPython版と生成済みC++版を比較します。PyPyと `c++` コマンドが必要です。HomebrewのBoostヘッダは自動検出します。他の場所にあるBoostはコンパイラから参照できるよう設定してください。

自分の問題に使う場合は `input.py` と `make_testcase.py` の入力生成関数を編集し、`run.py` で再変換してください。不一致時は入力を `failed_case.txt` に保存します。`output.cpp` には変換元のコードもコメントとして含まれます。

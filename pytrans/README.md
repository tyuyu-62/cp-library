# pytrans

競技プログラミング用のPythonサブセットをC++へ変換する、決定的なルールベースのツールです。変換処理はPython標準ライブラリだけで動作し、AIやネットワークを使用しません。Python全体の互換性を提供するものではありません。

## 必要な環境

- Python 3.10以降、または対応するPyPy。
- 生成C++のコンパイルには、C++20とGCC/Clang系の整数拡張（`__int128` など）に対応するコンパイラ。
- 既定のnative自動多倍長整数方式には、Boost Multiprecisionのヘッダ。

CLIのコンパイラ既定値は `c++` です。`--compiler` で変更できます。Boostをコンパイラから参照できるように設定してください。Homebrewの `/opt/homebrew/include/boost` が存在する場合、CLIはincludeパスを自動追加します。

## 使い方

`pytrans/` の親ディレクトリ（公開リポジトリのルート）から実行します。次の例では、同じ場所に `main.py` と比較用の入力 `sample.in` を用意します。

```sh
python3 -B -m pytrans check main.py --target native
python3 -B -m pytrans emit main.py --target native -o main.cpp
python3 -B -m pytrans build main.py --target native -o main
python3 -B -m pytrans diff main.py --target native --input sample.in --python pypy3
python3 -B -m pytrans explain main.py
```

`emit` はC++ソースを出力し、`build` は実行ファイルまで作成します。`diff` は指定した入力でPython版と生成版の出力を比較します。既定は出力のバイト一致、`--tokens` 指定時は空白区切りの比較です。`explain` は型解析の概要を表示します。

型注釈は任意です。リテラル、組み込み関数、コンテナ操作、関数呼び出し元から型を推論します。対応外の構文や解決できない型は変換エラーになります。変換に成功しても、対象の入力・値域で元のPythonと同じ結果になることを確認してください。

## 変換方式と整数

CLIで `--target` を省略すると従来方式の `cpp` が選ばれます。上の例では明示的に `native` を指定しています。

| 指定 | 整数の扱い |
| --- | --- |
| `--target native` | 通常は64bitで保持し、必要な演算でBoostの多倍長整数へ昇格する `Int` 方式。 |
| `--target native --native-int int64` | 固定64bit方式。Boost不要。通常の整数演算は中間結果も64bitに収まることが前提。 |
| `--target cpp` | 従来の変換方式。整数は符号付き64bit。 |
| `--target native --native-int int` / `"long long"` | 基本構文だけを扱う旧ネイティブ方式。 |
| `--target codon` | 対応範囲を検査してCodon用ソースを出力。ビルドには別途Codon環境が必要。 |

整数の契約、旧方式の制限、C++ライブラリ置換については [NATIVE.md](NATIVE.md) を参照してください。

## 対応範囲と制限

現在の主な対応範囲は、整数・浮動小数点数・文字列・リスト・tuple・deque・set・挿入順を保つdict、`defaultdict(int/list/set/dict/str/bool/float)`と引数なしlambda factory、`None`を含むOptional値、関数、再帰、デフォルト引数、`if`、`for`、`while`、ループの`else`、`range`、`enumerate`、`zip`、1段のlist/set/dict内包表記とジェネレータ式です。継承のない通常クラスについて、`__init__`、インスタンス属性、メソッド、`@staticmethod`、`__len__`、`__getitem__`、`__setitem__`、`__contains__`、`__str__`を変換できます。負の添字とstep付きslice・slice代入、`del`、固定長のsequence unpack、`print(*values)`、`reversed`、`all`、`any`、`pow`の剰余形式、`divmod`、`bit_length`、`bin`、`join`、`bisect`、`heapq`、`math`にも対応します。

局所関数は外側の値を引数に変換して扱います。単純なf-string・`str.format`の`{}`、`sorted`と`list.sort`の単一引数lambdaキー、`range(...)[::-1]`などのrangeスライス、tupleの反復、`bytes`の整数要素、整数用`array`の`i/I/H`型コード、`bytearray`のゼロ初期化と整数配列操作、`randint`・`randrange`・63bitまでの`getrandbits`、`math`の基本的な三角関数にも対応します。整数用`array`は生成C++で64bit整数のvectorとして扱うため、元の型コードの値域内で動作するコードを対象にします。動的な`tuple(sequence)`は添字参照用の配列として変換するため、tupleの不変性を利用するコードは対象外です。

競技プログラミング用テンプレートでよく使われる、`input = lambda: sys.stdin.readline().rstrip()`と、文字列バッファへ追加する独自`print`関数も認識します。この形式では生成C++から標準出力へ直接書き出し、末尾のバッファflushは省略します。文字列グリッドの転置で生じる`list[str]`と`list[list[char]]`の変換にも対応します。

継承、`@staticmethod`以外のクラス／メソッドデコレータ、`yield`を使うイテレータ、例外処理、`with`、一般のlambda、複数段の内包表記、ユーザー定義関数へのキーワード引数は未対応です。リスト引数の再束縛は、元のリストを変更しないと判定できる場合に限りコピーして扱います。f-stringの書式指定にも未対応です。従来の `--target cpp` ではPythonの整数が符号付き64bitになるため、その範囲に収まるコードを対象にします。新しいnative方式のIntについてはNATIVE.mdを参照してください。

`--target codon`は、対応範囲を検査した入力をCodon用ソースとして出力します。Codonのビルド環境がある場合は`build`で実行ファイルまで作成できます。

## 同梱サンプル

[contest/](contest/README.md) には整数配列の合計を求める小さなサンプルがあります。リポジトリのルートから次を実行します。

```sh
python3 -B pytrans/contest/run.py
pypy3 -B pytrans/contest/make_testcase.py
```

`run.py` は `contest/input.py` を読み、`contest/output.cpp` を上書きします。生成C++の先頭には変換元のコードをコメントとして含めます。比較スクリプトは、生成済みC++とPythonの出力をseed固定の100入力で確認します。詳細は[サンプルの説明](contest/README.md)を参照してください。

## Python API

リポジトリのルートからimportできます。

```python
from pytrans import transpile_native_cpp

source = "a, b = map(int, input().split())\nprint(a + b)\n"
cpp = transpile_native_cpp(source, "main.py")
fixed_cpp = transpile_native_cpp(source, "main.py", int_type="int64")
```

APIは生成C++を文字列で返します。ファイルへの保存やコンパイルは呼び出し側で行います。

## ライセンス

[CC0 1.0 Universal](LICENSE)。

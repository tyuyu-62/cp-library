# Competitive Programming Toolkit

競技プログラミング用のPyPyライブラリと、Pythonの対応サブセットをC++へ変換するツールです。

## 構成

| ディレクトリ | 内容 |
| --- | --- |
| [pypy/library/](pypy/library/) | PyPy向けのアルゴリズム・データ構造。必要なファイルを提出コードへコピーして利用します。 |
| [pytrans/](pytrans/) | PythonからC++へのルールベース変換ツール。使い方と対応範囲は[pytransのREADME](pytrans/README.md)を参照してください。 |

PyPyライブラリのAPI・入力の制約・計算量は、各ファイルのdocstringを参照してください。PyPy固有の機能を使う実装もあります。

## C++への変換

このREADMEがあるディレクトリで実行します。

```sh
python3 -B -m pytrans emit main.py --target native -o main.cpp
```

変換処理はPython標準ライブラリだけを使い、AIやネットワークを使用しません。生成C++はC++20を使用し、既定の自動多倍長整数方式にはBoostヘッダが必要です。固定64bit方式と未対応のPython構文については[ネイティブC++生成の説明](pytrans/NATIVE.md)を参照してください。

小さな入出力例と再生成・比較の手順は[変換サンプル](pytrans/contest/README.md)にあります。

## ライセンス

`pypy/` と `pytrans/` はともにCC0 1.0 Universalで公開します。

- [PyPyライブラリのライセンス](pypy/LICENSE)
- [pytransのライセンス](pytrans/LICENSE)

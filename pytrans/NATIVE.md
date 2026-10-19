# ネイティブC++生成

`--target native` / `transpile_native_cpp()` の既定は新しい `Int` 方式です。`contest/run.py` もこの方式を使用します。既存の型解析・IR・構文生成を再利用し、関数・クラス・コンテナ・内包表記・再帰・ループ・スライスなど従来の `cpp` 対応範囲を土台にしています。`--target cpp` は従来の方式のまま残しています。

以下のCLI例は `pytrans/` の親ディレクトリから実行します。

```sh
pypy3 -B -m pytrans emit main.py --target native -o main.cpp
pypy3 -B -m pytrans build main.py --target native --compiler c++ -o main
pypy3 -B -m pytrans diff main.py --target native --python pypy3 --input sample.in
```

```python
from pytrans import transpile_native_cpp
cpp = transpile_native_cpp(source, "main.py")
```

Python側は標準ライブラリだけで動作し、AIやネットワークを使用しません。生成C++はC++20を使用し、自動整数方式はBoost Multiprecisionのヘッダを使用します。固定64bit方式はBoost不要です。CLIは `/opt/homebrew/include/boost` が存在する場合にコンパイルのincludeパスを追加します。他の環境で自動整数方式を使う場合はBoostをコンパイラから参照できるように設定してください。

## 整数方式の選択

```python
from pytrans import transpile_native_cpp, transpile_native_full_cpp

cpp = transpile_native_full_cpp(source, "main.py", integer_mode="fixed64")
cpp = transpile_native_cpp(source, "main.py", int_type="int64")
```

```sh
pypy3 -B -m pytrans emit main.py --target native --native-int int64 -o main.cpp
```

固定64bitはnative_fixed64.pyで現在のネイティブ生成・IR・構文対応・元の名前・starred代入・helperの削減を共有します。基本構文だけの旧NativeEmitterとは別です。旧`--native-int int`/`long long`は互換のため保持します。Pythonの整数変数、関数引数、クラス属性、コンテナ要素には通常のlong long（i64 alias）を用い、IntとBoostへの自動昇格を生成しません。剰余式・回数計算には安全な計算のため一時的に128bitを使うことがあります。中間結果を含めて64bitを超える通常の演算は保証範囲外で、符号付きoverflowにより結果が不正・未定義となる可能性があります。範囲外の整数リテラルは変換時エラーです。負の床除算・剰余、整数のbit_length/bit_count、正確なisqrt、64bitに収まるbytes変換は扱います。最小値/-1の除算は範囲外エラーです。

登録済みC++ライブラリがIntを要求する場合、固定64bitではその置換を使わずPython定義から翻訳します。固定64bit互換の登録済み実装は引き続き使えます。浮動小数点の出力形式など、既存の対応制限は引き継ぎます。固定64bitで生成したC++の速度やACは、入力・値域を満たすケースで確認してください。

## 生成方針

- Pythonの変数名を保ちます。C++予約語・既存の標準名に衝突する名前は既存のエスケープ規則を使います。型の再束縛による別の記憶領域は `a_0` などにし、既存の同名変数との衝突を避けます。文字列リテラル内の名前は書き換えません。
- 非負の定数、非負を証明できる式、本文で書き換えられない非負rangeのループ変数による添字は `a.at(static_cast<size_t>(i))` へ直接生成します。`.at` により範囲外検査を保持します。符号不明の添字は負数補正を含む `py_at` を残します。
- 配列は `vector` / `deque`、セットは `set`、通常の追加・ソートなどは標準C++を使います。局所配列の独自ラッパーは新方式では使いません。整数の加減乗算・除算・剰余・比較・bit演算は `Int` の演算子へ生成します。2引数の整数min/maxは `std::min/max` にします。
- 必要なランタイム定義と依存先だけを出力します。使われないpow・slice・乱数などの `py_` 関数は含めません。型を使わない一部のコンテナoverloadも省きます。同じ名前の汎用overloadや必要なクラス内のメソッドはまとめて残る場合があります。
- list/dequeを右辺とする `t, *q = map(int, input().split())`、先頭／中間／末尾starred、空の残り要素、関数内・モジュール直下の代入に対応します。右辺を1回評価し、要素不足を検査して残りをvectorへコピーします。

## IntとBigInt

`int_runtime.py` が整数実装の主な交換箇所です。bytes変換のBoost接続はnative_full.pyにもあります。`Int` は通常 `long long` と空の所有ポインタを保持し、多倍長値が必要になった時だけ `boost::multiprecision::cpp_int` を確保します。コピーした多倍長値は独立した所有を持ち、共有した値を書き換えません。結果が64bitへ戻れば通常表現へ戻ります。

- 加減乗算は `__builtin_*_overflow` で判定し、失敗時に元のオペランドから多倍長で計算し直します。桁あふれしたC++の符号付き演算結果を使用しません。
- 最小整数の符号反転・絶対値・`// -1`、大きい入力・リテラル・累乗・左右shift・符号付きbit演算・負数の床除算と剰余を扱います。0除算・負のshiftは実行時エラーです。
- int64を超える入力も整数として読み込みます。辞書キー用hashと比較、十進出力、bit_length/bit_countに対応します。mod累乗やgcd/lcm、正確なisqrt、大きい整数のbytes変換もIntで動作します。
- 配列サイズ・添字や外部C++ライブラリの固定幅整数型への変換は範囲を確認し、収まらなければ例外にします。多倍長を暗黙に切り捨てません。
- 更新されない64bitのグローバル整数定数はconstexpr long longで保持します。値域を証明できる整数変数と式は下記の解析によりlong longへ特殊化し、それ以外はIntに残します。連続配列全体の昇格、自作BigIntへの交換は今後の範囲です。Intの検査と記憶量には費用があり、固定幅の従来方式より常に速いことは保証しません。

## 安全な64bit整数への特殊化

`integer_bounds.py` がスコープ内の全代入を集め、値の範囲を保守的に計算します。64bitの定数、組み込みコンテナのlen、範囲が分かるrange/enumerate、範囲を求められる整数演算が対象です。関数・メソッドのローカル変数とモジュール変数に適用し、引数・クラスの属性・コンテナ要素の整数型はIntを維持します。Pythonのintという型注釈や変数名が添字らしいことだけでは特殊化しません。64bit環境上のコンテナ長は符号付き64bitの範囲内として扱います。

入力による範囲不明の代入が1つでもあれば、その変数はIntに残します。更新によって範囲が広がり続ける累積値もIntです。global/nonlocal書き込み、型の再束縛、ループ内で巨大値に変更する変数を固定幅と決めつけません。整数式の各演算も範囲を確認し、結果を小さな剰余にする場合でも中間積が64bitを超えうるならIntで計算します。負数の床除算・剰余、最小値の符号反転・絶対値、64以上のshiftなどは必要に応じてIntに残します。

```python
def solve():
    a = list(map(int, input().split()))
    n = len(a)
    for i in range(n):
        j = i // 2
        print(a[j])
```

この例ではn/i/jとrangeの内部カウンタはlong long、aの要素はIntです。`range(int(input()))`の上限は不明なのでIntを維持します。rangeの端点・stepが64bitに収まると分かる場合は、ターゲット変数がIntでも内部カウンタを固定幅にできます。stepが±1以外の場合、回数を最初に128bitで求め、最後の加算を省略してオーバーフローを防ぎます。ループ本体のbreak/continueとターゲット変数の再代入に影響しません。enumerateの配列内offsetはsize_t、開始値と長さから安全性を証明できる場合はindexもlong longです。非負の固定幅添字はatを使い、符号不明の固定幅添字にはIntを作らないpy_at overloadを使います。

BigIntの演算実装とIntの保持・昇格・コピーの処理は別です。自作BigIntにより多倍長経路を改善できる可能性がありますが、64bitの表現判定やoverflow検査の費用は自動ではなくなりません。現在の外側のunique_ptrによる確保や値コピーも交換・測定の対象です。今回の変更ではBoostとIntの保持方法は変更しません。

## 事前に整備したC++ライブラリへの置換

`native_libraries.json` は初期状態で空です。PythonライブラリとC++ライブラリの対応は自動的に推測しません。コンテスト前にAIなどで最適化したC++を準備し、同じ公開API・所有権・整数の契約を満たすadapterを登録してください。整備済みの対応がなければ、普通のルールベース翻訳を行います。

Pythonの1つのclass/function定義を `definition_fingerprint()` へ渡すと、位置・通常のコメントを除いたASTのSHA-256が得られます。名前と指紋の両方が一致する定義だけC++へ置換します。Pythonを編集して指紋が変われば自動的に通常変換へ戻ります。C++本文自体は翻訳器で書き換えません。

```json
{
  "version": 1,
  "libraries": [
    {
      "kind": "class",
      "name": "Box",
      "python_sha256": "登録するPython定義の64桁の指紋",
      "cpp_path": "native-libraries/Box.hpp"
    }
  ]
}
```

`kind` はclass/function、パスはregistryからの相対パスです。引用符付きのローカルincludeを再帰的に展開し、未解決includeはエラーにします。システムincludeを残した単一のC++を生成するため、提出先に個人のヘッダは必要ありません。登録済みAPIの意味や速度の正しさは登録前に検証してください。固定幅で計算するC++ライブラリの内部overflowはIntの呼び出し側だけでは防げないため、内部もInt対応するか、制約を保証したadapterが必要です。

APIでは `transpile_native_cpp(..., registry=path_or_registry)` を指定できます。指定なしでは同じpytransディレクトリのJSONを読みます。

## 互換性と検証範囲

動的なPython全体の互換層ではありません。既存の型推論・alias・Unicode・例外処理などの制約は引き継ぎます。浮動小数点の出力表現や一部の標準ライブラリAPIはPython完全互換を保証しません。負の整数指数が整数結果として推論されたケースは実行時エラーにし、無言で1を返しません。

同梱サンプルの再生成とPython／C++の出力比較は[サンプルの説明](contest/README.md)を参照してください。

---

# 固定幅整数による旧ネイティブ方式

`--native-int int` / `"long long"` は、`native.py` の独立した AST 変換規則を使用します。共通の `TranspileError` を診断に使用します。生成する C++20 は標準ライブラリだけを使い、Python 互換ランタイム・`py_` 関数・動的な値のラッパーを含みません。既定の変換方式は従来の `cpp` のままです。

## 使用方法

`pytrans` の親ディレクトリから実行します。

```sh
python3 -m pytrans check main.py --target native --native-int "long long"
python3 -m pytrans emit main.py --target native --native-int int -o main.cpp
python3 -m pytrans build main.py --target native --native-int "long long" -o main
python3 -m pytrans diff main.py --target native --native-int "long long" --input sample.in --python pypy3
```

```python
from pytrans import transpile_native_cpp

cpp = transpile_native_cpp(source, "main.py", int_type="int")
```

`int_type` は `"int"` または `"long long"` で、この方式を選ぶには明示指定が必要です。`transpile_native_cpp()` と `--target native` の整数方式の既定は `Int` です。

## 変換例

```python
N = int(input())
for i in range(N):
    print(i + 1)
```

`--native-int int` で次の C++ を生成します。

```cpp
#include <iostream>
#include <limits>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int N;
    cin >> N;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    for (int i = 0; i < N; ++i) {
        cout << (i + 1) << "\n";
    }
    return 0;
}
```

入力行の残りを消費する `ignore` により、後続の `input()` を `getline` に変換できます。改行は `"\n"`、明示的な `print(..., flush=True)` は `cout.flush()` です。

## 現在の対応範囲

| Python | C++ |
| --- | --- |
| 整数・bool・文字列、代入・複合代入 | `int` / `long long`、`bool`、`string` |
| `if/elif/else`、比較、`and/or/not` | `if/else`、比較・短絡評価 |
| `for ... in range(stop)` | 通常のカウンタ for 文 |
| `range(start, stop, step)` | 正・負・変数 step のカウンタ for 文 |
| `break`、`continue`、for の `else` | `break`、`continue`、必要時だけ終了判定変数 |
| `n = int(input())` | `cin >> n` と行末の消費 |
| `s = input()` | `getline(cin, s)` |
| `a, b = map(int, input().split())` | `cin >> a >> b` と行末の消費 |
| `print(...)`、`sep/end/flush` | `cout << ...` と必要時だけ flush |
| 固定個数の代入、`a, b = b, a` | 型付き変数と必要時だけ一時変数 |
| `+ - * // % & \| ^`、単項 `+ - ~` | ネイティブ演算。`// %` は符号の違いを補正 |
| 条件式、`bool`、`abs`、文字列の `len` と連結 | 条件式、標準演算・標準ライブラリ |

基本の範囲だけを実装しています。関数・クラス・配列や辞書などのコンテナ・内包表記・import・while・浮動小数点数・例外処理・一般の関数呼び出しは未対応で、行番号付きのエラーにします。既存方式への自動フォールバックは行いません。

## 型と実行上の制約

- 変数の型はリテラル・入力・演算から推定します。同じ変数を異なる型へ変更するコードはエラーです。
- 未初期化の可能性がある変数の参照は拒否します。分岐の両側で代入された変数は分岐後にも使用できます。動的な range が必ず非空であることなど、一般的な条件の証明は行いません。
- range の引数は Python 同様、ループ開始時の値を使用します。ループ内で変更される境界式は事前に保存します。
- ループ変数を本文で書き換える場合やループ後にも参照する場合は、内部カウンタと Python の変数を分けます。基本の例では通常の `for (int i = ...)` を生成します。
- step が 0 の定数は変換エラー、実行時に 0 となった step は実行エラーです。step による終了時のカウンタのオーバーフローは回避します。
- 整数入力は選択した符号付き型に収まる ASCII の十進数を対象とします。1 回の `int(input())` は 1 行 1 個、unpack は 1 行に対象と同数の整数がある正常なコンテスト入力を前提にします。空行からの整数読み込み、EOF、余分な token などの Python の例外動作は再現しません。
- 文字列は UTF-8 として保持します。文字列の `len` は ASCII の範囲が対象です。Python の Unicode 文字数全般は再現しません。
- 全ての整数演算の途中結果が選択した型に収まるコードが対象です。Python の任意精度整数、最小整数に対する `abs` や `// -1` は保証しません。
- 負数の `//` と `%`、整数・文字列を返す `and/or` は必要な式だけ局所 lambda で意味を保ちます。通常の加減乗算・分岐・for・入出力は直接出力します。

## 固定64bitの追加最適化

固定個数の整数入力をstd::arrayで取得し、全constructorで同じ関数を使う変更不可のcallbackを直接呼び出す。中間演算まで64bitへ収まると証明した剰余だけ通常の%へする。符号不明の除算・剰余と未保証の広い積は従来の互換処理を維持する。

数値vectorの線形添字を使うstep=±1の単純なrangeループは、開始時に全アクセスの範囲を確認する。成功時はvector[]を使い、失敗時は負の添字と境界検査を保つ。配列の長さや境界を変更する処理、未知の関数呼び出し等には適用しない。変数名を保ち、必要なguardと通常のC++ループを生成する。自動Int方式へはこの最適化を適用しない。


## LIとstarred入力の特殊化

固定64bitでは、LIのような単純な整数行入力関数を本体のASTから認識し、固定個数分解をstack上のarrayへ直接読む。int(x)の変換と引数/整数定数の加減算を対象にし、副作用・再束縛・自由変数・filter・登録済みC++置換がある場合は従来処理を使う。

t,*q=LI()はprefix/suffixとqへ直接読み、元の全入力配列からのコピーを省く。ローカルqが他へ渡らず固定個数の整数分解にだけ使われる場合に限り容量を再利用する。毎回新規の場合も行の長さから容量を先に確保する。a,b,c=qは要素数を確認して分解し、不要なリスト全体コピーを作らない。別名・closure・例外観察などがある場合は容量再利用を使わない。一般のlistの別名/再代入に関する既存の互換性限界は残る。

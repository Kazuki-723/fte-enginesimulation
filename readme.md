# FTE enginesimulator 2.0.0.
過去に使っていたExcelデータのフローの改良及び，データ出力をおこなってくれるソフトです．

現在の最新版
- Stable 2.0.0. 

Devは他ブランチを参照．

現在，GUI版はFletのバージョン対応作業中．2.0.0.は動きません．

現状GUIの動作確認済みバージョンは，1.2.0まで，1.3.0も未確認ですが理論上動きます．

# 主要機能
- 要求された初期推力に対する初期条件の計算モード
- 初期設計に対するエンジンの燃焼設計

現状酸化剤は亜酸化窒素のみ，燃料はMMA，ABSの二種類から選択．
燃料形状は純粋な円筒のみ，そのほか形状に対する計算をdevで実装中．

ABSについては，AとBとSの質量分率を調整可能．
# Getting Started

## 要求環境
環境要求についてまとめる．
以下に記載しているのは開発者の環境となっている．
uvで管理しているので，そちらでも可．
CLIは適当な環境でも動作確認済み．

- windows11
- Python 3.11.15
- matplotlib 3.10.7
- numpy 2.3.4
- pandas 2.3.3
- scipy 1.16.3
- cea 3.3.4
- tqdm 4.70.1
- そのほか標準ライブラリ(mathなど)

このコードはpythonライブラリで完結しているので動作保証はしていませんが，Macでも動くと思われます．

## 実行まで
環境構築を行ったうえで，実際に実行するまでの手順をまとめる．

- プログラム側については，このリポジトリにあるもので動作する．
- inputdatasフォルダを用意し，N2Oの圧力-密度データをcsvで用意する．
    - 基本NISTの thermophysical properties から用意する．
- 動作用のjsonc fileを記述する，
    - sampleを参考に値を調整する．コメントのパース機能はつけてあるので，サンプル通りの書式のコメントはパースする．
- main.pyを実行する．

## main.pyの操作
実行後のコマンドライン操作についてまとめる．

- 実行すると，mode選択を聞かれる． 
    - 1だと初期条件の計算，2だと，時間発展の計算を行う．
    - 2のほうのjsoncには一部1の結果を張る部分があるので，ファイルを用意し一旦1を先に実行する．
- modeを選択すると，ファイル名の入力を聞かれる．
    - "filename.jsonc"を入力すると，データを読み込んでくれる．
- 問題がなければ出力に読み込み結果を記載し，計算を開始する．
- 終了するとresultを提示する．
- 初期条件モードは最後にjsonコピー用データが出てくるので，それを時間発展用のjsonにコピーし，ほかの必要条件を追記する．
- 再度実行し，今度は2を選択する．
- 同様にファイル名を聞かれるので，時間発展用のファイル名を入れる
- 問題がなければ出力に読み込み結果を記載し，計算を開始する．
- 終了すると，結果の概要と詳細のcsvファイル保存を開始する．
    - csvはファイル名を聞かれるので，"filename.csv"を入力する．

注意点
- ファイル名は拡張子まで入力する．
- 型チェックは行っているものの，負値のパースや単位系等のチェックはしていないので，入力の妥当性には注意する．


## 結果の確認
計算後のterminal出力と，resultviewer.pyで基本的な結果を確認できる．

terminal出力では，入力として用意した諸元値以外に
- トータルインパルス
- 燃焼後燃料内径
- 平均推力
- 比推力(燃焼時間平均)

を確認できる．

resultviewerでは，出力したcsvから各パラメータの時間推移を確認できる．

# Reference
全部書くと多いので，データベースおよび内容の被りが少ない主要なものを書いています．

詳細は各種コードのコメントアウトおよび，資料のほうを見てください．

NIST thermophysical properties https://webbook.nist.gov/chemistry/fluid/

George P Sutton and Oscar Biblarz. Rocket Propulsion Elements. John Wiley & Sons, 9th edition, 2016. ISBN 9781118753910.

Richard Nakka's Experimental Rocketry Web Site, https://www.nakka-rocketry.net/

Sanford Gordan and Bonnie J. McBride. Computer program for calculation of complex chemical equilibrium compositions and applications. I: analysis. 1994. NASA-RP-1311. URL: https://ntrs.nasa.gov/citations/19950013764.

Sanford Gordan and Bonnie J. McBride. Computer program for calculation of complex chemical equilibrium compositions and applications. II. users manual and program description. 1996. NASA-RP-1311. URL: https://ntrs.nasa.gov/citations/19960044559.

Rocket Engines – Introduction to Aerospace Flight Vehicles, https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/rocket-engines/

Karp, Ashley Chandler and Jens, Elizabeth Therese. Hybrid rocket propulsion design handbook, 2024 ISBN 9780128161999.
https://www.sciencedirect.com/book/monograph/9780128161999/hybrid-rocket-propulsion-design-handbook

# NASACEA license 
Modifications:

   Custom 4 species data added by FROM THE EARTH Tohoku University Student Rocket team, 2026.
   Original NASA CEA code is Apache License 2.0.
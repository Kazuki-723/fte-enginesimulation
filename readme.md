# FTE enginesimulator
過去に使っていたExcelデータのフローの改良及び，データ出力をおこなってくれるソフトのGUIVersionです．
※リリースしたCUIとGUIは計算式をそろえていますが，同一の結果が出るかは検証していません．

# Getting Started

## 要求環境
環境要求についてまとめる．
以下に記載しているのは開発者の環境となっている．
uvで管理して(もらって)いるので，そちらでも可．

- windows11
- Python 3.13.0
- flet 0.28.3
- matplotlib 3.10.7
- numpy 2.3.4
- pandas 2.3.3
- scipy 1.16.3

標準ライブラリについては，csv, re, base64, io, math, os, subprocessを使用している．
一部windows標準コマンドを実行しているためおそらくwindows以外は非対応(未確認)．

## 実行まで
環境構築を行ったうえで，実際に実行するまでの手順をまとめる．

- プログラム側については，このリポジトリにあるもので動作する．
- main-gui.pyと同階層にdataフォルダを作成し，内部にNASACEAを用意する．
- inputdatasフォルダを用意し，N2Oの圧力-密度データをcsvで用意する．
- main.guiを実行する．

# Help
ターミナルにエラーコードが出るものはそれで一度調べてみてください．
それ以外を大体記述しています．

両方
- NASACEA周りのエラー
    - 日本語が含まれていないディレクトリに全部おいて，main.pyのある階層から実行してください．
    - main.pyと同階層に\dataを作り，その中にFCEA2m.exe含めてすべてあるか確認してください．
    - NASACEA単体での動作を確認してください．
    - 単位整合含め入力として適切な値を入れているか確認してください．

GUI version
- UIのボタンが足りない，何も表示されない
    - flet側のバグ(かは知りませんが)でたまに何もない画面を起動することがあります．一度コードを再起動してください．
    - windows本体のディスプレイ設定から，画面のピクセル数，拡大率を調整してください．
        - 開発者がさぼってabsoluteで要素指定をしているので表示できる可能性があります．

- 結果が表示されない，動いているかわからない
    - vscode他のエディタ側のターミナルに逐次計算結果を表示するようにしています．
        - ターミナルに計算結果が流れていれば動いています．計算しきるまで待ってあげてください．
        - 表示されていない場合には，計算できていません．すべての入力に適切な値を入れているか確認してください．

- csvファイルの名前を変えたい
    - main_GUI.pyの280行目`filename = f"filename.csv"`で変更できます．(filename.csvが出力されます．)

CUI version
- 入力値の編集に関して
    - rocket_constants.pyに各パラメータの設定が要るので入力してください．全部floatです．

- csvファイルの名前を変えたい
    - rocket_simulation.pyの360行目`filename = f"filename.csv"`で変更できます．(filename.csvが出力されます．)



## GUIファイル以外について
main_gui.pyでGUIversionを回すために必要なのは，

- inputprograms
- inputdatas (NISTから密度曲線のcsvを持ってくるところ)
- data (CEA置いておくところ)

で，残りは，過去のCUIversion用のフォルダとなっている．

CUI versionは，
- データ系のフォルダ(GUIと共通)
- pythonfiles(endingsimu.pyは未使用)

で動作する．
FLXsimulation.pyはCUIversionをクラス化して分割する前のコード．アルゴリズム周りがかなり違っているので，使わない．

おおよそ，
FLXsimulation.pyが一番古いデータで，そこからフォルダ整理したのが，main.py-pythonfliesの部分
残りはおおよそ試行錯誤中の子たち．
# Version History
 2025/12/4 V1.0.0 リリース

# developer's Memo
実装予定機能，不審挙動についてまとめる．

- x,yのvalueを選択してプロットする機能
- UI関係
- 型判別
- 簡易再現法の比較検証機能
- フランジ回り

Fix
- UIのabsolute → relative実装
- result csvを出した後に計算ボタンが処理を受け付けない
- 値エラー時に実行ボタンを押した際の挙動の統一
    - ボタンにエラーを出して押せなくする．
- 有効数字まわり
-eqn.8の誤植
-ノズル周りにおける物性値の定義


# Reference
George P Sutton and Oscar Biblarz. Rocket Propulsion Elements. John Wiley & Sons, 9th edition, 2016. ISBN 9781118753910.

Richard Nakka's Experimental Rocketry Web Site, https://www.nakka-rocketry.net/

Sanford Gordan and Bonnie J. McBride. Computer program for calculation of complex chemical equilibrium compositions and applications. I: analysis. 1994. NASA-RP-1311. URL: https://ntrs.nasa.gov/citations/19950013764.

Sanford Gordan and Bonnie J. McBride. Computer program for calculation of complex chemical equilibrium compositions and applications. II. users manual and program description. 1996. NASA-RP-1311. URL: https://ntrs.nasa.gov/citations/19960044559.

Rocket Engines – Introduction to Aerospace Flight Vehicles, https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/rocket-engines/
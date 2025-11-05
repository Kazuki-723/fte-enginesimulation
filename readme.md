# FTE enginesimulator
過去に使っていたExcelデータのフローの改良及び，データ出力をおこなってくれるソフトのGUIVersionです．
※旧VersionのCUIとは一部計算式及びアルゴリズムが変更されているため微妙に違う結果が出る可能性があります．

# Getting Started

## 要求環境
環境要求についてまとめる．
以下に記載しているのは開発者の環境となている．
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
実行してみておかしかった時の確認点
デバッグ中なので後で作る．

## GUIファイル以外について
main_gui.pyでGUIversionを回すために必要なのは，

- inputprograms
- inputdatas (NISTから密度曲線のcsvを持ってくるところ)
- data (CEA置いておくところ)

で，残りは，過去のCUIversion用のフォルダとなっている．

おおよそ，
FLXsimulation.pyが一番古いデータで，そこからフォルダ整理したのが，main.py-pythonfliesの部分
残りはおおよそ試行錯誤中の子たち．
# Version History
リリース後追記予定

# developer's Memo
実装予定機能，不審挙動についてまとめる．

- 出力側に入力データのprint
- 出力にラベルを振る
- x,yのvalueを選択してプロットする機能
- UI関係

# Reference
George P Sutton and Oscar Biblarz. Rocket Propulsion Elements. John Wiley & Sons, 9th edition, 2016. ISBN 9781118753910.

Richard Nakka's Experimental Rocketry Web Site, https://www.nakka-rocketry.net/

Sanford Gordan and Bonnie J. McBride. Computer program for calculation of complex chemical equilibrium compositions and applications. I: analysis. 1994. NASA-RP-1311. URL: https://ntrs.nasa.gov/citations/19950013764.

Sanford Gordan and Bonnie J. McBride. Computer program for calculation of complex chemical equilibrium compositions and applications. II. users manual and program description. 1996. NASA-RP-1311. URL: https://ntrs.nasa.gov/citations/19960044559.
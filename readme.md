# simple wave calculation by BEM
時間発展する簡単な系の境界要素法による数値計算のコード.
先生にもらったコードを現代風に書き直したもの.

# 注意
このコードは時間軸の基底に区分一定基底を用いている. したがってそれほど正確にはならない.

また, 著作権を考慮してトラッキングしていないが, 先生にもらった`sphere.f`を用いてメッシュデータを作成しなければならない.
`mesh`ディレクトリを作成しその中に`sphere.f`を置き, 普通にコンパイルする. 生成したファイルを実行して整数を入力すると大きいほど細かいメッシュを作成できる.

# コンパイル
まず
```shell
gfortran -c -o utils utils.f90 -I/usr/local/lib -llapack -lblas
```

で依存モジュールをコンパイル. 次に
```shell
gfortran -g -O0 -o wave_simple_remake wave_simple_remake.f90 utils -I/usr/local/lib -llapack -lblas
```

でコンパイルできる. あとは`./wave_simple_remake`などで実行.

## 出力されるもの
* `mesh_gnup.dat`：扱っている境界を表す. `gnuplot -p plot_domain.gp`とすれば色々設定して出力し空間を確認できる（オリジナルでは法線も図示していたがまだ実装していない）.
* `fort.51~`：`50+時間ステップ数`のファイルに座標ごとの数値計算結果が出力される. `gnuplot -p -e "num=xx" plot_about_space.gp`とすると`num`に入れた番号のファイルの数値解と厳密解を重ねてプロットする.
* `fort.40, 41`：適当な点での波の時間変化が出力される. `gnuplot -p plot_about_time.gp`でプロットされる.

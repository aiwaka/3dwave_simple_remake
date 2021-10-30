# simple wave calculation by BEM
時間発展する簡単な系の境界要素法による数値計算のコード.
先生にもらったコードを現代風に書き直したもの.

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
* `mesh_gnup.dat`：扱っている境界を表す. gnuplotを起動し`splot 'mesh_gnup.dat' w l`で空間を確認できる（オリジナルでは法線も図示していたがまだ実装していない）.

# 実二次体の基本単数を列挙するプログラム

## 概要

実二次体の基本単数を求めるプログラムです．    

計算の過程で出てくる整数方程式の解法は，    
連分数展開によるペル方程式の解法を利用して実装しています．

ペル方程式とは，`x^2 - y^2m = ±1`の形の整数方程式のことで，    
√mの連分数展開を用いた解の構成法が知られています．    
右辺が±4のものもペル方程式と呼ばれ，基本単数を求める際にその場合の整数解も必要となりますが，    
この場合は(1+√m)/2の連分数展開を用いた解の構成法が有澤によって提案されています [2]．

理論上は，計算の過程で出てくる整数方程式は，単純な方針で解けるものですが，    
実際に計算を行うと，方針のわかりやすさに反して，実際に解くのはそう簡単でないことがわかります．    
本プログラムは，素朴な方法ではうまく計算できない場合があることを確認するためのプログラムです．

## プログラムの実行方法

以下の環境で動作確認をしています．    
* コンパイラ: g++ (Ubuntu 7.3.0-27ubuntu1~18.04) 7.3.0
* ライブラリ: GMP (gmpxx.h)

g++/GMPを使う場合は以下のコマンドで実行できます．    
実行には時間がかかるので，注意が必要です．
```
$ make
$ ./unit
```

実行すると，実二次体K=Q(√m) (2≦m≦200, mは平方因子をもたない) の基本単数(>1)を出力します．    
ただし，Qは有理数体，Q(α)はQにαを添加して得られる体を表します．    

GMPを使わない場合は，ヘッダファイルの
```cpp
#define GMP
```
をコメントアウトしてからコンパイルしてください．    
この場合，実行すると計算中にオーバフローが発生するため，一部の結果がエラーで表示されます．

g++を使わない場合は個別の環境に合わせてコンパイルして実行してください．


## 説明

このプログラムは，実二次体K=Q(√m) (mは平方因子をもたない) の基本単数(>1)を計算します．    
ただし，Qは有理数体，Q(α)はQに代数的数αを添加して得られる体を表します．    
基本単数の説明は省略します．

K=Q(√m)について，以下の方法でKの基本単数が得られることが知られています.    

DをKの判別式とします．T, Uに関する以下の整数方程式
<pre>
    T^2 - U^2D = ±4
</pre>
の整数解(T, U)を用いると，
<pre>
         T + U√D
    ε = ---------
            2
</pre>
はKの単数となります．特に最小正のU=uとそれに対応するT=tをとると、
<pre>
          t + u√D
    ε0 = ---------
             2
</pre>
がKの基本単数(>1)となります [1].

整数方程式 `T^2 - U^2D = ±4` の解法として，    
U=1,2,...の順に`±4 + U^2D`が平方数となる最小のUを求める単純な方法が考えられますが，    
より簡単に求まる場合は，別の方法で解いています．

計算の過程で64ビットに収まらない正整数を扱うこともあるので，GMPを用いて実装しました．    
ただし，ペル方程式の性質を利用した解法を導入した現在の版では，    
小規模な例であれば64ビット整数で計算できるようになっています．

以下に解法の方針をまとめます．

### m≡2, 3 (mod 4)の場合
この場合は，
<pre>
    S^2 - U^2m = ±1    (*)
</pre>
の解(S, U)，最小解(s, u)を用いて，
<pre>
    ε = S + U√m,    ε0 = s + u√m
</pre>
と書けます [1]．

(\*)の形の整数方程式はペル方程式と呼ばれ，    
√mの連分数展開を用いた最小解の構成法が知られています．    
このプログラムでは，この方法を用いて求めています．


### m≡1 (mod 4) のうち m = t^2 + 4 (tは正整数) と書ける場合

この場合は，
<pre>
    ε0 = t + √m
</pre>
と書けることが示せます [1]．

このプログラムでは，この方法を用いて求めています．


### m≡1 (mod 4) のうち m = t^2 + 4 (tは正整数) と書けない場合

この場合は，
<pre>
    T^2 - U^2m = ±4    (**)
</pre>
の解(T, U)，最小解(t, u)を用いて，
<pre>
    ε = (T + U√m)/2,    ε0 = (t + u√m)/2
</pre>
と書けます [1]．

(\*\*)の形の整数方程式もペル方程式と呼ばれますが，(\*)と全く同じようには解けません．    
ただし，(\*)の解法の拡張として，(1+√m)/2の連分数展開を用いた最小解の構成法が知られています [2]．    
このプログラムでは，この方法を用いて求めています．

## 二次体の基本単数

実二次体K=Q(√m) (mは平方因子をもたない) の基本単数(>1)は次の通りです．    
mはQ(√m)のm, N(ε0)は基本単数ε0のノルムを表します．    
overflowの列は，ペル方程式の性質を利用しないとオーバフローするものを「\*」で表しています．    
また，右辺が±4になるペル方程式の性質を利用しないとオーバフローするものを「\*\*」で表しています．

以下の表では，mが大きくても200程度のものしか載せていませんが，    
それでもいくつかのmで，TやUにあたる整数がかなり大きい値をとる場合があるとわかります．    
素朴な実装ではそのような場合に平方数判定が難しくなり，    
方針のわかりやすさに反して，実際に解くのはそう簡単でないことがわかります．

なお，ペル方程式の性質を利用した現在の版では，    
mが200未満であれば，すべてのQ(√m)の基本単数が64ビット整数で解くことができます (GMPは不要です)．

|m|N(ε0)|ε0|overflow|
|---:|---:|---:|:---:|
|2|-1|1+√2||
|3|1|2+√3||
|5|-1|(1+√5)/2||
|6|1|5+2√6||
|7|1|8+3√7||
|10|-1|3+√10||
|11|1|10+3√11||
|13|-1|(3+√13)/2||
|14|1|15+4√14||
|15|1|4+√15||
|17|-1|4+√17||
|19|1|170+39√19||
|21|1|(5+√21)/2||
|22|1|197+42√22||
|23|1|24+5√23||
|26|-1|5+√26||
|29|-1|(5+√29)/2||
|30|1|11+2√30||
|31|1|1520+273√31||
|33|1|23+4√33||
|34|1|35+6√34||
|35|1|6+√35||
|37|-1|6+√37||
|38|1|37+6√38||
|39|1|25+4√39||
|41|-1|32+5√41||
|42|1|13+2√42||
|43|1|3482+531√43||
|46|1|24335+3588√46||
|47|1|48+7√47||
|51|1|50+7√51||
|53|-1|(7+√53)/2||
|55|1|89+12√55||
|57|1|151+20√57||
|58|-1|99+13√58||
|59|1|530+69√59||
|61|-1|(39+5√61)/2||
|62|1|63+8√62||
|65|-1|8+√65||
|66|1|65+8√66||
|67|1|48842+5967√67|\*|
|69|1|(25+3√69)/2||
|70|1|251+30√70||
|71|1|3480+413√71||
|73|-1|1068+125√73||
|74|-1|43+5√74||
|77|1|(9+√77)/2||
|78|1|53+6√78||
|79|1|80+9√79||
|82|-1|9+√82||
|83|1|82+9√83||
|85|-1|(9+√85)/2||
|86|1|10405+1122√86||
|87|1|28+3√87||
|89|-1|500+53√89||
|91|1|1574+165√91||
|93|1|(29+3√93)/2||
|94|1|2143295+221064√94|\*|
|95|1|39+4√95||
|97|-1|5604+569√97||
|101|-1|10+√101||
|102|1|101+10√102||
|103|1|227528+22419√103|\*|
|105|1|41+4√105||
|106|-1|4005+389√106||
|107|1|962+93√107||
|109|-1|(261+25√109)/2||
|110|1|21+2√110||
|111|1|295+28√111||
|113|-1|776+73√113||
|114|1|1025+96√114||
|115|1|1126+105√115||
|118|1|306917+28254√118|\*|
|119|1|120+11√119||
|122|-1|11+√122||
|123|1|122+11√123||
|127|1|4730624+419775√127|\*|
|129|1|16855+1484√129||
|130|-1|57+5√130||
|131|1|10610+927√131||
|133|1|(173+15√133)/2||
|134|1|145925+12606√134|\*|
|137|-1|1744+149√137||
|138|1|47+4√138||
|139|1|77563250+6578829√139|\*|
|141|1|95+8√141||
|142|1|143+12√142||
|143|1|12+√143||
|145|-1|12+√145||
|146|1|145+12√146||
|149|-1|(61+5√149)/2||
|151|1|1728148040+140634693√151|\*|
|154|1|21295+1716√154||
|155|1|249+20√155||
|157|-1|(213+17√157)/2||
|158|1|7743+616√158||
|159|1|1324+105√159||
|161|1|11775+928√161||
|163|1|64080026+5019135√163|\*|
|165|1|(13+√165)/2||
|166|1|1700902565+132015642√166|\*|
|167|1|168+13√167||
|170|-1|13+√170||
|173|-1|(13+√173)/2||
|174|1|1451+110√174||
|177|1|62423+4692√177|\*\*|
|178|1|1601+120√178||
|179|1|4190210+313191√179|\*|
|181|-1|(1305+97√181)/2||
|182|1|27+2√182||
|183|1|487+36√183||
|185|-1|68+5√185||
|186|1|7501+550√186||
|187|1|1682+123√187||
|190|1|52021+3774√190|\*|
|191|1|8994000+650783√191|\*|
|193|-1|1764132+126985√193|\*\*|
|194|1|195+14√194||
|195|1|14+√195||
|197|-1|14+√197||
|199|1|16266196520+1153080099√199|\*|


## 参考文献
[1] 石田 信，数学全書5 代数的整数論，森北出版株式会社，東京，1985.    
[2] 有澤 健治，平方根の連分数とペル方程式 第3版，http://ar.nyx.link/cf/pell.pdf, 2018 (最終閲覧日: 2019/2/14).

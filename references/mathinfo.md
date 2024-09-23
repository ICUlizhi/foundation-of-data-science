# 1 高维空间

## 1.1 大数定律

我们直觉上能够感觉到，对于某项数据，当样本数量足够大的时候，其平均值就会越来越接近这项数据的期望。在统计学上用大数定理来描述，即

$Prob(\vert\frac{x_{1}+x_{2}+\cdots+x_{n}}{n}-E(x)\vert\geq\epsilon)\leq\frac{Var(x)}{n\epsilon^{2}}$ (1)

其中E(x)和Var(x)分别代表了x的数学期望和方差。为证明这一定律，我们使用两个不等式作为引理。

引理1.1（Markov＇s inequality）．设x是一个非负随机变量。对于∀a&gt;0，有

$P(x\geqa)\leq\frac{E(x)}{a}$ (2)

证明1.1．根据随机变量期望的定义可得

$E(x)=\int_{0}^{+\infty}xp(x)dx\int_{0}^{a}xp(x)dx+\int_{a}^{+\infty}xp(x)dx$

(3)

$\geq\int_{a}^{+\infty}xp(x)dx\geqa\int_{a}^{+\infty}p(x)dx=aP(x\geqa)$

因此有$P(x\geqa)\leq\frac{E(x)}{a}$。

推论1．$P(x\geqbE(x))\leq\frac{1}{b}$

马尔可夫不等式仅仅使用了数据的期望就控制住了数据“尾部”的概率。下面的切比雪夫不等式利用数据的方差给出了更强的约束。

引理1.2（Chebyshev＇s inequality）．设x是一个非负随机变量。对于∀c&gt;0，有

$P(\vertx-E(x)\vert\geqc)\leq\frac{Var(x)}{c^{2}}$ (4)

证明1.2．我们知道$P(\vertx-E(x)\vert\geqc)=P(\vertx-E(x)\vert^{2}\geqc^{2})$。不妨记$y=\vertx-E(x)\vert^{2}$，显然有E(y)=Var(x)，运用公式（2）即得：

$P(\vertx-E(x)\vert\geqc)=P(\vertx-E(x)\vert^{2}\geqc^{2})\leq\frac{E(y)}{c^{2}}=\frac{Var(x)}{c^{2}}$ (5)

此外，还有一些计算公式：

E(x+y)=E(x)+E(y)

Var(x-c)=Var(x) (6)

$Var(cx)=c^{2}Var(x)$

此外，如果x，y相互独立，那么有E(xy)=E(x)E(y)。因此对于两个独立的随机变量，有

$Var(x+y)=E(x+y)^{2}-E^{2}(x+y)=Var(x)+Var(y)$ (7)

定理1.3（大数定理）．设$x_{1},x_{2},\cdots,x_{n}$是随机变量x的n个独立样本。那么有

$Prob(\vert\frac{x_{1}+x_{2}+\cdots+x_{n}}{n}-E(x)\vert\geq\epsilon)\leq\frac{Var(x)}{n\epsilon^{2}}$ (8)

Proof.由于$E(\frac{x_{1}+x_{2}+\cdots+x_{n}}{n})=E(x)$，因此

$P(\vert\frac{\vertx_{1}+x_{2}+\cdots+x_{n}}{n}-E(x)\vert\geq\epsilon)=P(\frac{\vertx_{1}+x_{2}+\cdots+x_{n}}{n}-E(\frac{x_{1}+x_{2}+\cdots+x_{n}}{n})\vert\geq\epsilon$

$\leq\frac{Var(\frac{x_{1}+x_{2}+\cdots+x_{n}}{n})}{\epsilon^{2}}$

$=\frac{1}{n^{2}\epsilon^{2}}Var(x_{1}+x_{2}+\cdots+x_{n})$

$=\frac{1}{n^{2}\epsilon^{2}}(Var(x_{1})+Var(x_{2})+\cdots+Var(x_{n}))\\=\frac{Var(x)}{2}$ $=\frac{Var(x)}{n\epsilon^{2}}$

(9)

□

## 1.2 高维空间中的几何

现在我们来考虑$R^{d}$中的几何。高维空间中的几何与二维、三维有许多截然不同的性质。比如，设A是$R^{d}$中的一个单位球体，很容易得出d维空间中1-ε≤r≤1描述的单位球体表面厚度为ε的球壳的体积是$V_{\epsilon}=[1-(1-\epsilon)^{d}]volume(A)$，显然有

$\lim_{d\rightarrow\infty}V_{\epsilon}=volume(A)$ (10)

也就是说，高维球体的体积几乎全部集中在球的表面。由不等式$1-tx\leq(1-x)^{t}\leqe^{-tx}$，可以得出

$\frac{V_{1-\epsilon}}{V}=(1-\epsilon)^{d}\leqe^{-\epsilond}$ (11)

下文的符号规定：记V(d,r)和A(d,r)分别为$R^{d}$中半径为r的球体的“体积”和“表面积”，并在r=1时简写为V(d)和A(d).从量纲中我们不难得出$V(d,r)=V(d)r^{d},A(d,r)=A(d)r^{d-1}.$

那么单位球的体积显然应该定义为

$V(d)=\int_{x_{1}=-1}^{x_{1}=1}\int_{x_{2}=-\sqrt{1-x_{1}^{2}}}{1-x_{1}^{2}}\cdots\int_{-\sqrt{1-\sum_{i=1}^{4-1}x_{i}}}dx_{d$ (12)

在直角坐标系下积分形式较为复杂1。用类比的方法，求$R^{3}$中球体的体积，我们可以认为是用球面面积对半径进行积分，即

$V(3,R)=\int_{r=0}^{r=R}(4\pir^{2})dr=\frac{4\pi}{3}R^{3}=\int_{r=0}^{r=R}A(3,r)dr$

类比到高维情形，则有

$V(d,R)=\int_{r=0}^{r=R}A(d,R)dr=A(d)\int_{r=0}^{r=R}r^{d-1}dr=\frac{A(d)}{d}R^{d}$ (13)

因此，要求出d维单位球体的体积，就要计算其“面积”。利用正态分布中的结果$\int_{-\infty}^{+\infty}e^{-x^{2}}dx=$ $\sqrt{\pi},$，我们构造积分

1或许读者能够记得在数学分析第三册P297给出了从直角坐标到球坐标的积分过程。

$I=\int\cdots\int_{x\inR^{d}}e^{-\sumx_{i}^{2}}dx_{1}\cdotsdx_{d}=(\sqrt{\pi})d$ (14)

另一方面，如果我们作球坐标变换，因为$e^{-\sumx_{i}^{2}}=e^{-r^{2}},$

$I=\int_{0}^{+\infty}e^{-r^{2}}A(d,r)dr=A(d)\int_{0}^{+\infty}e^{-r^{2}}r^{d-1}=\frac{1}{2}A(d)\Gamma(\frac{d}{2})$ (15)

这样就确定出2

$A(d)=\frac{2(\sqrt{\pi})^{d}}{\Gamma(\frac{d}{2})}$

(16)

$V(d)=\frac{2(\sqrt{\pi})^{d}}{d\Gamma(\frac{d}{2})}$

#### 1.2.1 赤道附近的体积

另一个非常有趣的结果是，高维球体几乎所有的体积都集中在“赤道”附近，换句话说，如果指定向量v的方向为“北方”，那么绝大多数单位向量u满足$u\cdotv=O(1/\sqrt{d})$。换句话说，绝大多数的单位向量满足$\vertx_{1}\vert=O(1/\sqrt{d})$

定理1.4．对于∀c≥1,d≥3，满足$\vertx_{1}\vert\leq\frac{c}{\sqrt{d-1}}$的点集的体积至少占据球体的$1-\frac{2}{c}e^{-\frac{c^{2}}{2}}$。

Proof.根据对称性，我们只需要证明对于$x_{1}\geq0$的上半球体，至多有$\frac{2}{c}e^{-\frac{c^{2}}{2}}$的体积满足$x_{1}\geq$ $\frac{c}{\sqrt{d-1}}$。我们将半球的上半部分记为A，整个半球记为H，也就是说我们需要证明

$\frac{V(A)}{V(H)}\leq\frac{upperboundofV(A)}{V(H)}=\frac{2}{c}e^{-\frac{c^{2}}{2}}$ (17)

V(A)的精确表达式为

$V(A)=\int_{\frac{c}{\sqrt{d-1}}^{1}V(d-1)(1-x_{1}^{2})^{\frac{d-1}{2}}dx_{1}$ (18)

利用不等式$1-x\leqe^{x}$，我们有

$V(A)\leqV(d-1)\int_{\frac{c}{\sqrt{d-1}}^{\infty}e^{-\frac{d-1}{2}x_{1}^{2}}dx_{1}$ (19)

这个积分不容易计算，由于在$x\geq\frac{c}{\sqrt{d-1}}$时有$\frac{x\sqrt{d-1}}{c}\geq1,$，因此插入这一项，使得积分容易算出：

$V(A)\leqV(d-1)\frac{\sqrt{d-1}}{c}\int_{\frac{c}{\sqrt{a-1}}}^{\infty}x_{1}e^{-\frac{d-1}{2}x_{1}^{2}}dx_{1}=\frac{V(d-1)}{c\sqrt{d-1}}e^{$ (20)

虽然我们已经有了H的表达式，但是其表达式较为复杂，因此我们估计一个下界。为了能消去上面得出的$\sqrt{d-1}$，我们选取上半球中高度为$\frac{1}{\sqrt{d-1}}$的圆柱，从而有

$V(H)\leqV(d-1)(1-\frac{1}{d-1})^{\frac{d-1}{2}}\cdot\frac{1}{\sqrt{d-1}}\leq\frac{V(d-1)}{2\sqrt{d-1}}$ (21)

2关于Γ函数的定义和部分性质可参见附录一

综上就可以得出

$\frac{V(A)}{V(H)}\leq\frac{\frac{V(d-1)}{c\sqrt{d-1}}e^{-\frac{c^{2}}{2}}}{\frac{V(d-1)}{2\sqrt{d-1}}}=\frac{2}{c}e^{-\frac{c^{2}}{2}}$ □ (22)

-

从上面的分析中可以得知，任意选取两个点它们的内积大概率接近于0，也就是任选两个向量，它们都非常接近垂直。下面的定理可以更精确地告诉我们对于n个点的情形。

定理1.5．对于在单位球面上n个随机选取的点$x_{1},x_{2},\cdots,x_{n}$，有1-O(1/n)的可能性

1．对于所有的i都有$\vertx_{i}\vert\geq1-\frac{2\logn}{d}$

2．对于所有的i≠j都有$\vertx_{i}\cdotx_{j}\vert\leq\frac{\sqrt{6\logn}}{\sqrt{d-1}}$

Proof.第一条定理根据前面我们的结果，满足$\vertx_{i}\vert\leq1-\epsilon$的点的比例小于$e^{-\epsilond}$。因此有

$P(\vertx_{i}\vert\geq1-\frac{2\logn}{d})\leqe^{-(\frac{2\logn}{d})d}=\frac{1}{n^{2}}$ (23)

根据概率公式$P(V_{i=1}^{n}x_{i})\leq\sum_{i=1}^{n}P(x_{i})$，我们有

$P(-\cup_{i=1}^{n}(\vertx_{i}\vert\geq1-\frac{2\logn}{d}))\geq1-\frac{1}{n}$ (24)

对于第二条定理，由于共有$C_{n}^{2}=\frac{n(n-1)}{2}$ 种组合，我们前面已经证明一个单位向量与另一个单位向量内积大于 $\frac{c}{\sqrt{d-1}}$ 的概率至多是$\frac{2}{c}e^{-\frac{c^{2}}{2}}$，因此对于$c=\sqrt{6\logn}$，其每一对的概率为$l(e^{-\frac{6\log_{n}}{2}})=O(\frac{1}{n^{3}})$，因此$C_{n}^{2}$对中内积大于给定值的概率为O(1/n) -

### 1.3 随机投影和 Johnson-Lindenstrauss引理

在涉及高维数据的处理时，我们经常需要找到一个高维点的最近邻点。由于我们需要处理n个d维数据，而n，d通常都非常大，因此难以做到直接在原数据上进行搜索。一般来说，我们希望能够在O(logn),O(logd)的多项式时间内完成搜索，而在此之前的预处理可以花费n，d的多项式时间。通常的操作是将这些d维数据在尽量保证距离变化不大的情形下映射到维度较小的k维空间。使用Gaussian Annulus Theorem（高斯环定理），我们发现这样的映射确实存在，而且并不复杂。

定理1.6（Gaussian Annulus Theorem）．对于一个每个维度都是一个方差为1的d维球面高斯向量，对任意$\beta\leq\sqrt{d},$至多$3e^{-c\beta^{2}}$的概率其位置位于环$\sqrt{d}-\beta\leq\vertx\vert\leq\sqrt{d}+\beta之外,这里c是一$个正常数.

要使用的映射f:$R^{d}\rightarrowR^{k}$定义如下：随机选取$R^{d}$中k个坐标服从方差为1的高斯分布的Gaussian Vector$u_{1},u_{2},\cdots,u_{k}.$ 对于$\forallv\inR^{d},$定义

$f(v)=(u_{1}\cdotv,u_{2}\cdotv,\cdots,u_{k}\cdotv)$

我们将证明，会以高概率满足$\vertf(v)\vert\approx\sqrt{k}\vertv\vert.$因为映射满足$f(v_{2}-v_{1})=f(v_{2})-f(v_{1})$，故如果要估计$\vertv_{2}-v_{1}\vert$，只需要计算映射后的点在$R^{k}$中的距离即可，因为$\sqrt{k}$是一个已知的常数因子。距离会增大的原因是我们的$u_{i}$并不是单位向量。另外注意$u_{i}$并不是相互正交的。如果我们要求正交性，$u_{i}$之间将会失去独立性。

定理1.7．设v是Rd中的任意一个向量，f的定义如上所述。则存在c&gt;0使得∀ε∈(0,1)，有

$Prob(\vert\vertf(v)\vert-\sqrt{k}\vertv\vert\vert\geq\epsilon\sqrt{k}\vertv\vert)\leq3e^{-ckz^{2}}$ (25)

随机性来源于用于构造映射f的随机向量$u_{i}$的随机性。

Proof.由于v可被单位化，不妨认为|v|=1.因此映射后的每个坐标$u_{i}$·v实际上是d个Gaussian Vector的线性组合，均值为0，方差为

$Var(u_{i}\cdotv)=\sum_{j=1}^{d}v_{j}^{2}Var(u_{ij})=\sum_{j=1}^{d}v_{j}^{2}=1$

因此映射后的坐标就是在$R^{k}$中随机选取的单位Gaussian向量。再使用 Gaussian Annulus Theo-rem 即可证明。 -

课本上对于Gaussian Annulus Theorem的证明使用了一些结论来bound 随机变量的矩，对Chernoff Bound 及其证明过程较为熟悉的读者可以参照附录二对于 Johnson-Lindenstrauss 引理的证明3。

3笔者按：附录二的过程是孔老师在随机算法课上教授的，她认为两个证明过程我们会更喜欢用Gaussian Annulus Thm 的方法，但笔者觉得哪个方法都不是正常人能想出来的

### 2 奇异值分解

#### 2.1 奇异值分解简介

对于矩阵$A_{n\timesd}$，我们将其每一行视为一个d维的向量，奇异值分解的目标是寻找矩阵$U_{n\timesr_{1}}D_{r\timesr_{r}}V_{r\timesd}$使得

$A=UDV^{T}$ (26)

其中U，V是正交矩阵，D是对角元全为正数对角矩阵。其中V的列向量（也就是$V^{T}$的行向量）即为“最佳近似子空间”的标准正交基向量。而U的元素就是A在相应基向量的投影值。也就是

(27)

<!-- A = [ u _ { 1 } , u _ { 2 } \cdots u _ { r } ] d i a g [ v _ { 1 } , \sigma _ { 2 } , \cdots , \sigma _ { r } \} [ \begin{matrix} v _ { 1 } ^ { T } \\ \overline { v _ { r } ^ { T } \end{matrix} ]  -->
![](https://textin-image-store-1303028177.cos.ap-shanghai.myqcloud.com/external/ad2fb5922a4b4e43)

将上述写法展开，就是

$A=\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}$ (28)

其中的$u_{i}$是n×1的向量，$v_{i}$是d×1的向量。

就像将x投影到一组基向量$v_{1},v_{2},\cdots,v_{d}$中一样，对于SVD，前k个向量保证A的行向量在k维子空间中的投影平方和最大（距离平方和最小）。

我们已经了解过特征值分解。对于方阵A，满足Ax=λx的向量x称为特征向量，对应的λ称为特征值。如果A实对称，那么A必然有n个特征值，从而A就可以分解为

$A=P^{T}BP$ (29)

其中P是正交矩阵，$B=diag\{\lambda_{1},\lambda_{2},\cdots,\lambda_{n}\}$，P的行向量维对应的特征向量。

特征值分解对于A有着要求。首先必须是方阵，其次必须有n个特征向量。但是SVD对于矩阵没有要求，我们总可以找到这样的两个由单位正交向量组成的矩阵U，V和对角矩阵D。其中V的列向量称为右奇异向量，U的列向量称为左奇异向量。如果A是可逆的，那么

$A^{-1}=VD^{-1}U^{T}$ (30)

我们接下来将会看到，有下面的关系式成立：

$Av_{i}=\sigma_{i}u_{i}A^{T}u_{i}=\sigma_{i}v_{i}$ (31)

这一点不难验证，将公式27右乘$[v_{1}v_{2}\cdotsv_{r}]$可知

$A[v_{1}v_{2}\cdotsv_{r}]=[\sigma_{1}u_{1}\sigma_{2}u_{2}\cdots\sigma_{r}u_{r}]$ (32)

左乘$[u_{1}u_{2}\cdotsu_{r}]^{T}$可得

$A^{T}[u_{1}u_{2}\cdotsu_{r}]=[\sigma_{1}v_{1}\sigma_{2}v_{2}\cdots\sigma_{r}v_{r}]$ (33)

换言之，A作用于$v_{i}$会将其变换为$u_{i}$的$\sigma_{i}$倍，$A^{T}$作用于$u_{i}$会将其变换为$v_{i}$的$\sigma_{i}$倍。注$A^{T}Av_{i}=d_{i}^{2}v_{i}$，也就是说第i个奇异向量$v_{i}$是方阵$A^{T}A$的第i个特征向量。意到

我们选取奇异向量的标准是：其生成的子空间能够让矩阵A在其上的投影平方和最大。我们不妨考虑一维的情形。对于一个向量$a_{i}$和一条过原点的单位方向向量为v直线，用disti和proji 分别代表向量$a_{i}$到直线的距离和投影，我们有

$(disti)^{2}+(proji)^{2}=\vert\verta_{1}\vert\vert$ (34)

那么最小化距离平方和实际上就是最大化投影平方和。前者的表述类似于最小二乘法，后者代表在同维度的子空间中最大限度地保留了原数据的信息。可能你会觉得选取投影平方和作为数据留存度的指标略显武断，但是我们随后就会看到，“平方”拥有非常良好的性质，式34就是其中之

一。

#### 2.2 奇异向量

将n×d的矩阵A的行向量视为d维空间的点，考虑过原点的直线，其单位向量为v，对于$\foralla\inR^{d}$,v·a就是a在v上的投影，因此$\vert\vertAv\vert\vert^{2}$即为A中n个行向量在v上的投影平方和。我们的目的是使得||Av||最小，也就是找到

$v_{1}=\arg\max_{\vertv\vert=1}\vert\vertAv\vert\vert$ (35)

当然$v_{1}$可能不止一个（或者说必然不止一个），如果$v_{1}$是奇异向量，那么-v1也是。随便选择一个即可，后面我们都这样处理。$v_{1}$称为第一个奇异向量。而$\sigma_{1}(A)=\vert\vertAv_{1}\vert\vert$称为A的第一奇异值。因此有

$\sigma_{1}^{2}(A)=\sum_{i=1}^{n}(a_{i}\cdotv_{1})^{2}$ (36)

即为A的所有行向量在$v_{1}$上的投影平方和。

当然，数据未必会集中在一条线附近，可能是一个平面，换句话说，如果我们已经有了一个寻找$v_{1}$的算法，怎么去寻找到更高维的空间？

一种贪心的做法是，在v1的正交补空间用同样的方法寻找第二个奇异向量和相应的奇异值。也就是

$v_{2}=\arg\max_{\vertv\vert=1}\vert\vertAv\vert\vert$ (37)

同样的方法我们可以定义下面的奇异向量，直到

$v\botv_{1},\cdots,v_{r},\vert\vertAv\vert\vert=0\\\vertv\vert=1\end{matrix}\vert=1$ (38)

此时r=rank(A)

如何保证这种贪心法的正确性呢？我们给出下面的定理。

定理2.1．设A是一个n×d的矩阵，其奇异向量为$v_{1},v_{2}\cdots,v_{r}。$ 对于1≤k≤r,，令$V_{k}=$ $L(v_{1},v_{2}\cdots,v_{k})$，那么$V_{k}$就是A的k维最佳近似子空间。

Proof.当k=1时显然成立。对于k=2时，假设$W=L(w_{1},w_{2})$是一个二维的子空间，我们选取$w_{2}\botv_{1}。$方法是：如果$v_{1}\botW,$，则任选一个即可，否则，选取W中垂直于$v_{1}$投影的向量即可。这样，由$v_{1},v_{2}$的定义可知$\vert\vertAw_{1}\vert\vert^{2}\leq\vert\vertAv_{1}\vert\vert^{2}$,$\vert\vertAw_{2}\vert\vert^{2}\leq\vert\vertAv_{2}\vert\vert^{2}$，从而有

$\vert\vertAw_{1}\vert\vert^{2}+\vert\vertAw_{2}\vert\vert^{2}\leq\vert\vertAv_{1}\vert^{2}+\vert\vertAv_{1}\vert^{2}$ (39)

对于k维情形也可以用类似方法证明。 -

我们注意到n维向量Avi是一张记录了A的n个行向量在$v_{i}$的投影（带符号）的表格。假如我们认为$\sigma_{1}(A)=\vert\vertAv_{1}\vert\vert$就是矩阵A在v1方向上的一部分。要让这个认识符合我们的一般认识，那么所有这样部分的平方和应该等于"A的全部”。设$a_{j}$是A的第j个行向量，那么有

$\sum_{j=1}^{n}\verta_{j}\vert^{2}=\sum_{j=1}^{n}\sum_{i=1}^{r}(a_{j}\cdotv_{i})^{2}=\sum_{i=1}^{r}\sum_{j=1}^{r}(a_{j}\cdotv_{i})^$ (40)

而所有行向量模方之和实际就是A中所有元素的平方和。这就是我们前面提到的"A的全部”，称为Frobenius 范数，即

$\vert\vertA\vert\vertF=\sqrt{\frac{\sum_{j,k}^{a_{jk}^{2}}}$ (41)

我们前面已经证明了

$\sum_{i=1}^{r}\sigma_{i}^{2}(A)=\vert\vertA\vert\vertF$ (42)

向量$v_{1},v_{2},\cdots,v_{r}$被称为右奇异向量。用A与其作用后正交化，即

$u_{i}=\frac{1}{\sigma_{i}(A)}Av_{i}$ (43)

事实上，Mi就是最大化$\vert\vertu^{T}A\vert\vert$的向量，且彼此正交，称为左奇异向量。

#### 2.3 奇异值分解

设A是一个n×d的矩阵，其r个奇异向量为$v_{1},v_{2},\cdots,v_{n}$，对应的奇异值为$\sigma_{1},\sigma_{2},\cdots,\sigma_{r}.$我们定义其左奇异向量为$u_{i}=\frac{1}{\sigma_{i}}Av_{i}。$。那么$\sigma_{i}u_{i}v^{T}$是一个秩为1的矩阵，其行向量就是矩阵A的行向量“在$v_{i}$方向的部分”，也就是A的行向量在$v_{i}$方向的分量的向量坐标表示。我们将要证明A可以表示为一系列秩为1的矩阵之和，形如

$A=\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}$ (44)

从几何上来说，每个A中行向量代表的点被分解为它在$v_{i}$方向的分量之和。我们也会从代数角度证明这个结论。我们先从下面的引理出发。

引理2.2．如果A和B对于任意向量v都有Av=Bv，那么A=B.

证明2.1．选取标准正交基即可。

定理2.3．设A是一个n×d的矩阵，其r个右奇异向量为$v_{1},v_{2},\cdots,v_{n},$，对应的奇异值为$\sigma_{1},\sigma_{2},\cdots,\sigma_{r}$其左奇异向量为$u_{1},u_{2},\cdots,u_{r}。$ 那么

$A=\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}$ (45)

Proof.对于任意一个$v_{j},$，我们分别用A和$\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}$左乘，有

$\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}v_{j}=\sigma_{j}u_{j}=Av_{j}$ (46)

由于任意一个向量v都可以分解为一个与所有右奇异向量都垂直的向量和右奇异向量的线性组合，所以对于任意v都有$\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}v=\sigma_{j}u_{j}=Av$，由上面的引理即可得出结论 -

这样，我们就完成了对矩阵A的奇异值分解。如果奇异值互异，只需要规定其序关系即得到唯一的分解。如果有相同的奇异值，那么对应的奇异向量会张成子空间，此子空间的任意一组正交基都可以作为奇异向量。

#### 2.4 最佳k维近似

本节的标题是英语奇差的笔者直译的结果。

设A是一个n×d的矩阵，其奇异值分解为

$A=\sum_{i=1}^{r}\sigma_{i}u_{i}v_{i}^{T}$ (47)

对于任意1≤k≤r,，令

$A_{k}=\sum_{i=1}^{k}\sigma_{i}u_{i}v_{i}^{T}$ (48)

显然ran$(A_{k})=k$。我们随后会看到， $A_{k}$是对A的最佳的k维近似，且其误差可用F-范数估计。换言之，$v_{1},v_{2},\cdots,v_{k}$张成的k维空间是所有k维子空间中，使A的行向量距离平方和最小的那一个。为了证明这个结论，我们先给出下面的引理：

引理2.4．$A_{k}$的行向量是对应的A的行向量在$V_{k}$的投影。

证明2.2．设a是任意一个A的行向量，由于$v_{1},v_{2},\cdots,v_{k}$彼此正交，因此它在$V_{k}$的投影就是$\sum_{i=1}^{k}(a\cdotv_{i})v_{i}^{T}$。因此整个矩阵A的投影就是$\sum_{i=1}^{k}Av_{i}v_{i}^{T}.$.于是有

$\sum_{i=1}^{k}Av_{i}v_{i}^{T}=\sum_{i=1}^{k}\sigma_{i}u_{i}v_{i}^{T}=A_{k}$ (49)

定理2.5．对于任意秩不大于k的矩阵B，有

$\vert\vertA-A_{k}\vert\vert\leq\vert\vertA-B\vert\vert$ (50)

Proof.设B是所有秩不大于k的矩阵中使得$\vert\vertA-B\vert\vert^{2}$最小的那一个。设V是B的行向量张成的子空间。那么dimV≤k。由于B是使得||A-B||最小的矩阵，因此，B的行向量就是A相应行向量在V的投影。否则，用对应的投影代替B的行，这仍然保证了B的行向量包含在V中，但是这使得||A-B||减小了，与B的定义不符。

现在由于B的行向量就是A相应行向量在V的投影，而能使得距离平方和最小的子空间的投影构成的矩阵就是$A_{k},$，因此有$\vert\vertA-A_{k}\vert\vert\leq\vert\vertA-B\vert\vert.$-

-

如果我们要对若干数据x计算Ax，那么每次计算要进行nd次乘法，加上加法的时间，时间复杂度为O(nd)。但是如果我们用$A_{k}$来近似代替A，那么计算$A_{k}x$时，只需计算k(n+d)次乘法，时间复杂度为O(k(n+d))。这在k≪min{n,d}的时候能够显著提高计算效率。那么其误差是多少呢？由于x未知，我们应该给出一个对所有x都符合的估计，所以我们选取最大的$\vert(A_{k}-A)x\vert。$当然，在|x|没有约束的话，最大值也没有界限。因此我们限定|x|≤1。这样，我们就定义了A 的一种新的范数，即

$\vert\vertA\vert\vert_{2}=\max_{\vertx\vert\leq1}\vertAx\vert$ (51)

即2-范数。注意到它就是$\sigma_{1}(A)$。4

4本节的内容还有部分遗漏，比如为何左奇异向量相互正交，以及幂法（Power Method）求特征值等，这些在课本上都有较为完整的证明。

### 3 机器学习

### 3.1 误差估计

如果H是一个规则集，ε,δ&gt;0如果一个大小为

$n\geq\frac{1}{\epsilon}(\log\vertH\vert+\log(\frac{1}{\delta}))$

的训练集是依据分布D挑选的，那么有大于1-8的概率，任何hEH且真实错误$err_{D}(h)>\epsilon$都会使得训练错误$err_{S}(h)>0.$换句话说，有大于1-8的概率，任何h∈H且训练错误$err_{S}(h)=0$都会使真实错误$err_{D}(h)<\epsilon$

证明过程如下：设$h_{1},h_{2},\cdots$是H中真实错误率大于等于ε的概念，那么对于一个大小为n 的训练集S，设$A_{i}$是事件“$\cdotsh_{i}$的训练集错误率为0”，那么

$P(A_{i})\leq(1-\epsilon)^{n}$

那么根据联合概率的不等式，至少一个发生的概率

$P(\cup_{i}A_{i})\leq\vertH\vert(1-\epsilon)^{n}\leq\vertH\verte^{-n\epsilon}\leq\vertH\verte^{-1\log\vertn\vert-\log(1/\delta)}=\delta$

上述估计是建立在概念集中有训练错误率为0的概念的基础上的。如果表现最好的概念的错误率也只有5％呢？我们能否在足量的样本下，使得真实错误率也接近训练错误率呢？答案是可以的，要证明这个结论，首先需要一个二项分布的结论：

若随机变量X∼B(n,p)，那么对于任意α∈[0,1]有

$P(X/n>p+\alpha)\leqe^{-2na^{2}}$

$P(X/n<p+\alpha)\leqe^{-2na^{2}}$

利用上述结论很容易证明，如果

$n\geq\frac{1}{2\epsilon^{2}}(\log\vertH\vert+\log(2/\delta))$

就能保证以大于等于1-8的概率，所有h∈H都满足$\verterr_{S}(h)-errD(h)\vert\leq\epsilon.$

### 3.2 Online Learning

到目前为止我们考虑的都是批处理学习（batch learning）。也就是说，给定一批（batch）数据，以及一个训练样本S，你的目标是通过这个训练集，生成一个概念h使它在新的数据上也极少犯错。

我们现在将目光转向更困难的Online Learning的情形，这里我们不能假定数据是根据某一概率分布或是某一概率性的过程来选定的。

OL的过程是这样的：在每个时刻t=1,2,⋯，有两件事情发生：

首先，算法会获得一个随机的样例$x_{t}\in\chi$，并被要求给出一个预测lt作为它的标签。

然后，算法将被告知这个样例的真实标签$c^{*}(x_{t})$，如果$c^{*}(x_{t})\neql_{t}$，那么算法将会获得一个debuff。

#### 3.2.1 逻辑或的学习

这个学习算法的目标是尽可能少犯错误。我们举一个学习“逻辑或”的学习算法的例子。比如在垃圾邮件的判定中，假定“重要邮件”是指满足一系列“重要”要求中的至少一个的邮件，比如一共有d个条件，我们想知道这些条件中哪些是“重要”的。每个到来的邮件将被标记为｛0,1｝d 中的一个标记。从$h_{0}=x_{1}Vx_{2}V\cdotsVx_{d}$开始，每次如果犯错，就将犯错标签从“重要标签”中去掉，犯错至多不超过d次。

实际上，d次同时也是一个下界，我们可以证明，没有任何一种算法能够保证对任何一种样例输入的序列，犯错次数都小于d。只需要每次置第i个特征为1，然后总是告诉算法预测是错误的即可。

#### 3.2.2 折半算法 Halving Algorithm

如果我们不关心算法的运行时间，那么又一种简单的OL能够保证犯错次数不超过$\log_{2}(\vertH\vert)$。对每个给定样例，我们都询问所有H中的规则，并将更多人的选择作为输出，这样如果犯错一次，我们就能至少减少一半的规则。

当然这个算法在H中没有完美规则的时候无法给出出错率最低的那个，它可能一开始就把它剔除了。

#### 3.2.3 感知机算法

这一部分因为12月3日的课程回放只有一节，所以不知道讲了多少

##### 3.2.4 12.10 Random Weighted Majority

如果我们有一个大小为n的H，并且确定有一个完美规则，那么我们有$\log_{2}(n)$的折半算法。

如果并不存在完美规则呢？如何找到出错最少的规则（optimal expert）？经过T轮之后，我们选择出错最少的专家。

我们将会给每个专家一个初始权重为1，每次有一个样例输入，我们将会按照权重选择一个专家的回答作为这一轮的回答（但此时我们仍然知道其他专家的回答）。当我们知道正确答案之后，我们将把所有回答错误的专家的权重乘以1-ε.

由于我们是按照权重选择的，如果令 $\sumw_{i}^{t}$代表所有在第t轮中犯错的的专家的权重的和，wrong 令wi代表所有在第t轮中犯错的的专家的权重的和， $W_{t}=\sumw_{i}^{t}$代表这一轮所有专家的权重right 和，那么我们在这一轮选择错误的概率是

$F_{t}=\frac{\sum_{wrong}w_{i}^{t}}{\sumw_{i}^{t}}$

t轮过后，出错的次数的期望是

$E(wrongtimes)=\sum_{i=1}^{t}F_{i}$

每次更新的时候，我们将对犯错的权重更新，那么

$W_{t+1}=(\sum_{ri_{8}ht}^{t}+(1-\epsilon)\sum_{w_{0}ng}w_{i}^{t})=W_{t}(1-F_{t}\epsilon)$

那么经过T轮之后，假设最优的专家只犯错了$N_{Opt}$次，我们可以知道

$W_{T}\geq(1-\epsilon)^{N_{Opt}}$

（为什么采取如此粗略的估计，大概其他的如果次次犯错，权重下降得更快？）

初始权重为n，那么我们就有

$W_{T}=n\prod_{i=1}^{T}(1-F_{i}\epsilon)\geq(1-\epsilon)^{N_{opt}}$

两侧取对数，并利用log(1-x)≤-x，因此有

E(Total numo$f_{mistakes}=\sum_{i=1}^{T}F_{i}\leq\frac{1}{\epsilon}(\logn-N_{opt}\log(1-\epsilon))$

Batch Learning 比 Online Learning 更简单。

#### 3.2.5 Boosting 提升方法

提升方法是一种可以用来减小监督式学习中偏差的机器学习算法。面对的问题是迈可·肯斯（Michael Kearns）提出的：一组“弱学习者”的集合能否生成一个“强学习者”？弱学习者一般是指一个分类器，它的结果只比随机分类好一点点；强学习者指分类器的结果非常接近真值。

也就是说，弱学习者对于任何分布的样本，其给出的专家的错误率

$P(wrong)\leq\frac{1}{2}-\gamma$

其中$0<\gamma\leq\frac{1}{2}.$.一个强学习者则可以在足量样本（比如说，是$\frac{1}{\epsilon}$的多项式）下以极大概率到达任意小的错误率。以下我们将说明一个事实，即通过提升方法，一个在任何样本分布下都能工作的弱学习者可以提升为一个强学习者。简单地说，我们会给训练集的样本加权，每轮运行结束后我们将提升回答错误的样本的权重（乘以α&gt;1)，注意与RWM的区别是，这里是提升错误样本的权重，而非降低错误专家的权重。在足够多的轮数之后，我们将所有得到的专家$h_{1},h_{2},\cdots,h_{t}$投票决定作为最终的专家。

假如在第t轮运行结束之后，样本的总权值为$W_{t}$，那么在第t+1轮中，回答错误的样本所占的权重小于等于$\frac{1}{2}-\gamma$，从而更新后就有

$W_{t+1}\leqW_{t}((\frac{1}{2}-\gamma)\alpha+\frac{1}{2}+\gamma)$

考虑t轮运行之后，共得到$h_{1},h_{2},\cdots,h_{t}$这t个专家。设m是使得这t个专家投票之后仍然判断错误的样本个数，那么可以知道这m个样本至少令$\frac{t}{2}$个专家判断错误，取定$\alpha=\frac{\frac{1}{2}+\gamma}{\frac{1}{2}-\gamma}$，那么对于t轮运行后的总权值，我们有如下估计：

$m\alpha^{\frac{t}{2}}\leqW_{t}\leqn(1+2\gamma)^{t}$

解出

$m\leq(1-4\gamma^{2})^{\frac{1}{2}}\leqne^{-2t\gamma^{2}}$

上式用到了$1-x\leqe^{-x}$.当$t>\frac{\logn}{2\gamma^{2}}$时，m&lt;1.

本章的一个重要概念-VC dimension 因时间原因未能整理.

### 4 针对大量数据的算法：Streaming，Sketching and Sampling

### 4.1 导引

如果你有n个正数$a_{1},a_{2},\cdots,a_{n}$，并且都不超过m，如果让你以其大小为权重，随机选择一个数字作为输出，你会如何操作？直观的方法是，将这n个数字都存储起来，这样你需要O（nlogm）大小的空间。另一种方法是，存储两个值，sum代表当前输入的所有数字之和，num代表选择结果。初始时，$\sin=num=a_{1}$.当输入$a_{j}$时，更新$\sin<-\sinm+a_{j}$，并以$\frac{a_{j}}{\sin}$的概率将选择结果变为$a_{j}。$ 容易证明算法结束后，$a_{j}$被选中的概率是

$P(a_{j})=\frac{a_{j}}{\sum_{i=1}^{n}a_{i}}$

但是这个算法只需要O(logm+logn)的空间。

### 4.2 数据流中不同元素的个数

对于大量数据的流，统计其中不同元素的个数是很有意义的。考虑n个整数组成的序列$a_{1},a_{2},\cdots,a_{n}$ $1\leqa_{n}\leqm$，并且n，m都非常大。如何统计其中不同元素的个数？你可以用桶的方法在O（m）的空间下完成这一目标，也可以用全部存储的方法在O（nlogm）的空间下完成。但是无论如何，只要你想要精确求解其中的相异元素个数，你所需要的空间至少是O（m）（证明见188页）。但是我们可以用随机和近似的方法，在更小的空间代价之下获得一个较好的近似答案。

要引出我们的估计方法，我们不妨考虑在［0,1］内独立随机挑选s个实数。我们考虑这s个实数最小值的数学期望。由几何概型，最小值的概率分布函数为

$F(x)=P(\min\leqx)=1-(1-x)^{5}$

求导得到其概率密度为$f(x)=s(1-x)^{s-1}.$.因此其数学期望

$E(\min)=\int_{0}^{1}xf(x)dx=\frac{1}{s+1}$

类似的思想，在n，m非常大的时候，在{1,2,3,⋯,m}中独立随机挑选一些数字构成集合S，作出估计$\min\approx\frac{m}{\vertS\vert+1}$，所以$\vertS\vert=\frac{m}{m\in}-1。$

当然这样的估计是建立在数据相互独立的情况下。如果数据不相互独立，比如说，就挑选了{1,2,3,⋯,m}中最小的|S|个数，那么算法会给出非常糟糕的结果。将不相互独立的数据进行完全的随机映射需要大量的空间，但是我们将会看到，通过2-universal 哈希函数可以使得映射后的数据两两独立，这足以让我们完成精确度估计的关键步骤，而且只需要花费O（logm）的空间。

定义一个哈希函数集

H={h|h:{1,2,3,⋯,m}→{0,1,2,3,⋯,M-1}}}

我们说H是2-universal的或者两两独立（pairwise independent）的，是指如果∀x,yE {1,2,3,⋯,m}并且x≠y,h(x)和h(y)能够等概率地取到{0,1,2,3,⋯,M-1}中的值并且两两独立。换句话说，对∀w,z∈{0,1,2,3,⋯,M-1}，有

$P(h(x)=wandh(y)=z)=\frac{1}{M^{2}}$

我们给出一个2-universal的例子。令M是一个大于m的素数，对任意整数对(a,b)∈[0,M-1］2，定义函数

$h_{ab}(x)=ax+b$ mod M

存储这个函数只需要记录两个数字即可，空间为O(logM)。那么h(x)=wandh(y)=z当且仅当

$(\begin{matrix}x&1\\y&1\end{matrix})(\begin{matrix}a\\b\end{matrix})=(\begin{matrix}w\\z\end{matrix})modM$

当x≠y时，有

$(\begin{matrix}a\\b\end{matrix})=(\begin{matrix}x&1\\y&1\end{matrix})^{-1}(\begin{matrix}w\\z\end{matrix})$ mod M

因此(w,z)和(a,b)一一对应，也就是

$P(h(x)=wandh(y)=z)=\frac{1}{M^{2}}$

下面我们证明，有极大概率（大于2／3）使得

$\frac{m}{6s}\leq\min\leq\frac{6m}{s}$

成立。

首先我们估计$P(\min<\frac{m}{6s})$。这部分不需要两两独立，只需要联合概率公式就有

$P(\min<\frac{m}{6s})\leqS\times\frac{m}{6s}\times\frac{1}{m}=\frac{1}{6}$

另一个方向，估计$P(\min>\frac{6m}{s})$.此时我们不能用独立事件相乘的公式，因为我们无法保证全部独立，只能通过Hash方法保证两两独立。我们定义随机变量

$I_{i}=\{\begin{matrix}0,h(a_{i})\geq\frac{6m}{s}\\1,else\end{matrix}$

令$\gamma=\sumI_{i}$，那么Y=0就代表所有数字都大于$\frac{6m}{s}$，我们希望P(Y=0)越小越好。根据切比雪夫不等式，有

$P(Y=0)\leqP(\vertY-E(Y)\vert\geqE(Y))\leq\frac{Var(Y)}{E^{2}(Y)}$

这里$E(Y)=\sumE(I_{i})=s\times\frac{6m}{s}\times\frac{1}{m}=6$，重点在于Var(Y)的计算。

我们知道$Var(X_{1}\pmX_{2})=Var(X_{1})+Var(X_{2})\pm2Cov(X_{1},X_{2})$，其中

$Cov(X_{1},X_{2})=E[(X_{1}-E(X_{1}))(X_{2}-E(X_{2}))]=E(X_{1}X_{2})-E(X_{1})E(X_{2})$

是双线性函数。且按照定义有Cov(X,X)=Var(X).

那么

$Var(\sum_{i=1}^{s}I_{i})=\sum_{i=1}^{s}Var(I_{i})+\sum_{i\neqj}Cov(I_{i},I_{j})$

由于两两独立，因此协方差部分全部为0，也就是说此处仍然成立

$Var(\sum_{i=1}^{s}I_{i})=\sum_{i=1}^{s}Var($i)

而根据伯努利分布的方差有

$Var(\sum_{i=1}^{s}I_{i})\leqs\times\frac{6m}{sm}(1-\frac{6m}{sm})\leq6$

所以有

$P(Y=0)\leq\frac{1}{6}$

### 4.3 Majority Vote

如果有一串超长的字符序列（长度为n），其中有一个字符出现频率超过的一半，每次读取一个字符，如何确定这个超过一半的字符？

方法是类似“擂台”的方式，设置一个位置记录当前的胜者及其出现次数，每次到来一个字符，如果位置为空，则记下该字符并将次数设置为1；如果位置不为空且两字符相同，则次数＋1；如果两字符不同，则次数-1，若减到0则位置变为空。最后剩下的字符即为答案。此方法的正确性可由反证法证明。若最后留下的不是所求字符，由于每次“擂台”都会消耗两个字符，那么可得总字符数大于n，矛盾。

如果想求出字符频率超过$\frac{n}{k+1}$的所有字符，那么只需设置k个位置，下一个字符到来时，若还有空位，则放入；若没有空位且与所有位置内字符不同，则所有位置的次数-1。此算法正确性也可同样证明。

### 4.4 Sketching

在这一部分我们讨论两件事：查重与矩阵近似。本节标题译为“画素描，概述”，研究的问题就是如何提取数据的特征。在查重的时候，我们不可能将待查重文章与数据库内所有论文逐一比对，而是用一些各自的关键字句来代替待查重文章和数据库内文章，先进行简单的比对，若重复率高，再进行详细比对。同样，在处理超大型矩阵的乘法运算时，为了节省空间，我们需要牺牲一些精确度，用某些数字或向量代替原矩阵进行运算。（这里SVD是不行的，因为求大型矩阵的奇异值太慢了）

#### 4.4.1 查重

如何提取特征？我们最先想到的是，在两个部分里各自随机取样，但是这个方法相当糟糕，因为假如有一个部分是重复的，你必须在两个部分的取样中都抽到这一部分，这会使得概率变成原来的平方的量级。

孔老师在此用了两个班A，B的学生来举例。这两个班的学生可能有所交叉，我们如何估计两个班学生的重合度？比如用

$\frac{\vertA\capB\vert}{\vertA\cupB\vert}$

来估计。这在两个班人数都不大的时候很容易做到，但是如果人数非常庞大呢？我们的方法是，在AUB上进行随机排序（Random rank），获得一个排名。我们分别取出A中的后10名、B中的后10名、总体中的后10名构成集合F(A),F(B),，C。那么我们用

$\frac{\vertF(A)\capF(B)\capC\vert}{\vertC\vert}$

来作为度量指标。这里的思想是：“如果一个人在总体排名后10名，那么如果它在一个集合中，大概率也会是后10名”。当然10这个数字可以进行调整。

### 4.5 矩阵的 sketching

这一节的任务是估计两个超大型低秩矩阵的乘积AB，其中A，B的形状分别为n×k,k×m.若$A=(\alpha_{1}\alpha_{2}\cdots\alpha_{k}),B=(\beta_{1}^{T}\beta_{2}^{T}\cdots\beta_{k}^{T})^{T},$，那么

$AB=\sum_{i=1}^{k}\alpha_{i}\beta_{i}$

很自然地，我们会想到选择一些指标，用对应的列向量、行向量相乘再相加作为近似。那么应该如何选择呢？首先我们考虑选择一个指标的情形，那么我们就需要给出一个概率向量p=$(p_{1}p_{2}\cdotsp_{k})$，分别对应选择第k个指标的概率，如果选择的结果为i，那么我们将X=$\frac{i}{p_{i}}\alpha_{i}\beta_{i}$作为估计，其期望值为

$E(X)=\sum_{j=1}^{k}p_{j}\times\frac{1}{p_{j}}\alpha_{j}\beta_{j}=AB$

虽然上述估计方法保证了期望是正确的，但是不同的概率选择会影响方差，其结果的好坏也千差万别。一个可行的方案是等概率选择，即令$p_{i}=\frac{1}{k},i=1,2,\cdots,k.$.但是直觉告诉我们这并不可行，因为很显然不同的指标对于结果的贡献应该是不同的，似乎模长更大的贡献更大一些。我们用估计值与真实值的差的F-范数平方的期望来度量性能好坏，即

$E(\vert\vertX-AB\vert\vert_{F}^{2})$

令$X=(x_{ij})_{n\timesm,AB}=(c_{ij})_{n\timesm}$，则有$E(x_{ij})=c_{ij}$，因此

$E(\vert\vertX-AB\vert\vert_{2}^{2})=E(\sum_{i=1}^{n}\sum_{j=1}^{m}(x_{ij}-c_{ij})^{2})$

$=\sum_{i=1}^{n}\sum_{j=1}^{m}E(x_{ij}^{2}-2x_{ij}c_{ij}+c_{ij}^{2})$

$=E(\vert\vertX\vert\vert_{F}^{2})-E(\vert\vertAB\vert\vert_{F}^{2})$

也就是说要让E(∥X∥E)最小。

$E(\vert\vertX\vert\vert_{F}^{2})=\sum_{i=1}^{k}p_{i}E(\vert\vert\frac{\alpha_{i}\beta_{i}}{p_{i}}\vert\vert_{F}^{2})=\sum_{i=1}^{k}\frac{1}{p_{i}}\vert\alpha_{i}\vert^$

由于$\sump_{i}=1$，根据Cauchy不等式，当

$p_{i}=\frac{\vert\alpha_{i}\vert\vert\beta_{i}\vert}{\sum_{j=1}^{k}\vert\alpha_{j}\vert\vert\beta_{j}\vert}$

时，E(∥X∥2)最小，即概率正比于两个模长乘积。（插一句：孔老师在问大家应该怎么得到最小值的时候说了一句：“You are very good at it，right？ Otherwise you won＇t be here.＂）

如果可以选择多组那么按照上述概率选择即可。

这是两个矩阵乘法的情形。如果我们有一个高阶低秩的矩阵A，想对其自身进行估计，应该怎么做？根据我们对乘法的经验，由A=AI可以对AI进行上述估计，然而此处并不太可行，因为I的秩太高了，这样做会丢失大量信息。

另一个想法是，找一个矩阵P使得A≈AP，再进行估计。我们想这样操作：从A中独立地选择一些行、列（无需一一对应），构成矩阵R,C(R是A中线性无关的r个行向量组成的矩阵），寻找矩阵U使得CUR是A的一个估计。

将$A_{n\timesm}$看作是对m维向量的一个算子，则我们希望如果$x\inR^{m}$是重要的，则APx=Ax 若x不重要，那么APx=0

如何定义“重要”和“不重要”？我们将属于R的行向量空间S中的向量称为重要，若与之正交(Rx=0)则称为不重要。那么令

$P=R^{T}(RR^{T})^{-1}R$

容易看出如果x∈S⊥，那么APx=0。。如果x∈S，则x可被R中行向量线性表出，令$x=R^{T}y$，那么

$APx=AR^{T}(RR^{T})^{-1}RR^{T}y=AR^{T}y=Ax$

我们考察一下二者的相似度。这次采用2-范数，回忆

$\vert\vertA-AP\vert\vert_{2}=\max_{\vertx\vert=1}^{ax}\vert\vert(A-AP)x\vert\vert_{2}$

对于x∈S，上式恒等于0，因此只需考虑S⊥中的向量。因此

$\max_{\vertx\vert=1}\vert\vert(A-AP)x\vert\vert_{2}=\max_{\vertx\vert=1}\vert\vertAx\vert\vert_{2}$

范数最大与范数平方最大等价，又由Rx=0,，且$R^{T}R$可作为$A^{T}A$的近似，那么

$\vert\vertAx\vert\vert_{2}^{2}=\vert\vertx^{T}A^{T}Ax\vert\vert_{2}^{2}=\vert\vertx^{T}(A^{T}A-R^{T}R)x\vert\vert_{2}^{2}\approx0$

### 5 随机图

所谓随机图，就是边的存在是有一定概率的。本节主要讨论的是G(n,p)，也就是n个顶点、每条边出现的概率是p的随机图。当p=0时为零图，当p=1时为n阶完全图。

随机图的一个特点就是，图的很多性质（如连通性）都是在p大于某一阈值之后突然出现的，而在此阈值之前出现概率几乎为零。这被称为Threshold Property。比如，当p大于某个阈值之后，随机图几乎必然联通，等等。

另外需要注意的是，对于确定的n，其某种性质出现的概率当然可以计算，我们这里说的“突然出现”，是指在n→co的意义下的。后面我们会看到这一点。

我们先考虑一个简单的问题：对于随机图G(n,p)，其中构成三角形的个数的期望是多少？这个问题并不难回答，首先取定三个点，然后乘以相应概率即可。对每组三个点我们都定义一个随机变量$I_{\Delta}$，当这三个点之间三条边都存在时为1，反之为0。即

$E(\sum_{所有\Delta}I_{\Delta})=(_{3}^{n})p^{3}\simn^{3}p^{3}$

从直观上我们可以看出，当$p=o(\frac{1}{n})$时E→0(n→∞). 当$p\ggO(\frac{1}{n})$时E几乎不可能为0。下面我们来证明这一点，分别需要使用马尔可夫不等式和切比雪夫不等式。

低于阈值时，三角形几乎不会出现，这很容易证明。令$X=\sumI_{\Delta}$代表三角形个数，由马尔可夫不等式

$P(X>a)<\frac{E(X)}{a}$

由于E(X)→0，取一个小于1的a即可。

但是当p大于阈值的时候却并不能单纯地用期望来证明。孔老师举例如下：假如班级的平均分是99分，但是实际上这个班级的同学可以得到任意高的分数，那么并不能保证随便抽取一个同学，他的得分在99分附近。相反，完全有可能是一个同学得到了几千分，其他同学都很低。不过，如果还能给定这个班级的分数方差很小，那么这种断言就可信了。我们下面证明就是从方差入手。这里给出了一种感性的认识。

用我们前面用过的技术，有

$P(X=0)\leqP(\vertX-E(X)\vert\leqE(X))\leq\frac{Var(X)}{E^{2}(X)}$

以及

$Var(X)=\sum_{\Delta}Var(I_{\Delta})+\sum_{\Delta\neq\Delta^{\prime}}Cov(I_{\Delta},I_{\Delta^{\prime}})$

对于方差，有

$Var(I_{\Delta})=E(I_{\Delta}^{2})-E^{2}(I_{\Delta})\leqE(I_{\Delta})$

对于协方差，有

$Cov(I_{\Delta},I_{\Delta^{\prime}})=E(I_{\Delta}I_{\Delta^{\prime}})-E(I_{\Delta})E(I_{\Delta^{\prime}})\leqE(I_{\Delta}I_{\Delta^{\prime}})$

将不独立的形状画出来，是四个顶点五个边构成了两个相邻的三角形，出现概率约为$n^{4}p^{5}$。因此$Var(X)\simn^{3}p^{3}+n^{4}p^{5}$。所以

$P(X=0)\leq\frac{Var(X)}{E^{2}(X)}\sim\frac{n^{3}p^{3}+n^{4}p^{5}}{(n^{3}p^{3})^{2}}\rightarrow0(n\rightarrow\infty)$

我们用相同的思想处理一下图中出现4-clique（即4阶完全子图）的阈值。首先其数量的期望为$(_{4}^{n})p^{6},$，方差的阶为$n^{4}p^{6}+n^{5}p^{9}+n^{6}p^{11}$，同样的方法可以证明。

上述讨论可以推广到任意给定的图H上。假定图形H有v个定点，e条边，那么其出现个数的均值为$\Theta(n^{v}p^{e})$，我们自然会认为$p=n^{-v/e}$是一个可能的阈值。（注意到这里v／e是图H的平均度数倒数的两倍）

那么什么时候$p=n^{-v/e}$确实是图H出现的阈值呢？我们称一个图是平衡的，当且仅当这个图的平均度数大于等于任何其子图的平均度数。由此我们可以证明：

定理5.1．如果H是平衡的，那么$p=n^{-v/e}$是它出现的阈值。

对于一般的图H，定理仍然成立，但是此时的v／e必须是H所有子图的最小值，换言之取决于最“稠密”的子图。

最后我们来研究随机图直径问题。直径是指图中顶点最小路径长度的最大值。如果图很稀疏，那么直径会很大，反之很小。对于完全图， d=1。那么如何度量d&gt;2的概率临界值呢？

容易知道，对于顶点u，v，如果ヨw使得(u,w),(w,v)是边，那么d(u,v)≤2。对于每对顶点，有另外n-2个顶点可作为桥梁，而且它们自身不能直接相连，那么不能在两条边内到达的顶点数的期望是

$E(\sum_{pairs}I_{pairs})=(_{2}^{n})(1-p)(1-p^{2})^{n-2}$

## 6 附录一：Γ函数以及部分性质

$\Gamma(z)=\int_{0}^{+\infty}t^{z-1}e^{-t}dt(Re(z)>0)$

Γ函数有如下性质：

(1)Γ(1)=1

(2)Γ(z+1)=zΓ(z)

(3)Γ(n)=(n-1)!

(4)$\Gamma(z)\Gamma(1-z)=\frac{\pi}{\sin\piz}$

性质4的证明：

$\Gamma(z)\Gamma(1-z)=\iint_{t}\sum_{0}^{0}e^{-(t+s)}(\frac{t}{s})^{x}\frac{1}{t}dsdt$

换元令ξ=t+s,η=t/s，得到

$\Gamma(z)\Gamma(1-z)=\int_{0}^{+\infty}\frac{\eta^{z-1}}{1+\eta}d\eta=\frac{\pi}{\sin\piz}$

此积分可用复变函数求解。

在性质4中取$z=\frac{1}{2}$就有$\Gamma(\frac{1}{2})=\frac{\sqrt{\pi}}{2}$，这个结果也就是高斯积分

$\int_{0}^{+\infty}e^{-x^{2}}dx=\frac{\sqrt{\pi}}{2}$

从这个结果出发，结合递推公式，我们就有

$\Gamma(\frac{2n+1}{2})=\frac{(2n-1)!!}{2^{n}}\sqrt{\pi}$

在函数定义式中换元，令$t=r^{2}$可得

$\int_{0}^{\infty}r^{p}e^{-r^{2}}dr=\frac{1}{2}\Gamma(\frac{p+1}{2})$

我们在正文中求出了高维单位球体的表面积、体积表达式

$A(n)=\frac{2(\sqrt{\pi})^{n}}{\Gamma(\frac{n}{2})}$

$V(n)=\frac{2(\sqrt{\pi})^{n}}{n\Gamma(\frac{n}{2})}$

容易看出当n→∞时，n维球体的体积趋于零。画出函数图像如下：

<!-- 5 4 3 2 1 5 10 15 20  -->
![](https://textin-image-store-1303028177.cos.ap-shanghai.myqcloud.com/external/922611f7c438d9f1)

### Figure 1：函数$\frac{2(\sqrt{\pi})^{x}}{x\Gamma(\frac{x}{2})}$的图像

可以看出五维单位球体体积最大为$\frac{8\pi^{2}}{15}\approx5.264.$

## 7 附录二：Johnson-Lindenstrauss Lemma的另一种证明

定理7.1．对于任意的$R^{d}$中的n个点构成的集合X和∀ε∈(0,1)，存在一个映射$\varphi:R^{d}\rightarrowR^{k},$这里

$k=[\frac{4\lnn}{\epsilon^{2}/2-\epsilon^{2}/3}]\leq[\frac{24}{\epsilon^{2}}\lnn]$

使得∀u,v∈X,

$(1-\epsilon)\vert\vertu-v\vert\vert_{2}^{2}\leq\vert\vert\varphi(u)-\varphi(v)\vert\vert_{2}^{2}\leq(1+\epsilon)\vert\vertu-v\vert\vert_{2}^{2}$

Proof．对于这个结果有多种证明方法，我们采用［DG99］的证明方法。证明采用了概率方法，其构造的嵌入方式非常简单：仅仅将每个点映射到一个随机的k维超平面上。更准确地说，对于一个原空间中的点v，我们令它的象$\varphi(v)=\sqrt{\frac{d}{k}}v^{\prime}$，这里v代表v在这个超平面上的投影。

为了分析这个映射的性能，我们需要考虑当映射φ随机时，随机变量

$\frac{\vert\vert\varphi(u)-\varphi(v)\vert\vert_{2}^{2}}{\vert\vertu-v\vert\vert_{2}^{2}}$

的分布。我们不妨假定$\vert\vertu-v\vert\vert_{2}^{2}=1$，也就是u-v是单位向量。注意到$\parallel\varphi(u)-\varphi(v)\parallel_{2}^{2}$就是一个固定单位向量映射到一个随机超平面之后的模的平方，等价于一个随机单位向量映射到一个固定的超平面。因此，我们可以通过在d维单位球面上随机选取一个点，将其投影到前k个基向量构成的超平面上的方法来研究。首先生成一个随机向量$X=(X_{1},X_{2},\cdots,X_{d}),X_{i}\simN(0,1)$.i.d.随后将其单位化形成单位向量$Z=\frac{1}{\vert\vertX\vert\vert_{2}}(X_{1},X_{2},\cdots,X_{d}).$.随后将其投影，我们得到了向量Y=$\frac{1}{\vert\vertX\vert\vert_{2}}(X_{1},X_{2},\cdots,X_{k}).$

我们的目标就是分析随机变量

$L=\vert\vertY\vert\vert_{2}^{2}=\frac{X_{1}^{2}+\cdots+X_{k}^{2}}{X_{1}^{2}+\cdots+X_{d}^{2}}$

的分布。注意，根据对称性我们有$\mu=E[L]=\frac{k}{d}.$.这也正是我们要在映射前乘以系数$\sqrt{\frac{d}{k}}$的原因。这里，我们要用到的关键事实是一个如下的Chernoff型的边界估计：

## 命题7.2．对于如上定义的L和μ,，我们有

1.$Pr[L<(1-\epsilon)\mu]\leq\exp(-\frac{\epsilon^{2}k}{4})$

2.$Pr[L\geq(1+\epsilon)\mu]\leq\exp(-\frac{k}{2}(\frac{z^{2}}{2}-\frac{z^{2}}{3}))$

这个命题我们将稍后证明。我们首先完成原定理的证明。

根据命题，对于$k\geq[\frac{4\lnn}{\epsilon^{2}/2-\epsilon^{2}/3}]$，我们有

$Pr[\vertL-\mu\vert\geq\epsilon\mu]\leq2\exp(-2\lnn)=\frac{2}{n^{2}}$

因此这个映射“不好”的概率，也就是

Pr[|L-μ|≥εμfor any pair$u,v]\leq\frac{2}{n^{2}}(_{2}^{n})=(1-\frac{1}{n})$

因此映射是满足要求的概率大于$\frac{1}{n}$.所以必然存在满足要求的映射。我们也可以知道大约运行O（n）轮次我们就可以得到一个满足要求的映射。 -

下面我们要证明前面的命题。

Proof.同样使用证明Chernoff Bound的方法，我们有

$Pr[L\leq(1-\epsilon)\mu]=Pr[(\frac{X_{1}^{2}+\cdots+X_{k}^{2}}{X_{1}^{2}+\cdots+X_{d}^{2}})\leq(1-\epsilon)\frac{k}{d}]$

$=Pr[k(1-\epsilon)(X_{1}^{2}+\cdots+X_{d}^{2})-d(X_{1}^{2}+\cdots+X_{k}^{2}))\geq0]$

$=Pr[\exp\{t[k(1-\epsilon)(X_{1}^{2}+\cdots+X_{a}^{2})-d(X_{1}^{2}+\cdots+X_{k}^{2})]\}\geq1]$

$\leqE[\exp\{t[k(1-\epsilon)(X_{1}^{2}+\cdots+X_{a}^{2})-d(X_{1}^{2}+\cdots+X_{k}^{2})]\}]$

$=E[\exp\{tk(1-\epsilon)X_{1}^{2}\}]^{(d-k)}正$ $[\exp\{t(k(1-\epsilon)-d)X_{1}^{2}\}]^{k}$

$=(1-2tk(1-\epsilon))^{-(d-k)/2}(1-2t(k(1-\epsilon)-d))^{-k/2}$

最后一步的等式使用了如下的事实：如果X∼N(0,1)，那么当$s<\frac{1}{2}$时有$E[e^{sX^{2}}]=(1-2s)^{-\frac{1}{2}}.$这是显然的，因为

$IE[e^{sX^{2}}]=\int_{-\infty}^{+\infty}\frac{1}{\sqrt{2\pi}}e^{-(\frac{1}{2}-s)x^{2}}dx=(1-2s)^{-\frac{1}{2}}$

此广义积分当且仅当$s<\frac{1}{2}$时收敛。

余下要做的有两点：第一，取到上式的最小值；第二，证明取到最小值的t能够满足tk(1-ε)&lt;$\frac{1}{2}$和$t(k(1-\epsilon)-d)<\frac{1}{2},$，以便保证上式这一放缩是合理的。求导数可知极值点为

$t=\frac{\epsilon}{2(1-\epsilon)(d-k(1-\epsilon))}$

并且t满足上述两个条件。

$Pr[L\leq(1-\epsilon)\mu]\leq(1-\frac{k\epsilon}{d-k(1-\epsilon)})^{-(d-k)/2}(1+\frac{\epsilon}{1-\epsilon})^{-k/2}$

$=(1+\frac{k\epsilon}{d-k})^{(d-k)/2}(1-\epsilon)^{k/2}$

$<\exp(\frac{k\epsilon}{2})(1-\epsilon)^{k/2}[\since(1+\frac{x}{y})^{y}<e^{x}]$

$=\exp(\frac{k\epsilon}{2}+\frac{k}{2}\ln(1-\epsilon))$

$<\exp(-\frac{\epsilon^{2}k}{4})$ [by Taylor expansion: I$n(1-\epsilon)<(-\epsilon-\frac{\epsilon^{2}}{2})]$

另一个不等式可以用同样的方法进行证明。唯一的不同在于需要使用$\ln(1+\epsilon)<\epsilon-\frac{\epsilon^{2}}{2}+\frac{\epsilon^{3}}{3}$ -

## 参考文献

[DG99] S. DASGUPTA and A. GuPTA, "An elementary proof of the Johnson-Lindenstrauss Lemma,"Technical Report TR-99-006,International Computer Science Institute, Berkeley,CA,1999.


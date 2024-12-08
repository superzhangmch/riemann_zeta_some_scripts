# zeta(s) 解析延拓前的收敛情况

黎曼 zeta (或 ζ) 函数，即 $\zeta(s) = \sum_{n=1}^{\infty} \frac 1 {n^s}$。

#### 1. 在 s=a+ib, a > 1, $b \in \mathbb{R}$ 时

易证它是收敛的。因为 

$$\frac 1 {n^s} = \frac 1 {n^{a+ib}} = \frac 1 {n^a} \exp^{(\log n) (-ib)} = \frac 1 {n^a} (\cos(-b \log n)+i \sin(-b \log n)).$$

#### 2. 在 a <= 1 时

按一般说法，它发散，所以才需要“解析延拓”(可延拓到只有s=1上没法定义)。单就此级数 $\sum_{n=1}^{\infty} \frac 1 {n^s}$ 而论，为啥就发散？发散则意味着 $\sum_n \frac 1 {n^a} \cos(b \log n)$, $\sum_n \frac 1 {n^a} \sin(b \log n)$ 发散，而为啥后者发散？

关于 zeta(s) 在 0 < real(s) < 1 的发散：
- 有一个[解释](https://math.stackexchange.com/questions/2700579/how-to-prove-the-divergence-of-zetas), 未察看究竟。
- 看这里：https://www.zhihu.com/question/5653967704
- 另据[here](https://math.stackexchange.com/questions/4144986/how-to-show-sum-1-ns-doesnt-converge-for-0-leq-res-leq-1?noredirect=1&lq=1), 由[Apostol的书11.6节](https://dl.icdst.org/pdfs/files1/ebc2974176a03ab93756026a97b6d370.pdf)，
> **Theorem 11.8** if the series $\sum \frac {f(n)}{n^s}$ converges for $s = \sigma_0 + i t_0$, then it also
> converges for all s with $\sigma > \sigma_0$. If it diverges for $s = \sigma_0 + i t_0$, then it
> diverges for all s with $\sigma < \sigma_0$. 【zeta(1)不收敛，从而不可能有 real(s) < 1 使得 zeta(s) 收敛】

其证明截图如下：

![image](https://github.com/user-attachments/assets/7eb7be5d-bbfc-4ad7-bf19-647a81ed12ed)


可以[程序](https://github.com/superzhangmch/riemann_zeta_some_scripts/blob/main/zeta_diverge_for_real_le1.py)验证下，确实不太收敛的:
```
appr -8.24344142357     2.44850303408
appr 25.2110006318      12.3943927543
appr -21.9209636088     -87.4027792214
appr -157.170356008     234.336644821
appr 885.569923597      -127.297062687
should_be 0.143936427077     -0.722099743529
```

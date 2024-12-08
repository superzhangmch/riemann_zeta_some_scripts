# 欧拉-麦克劳林公式计算 zeta(s)


对解析延拓前的 zeta(s), 只有 real(s) > 1 时，才可按原始级数部分和逼近求值。

然而 $\sum_n 1/n^s$, real(s) > 1 收敛较慢，可借助[欧拉-麦克劳林求和公式](https://en.wikipedia.org/wiki/Euler%E2%80%93Maclaurin_formula)——该公式擅长的正是级数求和，要么计算不收敛级数的部分和，要么加速慢收敛级数。对应到 zeta 函数，该公式形式为：

$$\sum_{k=1}^{n} \frac{1}{k^s} \sim \zeta(s) - \frac{1}{(s-1)n^{s-1}} + \frac{1}{2n^s} - \sum_{i=1}^{\infty} \frac{B_{2i}}{(2i)!} \frac{(s+2i-2)!}{(s-1)!n^{s+2i-1}}$$

从而 

$$\zeta(s) \sim \sum_{k=1}^{n} \frac{1}{k^s} + \frac{1}{(s-1)n^{s-1}} - \frac{1}{2n^s} + \sum_{i=1}^{\infty} \frac{B_{2i}}{(2i)!} \frac{(s+2i-2)!}{(s-1)!n^{s+2i-1}}$$

python [代码](https://github.com/superzhangmch/riemann_zeta_some_scripts/blob/main/zeta_using_euler_maclaurin_formula.py)验证之：

output:

```
('xx', 1)
('xx', 2)
('xx', 3)
('xx', 4)
('xx', 5)
('xx', 6) # N=100，,只用了 6 项
appr 9.162109038032135e-11      -5.754777319825002e-10
('xx', 1)
('xx', 2) # N=1000, 只用了 2 项
appr 9.16337586146096e-11       -5.75611212710788e-10
should_be   9.16161231436e-11  -5.7548264844e-10
```

其实上面的欧拉麦克劳林展开本身就是 zeta(s) 的一种解析延拓，所以不只可以用来计算 s.real > 1 的情形。

所以这是选了 zeta 函数第一个零点位置，确实算得可以。

如果算 `calc(-1+0.0000001*1j)`, N=100时得出 `-0.08333333333332066-1.654211440989384e-08*1j`, 很接近 -1/12。这正是所谓自然数之和为 -1/12 的来源。

如果算 `calc(-1+0*1j)`, mpmath给出-1/12. 欧拉法给出NaN（符合预期）。

其他：欧拉麦克劳林展开为啥不只作用于 s.real > 1:

见 [here](http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html)提到了成立条件, 另外多处看到 euler-maclaurin其实正是zeta(s)的解析延拓，参[这里](https://math.stackexchange.com/questions/3159749/for-what-real-part-of-s-as-a-function-of-q-is-the-euler-maclaurin-formula-a), [这里](https://travorlzh.github.io/2020/12/26/bernoulli-polynomials.html)。


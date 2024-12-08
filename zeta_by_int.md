# zeta(s) 解析延拓到全空间


$\zeta(s) = \sum_{n=1}^{\infty} \frac 1 {n^s},\ \ real(s) > 1$ 可以解析延拓到除了 s=1 外的整个复空间(其中积分是围着正实轴： $+\infty \stackrel {实轴上方} {\rightarrow } 环绕原点0 \stackrel {实轴下方}{\rightarrow } +\infty$)：

$$zeta(s) = \frac {\Gamma(1-s)} {2 \pi i} \int_C \frac {(-z)^s} {(e^{z}-1) \cdot z} dz,\ \ \ s != 1$$

它在 real(s) <=1 一侧，除了 s=1 之外都有定义（在 real(s) > 1一侧，s 在整数点 gamma(1-s) 无定义，但这不属于延拓部分）。不管怎样推出的，实际[python 计算](https://github.com/superzhangmch/riemann_zeta_some_scripts/blob/main/zeta_by_int.py)验证下：
```
('calc', (0.5+14.134725141j), (-1.589580733999144e-06+3.1161849872000368e-06j))
('should_be', (0.5+14.134725141j), mpc(real='9.161612314364e-11', imag='-5.754826484396e-10')) # match
__
macll.py:99: RuntimeWarning: invalid value encountered in divide
  zeta_s = gamma(1-s) / (2*cmath.pi*1j) * x
('calc', 2, (nan+nanj))  # 无定义
('should_be', 2, mpf('1.644934066848')) 
__
('calc', 2.00001, (1.6581786099366211+3.665468959566902e-10j))
('should_be', 2.00001, mpf('1.644924691471')) # 基本match
__
('should_be', 2.000001, (1.777473800519907+3.9043379294062815e-09j))
('calc', 2.000001, mpf('1.644933129297')) # 不大match。大概只是误差导致
__
('calc', (1+1j), (0.5821575680493829-0.9268478056954552j))
('should_be', (1+1j), mpc(real='0.5821580597549', imag='-0.9268485643333')) # match
__
('calc', (-2-2j), (0.08635155276513681-0.020543236862247207j))
('should_be', (-2-2j), mpc(real='0.08638207303284', imag='-0.02053604281696')) # match
__
('calc', (2+2j), (0.8673519647118549-0.2751268457031929j))
('should_be', (2+2j), mpc(real='0.8673518296346', imag='-0.2751272388086')) # match
__
('calc', (-20+20j), (7.014837841351859e+22-1.2214195900731435e+23j)) 
('should_be', (-20+20j), mpc(real='295261595076.0', imag='167640361564.0')) # bad。大概是数值积分不准的缘故
```

两者大体是匹配的。在 real(s) > 1 一侧两者相等所以它确实是原始形式的解析延拓。

不过如果 s 模长稍为过大（几十之后），则算出的差别很大，乃数值积分不准的问题。

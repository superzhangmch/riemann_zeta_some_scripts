# zeta(s) 解析延拓到 0 < real(s) <= 1

$\zeta(s) = \sum_{n=1}^{\infty} \frac 1 {n^s},\ \  real(s) > 1$  可以解析延拓到 0 < real(s) <=1 （注意 s=1 时，仍没法定义）：

$$\zeta(s) = \frac {1} {1 - 2^{1-s}} \eta(s) = \frac 1 {1 - 2^{1-s}} \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n^s}.$$

为啥：η(s) := ∑交错 = ∑奇 - ∑偶 = (∑奇+∑偶) - 2∑偶， 而 ∑奇+∑偶=ζ(s), 且 $∑偶 = \sum_n \frac 1 {(2n)^s} = 2^{-s} \sum_n \frac 1 {n^s}=2^{-s} ζ(s)$，所以 $η(s) = ζ(s)-2 \cdot 2^{-s}ζ(s)$, 从而推出上式。

在 0 < real(s) <= 1 区间，ζ(s) 也差不多是交错的，为啥 η(s) 收敛，而 ζ(s) 发散？

根据[here](https://math.stackexchange.com/questions/2188438/proving-convergence-of-the-dirichlet-eta-function),  $n^s−(n+1)^s=s∫^{n+1}_n x^{−s−1}dx$, 令 s=a+bi， c=|s| 则 $|n^s−(n+1)^s| ≤ c ∫^{n+1}_n |x^{−s−1}|dx = c ∫^{n+1}_n |x^{−a−1}|dx ≤ c n^{−a−1} = c \frac 1 {n^{a+1}}$, 而 $\sum_n \frac 1 {n^{a+1}}$ 在 a > 0时是收敛的。

[程序](https://github.com/superzhangmch/riemann_zeta_some_scripts/blob/main/zeta_by_eta.py)验证下：
output:
```
('calc', (0.5+14.134725141j), (0.000202381631061822-5.8421068595907595e-05j))
('should_be', (0.5+14.134725141j), mpc(real='9.161612314364e-11', imag='-5.754826484396e-10')) # match
__
('calc', (0.5+1j), (0.14339954678105482-0.7222224653614551j))
('should_be', (0.5+1j), mpc(real='0.1439364270773', imag='-0.7220997435288')) # match
__
('calc', 2, 1.6449340668491672)
('should_be', 2, mpf('1.644934066848')) # match
__
('calc', 2.00001, 1.6449246914661195)
('should_be', 2.00001, mpf('1.644924691471')) # match
__
('calc', 2.000001, 1.6449331293019531)
('should_be', 2.000001, mpf('1.644933129297')) # match
__
('calc', (1+1j), (0.5821574820565929-0.9268490203508495j))
('should_be', (1+1j), mpc(real='0.5821580597549', imag='-0.9268485643333')) # match
__
('calc', (2+2j), (0.867351829635518-0.275127238807934j))
('should_be', (2+2j), mpc(real='0.8673518296346', imag='-0.2751272388086')) # match
__
('calc', (20+20j), (1.000000257966888-9.180469499612744e-07j))
('should_be', (20+20j), mpc(real='1.000000257962', imag='-9.180469499534e-7')) # match
_
```
在 real(s) > 0 上，都匹配。看来确实是其延拓。

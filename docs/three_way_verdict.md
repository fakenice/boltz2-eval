> **2026-11-20 补记**：下面三行 Boltz-2 的"正确/错序"标注按当时误认的方向
> （数值大＝预测结合强）写成，方向实际相反（官方文档：`affinity_pred_value` 是 log10(IC50) μM，
> lower = stronger binding）。逐组重判见 [`方向口径更正_20261120.md`](方向口径更正_20261120.md)：
> Boltz-2 应为 MKK7 正确、CDC2L1 与 BRAF 错序（对 1 错 2），Vina 的标注不受影响。正文按原文保留。

- **MKK7 (6yz4)**: true pKd 7.52 vs 5.28; Vina -9.20 vs -7.46 -> 正确; Boltz-2 -1.63 vs -0.16 -> **错序**
- **CDC2L1 (7ukz)**: true pKd 6.54 vs 5.89; Vina -8.00 vs -9.88 -> **错序**; Boltz-2 +0.53 vs -0.05 -> 正确
- **BRAF (8c7x)**: true pKd 9.73 vs 6.43; Vina -8.78 vs -11.43 -> **错序**; Boltz-2 +0.52 vs -1.12 -> 正确

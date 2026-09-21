# 光学后处理错误报告与源头修正交接

日期：2026-09-21。
项目：pristine and hydrogen-functionalized beryllene，Thesis Project 3 / 当前 Chapter 6。
用途：交给负责源项目的 AI，复核证据后修改原后处理库、notebooks、输出图件及受影响的文字。

## 1. 结论与当前状态

已确认 pristine α-/β-beryllene 六种光学量的两套图存在以下问题：

1. β 的有效厚度计算漏乘超胞高度，导致介电归一化及其全部派生光学量改变。
2. 吸收函数需要角频率，但调用时传入普通频率，遗漏一次 2π。
3. 损失谱统一使用约 2.1 的纵轴上限，裁掉了 α 的真实高峰及部分高能端数据。

前两项是已经进入旧论文图件的后处理实现错误；第三项是显示范围问题。
这些证据不等于证明原始 VASP 计算错误，也不等于完成了全部材料和物理结论的审核。

Thesis 仓库已从原始数组重新生成这十二张 PDF，并同步修改相关正文和 SI。
源项目仓库及其公开归档尚未修改。
Thesis 中其他多体系合图仍可能包含旧处理结果，不能因这十二张图已替换而视为全部修复。
本次交接要求从原库和实际调用入口修正，避免后续运行原 notebook 又生成旧结果。

## 2. 固定来源与复现证据

源仓库：<https://github.com/Photonico/H-Beryllene_20250718>。

本报告核验的固定提交：

```text
f863ab02ebf93d5804f5552967ea750571778335
```

入口与关键文件：

- [9.0_thesis_demo.ipynb](https://github.com/Photonico/H-Beryllene_20250718/blob/f863ab02ebf93d5804f5552967ea750571778335/9.0_thesis_demo.ipynb)：体系、厚度参数、绘图范围及两套 PDF 导出。
- [vmatplot/algorithms.py](https://github.com/Photonico/H-Beryllene_20250718/blob/f863ab02ebf93d5804f5552967ea750571778335/vmatplot/algorithms.py)：能量与频率换算。
- [vmatplot/linear_optical_properties.py](https://github.com/Photonico/H-Beryllene_20250718/blob/f863ab02ebf93d5804f5552967ea750571778335/vmatplot/linear_optical_properties.py)：吸收、n、k、R、L 及绘图调用。
- [vmatplot/dielectric_function.py](https://github.com/Photonico/H-Beryllene_20250718/blob/f863ab02ebf93d5804f5552967ea750571778335/vmatplot/dielectric_function.py)：介电数据读取与绘图。
- `6.0_dielectric_selection/a-Beryllene/` 和 `6.0_dielectric_selection/b-Beryllene/`：`vaspout.h5`、`POSCAR`、`CONTCAR`、`INCAR`、`KPOINTS`；实际参数另检查了 `vasprun.xml` 的 parameters 段。

注意归档目录的大小写。Notebook 中部分字符串写为 `a-beryllene` / `b-beryllene`，核验的归档文件路径是 `a-Beryllene` / `b-Beryllene`。在大小写敏感系统上应按真实目录处理。

读取的 HDF5 datasets：

```text
results/linear_response/energies_dielectric_function
results/linear_response/density_density_dielectric_function
```

后者的 `tensor[index, index, :, 0]` 与 `tensor[index, index, :, 1]` 分别为实部、虚部；本次图件的面内和面外对应 `index=0` 的 xx 与 `index=2` 的 zz。
不要未经比较改读另一个介电数据集。

四个关键输入文件的 SHA-256：

```text
6.0_dielectric_selection/a-Beryllene/CONTCAR
e36bcb9206cd198020d09178a6bc5a3bb02a73ff113dda84f70d7d217eccd5cd

6.0_dielectric_selection/a-Beryllene/vaspout.h5
50609bd82c41b81ba4c386da69e4f3bb219ddb8b66eafb4f80375c3bc76d7d19

6.0_dielectric_selection/b-Beryllene/CONTCAR
db7f6143175489b93ac9ca4caa14ac5a1e05903ea31080528dd09a69a0ef9a02

6.0_dielectric_selection/b-Beryllene/vaspout.h5
639db6a94ff71123d7758eb37014c7f1f9961a201feebdd32d395da2e9a4dbea
```

核验顺序是先复现旧处理，再修正：

1. 原归档 PDF 与原 Thesis 对应图的身份已经核对。
2. 按旧厚度、旧频率处理原 HDF5 数组。
3. 用旧 PDF 自己的刻度标定坐标，提取其矢量曲线顶点；没有通过拟合原数组来反推坐标。
4. 六种光学量共 28 条独立曲线、6,510 个可见保留顶点与旧处理吻合，最大页面坐标残差约 0.002595 pt。
5. 修正后重新导出两套图，56 条导出曲线共 33,516 个采样顶点通过比对，最大残差约 0.003245 pt，没有高峰被裁掉。

这些残差衡量 PDF 导出与数组的一致性，不是 DFT 的物理误差估计。
旧 PDF 使用折线简化，因此旧图的严格比对针对保留顶点；被裁掉的峰由原数组与坐标上限共同确认。

## 3. P1：β 厚度漏乘超胞高度

Notebook 的原注释和实际赋值为：

```python
# System factors (d_c, (z_max-z_min)*d_c+2*r_be) r_be = 1.98
factor_a = (40, 2*1.98)
factor_b = (40, 0.0473771641975619+2*1.98)
```

β 的 `CONTCAR` 使用 `Direct` 分数坐标，第三晶格矢量为 `(0, 0, 40)` Å。
两颗 Be 原子的分数 z 坐标为：

```text
0.1236885820987794
0.0763114179012234
```

故分数跨度约为 0.0473771641975560，原子位置的实际垂直跨度约为 1.89508656790224 Å。
原代码直接把分数跨度加到 `2*r_be`，实际使用厚度约 4.00737716419756 Å。
按原注释和原半径约定，应为：

```text
d_2D = fractional_z_span * cell_height + 2 * r_be
     = 0.0473771641975560 * 40 + 2 * 1.98
     = 5.85508656790224 Å
```

这里修复的是原模型的实现。`r_be=1.98 Å` 是采用的有效厚度约定，并非唯一物理厚度。
α 的有效厚度仍为 3.96 Å。
原 SI 写的 5.85 Å 也不等于程序实际使用的 4.007377 Å；修订 SI 已写为约 5.855 Å，计算保留未舍入值。

原项目采用以下光学体积归一化：

```text
eps1_2D = 1 + (d_SC / d_2D) * (eps1_SC - 1)
eps2_2D =     (d_SC / d_2D) * eps2_SC
```

β 修正后，`eps1-1` 和 `eps2` 均为旧值的约 0.6844266293 倍。
n、k、R、L 必须由新介电函数重新计算，不能把全部派生量直接乘同一因子。
归一化必须从原始超胞介电数据开始，避免对已经归一化的数组再次处理。

源头修正要求：

- 对这两份结构，从实际结构读取超胞高度和原子位置，执行已确认的厚度约定。
- 不把 β 的 5.855086... 硬编码为其他结构的厚度。
- 其余结构逐一检查坐标模式、结构身份、周期边界和有效厚度定义；跨边界的层不能直接盲用分数坐标 `max-min`。当前两份已核实结构没有该跨边界问题。
- 本次保持原研究对 xx、zz 采用的同一光学体积归一化，不擅自把 zz 换成另一种静电逆介电模型。

## 4. P2：吸收函数混用了普通频率与角频率

固定提交的 `vmatplot/algorithms.py:149–154`：

```python
def energy_to_frequency(energy_array):
    energy_array = np.array(energy_array)
    frequency = energy_array/h_ev
    return frequency
```

其中 `h_ev=4.135667662e-15 eV s`，返回的是普通频率 `f=E/h`。
`vmatplot/linear_optical_properties.py:33–35` 的实际吸收表达式可写成：

```python
coe = (np.sqrt(2) * frequency / c_nm) * np.sqrt(
    np.sqrt(eps1**2 + eps2**2) - eps1
)
```

由于 `k=sqrt((sqrt(eps1**2+eps2**2)-eps1)/2)`，此式等于 `2*frequency*k/c_nm`。
Thesis 原有吸收公式使用角频率：

```text
alpha = 2 * omega * k / c
omega = 2*pi*f = E/hbar
```

实际主绘图调用在同文件约 560–570 行，将 `energy_to_frequency(energy_real)` 的结果直接传入吸收分支，没有补上 2π。
所以这条实际调用链缺少一个 2π；仅凭变量名叫 frequency 或 hbar 不能判定，以上结论来自函数的实际运算和调用关系。

源头修正要求：

- 为吸收入口明确普通频率与角频率的接口语义，在一处完成转换，确保恰好一次 2π。
- 可以增加明确的 `energy_to_angular_frequency`，或在吸收调用处把 `energy_to_frequency(E)` 乘以 `2*pi`。
- 不把通用 `energy_to_frequency` 静默改成返回角频率，以免破坏其他调用者。
- 同时审查 `plot_linear_optical_property`、`plot_merged_linear_optical_property`、仍被使用的 backup 分支及直接调用 `current_lop` / `comp_absorption_coefficient` 的 notebooks。旧代码在约 300、567、768、794 行存在相关频率调用；以当前源码函数和调用关系定位，行号仅作固定提交参考。
- 保留单位关系：能量 eV，`h` 为 eV s，`c_nm=2.99792458e17 nm/s`，输出吸收为 `nm^-1`。若需要 cm^-1，另显式乘 `1e7` 并同步轴标，不能只换标签。
- `linear_optical_properties.py:28` 还把 `4.135667662e-15` 命名为 `hbar`，实际数值是 h。应核对用途后修正命名或常量，避免继续混淆；这条命名问题本身不是上述缺失因子的唯一证据。

只修复频率问题、固定介电函数时，吸收在每个非零点应变为旧值的 2π 倍。
实际 β 还受 P1 影响，因此最终 β 新旧吸收比不能要求处处为 2π。
α 的介电函数不变，可以作为这一检查的直接样本。

## 5. P3：损失函数图的纵轴裁剪

`9.0_thesis_demo.ipynb` 的 Energy-loss spectrum 绘图单元原设置：

```python
lim = 2
value_boundary = (-0.05*lim, 1.05*lim)
```

纵轴上限因此为 2.1。即使用旧处理，α 也已经存在：

- 面外主峰约 8.8781，位于约 10.6372 eV。
- 面外另一个约 2.7314 的局部峰。
- 面内接近 12 eV 的窗口末端约 2.2265；这是窗口边界值，不代表已经捕获完整高能峰。

这些数值没有因 β 厚度修正而产生；原来的显示范围隐藏了它们。

源头修正要求：

- 分别设置面内、面外损失轴范围，覆盖实际数组并留余量。
- 当前 Thesis 对这两份数据采用面内 `(-0.1, 2.5)`、面外 `(-0.1, 9.5)`。
- 不把该上限机械推广到所有材料；增加导出前检查，确认绘图区内全部曲线均在轴范围内。
- 吸收轴也要随新数值更新；当前两面板均为 `(-0.0125, 0.2625) nm^-1`。
- 检查其他已有图是否用统一硬上限裁掉峰，而不只检查本 notebook。

## 6. 数值验收参考

以下为同一固定提交、相同原始 density-density 数组、原有效厚度模型下的结果。
范围是 0–12 eV，不外推到不存在的 12 eV 样本；α 和 β 最后一个可用样本分别约为 11.9819406022 和 11.9902565340 eV。
峰值来自约 0.020 eV 网格，不应解释为比网格更精确的峰位预测。

吸收最大值，单位 nm^-1，顺序为旧值 → 修正值：

- α xx：0.021562 → 0.135476。
- α zz：0.027854 → 0.175010。
- β xx：0.031402 → 0.159298。
- β zz：0.046811 → 0.241453。

除吸收外，α 的介电函数、n、k、R、L 在本次修复中数值不变；PDF 轴范围和采样点保留方式可以变化。

β 的派生量最大值，均为无量纲，顺序为 xx、zz：

- n：5.891374 → 4.902274；5.598277 → 4.662807。
- k：3.836544 → 3.149900；4.871221 → 4.000426。
- R：0.607160 → 0.537446；0.811940 → 0.766017。
- 面内 R 的最大值位置由约 9.5040 eV 变为约 1.3835 eV；派生谱的峰高排序可以改变，不能要求所有峰位和相对高低都不变。

损失函数：

- β xx 可见区局部峰约 0.3443 → 0.4490。
- β xx 深紫外区域最强局部峰约 0.6534 → 1.3335；两者不位于完全相同的能量，不应当作固定能量点的倍率。
- β xx 修正后窗口末端约 1.506761，且仍在上升；不能把它称为已经捕获的完整峰。
- β zz 主峰约 0.530866 → 1.001770。
- α zz 主峰维持约 8.878109；完整显示后明显高于 β zz 主峰。

参考数值用于定位明显实现错误。验收应使用未舍入数组和同一套常量，不要求计算结果精确等于这里已经舍入的十进制字符串。

## 7. 已确认的十二个文件与待排查范围

每行依次为 Thesis Chapter 3 文件、Thesis Project 3 文件：

1. `figures_ch3/3.1_dielectric.pdf` ↔ `figures_proj3/fig3.13_dielectric.pdf`。
2. `figures_ch3/4.1_absorption.pdf` ↔ `figures_proj3/fig3.14a_absorption.pdf`。
3. `figures_ch3/4.2_refractive.pdf` ↔ `figures_proj3/fig3.14b_refractive.pdf`。
4. `figures_ch3/4.3_extinction.pdf` ↔ `figures_proj3/fig3.15_extinction.pdf`。
5. `figures_ch3/4.4_reflectivity.pdf` ↔ `figures_proj3/S3.16_reflectivity.pdf`。
6. `figures_ch3/4.5_energy-loss.pdf` ↔ `figures_proj3/S3.17_energy-loss.pdf`。

源项目对应输出目录为 `figures_for_thesis/` 和 `figures_for_publication/`，实际 notebook 保存语句可定位每个文件。
两套图的区别主要在字体和排版，不能分别维护互不一致的计算数组。

必须继续排查，但尚未逐图验证或修复的范围包括：

- bulk / cubic / hydrogenated 对比图中重复出现的 pristine α、β 曲线。
- 其他材料是否调用同一个缺少 2π 的吸收分支。
- 其他结构是否也用分数跨度直接加原子半径，或把 full cell height 写成纯 vacuum thickness。
- 收敛测试图是否复用了错误归一化；按实际脚本核对，不能只改最终结果图。
- 合图、SI、论文稿及缓存 CSV/NPY 等中是否仍保留旧数组。
- 其他项目仓库的 `vmatplot` 副本是否有同一问题。只有核对实际版本和调用后才能扩大确认范围。

源项目建议先检索：

```sh
rg -n 'factor_[A-Za-z0-9_]*|r_be|z_max|z_min|value_boundary' --glob '*.py' --glob '*.ipynb'
rg -n 'energy_to_frequency|comp_absorption_coefficient|current_lop|plot_merged_linear_optical_property' --glob '*.py' --glob '*.ipynb'
rg -n 'hbar|h_ev|hbar_ev|frequency_to_energy' vmatplot
```

Notebook 搜索命中后应检查实际 code cells、执行顺序及导出路径，不能只依据保存的旧输出判断当前代码正确。

## 8. 同步检查的物理解读与实际计算范围

这两个已核实光学输出记录：VASP 6.5.0，`LOPTICS=T`，`ISPIN=1`，`LSORBIT=F`，`ISMEAR=0`，`SIGMA=0.01 eV`，`CSHIFT=0.1 eV`。
输入请求 `NBANDS=128`，实际输出为总能带数 133；NBANDS 不能描述为未占据能带数。
实际 `vasprun.xml` parameters 段的 `WPLASMAI=0`，WPLASMA 九个分量均为零；原后处理没有另加 Drude。
不能因为 HDF5 中缺少 WPLASMAI 字段就凭默认值推断，这里已有实际 XML 输出证据。

本次重绘保留原 IPA 带间谱范围，没有加入 Drude、电子—空穴或微观局域场修正，也没有新增 DFT 计算。
若拟讨论完整金属低频响应，需另行确定带内响应与阻尼的依据。

需要随源结果重新检查的原有说法：

- β 面内“可见峰比深紫外峰更强”：旧数组和新数组均不支持该损失谱比较。
- β 面外损失峰比 α 更显著：完整显示 α 峰后不能保留。
- β 面外 extinction 向近紫外快速下降：最强峰实际约在 5.955 eV，应据曲线重述。
- β 低能介电实部异号区：在原采样网格上约由 1.303–2.306 eV 变为 1.323–2.226 eV；精确交界需插值并说明方法，不把这些采样端点当作严格零点。
- 所有依赖旧 ε1 零点、n=1 交点、绝对吸收值或峰高排序的解释，须从修正数组重新核对。
- ε1 负值、两分量异号、n<1 或任意 loss 峰，均不能单独证明有限二维薄片中的传播 plasmon 及其寿命；数值展宽也不能直接当作真实寿命。
- 归一化后的介电函数属于选定厚度约定下的 effective response；由其导出的半无限介质标量 R，不应无条件等同于真实有限薄片的测量反射率。

这部分用于保持结果与解释一致，不要求预设追加一整套 HSE06、GW/BSE 或经验 Drude 计算。

## 9. 源头修正与验收顺序

1. 记录源项目当前 HEAD、工作区已有修改及实际运行环境；与上述固定提交对照，不能假设当前分支仍完全相同。
2. 保留原始 VASP 文件及旧图的可追溯版本。先在单独输出目录复现这两套旧曲线，确认读取的数据集和处理入口。
3. 修改源库的频率语义与吸收调用，修复 notebook 的有效厚度生成，再检查 merged / direct / backup 等实际调用入口。
4. 从原始超胞介电数组重新生成 ε、n、k、absorption、R、L；同一份计算数组供两套版式使用。
5. 逐一识别其他调用同源代码的体系和图件，给出“已验证并重生成 / 仍待核验”的明确清单。
6. 更新相关正文、SI、图注、单位及源码说明；保留作者原有叙述和 TeX 习惯，只改错误及其必要连带内容。
7. 在实际论文显示尺寸下检查完整峰、坐标、图例和图注，再向 Thesis 端交付图件与修改记录。

建议验收检查：

- 结构读出的 α、β 有效厚度分别为 3.96 Å 和约 5.85508656790224 Å。
- 同一非零能量下，角频率与普通频率之比为 2π；固定介电函数时，新旧吸收之比为 2π。
- 用独立关系 `alpha=4*pi*k/lambda_vacuum`，其中 `lambda_vacuum=h*c/E`，交叉检查吸收，避免测试与实现重复同一个错误。
- 检查 `n**2-k**2=eps1`、`2*n*k=eps2`，以及 R、L 与其明确公式的一致性。
- β 的厚度变化应通过 ε 传播到所有派生量；不要统一线性缩放 n、k、R、L。
- 全绘图区检查上下轴限，包含 α 面外约 8.8781 的 loss 峰。
- 两套导出图使用相同数值和物理坐标；从导出的 PDF 或明确数值文件抽查，而不只检查保存命令成功。
- 修正前后原始输入文件的哈希不变；明确是否保留或关闭绘图折线简化。

交回 Thesis 端时，请提供实际源代码修改、版本/commit、重生成文件清单、旧新数值核验结果，以及仍未解决的范围。
不要只返回重新命名的 PDF 或“全部已修复”的概括。

## 10. 可用的 Thesis 端复现材料

Thesis 仓库：<https://github.com/Photonico/PhD_thesis_20251216>。
以下文件已在当前本地工作区核对；本报告没有核验远程 HEAD 是否包含同一版本，接收方应核对版本或随报告取得所需文件。

- `figures_ch3/regenerate_optics.py`：从四个已固定哈希的原始文件生成六对修正图；常量和哈希内嵌，无需 `agents/`。
- `figures_ch3/verify_optics.py`：以 PDF 刻度和矢量路径校准并检查导出曲线。
- `figures_ch3/README.md`：运行方法和范围说明。
- `agents/ch3_verification/source-and-numeric-audit.json`：修正前来源、旧图复现、计算参数及旧新比较。
- `agents/ch3_verification/corrected-export-audit.json`：新图文件哈希、坐标范围和数值摘要。
- `agents/ch3_verification/corrected-pdf-vs-h5-numeric-audit.json`：新 PDF 与数组逐曲线比对。
- `agents/ch3_verification/a-vasprun-parameters-excerpt.xml` 和 `b-vasprun-parameters-excerpt.xml`：实际输出参数摘录。
- `agents/ai_application.md`：Thesis 已执行范围和剩余工作。

`agents/` 被 Git 忽略；需要时应随本报告单独传递证据文件。只收到本报告时，也可以依照固定提交、输入路径、哈希、原代码摘录和参考数值独立复核。
Thesis 的重画脚本是目前的独立复现参考，不能替代对源库和所有实际调用者的修复。

## 11. 额外代码审查线索：影响范围尚未追踪

本次重新阅读同一固定提交时，还看到 `vmatplot/algorithms.py:156–161`：

```python
def frequency_to_energy(frequency_array):
    frequency_array = np.array(frequency_array)
    energy = frequency_array/h_ev
    return energy
```

若此函数按名称接收普通频率，应使用 `E=h*f`，而不是除以 h。
这段实现与预期换算不符，但尚未确认哪些项目输出实际调用它；它不是本报告已经证明影响六对图的缺失 2π 调用链。
请单独追踪调用并检查 `frequency_to_energy(energy_to_frequency(E))` 的往返关系，再登记实际影响范围，不把它与 P2 的已确认影响混为一谈。


## 12. 图件更新与提取清单

本节按文件名列出提取清单，便于源项目修复后交回图件。
文件名中的 `3.13`、`S3.17` 等是原项目命名，不是当前 Thesis 自动生成的 Figure 编号；提取时以完整路径和文件名为准。
下面分为 A、B、C 三组，三组的核验状态不同。

### 12.1 A 组：已确认必须从修正后的源程序导出的十二张图

这是本报告已经逐曲线验证的六对图。Thesis 工作区已有修正版；源项目修改完成后，仍需从修复后的原库和 notebook 重新生成并核验，保证以后不会重新产生旧图。
两套都要交付，不能只交 publication 版或只交 thesis 版。

源项目的十二个导出文件，以下路径相对于 `H-Beryllene_20250718` 仓库根目录；可直接按此清单提取：

```text
figures_for_thesis/3.1_dielectric.pdf
figures_for_publication/fig3.13_dielectric.pdf
figures_for_thesis/4.1_absorption.pdf
figures_for_publication/fig3.14a_absorption.pdf
figures_for_thesis/4.2_refractive.pdf
figures_for_publication/fig3.14b_refractive.pdf
figures_for_thesis/4.3_extinction.pdf
figures_for_publication/fig3.15_extinction.pdf
figures_for_thesis/4.4_reflectivity.pdf
figures_for_publication/S3.16_reflectivity.pdf
figures_for_thesis/4.5_energy-loss.pdf
figures_for_publication/S3.17_energy-loss.pdf
```

对应 Thesis 的十二个目标文件，以下路径相对于 Thesis 仓库根目录：

```text
figures_ch3/3.1_dielectric.pdf
figures_proj3/fig3.13_dielectric.pdf
figures_ch3/4.1_absorption.pdf
figures_proj3/fig3.14a_absorption.pdf
figures_ch3/4.2_refractive.pdf
figures_proj3/fig3.14b_refractive.pdf
figures_ch3/4.3_extinction.pdf
figures_proj3/fig3.15_extinction.pdf
figures_ch3/4.4_reflectivity.pdf
figures_proj3/S3.16_reflectivity.pdf
figures_ch3/4.5_energy-loss.pdf
figures_proj3/S3.17_energy-loss.pdf
```

复制映射保持文件名不变：`figures_for_thesis/` 对应 `figures_ch3/`，`figures_for_publication/` 对应 `figures_proj3/`。
上面两段只是源端与目标端的对应关系，共十二张图，不是二十四张。

各对图需要包含的修正：

1. `3.1_dielectric.pdf` / `fig3.13_dielectric.pdf`：修正 β 厚度归一化后的 ε1、ε2，xx 与 zz 同步。
2. `4.1_absorption.pdf` / `fig3.14a_absorption.pdf`：α、β 均修正角频率因子，β 另修正厚度；轴标为 nm^-1，完整显示最高约 0.241453 的吸收峰。
3. `4.2_refractive.pdf` / `fig3.14b_refractive.pdf`：从修正后的 ε 重新计算 β 的 n。
4. `4.3_extinction.pdf` / `fig3.15_extinction.pdf`：从修正后的 ε 重新计算 β 的 k。
5. `4.4_reflectivity.pdf` / `S3.16_reflectivity.pdf`：从修正后的 ε、n、k 重新计算 β 的 R。
6. `4.5_energy-loss.pdf` / `S3.17_energy-loss.pdf`：重新计算 β 的 L，并恢复 α 被裁掉的高峰；面内、面外显示范围分别覆盖约 2.5 和 9.5。

`fig3.14a_absorption.pdf` 与 `fig3.14b_refractive.pdf` 虽然在 Thesis 中排为同一个 figure，仍是两个独立文件，提取时不能漏掉其中一个。
`S3.16_reflectivity.pdf` 与 `S3.17_energy-loss.pdf` 位于 SI，但也是 A 组必需文件。

### 12.2 B 组：必须追踪并同步核对的十二张多体系对比图

以下文件在当前 Thesis 中实际被引用，且图注说明包含 pristine α、β 或与氢化体系的对比。
它们尚未完成本报告 A 组那样的逐曲线核验，也尚未在 Thesis 端重生成。
应交给源项目 AI 检查生成入口；凡仍使用旧厚度、旧吸收函数或旧缓存数组的曲线，均须修正并重导出整个文件。
不要只改 A 组，留下同一材料在合图中使用另一套数值。

第一组：bulk / pristine / cubic trilayer 对比。
前三个文件分别为介电函数、吸收、损失；后三个分别为折射率、反射率、消光系数：

```text
figures_proj3/fig3.16_dielec.pdf
figures_proj3/fig3.17a_abs.pdf
figures_proj3/fig3.17b_energy-loss.pdf
figures_proj3/S3.18a_refractive.pdf
figures_proj3/S3.18b_reflectivity.pdf
figures_proj3/S3.18c_extinction.pdf
```

第二组：pristine / hydrogen-functionalized 对比。
同样依次为介电函数、吸收、损失、折射率、反射率、消光系数：

```text
figures_proj3/fig3.18_dielec_H.pdf
figures_proj3/fig3.19a_abs_H.pdf
figures_proj3/fig3.19b_energy-loss_H.pdf
figures_proj3/S3.19a_refractive_H.pdf
figures_proj3/S3.19b_reflectivity_H.pdf
figures_proj3/S3.19c_extinction_H.pdf
```

B 组列出的是 Thesis 的准确目标路径。源项目的实际 notebook、导出目录和文件名应通过保存语句确认，不能未经核实把它们全部假定为 `figures_for_publication/` 下的同名文件。
这些合图含 xx、yy、zz 分量；本次 A 组核验只有 xx、zz。处理合图时须读取实际 yy 数据，不能直接用 xx 代替 yy。
其他结构的厚度与计算设置也要分别核实，不能套用 pristine β 的 5.855086... Å。

### 12.3 C 组：十张收敛图，先核对归一化，再决定是否更新

这些图的文件在当前 Paper 和 Thesis 中存在，均须检查是否复用了错误厚度或旧后处理。
它们是介电函数收敛图，不能仅因吸收缺少 2π 就判定必须重画。

```text
figures_proj3/S3.3_dielec_hcp.pdf
figures_proj3/S3.4_dielec_hcp.pdf
figures_proj3/S3.5_dielec_alpha.pdf
figures_proj3/S3.6_dielec_alpha.pdf
figures_proj3/S3.7_dielec_beta.pdf
figures_proj3/S3.8_dielec_beta.pdf
figures_proj3/S3.9_dielec_bcc.pdf
figures_proj3/S3.10_dielec_bcc.pdf
figures_proj3/S3.11_dielec_cubic.pdf
figures_proj3/S3.12_dielec_cubic.pdf
```

每相两张分别为 NBANDS 和 k-point 收敛图：hcp、α、β、bcc、cubic。
其中 `S3.7_dielec_beta.pdf`、`S3.8_dielec_beta.pdf` 优先检查是否使用了旧 β 厚度。
如果使用原始超胞介电函数且没有错误归一化，应如实记录无需修改，而不是为了统一文件时间全部重导出。

### 12.4 提取与交付时的最小清单

- A 组十二张：源程序修复后必须重生成、核验并交付；同种物理量两套版式的数组与物理坐标应一致。
- B 组十二张：逐图追踪；确认复用错误处理的文件须更新，未核实的文件不得标记为已修复。
- C 组十张：逐图核对归一化；仅对实际受影响者更新，并给出判断依据。
- 总计三十四个需要提取或核对的现有文件，不代表已经确认三十四张全都有错误。
- 交付保留原文件名、目录映射和两套版式，并附更新清单与各图状态，避免用旧图覆盖 Thesis 中已经修正的 A 组文件。
- 本节只补充文件清单；没有据此修改 B、C 组图件，也没有改动源项目。
- 请你更新需要更新的图像, 然后这些需要提取的图像复制到一个新的文件夹中, 名字叫做 exported_figures。

# NIST SP800-22 随机性检测（Python）

这个项目提供一套 Python 版的 NIST SP800-22 随机性检测（15 项），用于对二进制序列做统计检验。

主要特点：

- 支持“单条序列详细结果”和“批量序列汇总表”两种使用方式
- 批量模式直接输出每个测试的 `通过/失败/跳过`、通过率，以及 p 值的简单统计（最小/平均/最大）
- 支持按序列长度给一套默认推荐参数（偏实用），并且会提示哪些参数属于“能跑但不一定符合标准”的情况

> 说明：SP800-22 是一套统计测试，不是“算一个分数就能证明随机”。更可靠的结论通常需要更长的序列、更多样本，再结合测试套件的整体表现一起看。

---

## 目录结构

- `opennist/`：测试代码都在这里
  - `sp800_22_suite.py`：统一入口（推荐从这里用）
  - `bitio.py`：读取输入数据（支持文件/目录，空行分段）
  - `*.py`：15 个测试实现（按功能拆成多个模块）

---

## 环境要求

- Python 3
- `numpy`（DFT 测试用）
- `scipy`（可选，建议装）：`gamma_functions.py` 会优先用 `scipy.special`
- 没有 `scipy` 的话也可以装 `mpmath`（可选）

安装示例：

```bash
pip install numpy
pip install scipy
```

---

## 包含的测试（15 项）

批量模式输出的 15 项测试名称如下：

- monobit
- frequency_within_block
- runs
- dft
- cumulative_sums
- approximate_entropy
- serial
- binary_matrix_rank
- longest_run_ones
- maurers_universal
- linear_complexity
- non_overlapping_template
- overlapping_template
- random_excursion
- random_excursion_variant

---

## 输入数据格式

### 1）bit 串（推荐）

一个文件里可以放多条序列，用空行分开，例如：

```
0101010101...（一整条）

1110001010...（第二条）
```

程序会把每一段当成一条序列来测。

### 2）目录

`key_path` 也可以给目录，会递归（可选）读取目录下的文件，每个文件同样按空行分段。

---

## 最常用：批量汇总（推荐）

用于对一批序列做汇总统计，例如：1000 条、每条 256 位，想看每项测试通过了多少。

```python
from opennist.sp800_22_suite import run_batch_recommended

summary = run_batch_recommended(
    r"E:\\data\\keys.txt",
    256,                  # 每条序列的长度（用来推荐参数）
    alpha=0.01,
    input_format="bit",   # bit / hex / auto
    print_results=True,
    max_items=1000
)
```

输出是一张表，包含：

- `通过/失败/跳过`
- `有效 = 通过 + 失败`（跳过不进分母）
- `通过率 = 通过 / 有效`
- `P最小 / P平均 / P最大`（把这个测试在所有有效样本上的 p 值做个简单统计）
- `备注`（主要放跳过原因）

表头会打印：

- `alpha=0.01（p >= alpha 算通过）`
- 本次运行用到的关键参数
- 如果用了“降级参数”，会额外打印“结果可能不符合标准”的提示

---

## 看单条/逐条细节（调试用）

用于定位某一条序列为什么 fail/skip，查看每项测试的具体 p 值：

```python
from opennist.sp800_22_suite import run_detail

rows = run_detail(
    r"E:\\data\\keys.txt",
    alpha=0.01,
    input_format="bit",
    print_results=True,
    max_items=1
)
```

---

## 自己手动配参数（可选）

推荐参数只是一套默认值，需要时可以手动覆盖：

```python
from opennist.sp800_22_suite import run_batch

summary = run_batch(
    r"E:\\data\\keys.txt",
    alpha=0.01,
    input_format="bit",
    max_items=1000,
    print_results=True,
    params={
        "serial_m": 4,                 # serial 的 m
        "freq_block_size": 32,         # block frequency 的块大小
        "nonoverlap_group": 0,         # non-overlapping 模板组
        "nonoverlap_index": 0,         # non-overlapping 模板下标
        "random_excursion_min_J": 50,  # 低于 500 属于“能跑但不标准”
    }
)
```

---

## 常见问题

### 1）为什么很多项显示“全跳过”？

有些测试对长度要求很高（比如 Maurer、Linear Complexity、Overlapping Template），当序列长度只有 256 位时，按标准条件通常做不了，就会显示“全跳过”。

### 2）`P最小/平均/最大` 到底怎么看？

批量模式不适合把 1000 个 p 值全打印出来，所以用三个数做个概览。真正的判定还是看 `通过/失败/跳过` 和 `通过率`。

---

## 入口文件

- `opennist/sp800_22_suite.py`
  - `run_batch_recommended(key_path, sequence_length, ...)`：批量 + 推荐参数（最常用）
  - `run_batch(key_path, ...)`：批量 + 手动参数
  - `run_detail(key_path, ...)`：逐条明细

# ChemSheet

> 批量化学品信息采集：输入名称/CAS，汇总 PubChem、EPA CompTox、RDKit、MolGpKa，输出一份 Excel。

把一张**物质清单**（英文名 / CAS）变成一份结构化 Excel：结构、分子描述符、理化性质、环境归趋、可选 pKa，一次查完。

数据来自 **PubChem**、**CompTox**、**RDKit**，需要时再加 **MolGpKa**。

> 适合环境科学、毒理学、化学品筛查等场景的批量物质信息整理。

---

## 功能

- 输入 `.xlsx` / `.csv` / `.txt`，自动识别分隔符（Tab / 分号 / 逗号）和编码（UTF-8 / GBK）
- 名称和 CAS 至少填一个即可；名称含逗号、空格不受影响
- 自动判断盐类：盐进 `salt_result`，结构描述符基于**最大有机片段**（去盐 + 中和）
- 一致性检查：PubChem 与 CompTox 的 InChIKey / 母体 InChIKey 对比，结果写入 `qc`
- 输出 5 个 sheet：`results`、`salt_result`、`pka_sites`（可选）、`qc`、`meta`
- xlsx 里直接内嵌结构图（`photo`）和 pKa 位点高亮图
- 查询失败自动重试；CompTox 证书异常时自动降级；输出被占用时自动改存

---

## 下载即用（Windows）

到 [Releases](https://github.com/CyMing-MX/ChemSheet_Auto/releases/latest) 下载 **`chemsheet-portable-win64.zip`**（约 54 MB）：

1. 解压到普通目录（例如 `D:\chemsheet`）；
2. 双击 `start.cmd`；
3. 把输入文件拖进窗口，按回车；
4. 输出名直接回车，等待完成。

自带 Python 3.12 + RDKit + 全部依赖，**不需要安装 Python**。

> 支持 Windows 10 / 11 64 位。macOS / Linux 请从源码运行（见下）。

---

## 从源码运行

需要 Python **3.10+**：

```bash
git clone https://github.com/CyMing-MX/ChemSheet_Auto.git
cd chemsheet
pip install -r requirements.txt
python workflow.py --input tests/sample_input.csv --out result.xlsx
```

Windows 也可以双击 `install_requirements.cmd` 安装依赖，再运行 `start.cmd`。

依赖：`rdkit`、`numpy`、`openpyxl`、`Pillow`。

---

## 命令行

```bash
python workflow.py --input 我的物质.csv --out 结果.xlsx
python workflow.py --input 我的物质.csv --out 结果.xlsx --with-molgpka
python workflow.py --input 我的物质.csv --out 结果.xlsx --with-molgpka --no-images
python workflow.py --input 我的物质.csv --out 结果.xlsx --optional
python workflow.py --input 我的物质.csv --out 结果.xlsx --format csv
```

| 参数 | 说明 |
|---|---|
| `--input` / `--out` | 输入文件 / 输出文件 |
| `--format xlsx\|csv` | 默认 xlsx |
| `--with-molgpka` | 采集 MolGpKa 位点 pKa（默认关闭） |
| `--no-images` | 采集 pKa 时不嵌入结构图 |
| `--optional` | 追加 9 个可选 CompTox 属性 |
| `--insecure` | 强制跳过 TLS 校验（默认自动检测，一般不用加） |
| `--resume` | 复用 `.cache` 中已完成的行 |

---

## 输入格式

第一行是表头，至少包含 `name` 或 `cas` 之一：

| name | cas |
|---|---|
| aspirin | 50-78-2 |
| caffeine | |
| | 64-17-5 |

- 表头大小写、空格、下划线不敏感；支持 `englishname`、`casrn`、`cas_no` 等别名
- 无表头时：2 列按 `name, cas`；1 列按 CAS 校验位判断
- 名称留空时，输出名按 PubChem Title → CompTox preferredName → CAS 回填

---

## 输出

### results（非盐）

默认 48 列（含 `photo`）；`--with-molgpka` 后 50 列；`--optional` 再加 9 列。

| 分类 | 列 |
|---|---|
| 标识 | `name, cas, iupac_name, smiles, photo` |
| 分子 | `molecular_formula, molecular_weight, exact_mass, heavy_atom_count, formal_charge` |
| 脂溶性 | `logP_pubchem, logP_comptox, logD_7.4, logD_5.5, logKoa` |
| pKa | `pka_acidic_comptox, pka_basic_comptox`（+ `pka_acid_min_molgpka, pka_base_max_molgpka`） |
| 氢键 | `hba_pubchem, hbd_pubchem, lipinski_hba, lipinski_hbd` |
| 极性与尺寸 | `tpsa, rotatable_bond_count, labute_surface_area, fraction_csp3, atom_count, heteroatom_count` |
| 环 | `ring_count` 等 8 列 |
| 溶解与挥发 | `water_solubility, vapor_pressure, henry_law_constant` |
| 热物性 | `boiling_point, melting_point, density, surface_tension` |
| 环境归趋 | `biodegradation_half_life` 等 6 列 |

### salt_result（盐）

默认 29 列（含 `photo`）；`--with-molgpka` 后增加 `parent_formula`、`parent_mw` 和两列 pKa，共 33 列。

### pka_sites（仅 `--with-molgpka`）

`name, cas, site_count, site1_no, site1_type, site1_pka, site1_atom_index, site1_image, ...`

没有可电离位点的物质也会保留一行（`site_count=0`）。

### qc / meta

`qc` 记录每个物质的解析状态、InChIKey 对比、盐/无机判断、MolGpKa 状态和错误信息；
`meta` 记录运行参数、数据源版本和单位。

---

## 数据源与限制

| 来源 | 用途 | 说明 |
|---|---|---|
| PubChem PUG REST | CID、SMILES、分子式、XLogP3、TPSA、氢键计数、CAS 同义词 | 无水印公开接口，偶发断连自动重试 |
| CompTox Dashboard API | logD、溶解度、蒸气压、沸点/熔点、环境归趋 | 证书过期时程序自动跳过校验 |
| RDKit | 去盐/中和、环计数、Lipinski、结构图 | 本地计算 |
| MolGpKa | 逐位点 pKa | 非官方公开接口，无 SLA；盐用母体片段 |

注意事项：

- CompTox 属性只在 `identity_status=OK` 且非盐时填充；
- 盐的 `molecular_formula` / `molecular_weight` 是盐整体，`parent_formula` / `parent_mw` 是中和后母体；
- MolGpKa 对季铵盐等永久电荷物质可能返回 0 位点，也可能给出可疑值，请谨慎使用；
- 公开数据仅供参考，正式使用前请核对原始来源。

---

## 目录结构

```
chemsheet/
  workflow.py            主流程
  start.py               新手向导
  start.cmd / workflow.cmd
  plugins/
    pubchem_plugin.py
    rdkit_plugin.py
    comptox_plugin.py
    molgpka_plugin.py
    _common.py
  vendor/
    pubchem_api/pubchem_client.py
    comptox_api/comptox_client.py
  tests/
    sample_input.csv
    sample_output.xlsx
  README.md
  使用说明.md
  requirements.txt
```

---

## 开发

- 编排器把任务拆成 4 个插件子进程，通过临时 CSV/JSON 交换数据；
- 每个插件可独立运行：`python plugins/pubchem_plugin.py --in rows.csv --out out.json`
- 缓存：`.cache/`（按输入文件哈希），`--resume` 可断点续跑；
- 环境变量：`ENV_PIPELINE_PYTHON`、`RDKIT_LIBS`、`PUBCHEM_API_DIR`、`COMPTOX_API_DIR`。

---

## 许可

[MIT](LICENSE)

数据分别归 PubChem / EPA CompTox / MolGpKa 所有，使用时请遵守各自的使用条款。

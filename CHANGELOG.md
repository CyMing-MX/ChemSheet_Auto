# Changelog

## 1.0.0 — 2026-10-30

首个公开版本。

- 输入支持 `.xlsx` / `.csv` / `.txt`，自动识别分隔符与 UTF-8/GBK 编码
- PubChem + CompTox + RDKit 合并输出，可选 MolGpKa 逐位点 pKa
- 盐类自动去盐/中和，母体与整体分别记录
- `qc` 表做 PubChem / CompTox 结构一致性检查
- xlsx 内嵌结构图与 pKa 位点高亮图
- 输出被占用或目录不可写时自动改存并提示
- 提供 Windows 便携版（Release 附件）

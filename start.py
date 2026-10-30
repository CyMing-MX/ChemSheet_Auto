# -*- coding: utf-8 -*-
"""新手向导：问两个问题（输入文件、输出文件），然后调用 workflow.py。"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

HELP = """\
============================================================
环境物质信息采集工具 - 帮助
============================================================

双击 start.cmd 后，只需要回答两个问题：
  1. 输入文件在哪里（把文件拖进窗口即可）
    支持 .xlsx / .csv / Tab 或逗号分隔的 .txt；名称里的逗号、空格都没问题
  2. 输出文件叫什么（直接回车 = 输入文件旁边的 result.xlsx）

其余都用默认设置：不采集 MolGpKa、不输出可选属性、证书自动处理。

如果你想改这些默认值，请用命令行（在文件夹地址栏输入 cmd 回车）：

  workflow.cmd --input 输入.csv --out 结果.xlsx
  workflow.cmd --input 输入.csv --out 结果.xlsx --with-molgpka
  workflow.cmd --input 输入.csv --out 结果.xlsx --with-molgpka --no-images
  workflow.cmd --input 输入.csv --out 结果.xlsx --optional

参数说明：
  --with-molgpka   采集 MolGpKa 位点 pKa（默认不采集）
  --no-images      采集 pKa 时不嵌入结构图
  --optional       额外输出 9 个不常用的 CompTox 属性
  --insecure       强制跳过证书校验（默认自动检测，一般不用加）
  --resume         中断后继续（复用已完成的行）

输出说明：
  results       普通物质
  salt_result   盐类
  pka_sites     pKa 位点（采集 MolGpKa 时才有）
  qc            查不到或对不上的记录
  meta          运行信息
"""


def ask(prompt, default=""):
    try:
        value = input(prompt)
    except EOFError:
        return default
    value = value.strip().strip('"').strip("'")
    return value or default


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help", "/?"):
        print(HELP)
        return 0

    print("=" * 60)
    print("            环境物质信息采集工具 - 新手向导")
    print("=" * 60)
    print()
    print("第一步：把输入文件（.xlsx / .csv / .txt）直接拖到这个窗口，然后按回车")
    print("        名称里有逗号也没关系，xlsx 或 Tab 分隔最稳")
    print()
    raw = ask("输入文件：")
    if not raw:
        print("没有输入文件，退出。")
        return 1
    src = Path(raw)
    if not src.exists():
        print(f"找不到文件：{src}")
        print("请重新运行，把文件拖进窗口后再按回车。")
        return 1

    sys.path.insert(0, str(HERE))
    try:
        from workflow import load_input
        parsed, _ = load_input(str(src))
    except Exception as exc:  # noqa: BLE001
        print(f"读取输入文件失败：{exc}")
        return 1
    if not parsed:
        print("文件里没有读到任何物质（请检查表头和内容）。")
        return 1
    print()
    print(f"读取到 {len(parsed)} 个物质：")
    for item in parsed[:5]:
        label = item["name"] or "（无名称）"
        print(f"  - {label}  {item['cas']}")
    if len(parsed) > 5:
        print(f"  ... 其余 {len(parsed) - 5} 个")
    print("如果数量或名称不对，请先检查输入文件再继续。")

    default_out = src.with_name("result.xlsx")
    print()
    print("第二步：输出文件名")
    print(f"        直接回车 = {default_out}")
    out = Path(ask("输出文件：", str(default_out)))

    cmd = [PY, str(HERE / "workflow.py"), "--input", str(src), "--out", str(out)]

    print()
    print("-" * 60)
    print("开始查询，请稍等（每个物质大约 3-5 秒）...")
    print("-" * 60)
    print()
    sys.stdout.flush()

    code = subprocess.call(cmd)
    print()
    if code == 0:
        print(f"运行完成！输出文件：{out}")
    else:
        print("运行结束，但返回了错误。请把上面的提示截图发给管理员。")
    print()
    print("按任意键退出...")
    try:
        input()
    except EOFError:
        pass
    return code


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n已取消。")
        raise SystemExit(1)

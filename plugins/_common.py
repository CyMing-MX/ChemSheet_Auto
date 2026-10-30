# -*- coding: utf-8 -*-
"""插件共用工具：路径、CSV/JSON 读写、TLS 开关。"""
import csv
import json
import os
import ssl
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

# 客户端默认用仓库内 vendor 目录，可用环境变量覆盖
PUBCHEM_DIR = Path(
    os.environ.get("PUBCHEM_API_DIR") or BASE_DIR / "vendor" / "pubchem_api"
)
COMPTOX_DIR = Path(
    os.environ.get("COMPTOX_API_DIR") or BASE_DIR / "vendor" / "comptox_api"
)

# RDKit 通过 pip 安装时不需要额外路径；非 pip 安装可用 RDKIT_LIBS 指定
RDKIT_LIBS = Path(os.environ["RDKIT_LIBS"]) if os.environ.get("RDKIT_LIBS") else None


def enable_rdkit(path=None):
    lib = path or RDKIT_LIBS
    if lib and Path(lib).exists() and str(lib) not in sys.path:
        sys.path.insert(0, str(lib))


def allow_insecure(enabled: bool):
    if enabled:
        ssl._create_default_https_context = ssl._create_unverified_context


def read_rows(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)


def read_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def fail(msg, code=1):
    print(f"错误：{msg}", file=sys.stderr)
    raise SystemExit(code)

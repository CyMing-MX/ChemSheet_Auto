@echo off
setlocal
set "PY=%~dp0runtime\python.exe"
if not exist "%PY%" set "PY=%ENV_PIPELINE_PYTHON%"
if not defined PY set "PY=python"
"%PY%" "%~dp0start.py" %*

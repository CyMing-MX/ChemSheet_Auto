@echo off
chcp 65001 >nul
setlocal
set "PY=%ENV_PIPELINE_PYTHON%"
if not defined PY set "PY=python"
echo Using Python: %PY%
"%PY%" -m pip install -r "%~dp0requirements.txt"
echo.
echo Done. Run start.cmd or workflow.cmd.
pause

@echo off
setlocal enabledelayedexpansion

:: Tu dong chuyen thu muc lam viec ve goc repo (parent directory cua thu muc scripts)
cd /d "%~dp0.."
set PYTHONPATH=%CD%\src;%PYTHONPATH%

echo =====================================================================
echo CRISPR-Cas12a Explainable AI Project
echo Step: Model Validation & Literature Benchmarking
echo Author: Minh Tran (Independent Researcher, Perth, Western Australia, Australia)
echo =====================================================================
echo [*] Thu muc lam viec: %CD%
echo.

:: ====================================================================
:: GIAI DOAN 1: Uu tien tim moi truong Python 3 da cai san pandas & scipy
:: ====================================================================
set PYTHON_CMD=

if exist ".venv\Scripts\python.exe" (
    .venv\Scripts\python.exe -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
    if !ERRORLEVEL! equ 0 (
        set PYTHON_CMD=.venv\Scripts\python.exe
        goto :found_python
    )
)
if exist "venv\Scripts\python.exe" (
    venv\Scripts\python.exe -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
    if !ERRORLEVEL! equ 0 (
        set PYTHON_CMD=venv\Scripts\python.exe
        goto :found_python
    )
)

py -3 -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=py -3
    goto :found_python
)

python3 -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=python3
    goto :found_python
)

python -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=python
    goto :found_python
)

py -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=py
    goto :found_python
)

for %%P in (
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python39\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python38\python.exe"
    "C:\Program Files\Python312\python.exe"
    "C:\Program Files\Python311\python.exe"
    "C:\Program Files\Python310\python.exe"
    "C:\Program Files\Python39\python.exe"
    "C:\Program Files\Python38\python.exe"
    "C:\ProgramData\anaconda3\python.exe"
    "C:\ProgramData\miniconda3\python.exe"
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
) do (
    if exist "%%~P" (
        "%%~P" -c "import sys, pandas; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
        if !ERRORLEVEL! equ 0 (
            set PYTHON_CMD="%%~P"
            goto :found_python
        )
    )
)

:: ====================================================================
:: GIAI DOAN 2: Neu khong co moi truong nao co san pandas,
:: tim bat ky trinh thong dich Python 3.6+ hop le nao de tu dong cai dat
:: ====================================================================
if exist ".venv\Scripts\python.exe" (
    .venv\Scripts\python.exe -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
    if !ERRORLEVEL! equ 0 (
        set PYTHON_CMD=.venv\Scripts\python.exe
        goto :install_deps
    )
)
if exist "venv\Scripts\python.exe" (
    venv\Scripts\python.exe -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
    if !ERRORLEVEL! equ 0 (
        set PYTHON_CMD=venv\Scripts\python.exe
        goto :install_deps
    )
)

py -3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=py -3
    goto :install_deps
)

python3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=python3
    goto :install_deps
)

python -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=python
    goto :install_deps
)

py -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    set PYTHON_CMD=py
    goto :install_deps
)

for %%P in (
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python39\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python38\python.exe"
    "C:\Program Files\Python312\python.exe"
    "C:\Program Files\Python311\python.exe"
    "C:\Program Files\Python310\python.exe"
    "C:\Program Files\Python39\python.exe"
    "C:\Program Files\Python38\python.exe"
    "C:\ProgramData\anaconda3\python.exe"
    "C:\ProgramData\miniconda3\python.exe"
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
) do (
    if exist "%%~P" (
        "%%~P" -c "import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)" >nul 2>&1
        if !ERRORLEVEL! equ 0 (
            set PYTHON_CMD="%%~P"
            goto :install_deps
        )
    )
)

echo.
echo [!] LOI NGHIEP TRONG: Khong tim thay Python 3.6+ tren he thong!
echo     Lenh 'python' hien tai dang goi phien ban cu (Python 2 hoac khong ho tro f-string).
echo     Vui long:
echo       1. Cai dat Python 3.8+ tu https://www.python.org/
echo       2. Hoac tao va kich hoat virtual environment:
echo          py -3 -m venv .venv
echo          .venv\Scripts\pip install -r requirements.txt
echo.
pause
exit /b 1

:install_deps
echo [*] Phat hien trinh thong dich Python 3: %PYTHON_CMD%
echo [!] Moi truong nay chua cai dat cac thu vien can thiet (chua co pandas, numpy, scikit-learn,...).
echo [*] Dang tien hanh tu dong cai dat thu vien tu requirements.txt...
echo.
%PYTHON_CMD% -m pip install -r requirements.txt
if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [!] Khong the tu dong cai dat thu vien qua pip.
    echo     Vui long mo Terminal/Command Prompt va chay thu cong lenh:
    echo     %PYTHON_CMD% -m pip install -r requirements.txt
    echo.
    pause
    exit /b 1
)
echo.
echo [OK] Da cai dat xong cac goi phu thuoc thanh cong!
echo.
goto :run_steps

:found_python
echo [*] Su dung trinh thong dich: %PYTHON_CMD%
echo.

:run_steps

%PYTHON_CMD% src\05_model_validation_benchmark.py

echo.
echo =====================================================================
echo Validation complete! Check:
echo   - manuscript\figures\Figure2_Model_Benchmark_Validation.png
echo   - data\results\benchmark_metrics_summary.json
echo   - data\results\model_comparison_benchmark.csv
echo =====================================================================
pause

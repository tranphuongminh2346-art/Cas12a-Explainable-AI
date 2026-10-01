@echo off
setlocal enabledelayedexpansion

:: Tu dong chuyen thu muc lam viec ve goc repo (parent directory cua thu muc scripts)
cd /d "%~dp0.."
set PYTHONPATH=%CD%\src;%PYTHONPATH%

echo ====================================================================
echo        CRISPR-Cas12a Explainable AI Machine Learning Pipeline
echo        Author: Minh Tran (Independent Researcher, Perth, Western Australia, Australia)
echo ====================================================================
echo [*] Thu muc lam viec: %CD%

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

echo [*] BUOC 1: Tien xu ly du lieu va Data Quality Control...
%PYTHON_CMD% src\01_preprocess_data.py
if %ERRORLEVEL% NEQ 0 (
    echo [!] Loi o Buoc 1. Vui long kiem tra lai file du lieu raw.
    pause
    exit /b %ERRORLEVEL%
)
echo.

echo [*] BUOC 1b: Kham pha du lieu (Exploratory Data Analysis - EDA) va tao Figure 1...
%PYTHON_CMD% src\01b_exploratory_data_analysis.py
if %ERRORLEVEL% NEQ 0 (
    echo [!] Loi o Buoc 1b.
    pause
    exit /b %ERRORLEVEL%
)
echo.

echo [*] BUOC 2: Trich xuat dac trung Sinh hoc va Nhiet dong hoc (Thermodynamics)...
%PYTHON_CMD% src\02_extract_features.py
if %ERRORLEVEL% NEQ 0 (
    echo [!] Loi o Buoc 2.
    pause
    exit /b %ERRORLEVEL%
)
echo.

echo [*] BUOC 3: Kiem tra mo hinh va du lieu du doan 10-Fold Cross-Validation...
if "%1"=="--retrain" (
    echo [*] Che do bat buoc huan luyen lai (--retrain duoc kich hoat)...
    %PYTHON_CMD% src\03_train_models.py
    if !ERRORLEVEL! NEQ 0 (
        echo [!] Loi o Buoc 3 khi huan luyen model.
        pause
        exit /b !ERRORLEVEL!
    )
) else (
    if exist "data\results\cas12a_predictions_cv.csv" (
        echo [+] Phat hien ket qua 10-Fold CV da duoc tinh toan san va kiem toan chat che:
        echo     data\results\cas12a_predictions_cv.csv (N = 11,365 canonical targets).
        echo [*] Bo qua buoc huan luyen local de tiet kiem thoi gian va tai nguyen CPU/RAM.
        echo     (Neu muon huan luyen lai toan bo tu dau, khuyen nghi dung Google Colab:
        echo      notebooks\Cas12a_Complete_Master_Pipeline.ipynb, hoac chay: scripts\run_pipeline.bat --retrain).
    ) else (
        echo [*] Khong tim thay ket qua co san, tien hanh huan luyen mo hinh 10-Fold CV...
        %PYTHON_CMD% src\03_train_models.py
        if !ERRORLEVEL! NEQ 0 (
            echo [!] Loi o Buoc 3 khi huan luyen model.
            pause
            exit /b !ERRORLEVEL!
        )
    )
)
echo.

echo [*] BUOC 4: Kiem tra do thi SHAP va giai thich co che sinh hoc...
if "%1"=="--retrain" (
    echo [*] Che do ve lai toan bo do thi SHAP...
    %PYTHON_CMD% src\04_plot_shap_figures.py
    if !ERRORLEVEL! NEQ 0 (
        echo [!] Loi o Buoc 4 khi ve do thi SHAP.
        pause
        exit /b !ERRORLEVEL!
    )
) else (
    if exist "manuscript\figures\Figure3A_SHAP_Summary_Beeswarm.png" (
        echo [+] Phat hien cac hinh ve TreeSHAP va Epistasis 300 DPI da co san:
        echo     - Figure 3A: TreeSHAP Global Beeswarm Summary
        echo     - Figure 3B: PAM x Seed Epistatic Interaction Matrix
        echo     - Figure 4:  Thermodynamic Interplay Landscape
        echo     - Figure S2: Outlier Error Analysis
        echo [*] Giu nguyen hinh ve chuan xuat ban da qua kiem dinh.
    ) else (
        echo [*] Tien hanh ve do thi SHAP...
        %PYTHON_CMD% src\04_plot_shap_figures.py
        if !ERRORLEVEL! NEQ 0 (
            echo [!] Loi o Buoc 4 khi ve do thi SHAP.
            pause
            exit /b !ERRORLEVEL!
        )
    )
)
echo.

echo [*] BUOC 5: Danh gia Benchmark va Tao Figure 2 So Sanh Literature (300 DPI)...
%PYTHON_CMD% src\05_model_validation_benchmark.py
if %ERRORLEVEL% NEQ 0 (
    echo [!] Loi o Buoc 5.
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ====================================================================
echo [OK] TOAN BO PIPELINE DA CHAY XONG THANH CONG!
echo     - Du lieu sach:  data\processed\clean_kim_2018_targets.csv
echo     - Ma tran dac trung: data\processed\cas12a_features_master.csv
echo     - Ket qua danh gia:  data\results\cas12a_predictions_cv.csv
echo     - Hinh ve bai bao:   manuscript\figures\
echo ====================================================================
pause

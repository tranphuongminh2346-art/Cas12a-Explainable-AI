@echo off
title Don dep va sap xep thu muc du an Cas12a-XAI
echo =====================================================================
echo TIEN HANH DON DEP VA SAP XEP LAI TOAN BO THU MUC DU AN...
echo Author: Minh Tran (UWA / Independent Researcher)
echo =====================================================================
echo.

:: 1. Tao thu muc scripts neu chua co
if not exist "scripts" mkdir "scripts"

:: 2. Di chuyen cac script tien ich vao thu muc scripts\
echo [*] Dang gom cac file script vao thu muc scripts\...
if exist "download_skills.bat" move /Y "download_skills.bat" "scripts\" >nul
if exist "download_skills.py" move /Y "download_skills.py" "scripts\" >nul
if exist "download_skills.ps1" move /Y "download_skills.ps1" "scripts\" >nul
if exist "download_skills_pure.ps1" move /Y "download_skills_pure.ps1" "scripts\" >nul
if exist "run_pipeline.bat" move /Y "run_pipeline.bat" "scripts\" >nul
if exist "run_validation.bat" move /Y "run_validation.bat" "scripts\" >nul
if exist "copy_figures_to_manuscript.bat" move /Y "copy_figures_to_manuscript.bat" "scripts\" >nul
if exist "organize_files.bat" move /Y "organize_files.bat" "scripts\" >nul
echo [v] Da don dep cac file script vao scripts\

:: 3. Sap xep cac file ket qua vao data\results\
echo [*] Dang chuyen cac file ket qua phan tich vao data\results\...
if exist "outlier_error_analysis.csv" (
    move /Y "outlier_error_analysis.csv" "data\results\" >nul
    echo [v] Da chuyen: outlier_error_analysis.csv -^> data\results\
)
if exist "pam_stratified_validation.csv" (
    move /Y "pam_stratified_validation.csv" "data\results\" >nul
    echo [v] Da chuyen: pam_stratified_validation.csv -^> data\results\
)
if exist "cas12a_predictions_cv.csv" (
    if exist "data\results\cas12a_predictions_cv.csv" (
        del /F /Q "cas12a_predictions_cv.csv"
        echo [v] Da xoa file trung lap cas12a_predictions_cv.csv o thu muc goc (da co san trong data\results\)
    ) else (
        move /Y "cas12a_predictions_cv.csv" "data\results\" >nul
        echo [v] Da chuyen: cas12a_predictions_cv.csv -^> data\results\
    )
)

:: 4. Xoa cac file anh trung lap o thu muc goc (neu con)
if exist "Figure3A_SHAP_Summary_Beeswarm.png" del /F /Q "Figure3A_SHAP_Summary_Beeswarm.png"
if exist "Figure4_Thermodynamic_Interplay.png" del /F /Q "Figure4_Thermodynamic_Interplay.png"

echo.
echo =====================================================================
echo [HOAN TAT] Thu muc du an da duoc sap xep 100%% chuan muc!
echo Thu muc goc hien chi chua cac file chuan cua GitHub repository:
echo   - README.md, LICENSE, setup.py, pyproject.toml, requirements.txt
echo   - manuscript\ (Ban thao, figures 300 DPI, supplementary S1-S5)
echo   - data\results\ (Du lieu ket qua phan tich & out-of-fold)
echo   - notebooks\ (3 Jupyter Notebooks chay Colab)
echo   - src\ (Ma nguon goi cas12a_xai)
echo   - scripts\ (Cac script tien ich va installer)
echo =====================================================================
pause

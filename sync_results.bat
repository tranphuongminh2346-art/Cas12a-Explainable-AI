@echo off
echo =====================================================================
echo Syncing 10-Fold Results into Manuscript and Data Folders...
echo =====================================================================

if not exist "cas12a_10fold_results_package" (
    echo [ERROR] cas12a_10fold_results_package folder not found!
    pause
    exit /b 1
)

echo [+] Copying high-resolution 300 DPI figures to manuscript\figures\...
copy /Y "cas12a_10fold_results_package\figures\Figure2_Model_Benchmark_Validation.png" "manuscript\figures\"
copy /Y "cas12a_10fold_results_package\figures\Figure3A_SHAP_Summary_Beeswarm.png" "manuscript\figures\"
copy /Y "cas12a_10fold_results_package\figures\Figure3B_SHAP_Interaction_Matrix.png" "manuscript\figures\"
copy /Y "cas12a_10fold_results_package\figures\Figure_S2_Outlier_Error_Analysis.png" "manuscript\figures\"

echo [+] Copying predictions and validation tables to data\results\...
copy /Y "cas12a_10fold_results_package\results\cas12a_predictions_10fold.csv" "data\results\"
copy /Y "cas12a_10fold_results_package\results\cas12a_predictions_10fold.csv" "data\results\cas12a_predictions_cv.csv"
copy /Y "cas12a_10fold_results_package\results\pam_stratified_validation.csv" "data\results\"
copy /Y "cas12a_10fold_results_package\results\outlier_error_analysis.csv" "data\results\"

echo [+] Copying Supplementary Tables S3 and S5 to manuscript\supplementary\...
copy /Y "cas12a_10fold_results_package\results\pam_stratified_validation.csv" "manuscript\supplementary\Supplementary_Table_S3_PAM_Stratified_Benchmark.csv"
copy /Y "cas12a_10fold_results_package\results\Supplementary_Table_S5_Ablation_Study.csv" "manuscript\supplementary\"

echo.
echo =====================================================================
echo [SUCCESS] All 10-fold figures, tables, and predictions synchronized!
echo =====================================================================
pause

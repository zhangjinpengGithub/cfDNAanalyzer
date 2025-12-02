"""
This script performs single-modality machine learning analysis for classification.

1. Parse arguments (input folder, FS method, classifier list, CV type, etc.)
2. For each input file in the input_dir:
    a. Load data using load_data()
    b. Build feature selection combinations (filter, wrapper, embedded, hybrid)
    c. For each FS combination:
        d. For each classifier:
            i. Initialize CV splitter
            ii. Record runtime and memory usage
            iii. For each fold:
                - Perform feature selection on training set
                - Apply same features on test set
                - Train classifier and predict probabilities
            iv. Collect prediction results
            v. Evaluate classification metrics (accuracy, F1, AUC, etc.)
            vi. Save confusion matrix and ROC (if binary classification)
3. Output prediction results and metrics into modality-specific folders under output_dir

Usage example:

nohup python run_single_modality.py \
    --modality single \
    --input_dir /home/zky/projects/202312_cfDNA_integrated_tools/2504DA-npy/testdata \
    --DA_output_dir /home/zky/projects/202312_cfDNA_integrated_tools/2504DA-npy/testdata_output \
    --modality single \
    --classNum 2 \
    --cvSingle LOO \
    --classifierSingle KNN SVM \
    --filterMethod IG CHI FS FCBF CC LVF MAD DR MI RLF SURF MSURF \
    --filterFrac 0.023618328 \
    --wrapperMethod BOR \
    --wrapperFrac 0.023618328 \
    --embeddedMethod LASSO RIDGE ELASTICNET RF \
    --embeddedFrac 0.023618328 \
    --hybridType FE \
    --hybridMethod1 IG \
    --hybridFrac1 0.023618328 \
    --hybridMethod2 LASSO \
    --hybridFrac2 0.023618328 > tst_run_single_modality.log &

nohup python run_single_modality.py \
    --modality single \
    --input_dir /home/zjp/projects/202402_cfDNAIntegratedTool/250402_downstream_result/two-class/01preprocessing/0317_format_var_std/CNA \
    --classNum 2 \
    --cvSingle KFold \
    --nsplitSingle 5 \
    --classifierSingle SVM KNN \
    --filterMethod IG CHI \
    --filterFrac 0.2  \
    --explain SHAP perm  \
    --DA_output_dir two_class_output \
    > run_single_modality11.log 2>&1 &
"""

import os
import pandas as pd
import numpy as np
import time
import tracemalloc
import psutil
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from config_parser import get_args
from methods.machine_learning.train_test_split import get_cv
from methods.machine_learning.classifiers import get_classifier
from methods.machine_learning.generate_metrics import evaluate_classification
from methods.machine_learning.model_explain import compute_model_explanations
from run_feature_selection import run_feature_selection, infer_fs_type, build_fs_combinations_from_args
import glob

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"


def load_single_modalities_csv(input_dir):
    files = sorted(glob.glob(os.path.join(input_dir, '*.csv')))
    datasets = {}
    for f in files:
        name = os.path.splitext(os.path.basename(f))[0]
        df = pd.read_csv(f)
        df.columns = df.columns.str.lower()
        X = df.drop(columns=["label", "sample"], errors="ignore")
        X.columns = X.columns.str.replace('.csv', '', regex=False)
        y = df["label"]
        y = y.astype('category').cat.codes
        sample_ids = df["sample"].values if "sample" in df.columns else df.index.to_numpy()
        datasets[name] = (X, y, sample_ids)
    return datasets


def load_single_modalities_npy(input_dir):
    files = sorted(glob.glob(os.path.join(input_dir, '*.npz')))
    datasets = {}
    for f in files:
        name = os.path.splitext(os.path.basename(f))[0]
        npz = np.load(f, allow_pickle=True)
        X = pd.DataFrame(npz["data"])
        y = pd.Series(npz["label"])
        sample_ids = npz["sample_ids"]
        datasets[name] = (X, y, sample_ids)
    return datasets


def smart_load_modalities(input_dir):
    npy_files = glob.glob(os.path.join(input_dir, '*.npz'))
    if npy_files:
        print(f"[INFO] Detected cached .npz files in {input_dir}, using them.")
        return load_single_modalities_npy(input_dir)
    else:
        print(f"[INFO] No cached .npz found, loading from .csv files.")
        return load_single_modalities_csv(input_dir)


def process_fold(fold_idx, train_idx, test_idx, X, y, sample_ids, args, clf_name, output_dir_combo):
    clf = get_classifier(clf_name)
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    fs_type = infer_fs_type(args)
    idx_list = run_feature_selection(
        X_train, y_train, args,
        fs_type=fs_type,
        verbose=False,
        fs_tag=f"{clf_name}_fold{fold_idx}"
    )
    #if args.filterFrac and isinstance(args.filterFrac, float) and args.filterFrac >= 1.0:
    #    idx_list = list(range(X_train.shape[1]))
    #elif len(idx_list) == 0:
    #    print(f"[WARN] FS returned empty list for {clf_name}_fold{fold_idx}, fallback to all features")
    #    idx_list = list(range(X_train.shape[1]))
    X_train_fs = X_train.iloc[:, idx_list]
    X_test_fs = X_test.iloc[:, idx_list]
    clf.fit(X_train_fs, y_train)

    # ---- Model explain ----
    if getattr(args, "explain", None):
        global_model_explain_dir = os.path.join(
            os.path.dirname(output_dir_combo), "model_explain", os.path.basename(output_dir_combo)
        )
        model_explain_dir = os.path.join(global_model_explain_dir, clf_name, f"fold_{fold_idx}")
        os.makedirs(model_explain_dir, exist_ok=True)
        compute_model_explanations(
            clf,
            X_train_fs,
            y_train,
            X_test_fs,
            y_test,
            output_dir=model_explain_dir,
            explain_methods=args.explain
        )

    if hasattr(clf, "predict_proba"):
        y_probs = clf.predict_proba(X_test_fs)
    else:
        pred_labels = clf.predict(X_test_fs)
        y_probs = []
        for label in pred_labels:
            prob = [0.0] * args.classNum
            prob[int(label)] = 1.0
            y_probs.append(prob)

    rows = []
    for i, idx in enumerate(test_idx):
        sample_id = sample_ids[idx]
        true_label = y.iloc[idx]
        prob = y_probs[i]
        row = {'SampleID': sample_id, 'TrueLabel': true_label, 'Classifier': clf_name}
        for cls in range(args.classNum):
            row[f'Prob_Class{cls}'] = prob[cls]
        rows.append(row)
    return rows


def run_single_modality_analysis(datasets, args):
    fs_combinations = build_fs_combinations_from_args(args)
    if not fs_combinations:
        raise ValueError("No valid FS parameters provided.")
    for name, (X, y, sample_ids) in datasets.items():
        print(f"\n===== Processing single modality: {name} =====")     
        cv = get_cv(method=args.cvSingle, nsplit=args.nsplitSingle, y=y, test_ratio=args.cvSingle_test_ratio)
        output_base_dir = os.path.join(args.DA_output_dir, name)
        os.makedirs(output_base_dir, exist_ok=True)

        all_predictions, all_metrics = [], []

        for fs_combo in fs_combinations:
            fs_label = fs_combo.get("fs_label", "default")
            output_dir_combo = os.path.join(output_base_dir, fs_label)

            for clf_name in args.classifierSingle:
                print(f"\n>>> Running classifier: {clf_name} with FS combo: {fs_label}")
                start_time = time.time()
                tracemalloc.start()
                process = psutil.Process()

                fold_results = []
                with ThreadPoolExecutor() as executor:
                    futures = [executor.submit(process_fold, i, tr, te, X, y, sample_ids, args, clf_name, output_dir_combo)
                               for i, (tr, te) in enumerate(cv.split(X, y))]
                    for fut in tqdm(as_completed(futures), total=len(futures)):
                        fold_results.extend(fut.result())

                total_time = round(time.time() - start_time, 4)
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                peak_mem = round(peak / 1024 / 1024, 4)

                prediction_df = pd.DataFrame(fold_results)
                prediction_df["Classifier"] = clf_name
                prediction_df["FS_Combination"] = fs_label
                #all_predictions.append(prediction_df)
                
                prob_cols = [col for col in prediction_df.columns if col.startswith("Prob_Class")]
                for col in prob_cols:
                    if col in prediction_df.columns:
                        prediction_df[col] = prediction_df[col].round(6)
                # ordered_cols = ['SampleID', 'TrueLabel', 'Classifier', 'FS_Combination']
                ordered_cols = ['SampleID', 'TrueLabel', 'FS_Combination', 'Classifier']
                
                prediction_df = prediction_df[ordered_cols + prob_cols]
                all_predictions.append(prediction_df)

                metrics = evaluate_classification(prediction_df, class_num=args.classNum)
                metrics["Classifier"] = clf_name
                metrics["FS_Combination"] = fs_label
                metrics["TotalTime_sec"] = total_time
                # metrics["FinalMemory_MB"] = mem_used
                metrics["PeakMemory_MB"] = peak_mem
                all_metrics.append(metrics)

        # === Save metrics & probabilities ===
        final_pred_df = pd.concat(all_predictions, axis=0)
        final_pred_path = os.path.join(output_base_dir, "single_modality_probabilities.csv")
        final_pred_df.to_csv(final_pred_path, index=False)
        print(f"Saved prediction results to: {final_pred_path}")

        metrics_df = pd.DataFrame(all_metrics)
        # metric_order = ['Classifier', 'FS_Combination', 'TotalTime_sec', 'FinalMemory_MB', 'PeakMemory_MB']
        metric_order = ['FS_Combination', 'Classifier', 'TotalTime_sec', 'PeakMemory_MB']
        
        other_metrics = [col for col in metrics_df.columns if col not in metric_order]
        metrics_df = metrics_df[metric_order + other_metrics]
        metrics_path = os.path.join(output_base_dir, "single_modality_metrics.csv")
        metrics_df.to_csv(metrics_path, index=False)
        print(f"Saved performance summary to: {metrics_path}")


        # === Merge model_explain summaries ===
        explain_base = os.path.join(output_base_dir, "model_explain")
        merged_rows = []

        for fs_combo in fs_combinations:
            fs_label = fs_combo.get("fs_label", "default")
            fs_dir = os.path.join(explain_base, fs_label)
            if not os.path.exists(fs_dir):
                continue

            for clf_name in args.classifierSingle:
                clf_dir = os.path.join(fs_dir, clf_name)
                fold_dirs = sorted(glob.glob(os.path.join(clf_dir, "fold_*")))
                if not fold_dirs:
                    continue

                shap_fold_tables = []
                perm_fold_tables = []

                for fold_path in fold_dirs:
                    fold_name = os.path.basename(fold_path)
                    fold_idx = fold_name.replace("fold_", "")

                    shap_file = os.path.join(fold_path, "shap_values.csv")
                    perm_file = os.path.join(fold_path, "permutation_importance.csv")

                    # ---- SHAP ----
                    if os.path.exists(shap_file):
                        df = pd.read_csv(shap_file)

                        if "Feature" in df.columns:
                            # long format
                            df = df.set_index("Feature").iloc[:, 0]
                        else:
                            # wide format: columns = features, rows = samples
                            df = df.mean(axis=0)
                            df.index.name = "Feature"

                        df = df.rename(f"SHAP_Fold{fold_idx}")
                        shap_fold_tables.append(df)

                    # ---- PERM ----
                    if os.path.exists(perm_file):
                        df = pd.read_csv(perm_file)

                        if "Feature" in df.columns and "ImportanceMean" in df.columns:
                            df = df.set_index("Feature")["ImportanceMean"]
                        else:
                            raise ValueError(f"Permutation file format invalid: {perm_file}")

                        df = df.rename(f"Perm_Fold{fold_idx}")
                        perm_fold_tables.append(df)

                # === Aggregate SHAP ===
                if shap_fold_tables:
                    shap_all = pd.concat(shap_fold_tables, axis=1)
                    shap_mean = shap_all.abs().mean(axis=1).rename("MeanAbsSHAP")
                    shap_std = shap_all.abs().std(axis=1).rename("MeanStdSHAP")
                    shap_summary = pd.concat([shap_mean, shap_std, shap_all], axis=1).reset_index()
                else:
                    shap_summary = pd.DataFrame(columns=["Feature", "MeanAbsSHAP", "MeanStdSHAP"])

                # === Aggregate PERM ===
                if perm_fold_tables:
                    perm_all = pd.concat(perm_fold_tables, axis=1)
                    perm_mean = perm_all.abs().mean(axis=1).rename("MeanAbsPermImportance")
                    perm_std = perm_all.abs().std(axis=1).rename("MeanStdPermImportance")
                    perm_summary = pd.concat([perm_mean, perm_std, perm_all], axis=1).reset_index()
                else:
                    perm_summary = pd.DataFrame(columns=["Feature", "MeanAbsPermImportance", "MeanStdPermImportance"])

                # === Merge SHAP + Perm ===
                merged = pd.merge(shap_summary, perm_summary, on="Feature", how="outer")
                merged.insert(0, "Classifier", clf_name)
                merged.insert(0, "FS_Combination", fs_label)

                # === reorder columns ===
                fold_cols_shap = sorted([c for c in merged.columns if c.startswith("SHAP_Fold")],
                                        key=lambda x: int(x.replace("SHAP_Fold", "")))
                fold_cols_perm = sorted([c for c in merged.columns if c.startswith("Perm_Fold")],
                                        key=lambda x: int(x.replace("Perm_Fold", "")))

                base_cols = [
                    "FS_Combination", "Classifier", "Feature",
                    "MeanAbsSHAP", "MeanStdSHAP",
                    "MeanAbsPermImportance", "MeanStdPermImportance"
                ]

                merged = merged[base_cols + fold_cols_shap + fold_cols_perm]
                merged_rows.append(merged)

        # === Save and cleanup ===
        if merged_rows:
            final_df = pd.concat(merged_rows, ignore_index=True)
            out_path = os.path.join(output_base_dir, "single_modality_model_explain.csv")
            final_df.to_csv(out_path, index=False)
            print(f"Saved merged explain summary to {out_path}")

            # Clean up
            try:
                shutil.rmtree(explain_base)
                print(f"Removed intermediate explain files under {explain_base}")
            except Exception as e:
                print(f"[WARN] Cleanup failed: {e}")

def main():
    args = get_args()
    os.makedirs(args.DA_output_dir, exist_ok=True)
    datasets = smart_load_modalities(args.input_dir)
    run_single_modality_analysis(datasets, args)


if __name__ == '__main__':
    main()
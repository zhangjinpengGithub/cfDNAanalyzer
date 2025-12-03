#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_single_independent_validation.py

Architecture
------------
This script implements single-modality training + (optional) external evaluation,
with special support for an Independent mode:

1. Internal training / CV:
   - Load single-modality datasets from args.input_dir.
   - For each dataset (file), run classification with:
       * A chosen CV scheme (KFold / LOO / Independent) via get_cv.
       * One or more feature-selection (FS) combinations from
         build_fs_combinations_from_args(args).
       * One or more classifiers from args.classifierSingle.
   - In Independent mode (cvSingle == "Independent"):
       * get_cv is expected to produce a single fold whose test_idx
         is None or empty.
       * The script then only trains on the full dataset and does not
         generate internal predictions/metrics.
       * The trained model plus FS info are saved as final_model.joblib under:
           <DA_output_dir>/<modality>/<FS_Combination>/models/<Classifier>/.

2. External evaluation (optional):
   - If args.external_input_dir is provided and contains matching single-modality
     data (same feature names), the script:
       * Loads the previously saved Independent-mode models.
       * Checks feature consistency between external data and trained model
         (missing features -> error; extra features -> warning, ignored).
       * Aligns external features to the exact order used in training.
       * Applies the same feature-selection indices and runs predictions
         on the external data.
       * Uses evaluate_classification (same as internal logic) to compute
         classification metrics for the external set.
       * For binary classification (classNum == 2), it also outputs TN/FP/FN/TP
         based on Prob_Class1 >= 0.5.
   - External prediction probabilities are saved to:
       <DA_output_dir>/<modality>/external_eval_probabilities.csv
   - External metrics are saved to:
       <DA_output_dir>/<modality>/external_eval_metrics.csv

Input formats
-------------
The script supports:
.csv:
   - Columns:
       - label (required): class label
       - sample (optional): sample ID; if absent, row index is used.
       - All other columns are treated as feature columns.
   - The label column is converted to categorical codes
     (assuming train and external share the same label set).

Key arguments (from config_parser.get_args)
-------------------------------------------
- input_dir:            Folder with single-modality training data (.csv or .npz).
- external_input_dir:   (Optional) Folder with external single-modality data.
- DA_output_dir:        Output directory root.
- classNum:             Number of classes.
- cvSingle:             CV scheme, e.g. "LOO", "KFold", or "Independent".
- nsplitSingle:         Number of splits for single-modality CV.
- cvSingle_test_ratio:  Test ratio for CV (if applicable).
- classifierSingle:     List of classifier names for single-modality models.
- FS-related arguments:   Parsed and expanded by:
      build_fs_combinations_from_args(args)
- explain:              List of explanation methods for
      compute_model_explanations (e.g. ["shap", "perm"]).
- save_model:           Whether to save trained models (used in Independent mode).

Output structure (internal)
---------------------------
For each modality <name> from input_dir, internal results are stored under:
    <DA_output_dir>/<name>/

- single_modality_probabilities.csv:
    Aggregated per-sample probabilities for all classifiers and FS combos
    (only if internal predictions exist, i.e. not pure Independent mode).

- single_modality_metrics.csv:
    Aggregated internal metrics across classifiers and FS combos
    (if internal predictions exist).

- single_modality_model_explain.csv:
    Aggregated SHAP/permutation importance summaries across folds,
    classifiers, and FS combos.

Intermediate explanation files under <DA_output_dir>/<name>/model_explain/
are merged and then deleted.

Example
-------
nohup python run_single_independent_validation.py \
  --modality single \
  --input_dir ./test_data_for_independent_validation/new_two \
  --classifierSingle KNN SVM \
  --cvSingle Independent \
  --filterMethod IG \
  --filterFrac 0.99999 \
  --DA_output_dir /home/keyao/251014_tst_cfdna/Independent_validation_1202/tst3 > independent_single.out &

Notes
-----
- This script is designed so that Independent mode handles:
    * Full-data training with no internal predictions.
    * Saving final models to be used later for external evaluation.
- External evaluation is only run if external_input_dir is provided.
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
from run_feature_selection import (
    run_feature_selection,
    infer_fs_type,
    build_fs_combinations_from_args,
)
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix
import glob
import joblib

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"


def load_single_modalities_csv(input_dir):
    """Load single-modality datasets from CSV files in input_dir.

    For each .csv file:
      - Convert column names to lower case.
      - Use all columns except 'label' and 'sample' as features.
      - Convert label column to categorical codes.
      - Use 'sample' as sample IDs if present, otherwise row index.
    """
    files = sorted(glob.glob(os.path.join(input_dir, '*.csv')))
    datasets = {}
    for f in files:
        name = os.path.splitext(os.path.basename(f))[0]
        df = pd.read_csv(f)
        df.columns = df.columns.str.lower()
        X = df.drop(columns=["label", "sample"], errors="ignore")
        X.columns = X.columns.str.replace('.csv', '', regex=False)
        y = df["label"]
        # still use category -> codes, assuming train and external share the same label set
        y = y.astype('category').cat.codes
        sample_ids = df["sample"].values if "sample" in df.columns else df.index.to_numpy()
        datasets[name] = (X, y, sample_ids)
    return datasets


def load_single_modalities_npy(input_dir):
    """Load single-modality datasets from .npz files in input_dir."""
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
    """Automatically choose between .npz and .csv loaders based on presence of .npz files."""
    npy_files = glob.glob(os.path.join(input_dir, '*.npz'))
    if npy_files:
        print(f"[INFO] Detected cached .npz files in {input_dir}, using them.")
        return load_single_modalities_npy(input_dir)
    else:
        print(f"[INFO] No cached .npz found, loading from .csv files.")
        return load_single_modalities_csv(input_dir)


def process_fold(fold_idx, train_idx, test_idx, X, y, sample_ids, args, clf_name, output_dir_combo):
    """Train and evaluate one CV fold (or full data in Independent mode).

    - If test_idx is None or empty:
        * Perform feature selection and fit on all training samples.
        * In Independent mode, save the final model to disk.
        * Return an empty list (no per-sample predictions).
    - Otherwise:
        * Perform feature selection on training set.
        * Fit classifier and predict probabilities for the test set.
        * Return a list of prediction rows for this fold.
    """
    clf = get_classifier(clf_name)

    # ---- Training part ----
    X_train = X.iloc[train_idx]
    y_train = y.iloc[train_idx]

    # If test_idx is empty: full-data training only (e.g. Independent mode)
    if test_idx is None or len(test_idx) == 0:
        fs_type = infer_fs_type(args)
        idx_list = run_feature_selection(
            X_train,
            y_train,
            args,
            fs_type=fs_type,
            verbose=False,
            fs_tag=f"{clf_name}_fold{fold_idx}",
        )

        X_train_fs = X_train.iloc[:, idx_list]
        clf.fit(X_train_fs, y_train)

        # In Independent mode, save the full-data final model
        if str(getattr(args, "cvSingle", "")).upper() == "INDEPENDENT" and getattr(args, "save_model", True):
            model_dir = os.path.join(output_dir_combo, "models", clf_name)
            os.makedirs(model_dir, exist_ok=True)
            model_bundle = {
                "model": clf,
                "selected_feature_indices": idx_list,
                "feature_names": X.columns.tolist(),
                "classNum": args.classNum,
            }
            model_path = os.path.join(model_dir, "final_model.joblib")
            joblib.dump(model_bundle, model_path)
            print(f"[INFO] Saved trained model to: {model_path}")

        # No internal predictions in this case
        return []

    # If test_idx is not empty, perform normal train/test within this fold
    X_test = X.iloc[test_idx]
    y_test = y.iloc[test_idx]

    fs_type = infer_fs_type(args)
    idx_list = run_feature_selection(
        X_train,
        y_train,
        args,
        fs_type=fs_type,
        verbose=False,
        fs_tag=f"{clf_name}_fold{fold_idx}",
    )

    X_train_fs = X_train.iloc[:, idx_list]
    X_test_fs = X_test.iloc[:, idx_list] if X_test is not None else None

    clf.fit(X_train_fs, y_train)

    # ---- Optional: model explanations ----
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
            explain_methods=args.explain,
        )

    # ---- Predict probabilities ----
    if hasattr(clf, "predict_proba"):
        y_probs = clf.predict_proba(X_test_fs)
    else:
        pred_labels = clf.predict(X_test_fs)
        y_probs = []
        for label in pred_labels:
            prob = [0.0] * args.classNum
            prob[int(label)] = 1.0
            y_probs.append(prob)

    # Build per-sample rows for this fold
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
    """Run internal training/CV for each single-modality dataset."""
    fs_combinations = build_fs_combinations_from_args(args)
    if not fs_combinations:
        raise ValueError("No valid FS parameters provided.")

    for name, (X, y, sample_ids) in datasets.items():
        print(f"\n===== Processing single modality: {name} =====")
        cv = get_cv(method=args.cvSingle, nsplit=args.nsplitSingle, y=y, test_ratio=args.cvSingle_test_ratio)
        output_base_dir = os.path.join(args.DA_output_dir, name)
        os.makedirs(output_base_dir, exist_ok=True)

        all_predictions, all_metrics = [], []
        is_independent = str(args.cvSingle).upper() == "INDEPENDENT"

        for fs_combo in fs_combinations:
            fs_label = fs_combo.get("fs_label", "default")
            output_dir_combo = os.path.join(output_base_dir, fs_label)

            for clf_name in args.classifierSingle:
                print(f"\n>>> Running classifier: {clf_name} with FS combo: {fs_label}")
                start_time = time.time()
                tracemalloc.start()
                process = psutil.Process()  # currently only used to trigger tracemalloc; can be extended

                fold_results = []
                with ThreadPoolExecutor() as executor:
                    futures = [
                        executor.submit(
                            process_fold,
                            i,
                            tr,
                            te,
                            X,
                            y,
                            sample_ids,
                            args,
                            clf_name,
                            output_dir_combo,
                        )
                        for i, (tr, te) in enumerate(cv.split(X, y))
                    ]
                    for fut in tqdm(as_completed(futures), total=len(futures)):
                        rows = fut.result()
                        # In Independent mode, rows may be an empty list (no internal predictions)
                        if rows:
                            fold_results.extend(rows)

                total_time = round(time.time() - start_time, 4)
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                peak_mem = round(peak / 1024 / 1024, 4)

                # No fold_results: likely Independent mode (no internal predictions)
                if not fold_results:
                    if is_independent:
                        print(
                            f"[INFO] cvSingle=Independent: only trained final model for "
                            f"{clf_name}, FS={fs_label}; no internal predictions/metrics."
                        )
                    else:
                        print(
                            f"[WARN] No fold results produced for {clf_name}, FS={fs_label}; "
                            "skip internal predictions/metrics."
                        )
                    # Skip prediction_df/metrics for this combination
                    continue

                prediction_df = pd.DataFrame(fold_results)
                prediction_df["Classifier"] = clf_name
                prediction_df["FS_Combination"] = fs_label

                prob_cols = [col for col in prediction_df.columns if col.startswith("Prob_Class")]
                for col in prob_cols:
                    if col in prediction_df.columns:
                        prediction_df[col] = prediction_df[col].round(6)

                ordered_cols = ['SampleID', 'TrueLabel', 'FS_Combination', 'Classifier']
                prediction_df = prediction_df[ordered_cols + prob_cols]
                all_predictions.append(prediction_df)

                metrics = evaluate_classification(prediction_df, class_num=args.classNum)
                metrics["Classifier"] = clf_name
                metrics["FS_Combination"] = fs_label
                metrics["TotalTime_sec"] = total_time
                metrics["PeakMemory_MB"] = peak_mem
                all_metrics.append(metrics)

        # === Save internal metrics & probabilities ===
        if all_predictions:
            final_pred_df = pd.concat(all_predictions, axis=0)
            final_pred_path = os.path.join(output_base_dir, "single_modality_probabilities.csv")
            final_pred_df.to_csv(final_pred_path, index=False)
            print(f"Saved prediction results to: {final_pred_path}")
        else:
            print("[INFO] No internal prediction results to save (maybe Independent mode only).")

        if all_metrics:
            metrics_df = pd.DataFrame(all_metrics)
            metric_order = ['FS_Combination', 'Classifier', 'TotalTime_sec', 'PeakMemory_MB']
            other_metrics = [col for col in metrics_df.columns if col not in metric_order]
            metrics_df = metrics_df[metric_order + other_metrics]
            metrics_path = os.path.join(output_base_dir, "single_modality_metrics.csv")
            metrics_df.to_csv(metrics_path, index=False)
            print(f"Saved performance summary to: {metrics_path}")
        else:
            print("[INFO] No internal metrics to save (maybe Independent mode only).")

        # === Merge model_explain summaries across folds/classifiers/FS combos ===
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

                    # ---- Permutation importance ----
                    if os.path.exists(perm_file):
                        df = pd.read_csv(perm_file)

                        if "Feature" in df.columns and "ImportanceMean" in df.columns:
                            df = df.set_index("Feature")["ImportanceMean"]
                        else:
                            raise ValueError(f"Permutation file format invalid: {perm_file}")

                        df = df.rename(f"Perm_Fold{fold_idx}")
                        perm_fold_tables.append(df)

                # Aggregate SHAP across folds
                if shap_fold_tables:
                    shap_all = pd.concat(shap_fold_tables, axis=1)
                    shap_mean = shap_all.abs().mean(axis=1).rename("MeanAbsSHAP")
                    shap_std = shap_all.abs().std(axis=1).rename("MeanStdSHAP")
                    shap_summary = pd.concat([shap_mean, shap_std, shap_all], axis=1).reset_index()
                else:
                    shap_summary = pd.DataFrame(columns=["Feature", "MeanAbsSHAP", "MeanStdSHAP"])

                # Aggregate permutation importance across folds
                if perm_fold_tables:
                    perm_all = pd.concat(perm_fold_tables, axis=1)
                    perm_mean = perm_all.abs().mean(axis=1).rename("MeanAbsPermImportance")
                    perm_std = perm_all.abs().std(axis=1).rename("MeanStdPermImportance")
                    perm_summary = pd.concat([perm_mean, perm_std, perm_all], axis=1).reset_index()
                else:
                    perm_summary = pd.DataFrame(columns=["Feature", "MeanAbsPermImportance", "MeanStdPermImportance"])

                # Merge SHAP + Permutation summaries
                merged = pd.merge(shap_summary, perm_summary, on="Feature", how="outer")
                merged.insert(0, "Classifier", clf_name)
                merged.insert(0, "FS_Combination", fs_label)

                # Reorder columns
                fold_cols_shap = sorted(
                    [c for c in merged.columns if c.startswith("SHAP_Fold")],
                    key=lambda x: int(x.replace("SHAP_Fold", "")) if x.replace("SHAP_Fold", "").isdigit() else 0,
                )
                fold_cols_perm = sorted(
                    [c for c in merged.columns if c.startswith("Perm_Fold")],
                    key=lambda x: int(x.replace("Perm_Fold", "")) if x.replace("Perm_Fold", "").isdigit() else 0,
                )

                base_cols = [
                    "FS_Combination",
                    "Classifier",
                    "Feature",
                    "MeanAbsSHAP",
                    "MeanStdSHAP",
                    "MeanAbsPermImportance",
                    "MeanStdPermImportance",
                ]

                merged = merged[base_cols + fold_cols_shap + fold_cols_perm]
                merged_rows.append(merged)

        # Save merged explanation summary and clean up intermediates
        if merged_rows:
            final_df = pd.concat(merged_rows, ignore_index=True)
            out_path = os.path.join(output_base_dir, "single_modality_model_explain.csv")
            final_df.to_csv(out_path, index=False)
            print(f"Saved merged explain summary to {out_path}")

            try:
                shutil.rmtree(explain_base)
                print(f"Removed intermediate explain files under {explain_base}")
            except Exception as e:
                print(f"[WARN] Cleanup failed: {e}")


def evaluate_on_external(external_datasets, args):
    """
    Use already trained Independent-mode models (final_model.joblib) to evaluate
    external datasets.

    Behavior:
    --------
    - For each modality and FS combination:
        * Load the corresponding final model.
        * Validate that external features match the model's feature set
          (missing -> error, extra -> warning and ignored).
        * Align external features to the same order as used during training.
        * Apply stored feature-selection indices.
        * Predict probabilities on the external samples.
    - Use evaluate_classification to compute all metrics in a way that is
      fully consistent with internal evaluation (binary/multi-class).
    - For binary classification (classNum == 2), additionally compute:
        TN, FP, FN, TP using a confusion matrix derived from:
          Prob_Class1 >= 0.5 as positive class prediction.

    Outputs:
    --------
    - For each modality <name>:
        * <DA_output_dir>/<name>/external_eval_probabilities.csv
        * <DA_output_dir>/<name>/external_eval_metrics.csv
    """
    fs_combinations = build_fs_combinations_from_args(args)
    if not fs_combinations:
        raise ValueError("No valid FS parameters provided for external evaluation.")

    for name, (X_ext, y_ext, sample_ids_ext) in external_datasets.items():
        print(f"\n===== External evaluation on modality: {name} =====")
        output_base_dir = os.path.join(args.DA_output_dir, name)

        all_predictions = []
        all_metrics = []

        for fs_combo in fs_combinations:
            fs_label = fs_combo.get("fs_label", "default")

            for clf_name in args.classifierSingle:
                model_dir = os.path.join(output_base_dir, fs_label, "models", clf_name)
                model_path = os.path.join(model_dir, "final_model.joblib")
                if not os.path.exists(model_path):
                    print(f"[WARN] Model not found for modality={name}, FS={fs_label}, clf={clf_name}: {model_path}")
                    continue

                print(f"[INFO] Loading model from: {model_path}")
                bundle = joblib.load(model_path)
                clf = bundle["model"]
                idx_list = bundle["selected_feature_indices"]
                feature_names = bundle["feature_names"]
                class_num_trained = bundle.get("classNum", args.classNum)

                # ---- Feature consistency checks: external vs model ----
                expected_features = set(feature_names)
                actual_features = set(X_ext.columns)

                # 1) Missing features -> error
                missing_cols = list(expected_features - actual_features)
                if missing_cols:
                    raise ValueError(
                        f"External dataset '{name}' is missing columns required by model "
                        f"(FS={fs_label}, clf={clf_name}): {missing_cols}"
                    )

                # 2) Extra features -> warning but ignored
                extra_cols = list(actual_features - expected_features)
                if extra_cols:
                    print(
                        f"[WARN] External dataset '{name}' has extra columns not used by the model "
                        f"(FS={fs_label}, clf={clf_name}). They will be ignored: {extra_cols}"
                    )

                # 3) Align columns strictly to the model's feature_names
                X_ext_aligned = X_ext[feature_names]

                # Sanity check: column order must match exactly
                if list(X_ext_aligned.columns) != list(feature_names):
                    raise ValueError(
                        f"Column order mismatch after alignment for modality '{name}', "
                        f"FS={fs_label}, clf={clf_name}."
                    )

                # 4) Apply the stored feature-selection indices
                X_ext_fs = X_ext_aligned.iloc[:, idx_list]

                # ---- Predict probabilities on external data ----
                if hasattr(clf, "predict_proba"):
                    y_probs = clf.predict_proba(X_ext_fs)
                else:
                    # No predict_proba: build one-hot probabilities from hard predictions
                    pred_labels = clf.predict(X_ext_fs)
                    y_probs = []
                    for label in pred_labels:
                        prob = [0.0] * class_num_trained
                        prob[int(label)] = 1.0
                        y_probs.append(prob)
                    y_probs = np.array(y_probs)

                # Predicted labels (for reference)
                y_pred = np.argmax(y_probs, axis=1)

                # Build probability detail rows (align with internal single_modality_probabilities)
                rows = []
                for i, sample_id in enumerate(sample_ids_ext):
                    row = {
                        "SampleID": sample_id,
                        "TrueLabel": int(y_ext.iloc[i]),
                        "FS_Combination": fs_label,
                        "Classifier": clf_name,
                    }
                    for cls in range(class_num_trained):
                        row[f"Prob_Class{cls}"] = float(y_probs[i, cls])
                    rows.append(row)

                pred_df = pd.DataFrame(rows)

                # Round probability columns
                prob_cols = [c for c in pred_df.columns if c.startswith("Prob_Class")]
                for col in prob_cols:
                    pred_df[col] = pred_df[col].round(6)

                # Reorder columns for readability
                ordered_cols = ['SampleID', 'TrueLabel', 'FS_Combination', 'Classifier']
                pred_df = pred_df[ordered_cols + prob_cols]

                all_predictions.append(pred_df)

                # Compute metrics using the same logic as internal evaluation
                metrics = evaluate_classification(pred_df, class_num=args.classNum)

                # ---- For binary classification, compute TN/FP/FN/TP ----
                tn = fp = fn = tp = np.nan
                if args.classNum == 2:
                    y_true_bin = pred_df['TrueLabel'].values
                    # Use Prob_Class1 >= 0.5 as positive threshold (consistent with evaluate_classification)
                    if 'Prob_Class1' in pred_df.columns:
                        y_prob_pos = pred_df['Prob_Class1'].values
                        y_pred_bin = (y_prob_pos >= 0.5).astype(int)

                        cm = confusion_matrix(y_true_bin, y_pred_bin, labels=[0, 1])
                        if cm.shape == (2, 2):
                            tn, fp, fn, tp = cm.ravel()

                metrics["TN"] = None if np.isnan(tn) else int(tn)
                metrics["FP"] = None if np.isnan(fp) else int(fp)
                metrics["FN"] = None if np.isnan(fn) else int(fn)
                metrics["TP"] = None if np.isnan(tp) else int(tp)

                # Add external evaluation metadata
                metrics["Modality"] = name
                metrics["FS_Combination"] = fs_label
                metrics["Classifier"] = clf_name

                all_metrics.append(metrics)

        # === Save external prediction probabilities ===
        if all_predictions:
            external_pred_df = pd.concat(all_predictions, axis=0)
            prob_path = os.path.join(output_base_dir, "external_eval_probabilities.csv")
            external_pred_df.to_csv(prob_path, index=False)
            print(f"[INFO] Saved external prediction probabilities to: {prob_path}")

        # === Save external evaluation metrics ===
        if all_metrics:
            metrics_df = pd.DataFrame(all_metrics)

            preferred_order = [
                "Modality",
                "FS_Combination",
                "Classifier",
                "accuracy",
                "precision",
                "recall",
                "sensitivity",
                "specificity",
                "f1",
                "auc",
                "TN",
                "FP",
                "FN",
                "TP",
            ]
            existing_cols = [c for c in preferred_order if c in metrics_df.columns]
            other_cols = [c for c in metrics_df.columns if c not in existing_cols]
            metrics_df = metrics_df[existing_cols + other_cols]

            metrics_path = os.path.join(output_base_dir, "external_eval_metrics.csv")
            metrics_df.to_csv(metrics_path, index=False)
            print(f"[INFO] Saved external evaluation metrics to: {metrics_path}")


def main():
    args = get_args()
    os.makedirs(args.DA_output_dir, exist_ok=True)

    # 1) Internal training + CV (in Independent mode, this also saves final_model.joblib)
    datasets = smart_load_modalities(args.input_dir)
    run_single_modality_analysis(datasets, args)

    # 2) If external_input_dir is provided, evaluate external datasets
    #    using the already saved models (Independent mode).
    external_input_dir = getattr(args, "external_input_dir", None)
    if external_input_dir:
        if not os.path.isdir(external_input_dir):
            raise ValueError(f"external_input_dir is not a directory: {external_input_dir}")
        external_datasets = smart_load_modalities(external_input_dir)
        evaluate_on_external(external_datasets, args)


if __name__ == '__main__':
    main()
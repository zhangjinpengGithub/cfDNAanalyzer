#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_multi_independent_validation.py

Architecture
------------
This script implements multi-modal external (independent) validation.

- Training set:         args.input_dir
- Independent test set: args.external_input_dir
- It performs one single full-data training + external validation.
  There is no internal KFold / LOO cross-validation in this script.

Workflow
-------------------
1. Load multi-modal training data from input_dir.
2. Load multi-modal external validation data from external_input_dir.
3. For each modality:
   - Check that train and external data have the same feature set.
   - If feature order differs, reorder external features to match train.
4. Concatenate training and external samples for each modality to form
   a combined matrix:
       [ X_train
         X_external ]
   and similarly combine labels and sample IDs.
5. Construct an ExternalCV object that mimics a scikit-learn CV splitter:
   - train_idx = [0 ... n_train-1]
   - test_idx  = [n_train ... n_train+n_ext-1]
6. For each fusion type in args.fusion_type:
   - For each feature-selection (FS) combination (from build_fs_combinations_from_args):
     - For each classifier in args.classifierMulti:
       - Run one of:
         * Concatenation fusion: run_concatenation_fusion
         * Model-based fusion:   run_voting_fusion / run_stacking_fusion
         * Transformation-based: run_trans_fusion
       - All feature selection and model training are done only on the training part (internal samples).
       - The predictors are then applied to the external part (test samples).
7. The predictions and evaluation metrics are saved under:
   args.DA_output_dir/external_validation/...
   specifically:
   external_validation/<fusion_type>/<fs_label>/<method>/.

Input formats
-------------
The script supports:
.csv:
   - Columns:
       - label (required): class label
       - sample (optional): sample ID. If missing, row index is used.
       - All other columns are treated as feature columns.
   - File name (without extension) is treated as the modality name.

The training and external directories must contain the same set of files
with the same names (one per modality), and the feature names must be identical.

Key arguments (from config_parser.get_args)
-------------------------------------------
- input_dir:           Folder with training multi-modal data.
- external_input_dir:  Folder with external validation multi-modal data.
- DA_output_dir:       Output directory root.
- classNum:            Number of classes.
- fusion_type:         List of fusion types: e.g. ["concat", "model", "trans"].
- classifierMulti:     List of classifier names for multi-modal models.
- model_method:        List of model-based fusion schemes (e.g. ["average", "weighted", "stack"]).
- trans_method:        List of transformation-based methods (e.g. ["pca", "rbf", "snf"]).
- Feature-selection related arguments are parsed by:
    build_fs_combinations_from_args(args)

Output structure
----------------
Under args.DA_output_dir/external_validation/:

For each fusion_type (e.g. "concat", "model", "trans"):
    external_validation/
        └── <fusion_type>/
            └── <fs_label>/
                ├── <method>/  # e.g. "concat", "average", "stack", "pca", ...
                │    ├── <Classifier>_predictions.csv
                │    ├── <Classifier>_metrics.csv
                │    ├── all_classifiers_predictions.csv
                │    └── all_classifiers_metrics.csv
                └── all_methods_all_classifiers_predictions.csv
                └── all_methods_all_classifiers_metrics.csv

Each *_predictions.csv contains:
- SampleID, TrueLabel, FusionType, FusionMethod, FS_Combination, Classifier,
  and probability columns (Prob_Class0, Prob_Class1, ...).

Each *_metrics.csv contains:
- Standard classification metrics from evaluate_classification
  plus runtime and memory usage.

Example usage
-------------
Assuming your main config parser supports the following arguments:

nohup python run_multi_independent_only.py \
  --modality multi \
  --input_dir ./train_data/raw_two/raw_no_stand/std_train \
  --external_input_dir ./test_data/new_two/std_test \
  --fusion_type concat model trans \
  --model_method average weighted stack \
  --trans_method pca rbf snf \
  --classifierMulti KNN SVM \
  --cvMulti Independent \
  --filterMethod IG \
  --filterFrac 0.99999 \
  --DA_output_dir ./Independent_validation_1202/tst4 > independent_multi.out

Notes
-----
- This script ignores args.cvMulti and always uses a single split:
  all training samples vs. all external samples.
- It is designed to be parallel in spirit to your other multi-modal scripts
  but specialized only for independent external validation.
"""

import os
import glob
import numpy as np
import pandas as pd
import time
import tracemalloc
import psutil
from tqdm import tqdm

from config_parser import get_args
from methods.machine_learning.generate_metrics import evaluate_classification
from run_feature_selection import build_fs_combinations_from_args
from methods.machine_learning.model_based_fusion import run_voting_fusion, run_stacking_fusion
from methods.machine_learning.concatenation_fusion import run_concatenation_fusion
from methods.machine_learning.transformation_based_fusion import run_trans_fusion
# from methods.machine_learning.train_test_split import get_cv

# import joblib

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"


# -------- Data loading --------
def load_modalities_smart(input_dir):
    """Load multi-modal data from input_dir, preferring .npz over .csv.

    Returns
    -------
    modality_names : list of str
        List of modality names (file names without extension).
    modality_data : list of dict
        Each dict contains:
            - "X" : pd.DataFrame of features
            - "y" : pd.Series of labels
            - "sample_ids" : array-like of sample IDs
    """
    npz_files = sorted(glob.glob(os.path.join(input_dir, '*.npz')))
    csv_files = sorted(glob.glob(os.path.join(input_dir, '*.csv')))

    modality_names = []
    modality_data = []

    if npz_files:
        print(f"[INFO] Loading .npz files from {input_dir}")
        for npz_file in npz_files:
            name = os.path.splitext(os.path.basename(npz_file))[0]
            with np.load(npz_file, allow_pickle=True) as data:
                X = pd.DataFrame(data["data"])
                y = pd.Series(data["label"])
                sample_ids = data["sample_ids"]
            modality_names.append(name)
            modality_data.append({"X": X, "y": y, "sample_ids": sample_ids})
    else:
        print(f"[INFO] Loading .csv files from {input_dir}")
        for csv_file in csv_files:
            name = os.path.splitext(os.path.basename(csv_file))[0]
            df = pd.read_csv(csv_file)
            df.columns = df.columns.str.lower()

            y = df["label"]
            if "sample" in df.columns:
                sample_ids = df["sample"].tolist()
            else:
                sample_ids = df.index.tolist()

            X = df.drop(columns=["label", "sample"], errors="ignore")

            modality_names.append(name)
            modality_data.append({"X": X, "y": y, "sample_ids": sample_ids})

    return modality_names, modality_data


# -------- Utilities: merge and save results --------
def merge_classifier_outputs(method_dir, classifier_list):
    """Within one fusion method directory, merge prediction/metrics across classifiers."""
    all_preds, all_metrics = [], []
    for clf_name in classifier_list:
        pred_path = os.path.join(method_dir, f"{clf_name}_predictions.csv")
        metrics_path = os.path.join(method_dir, f"{clf_name}_metrics.csv")
        if os.path.exists(pred_path):
            all_preds.append(pd.read_csv(pred_path))
        if os.path.exists(metrics_path):
            all_metrics.append(pd.read_csv(metrics_path))
    if all_preds:
        pd.concat(all_preds).to_csv(
            os.path.join(method_dir, "all_classifiers_predictions.csv"), index=False
        )
    if all_metrics:
        pd.concat(all_metrics).to_csv(
            os.path.join(method_dir, "all_classifiers_metrics.csv"), index=False
        )


def merge_all_methods_outputs(base_dir, method_list, classifier_list):
    """Within one FS combination directory, merge results across fusion methods."""
    all_preds, all_metrics = [], []
    for method in method_list:
        method_dir = os.path.join(base_dir, method)
        pred_file = os.path.join(method_dir, "all_classifiers_predictions.csv")
        metric_file = os.path.join(method_dir, "all_classifiers_metrics.csv")
        if os.path.exists(pred_file):
            all_preds.append(pd.read_csv(pred_file))
        if os.path.exists(metric_file):
            all_metrics.append(pd.read_csv(metric_file))

    if all_preds:
        combined_preds = pd.concat(all_preds, ignore_index=True)
        combined_preds = combined_preds.drop_duplicates()
        combined_preds.to_csv(
            os.path.join(base_dir, "all_methods_all_classifiers_predictions.csv"),
            index=False,
        )
    if all_metrics:
        combined_metrics = pd.concat(all_metrics, ignore_index=True)
        combined_metrics = combined_metrics.drop_duplicates()
        combined_metrics.to_csv(
            os.path.join(base_dir, "all_methods_all_classifiers_metrics.csv"),
            index=False,
        )


def save_outputs(
    prediction_df,
    clf_name,
    fs_label,
    method_dir,
    args,
    fusion_type=None,
    fusion_method=None,
    elapsed_time=None,
    peak_memory=None,
):
    """Save predictions and metrics for one (fusion_type, fusion_method, FS, classifier) combination."""
    process = psutil.Process()
    mem_used = round(process.memory_info().rss / 1024 / 1024, 4)
    total_time = elapsed_time if elapsed_time is not None else None
    peak_mem = peak_memory if peak_memory is not None else None

    prediction_df["Classifier"] = clf_name
    prediction_df["FS_Combination"] = fs_label
    prediction_df["FusionType"] = fusion_type
    prediction_df["FusionMethod"] = fusion_method

    prob_cols = [f"Prob_Class{i}" for i in range(args.classNum)]
    for col in prob_cols:
        if col in prediction_df.columns:
            prediction_df[col] = prediction_df[col].round(6)

    pred_order = [
        "SampleID",
        "TrueLabel",
        "FusionType",
        "FusionMethod",
        "FS_Combination",
        "Classifier",
    ] + prob_cols

    available_cols = [c for c in pred_order if c in prediction_df.columns]
    other_cols = [c for c in prediction_df.columns if c not in available_cols]
    prediction_df = prediction_df[available_cols + other_cols]

    metrics = evaluate_classification(prediction_df, class_num=args.classNum)
    metrics["Classifier"] = clf_name
    metrics["FS_Combination"] = fs_label
    metrics["FusionType"] = fusion_type
    metrics["FusionMethod"] = fusion_method
    metrics["TotalTime_sec"] = total_time
    metrics["FinalMemory_MB"] = mem_used
    metrics["PeakMemory_MB"] = peak_mem

    metrics_df = pd.DataFrame([metrics])
    metric_order = ["FusionType", "FusionMethod", "FS_Combination", "Classifier"]
    other_metrics = [c for c in metrics_df.columns if c not in metric_order]
    metrics_df = metrics_df[metric_order + other_metrics]

    os.makedirs(method_dir, exist_ok=True)
    prediction_df.to_csv(
        os.path.join(method_dir, f"{clf_name}_predictions.csv"), index=False
    )
    metrics_df.to_csv(
        os.path.join(method_dir, f"{clf_name}_metrics.csv"), index=False
    )

    print(f"[INFO] Saved results to: {method_dir}")
    print(f"[INFO] Time = {total_time}s | Memory = {mem_used}MB | Peak = {peak_mem}MB")


# -------- Simple CV wrapper: ExternalCV --------
class ExternalCV:
    """
    A simple CV-like object for external validation, mimicking the sklearn interface.

    - Training indices: [0 .. n_train-1]
    - Test indices:     [n_train .. n_total-1]
    """

    def __init__(self, n_train, n_total):
        self.n_train = n_train
        self.n_total = n_total
        self.train_idx = np.arange(n_train)
        self.test_idx = np.arange(n_train, n_total)

    def split(self, X, y=None, groups=None):
        # X, y, groups are unused; they are only here for interface compatibility.
        yield self.train_idx, self.test_idx

    def get_n_splits(self, X=None, y=None, groups=None):
        return 1


# -------- External Independent validation logic --------
def run_external_validation(
    train_modality_names,
    train_modality_data,
    ext_modality_names,
    ext_modality_data,
    args,
):
    """
    Use input_dir as the training set and external_input_dir as the
    independent validation set, and perform a single external evaluation.

    External results are saved under:
        DA_output_dir/external_validation/<fusion_type>/<fs_label>/<method>/.
    """
    print("\n================ External Validation (Independent) ================")

    # 1. Check consistency of modality counts and names
    if len(train_modality_names) != len(ext_modality_names):
        raise ValueError("Train and external modalities have different number of files.")
    for tn, en in zip(train_modality_names, ext_modality_names):
        if tn != en:
            raise ValueError(f"Modality name mismatch: train={tn}, external={en}")

    # 2. Build combined Xs / y / sample_ids across train + external
    Xs_combined = []
    ys_train, ys_ext = None, None
    sample_ids_train, sample_ids_ext = None, None

    for (train_m, ext_m) in zip(train_modality_data, ext_modality_data):
        X_train = train_m["X"].copy()
        X_ext = ext_m["X"].copy()

        # Ensure the feature sets are identical; if only order differs, align to training order
        if list(X_train.columns) != list(X_ext.columns):
            if set(X_train.columns) != set(X_ext.columns):
                raise ValueError(
                    "Feature columns between train and external data are not identical."
                )
            X_ext = X_ext[X_train.columns]

        X_combined = pd.concat([X_train, X_ext], axis=0, ignore_index=True)
        Xs_combined.append(X_combined)

        if ys_train is None:
            ys_train = train_m["y"].reset_index(drop=True)
            ys_ext = ext_m["y"].reset_index(drop=True)
            sample_ids_train = list(train_m["sample_ids"])
            sample_ids_ext = list(ext_m["sample_ids"])

    ys_combined = pd.concat([ys_train, ys_ext], axis=0, ignore_index=True)
    sample_ids_combined = sample_ids_train + sample_ids_ext

    n_train = len(ys_train)
    n_ext = len(ys_ext)
    n_total = n_train + n_ext

    print(f"[INFO] #Train samples = {n_train}, #External samples = {n_ext}")

    # ExternalCV: all training samples for train, all external samples for test
    cv_external = ExternalCV(n_train=n_train, n_total=n_total)

    fs_combinations = build_fs_combinations_from_args(args)
    if not fs_combinations:
        raise ValueError("No valid FS parameters provided for external validation.")

    ext_root_dir = os.path.join(args.DA_output_dir, "external_validation")
    os.makedirs(ext_root_dir, exist_ok=True)

    # If the user set cvMulti to something else in the config, warn and ignore it
    if str(getattr(args, "cvMulti", "Independent")).upper() != "INDEPENDENT":
        print(
            "[WARN] This script only implements Independent-style external validation."
        )
        print("[WARN] Ignoring args.cvMulti and using ExternalCV (train vs external).")

    for fusion_type in args.fusion_type:
        output_base_dir = os.path.join(ext_root_dir, fusion_type)
        os.makedirs(output_base_dir, exist_ok=True)

        print(f"\n[External] Fusing modalities using: {fusion_type}")
        for fs_combo in fs_combinations:
            # Apply FS combo parameters into args
            for key in fs_combo:
                setattr(args, key, fs_combo[key])
            fs_label = fs_combo.get("fs_label", "default")
            method_names = []

            for clf_name in args.classifierMulti:
                Xs = Xs_combined
                ys = ys_combined
                sample_ids = sample_ids_combined

                # -------- Concatenation fusion --------
                if fusion_type == "concat":
                    method_name = "concat"
                    method_dir = os.path.join(output_base_dir, fs_label, method_name)
                    os.makedirs(method_dir, exist_ok=True)
                    method_names.append(method_name)

                    print(
                        f"\n>>> [External Concat] Classifier: {clf_name}, FS: {fs_label}"
                    )
                    start_time = time.time()
                    tracemalloc.start()

                    prediction_df = run_concatenation_fusion(
                        Xs,
                        clf_name,
                        cv_external,
                        ys,
                        sample_ids,
                        args,
                        fs_label,
                        use_thread=True,
                    )

                    end_time = time.time()
                    current, peak = tracemalloc.get_traced_memory()
                    tracemalloc.stop()

                    elapsed_time = round(end_time - start_time, 4)
                    peak_mem = round(peak / 1024 / 1024, 4)

                    # prediction_df should contain only external test samples
                    # according to ExternalCV's test indices.
                    save_outputs(
                        prediction_df,
                        clf_name,
                        fs_label,
                        method_dir,
                        args,
                        fusion_type="concat",
                        fusion_method="concat",
                        elapsed_time=elapsed_time,
                        peak_memory=peak_mem,
                    )

                # -------- Model-based fusion --------
                elif fusion_type == "model":
                    for model_method in args.model_method:
                        method_dir = os.path.join(output_base_dir, fs_label, model_method)
                        os.makedirs(method_dir, exist_ok=True)
                        if model_method not in method_names:
                            method_names.append(model_method)

                        print(
                            f"\n>>> [External Model-{model_method}] Classifier: {clf_name}, FS: {fs_label}"
                        )
                        start_time = time.time()
                        tracemalloc.start()

                        if model_method in ["average", "weighted", "majority"]:
                            prediction_df = run_voting_fusion(
                                Xs,
                                clf_name,
                                cv_external,
                                ys,
                                sample_ids,
                                args,
                                fs_label,
                                method=model_method,
                                use_thread=True,
                            )
                        elif model_method == "stack":
                            prediction_df = run_stacking_fusion(
                                Xs,
                                clf_name,
                                cv_external,
                                ys,
                                sample_ids,
                                args,
                                fs_label,
                                use_thread=True,
                            )
                        else:
                            raise ValueError(
                                "Unsupported model fusion method for external validation."
                            )

                        end_time = time.time()
                        current, peak = tracemalloc.get_traced_memory()
                        tracemalloc.stop()

                        elapsed_time = round(end_time - start_time, 4)
                        peak_mem = round(peak / 1024 / 1024, 4)

                        save_outputs(
                            prediction_df,
                            clf_name,
                            fs_label,
                            method_dir,
                            args,
                            fusion_type="model",
                            fusion_method=model_method,
                            elapsed_time=elapsed_time,
                            peak_memory=peak_mem,
                        )

                # -------- Transformation-based fusion --------
                elif fusion_type == "trans":
                    for trans_method in args.trans_method:
                        method_dir = os.path.join(output_base_dir, fs_label, trans_method)
                        os.makedirs(method_dir, exist_ok=True)
                        if trans_method not in method_names:
                            method_names.append(trans_method)

                        print(
                            f"\n>>> [External Trans-{trans_method}] Classifier: {clf_name}, FS: {fs_label}"
                        )
                        start_time = time.time()
                        tracemalloc.start()

                        prediction_df = run_trans_fusion(
                            Xs,
                            clf_name,
                            cv_external,
                            ys,
                            sample_ids,
                            args,
                            fs_label,
                            method=trans_method,
                            use_thread=True,
                        )

                        end_time = time.time()
                        current, peak = tracemalloc.get_traced_memory()
                        tracemalloc.stop()

                        elapsed_time = round(end_time - start_time, 4)
                        peak_mem = round(peak / 1024 / 1024, 4)

                        save_outputs(
                            prediction_df,
                            clf_name,
                            fs_label,
                            method_dir,
                            args,
                            fusion_type="trans",
                            fusion_method=trans_method,
                            elapsed_time=elapsed_time,
                            peak_memory=peak_mem,
                        )

                else:
                    raise ValueError(
                        f"Unsupported fusion type for external validation: {fusion_type}"
                    )

            # Merge results across classifiers for each method under this FS combination
            for method in method_names:
                method_dir = os.path.join(output_base_dir, fs_label, method)
                merge_classifier_outputs(method_dir, args.classifierMulti)

            # Merge results across methods for this FS combination
            merge_all_methods_outputs(
                os.path.join(output_base_dir, fs_label), method_names, args.classifierMulti
            )

    print("\n[INFO] External validation finished. Results are under 'external_validation' folder.")


def main():
    args = get_args()
    os.makedirs(args.DA_output_dir, exist_ok=True)

    # 1) Load training data (from input_dir)
    train_names, train_data = load_modalities_smart(args.input_dir)
    print(f"[INFO] Loaded training modalities: {train_names}")

    # 2) Require external_input_dir as the independent validation set
    external_input_dir = getattr(args, "external_input_dir", None)
    if external_input_dir is None or not os.path.isdir(external_input_dir):
        raise ValueError(
            "This Independent script requires a valid external_input_dir for validation."
        )

    ext_names, ext_data = load_modalities_smart(external_input_dir)
    print(f"[INFO] Loaded external modalities: {ext_names}")

    # 3) Run Independent external validation
    run_external_validation(train_names, train_data, ext_names, ext_data, args)


if __name__ == "__main__":
    main()
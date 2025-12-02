"""
This module provides optional model interpretation utilities:
    - SHAP values
    - Permutation importance
"""
import os
import warnings
import numpy as np
import pandas as pd
import shap
from sklearn.inspection import permutation_importance


def compute_model_explanations(clf, X_train, y_train, X_test, y_test=None,
                               output_dir=None, explain_methods=None, random_state=42):
    """
    Compute and save model interpretation results.

    Parameters
    ----------
    clf : sklearn-like trained classifier
        The fitted model.
    X_train : pandas.DataFrame
        Training features.
    y_train : pandas.Series or np.ndarray
        Training labels.
    X_test : pandas.DataFrame
        Test features (used for evaluation and explanation).
    y_test : pandas.Series or np.ndarray, optional
        Test labels (required for permutation importance).
    output_dir : str
        Directory where explanation results will be saved.
    explain_methods : list
        List of explanation types, e.g. ['SHAP', 'perm'].
    random_state : int
        Random seed for reproducibility.
    """

    if not explain_methods or output_dir is None:
        return

    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------
    # Ensure consistent feature names
    # -------------------------------
    # Convert to DataFrame with feature names compatible with the trained model
    if hasattr(clf, "feature_names_in_"):
        X_train = pd.DataFrame(X_train, columns=clf.feature_names_in_)
        X_test = pd.DataFrame(X_test, columns=clf.feature_names_in_)
    else:
        X_train = pd.DataFrame(X_train, columns=getattr(X_train, "columns", None))
        X_test = pd.DataFrame(X_test, columns=getattr(X_test, "columns", None))

    # Silence deprecation and irrelevant sklearn warnings
    warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
    warnings.filterwarnings("ignore", message="X does not have valid feature names")

    # -------------------------------
    # 1. Permutation Importance
    # -------------------------------
    if 'perm' in [m.lower() for m in explain_methods]:
        print("[INFO] Computing permutation importance...")
        try:
            if y_test is None:
                print("[WARN] y_test not provided. Skipping permutation importance.")
            else:
                perm_result = permutation_importance(
                    clf,
                    X_test,
                    y_test,
                    n_repeats=5,
                    random_state=random_state,
                    n_jobs=-1
                )
                perm_df = pd.DataFrame({
                    'Feature': X_test.columns,
                    'ImportanceMean': perm_result.importances_mean,
                    'ImportanceStd': perm_result.importances_std
                }).sort_values('ImportanceMean', ascending=False)
                perm_path = os.path.join(output_dir, "permutation_importance.csv")
                perm_df.to_csv(perm_path, index=False)
                print(f"[INFO] Saved permutation importance to {perm_path}")
        except Exception as e:
            print(f"[WARN] Permutation importance failed: {e}")

    # -------------------------------
    # 2. SHAP Values
    # -------------------------------
    if 'shap' in [m.lower() for m in explain_methods]:
        print("[INFO] Computing SHAP values...")
        try:
            # Decide which prediction function to use
            if hasattr(clf, "predict_proba"):
                predict_func = clf.predict_proba
            elif hasattr(clf, "decision_function"):
                predict_func = clf.decision_function
            else:
                predict_func = clf.predict

            # Pick appropriate SHAP explainer
            model_name = str(type(clf)).lower()
            if "tree" in model_name:
                explainer = shap.TreeExplainer(clf)
            elif "linear" in model_name:
                explainer = shap.LinearExplainer(clf, X_train)
            else:
                # KernelExplainer fallback for SVM / KNN / etc.
                explainer = shap.KernelExplainer(predict_func, shap.sample(X_train, 50, random_state))

            # Compute SHAP values on a sample of the test set for speed
            shap_values = explainer.shap_values(shap.sample(X_test, 100, random_state))
            if isinstance(shap_values, list):
                shap_df = pd.DataFrame(shap_values[0], columns=X_test.columns)
            else:
                shap_df = pd.DataFrame(shap_values, columns=X_test.columns)

            shap_path = os.path.join(output_dir, "shap_values.csv")
            shap_df.to_csv(shap_path, index=False)
            print(f"[INFO] Saved SHAP values to {shap_path}")

        except Exception as e:
            print(f"[WARN] SHAP computation failed: {e}")

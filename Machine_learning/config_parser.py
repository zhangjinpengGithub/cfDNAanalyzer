import argparse

def get_args():
    parser = argparse.ArgumentParser(description="cfDNAanalyzer downstream analysis")

    # Basic control
    parser.add_argument('--noDA', action='store_true', help='Skip all the downstream analysis.')
    parser.add_argument('--modality', type=str, choices=['single', 'multi'], default='single', help='Analysis modality: single or multi')
    parser.add_argument('--input_dir', type=str, required=True, help='Input directory containing CSV files for single or multi modality')

    # Feature Selection: Filter
    parser.add_argument('--filterMethod', nargs='+', choices=['IG', 'CHI', 'FS', 'FCBF', 'PI', 'CC', 'LVF', 'MAD', 'DR', 'MI', 'RLF', 'SURF', 'MSURF', 'TRF'], default=None, help='Filter methods to use. Options: IG, CHI, FS, FCBF, PI, CC, LVF, MAD, DR, MI, RLF, SURF, MSURF, TRF')
    parser.add_argument('--filterFrac', type=float, default=0.2,help='Fraction of top features to retain based on importance ranking; if N ≥ 1 keep the top N features, if 0 < N < 1 keep the top N fraction of features (e.g. 0.1→10%).')
    # Feature Selection: Wrapper
    parser.add_argument('--wrapperMethod', nargs='+', choices=['BOR', 'FFS', 'BFS', 'EFS', 'RFS'], default=None, help='Wrapper methods to use. Options: LASSO, RIDGE, ELASTICNET, RF')
    parser.add_argument('--wrapperFrac', type=float, default=0.2, help='Fraction of top features to retain based on importance ranking; if N ≥ 1 keep the top N features, if 0 < N < 1 keep the top N fraction of features (e.g. 0.1→10%).')

    # Feature Selection: Embedded
    parser.add_argument('--embeddedMethod', nargs='+', choices=['LASSO', 'RIDGE', 'ELASTICNET', 'RF'], default=None, help='Embedded methods to use. Options: BOR, FFS, BFS, EFS, RFS')
    parser.add_argument('--embeddedFrac', type=float, default=0.2, help='Fraction of top features to retain based on importance ranking; if N ≥ 1 keep the top N features, if 0 < N < 1 keep the top N fraction of features (e.g. 0.1→10%).')

    # Feature Selection: Hybrid
    parser.add_argument('--hybridType', type=str, choices=['FE', 'FW'], default=None, help='Hybrid method type. Options: FE, FW')
    parser.add_argument('--hybridMethod1', type=str, default=None, help='Hybrid method 1')
    parser.add_argument('--hybridMethod2', type=str, default=None, help='Hybrid method 2')
    parser.add_argument('--hybridFrac1', type=float, default=0.2, help='Fraction of top features to retain in hybrid stage 1 based on importance ranking; if N ≥ 1 keep the top N features, if 0 < N < 1 keep the top N fraction (e.g. 0.1→10%).')
    parser.add_argument('--hybridFrac2', type=float, default=0.2, help='Fraction of top features to retain in hybrid stage 2 based on importance ranking; if N ≥ 1 keep the top N features, if 0 < N < 1 keep the top N fraction (e.g. 0.1→10%).')

    # ML Task Settings
    parser.add_argument('--classNum', type=int, default=2, help='Number of classification categories (2 for binary, >=3 for multi-class)')
    
    parser.add_argument('--cvSingle', type=str, choices=['LOO', 'KFold', 'Single', 'Independent'], default='LOO', help="CV method for single modality. Options: LOO, KFold, Single, Independent")
    parser.add_argument('--cvSingle_test_ratio', type=float, default=0.2, help="Test set ratio for Single train/test split mode")
    parser.add_argument('--nsplitSingle', type=int, default=5, help='Number of folds for KFold in single modality')
    parser.add_argument('--classifierSingle', nargs='+', choices=['KNN', 'SVM', 'RandomForest', 'GaussianNB', 'LogisticRegression', 'XGB', 'GBM'], default=None, help='Classifiers for single modality. Options: KNN, SVM, RandomForest, GaussianNB, LogisticRegression, XGB')

    parser.add_argument('--cvMulti', type=str, choices=['LOO', 'KFold', 'Single', 'Independent'], default='LOO', help='CV method for multi modality. Options: LOO, KFold, Single')
    parser.add_argument('--cvMulti_test_ratio', type=float, default=0.2, help="Test set ratio for Single train/test split mode")
    parser.add_argument('--nsplitMulti', type=int, default=5, help='Number of folds for KFold in multi modality')
    parser.add_argument('--classifierMulti', nargs='+', choices=['KNN', 'SVM', 'RandomForest', 'GaussianNB', 'LogisticRegression', 'XGB', 'GBM'], default=None, help='Classifiers for multi modality. Options: KNN, SVM, RandomForest, GaussianNB, LogisticRegression, XGB')
    parser.add_argument('--fusion_type', type=str, nargs='+', choices=['concat', 'model', 'trans'], default=None, help='Fusion strategy for multi-modality. Options: concat, model, trans (support multiple)')
    parser.add_argument('--model_method', type=str, nargs='+', choices=['average', 'weighted', 'majority', 'stack'], default=None, help='Model fusion methods: average, weighted, majority, stack (for fusion_type=model)')
    parser.add_argument('--trans_method', type=str, nargs='+', choices=['pca', 'linear', 'polynomial', 'rbf', 'sigmoid', 'snf'], default=None, help='Transformation methods for trans fusion: pca, linear, polynomial, rbf, sigmoid, snf')
    parser.add_argument('--external_input_dir', type=str, required=True)
    parser.add_argument('--explain', nargs='*', default=[], choices=['SHAP', 'perm'], help='Optional: compute SHAP values or permutation importance (time-consuming)')
    
    # Output control
    parser.add_argument('--DA_output_dir', type=str, default="Machine_Learning", help="Base directory for output files. Default is 'Machine_Learning'")

    args = parser.parse_args()

    print("\n===== Parameters Summary =====")
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")
    print("================================\n")

    return args

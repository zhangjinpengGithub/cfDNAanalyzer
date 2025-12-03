from sklearn.model_selection import LeaveOneOut, StratifiedKFold, train_test_split
import numpy as np

def get_cv(method, nsplit, y, test_ratio=0.2):
    """
    :param method: 'LOO' or 'KFold'
    :param nsplit: KFold's fold
    :param y: Labels for stratify
    """
    if method.upper() == 'LOO':
        return LeaveOneOut()
    elif method.upper() == 'KFOLD':
        return StratifiedKFold(n_splits=nsplit, shuffle=True, random_state=42)
    elif method.upper() == 'SINGLE':
        # Create a dummy CV object that returns ONE fold
        class SingleSplit:
            def split(self, X, y):
                train_idx, test_idx = train_test_split(
                    np.arange(len(y)),
                    test_size=test_ratio,
                    random_state=42,
                    stratify=y
                )
                yield train_idx, test_idx

        return SingleSplit()
    elif method.upper() == 'INDEPENDENT':
        class IndependentSplit:
            def split(self, X, y):
                idx = np.arange(len(y))
                yield idx, idx

            def get_n_splits(self, X=None, y=None, groups=None):
                return 1

        return IndependentSplit()

    else:
        raise ValueError("Supported CV: 'LOO', 'KFold', 'Single', 'Independent'")

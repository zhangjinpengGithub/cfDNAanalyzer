from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from xgboost import XGBClassifier


def get_classifier(name):
    # name = name.upper()
    classifiers = {
        'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
        'KNN': KNeighborsClassifier(n_neighbors=5),
        'GaussianNB': GaussianNB(),
        'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
        'SVM': svm.SVC(probability=True, kernel='rbf', random_state=42),
        'XGB': XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42),
        'GBM': GradientBoostingClassifier(
            n_estimators=150,           # 对应 R 中的 n.trees = 150
            max_depth=3,                # 对应 R 中的 interaction.depth = 3
            learning_rate=0.1,          # 对应 R 中的 shrinkage = 0.1
            min_samples_split=20,       # 近似对应 n.minobsinnode = 10
            min_samples_leaf=10,        # 直接对应 n.minobsinnode = 10
            random_state=42,
            subsample=1.0,              # 默认值，对应 R 中的 bag.fraction = 0.5（但你的R代码中没指定）
            verbose=0                   # 对应 R 中的 verbose = FALSE
        )
    }
    if name not in classifiers:
        raise ValueError(f"Unsupported classifier: {name}")
    return classifiers[name]

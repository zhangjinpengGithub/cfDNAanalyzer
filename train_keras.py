#!/usr/bin/env python3
"""
AutoKeras模型训练脚本
用于自动机器学习建模，支持分类任务
用法: python train_model.py --data 特征数据.csv --train 训练集分组.txt --valid 验证集分组.txt --model 模型名称
"""

import sys
import os
import json
import argparse
import logging
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import LabelEncoder, StandardScaler, MinMaxScaler
import autokeras as ak
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score

# 设置随机种子确保可重复性
def set_seed(seed=42):
    np.random.seed(seed)
    tf.random.set_seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    os.environ['TF_DETERMINISTIC_OPS'] = '1'

set_seed(42)

# 配置日志
def setup_logging(model_name):
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"{model_name}_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

# 数据加载和预处理函数
def load_and_preprocess(data_path, train_info_path, valid_info_path, 
                       scaler_type='minmax', fit_scaler=True):
    """
    加载和预处理数据
    """
    logger = logging.getLogger(__name__)
    
    # 检查文件是否存在
    for path in [data_path, train_info_path, valid_info_path]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"文件不存在: {path}")
    
    logger.info(f"加载特征数据: {data_path}")
    df = pd.read_csv(data_path)
    logger.info(f"原始数据形状: {df.shape}")
    
    # 处理缺失值
    df_clean = df.dropna(axis=1)
    dropped_cols = set(df.columns) - set(df_clean.columns)
    if dropped_cols:
        logger.info(f"删除了 {len(dropped_cols)} 个包含缺失值的列")

    # 进一步处理无穷大和极端异常值，避免 scaler 报错
    numeric_cols = df_clean.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 0:
        # 将 inf / -inf 替换为 NaN
        df_clean[numeric_cols] = df_clean[numeric_cols].replace([np.inf, -np.inf], np.nan)
        # 再次删除包含 NaN 的列（与前面策略一致）
        before_cols = set(df.columns)
        df_clean = df_clean.dropna(axis=1)
        after_cols = set(df_clean.columns)
        dropped_inf_cols = before_cols - after_cols - dropped_cols
        if dropped_inf_cols:
            logger.info(f"因包含 inf/-inf 或异常值额外删除了 {len(dropped_inf_cols)} 个数值列")
    
    # 加载分组信息
    logger.info(f"加载训练分组: {train_info_path}")
    train_info = pd.read_csv(train_info_path, sep='\t')
    valid_info = pd.read_csv(valid_info_path, sep='\t')
    
    # 确保必要的列存在
    required_cols = ['SampleID', 'Response']
    for info, name in [(train_info, '训练集'), (valid_info, '验证集')]:
        if not all(col in info.columns for col in required_cols):
            raise ValueError(f"{name}分组文件必须包含 {required_cols} 列")
    
    train_info = train_info[required_cols]
    valid_info = valid_info[required_cols]
    
    # 合并特征和标签
    train_df = pd.merge(df_clean, train_info, on='SampleID', how='inner')
    valid_df = pd.merge(df_clean, valid_info, on='SampleID', how='inner')
    
    logger.info(f"训练集样本数: {len(train_df)}")
    logger.info(f"验证集样本数: {len(valid_df)}")
    
    # 检查标签分布
    logger.info("训练集标签分布:")
    logger.info(train_df['Response'].value_counts().to_dict())
    logger.info("验证集标签分布:")
    logger.info(valid_df['Response'].value_counts().to_dict())
    
    # 分离特征和标签
    X_train = train_df.drop(['SampleID', 'Response'], axis=1)
    y_train = train_df['Response']
    sample_ids_train = train_df['SampleID'].tolist()
    
    X_valid = valid_df.drop(['SampleID', 'Response'], axis=1)
    y_valid = valid_df['Response']
    sample_ids_valid = valid_df['SampleID'].tolist()
    
    # 标签编码
    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)
    y_valid_encoded = label_encoder.transform(y_valid)
    
    logger.info(f"标签编码映射: {dict(zip(label_encoder.classes_, label_encoder.transform(label_encoder.classes_)))}")
    
    # 特征标准化
    if scaler_type.lower() == 'minmax':
        scaler = MinMaxScaler()
    elif scaler_type.lower() == 'standard':
        scaler = StandardScaler()
    else:
        raise ValueError(f"不支持的scaler类型: {scaler_type}")
    
    if fit_scaler:
        X_train_scaled = scaler.fit_transform(X_train)
        logger.info("拟合并转换训练集特征")
    else:
        X_train_scaled = scaler.transform(X_train)
    
    X_valid_scaled = scaler.transform(X_valid)
    logger.info("转换验证集特征")
    
    feature_names = X_train.columns.tolist()
    
    return (X_train_scaled, y_train_encoded, X_valid_scaled, y_valid_encoded,
           sample_ids_train, sample_ids_valid, label_encoder, scaler, feature_names)

# 训练AutoKeras模型
def train_autokeras_model(X_train, y_train, X_valid, y_valid, 
                         model_name, num_classes=2, max_trials=20, epochs=100):
    """
    训练AutoKeras模型
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"开始AutoKeras模型训练，最大尝试次数: {max_trials}, 最大轮数: {epochs}")
    logger.info(f"输入形状: {X_train.shape}, 类别数: {num_classes}")
    
    # 定义回调函数
    callbacks = [
        EarlyStopping(
            monitor='val_loss',
            patience=20,  # 增加耐心值
            restore_best_weights=True,
            verbose=1,
            mode='min'
        ),
        ReduceLROnPlateau(
            monitor='val_loss',
            factor=0.5,  # 学习率减少因子
            patience=10,  # 等待轮数
            min_lr=1e-6,  # 最小学习率
            verbose=1
        )
    ]
    
    # 初始化AutoKeras模型
    # 重要参数调整：
    # - max_trials: 搜索的不同模型架构数量
    # - overwrite: 是否覆盖之前的搜索结果
    # - tuner: 优化器选择（可选：'random', 'bayesian', 'hyperband', 'greedy'）
    # - seed: 随机种子
    model = ak.StructuredDataClassifier(
        max_trials=max_trials,
        overwrite=True,
        project_name=model_name + "_autokeras_project",
        directory="./autokeras_trials",  # 指定保存搜索结果的目录
        tuner='bayesian',  # 使用贝叶斯优化进行超参数调优
        seed=42,
        loss='categorical_crossentropy' if num_classes > 2 else 'binary_crossentropy',
        metrics=['accuracy'],
        objective='val_accuracy',
        max_model_size=None,  # 不限制模型大小
    )
    
    # 训练模型
    # 重要参数：
    # - validation_split: 从训练集中划分验证集的比例
    # - batch_size: 自动确定
    # - verbose: 日志详细程度
    history = model.fit(
        X_train,
        y_train,
        validation_data=(X_valid, y_valid),  # 使用独立的验证集
        epochs=epochs,
        callbacks=callbacks,
        verbose=2,  # 每个epoch显示进度条
        batch_size=min(32, len(X_train))  # 动态设置batch_size
    )
    
    # 评估模型
    logger.info("评估模型性能...")
    test_loss, test_accuracy = model.evaluate(X_valid, y_valid, verbose=0)
    logger.info(f"验证集损失: {test_loss:.4f}")
    logger.info(f"验证集准确率: {test_accuracy:.4f}")
    
    # 获取最佳模型
    best_model = model.export_model()
    
    return model, best_model, history, test_loss, test_accuracy

# 评估模型性能
def evaluate_model(model, X_test, y_test, label_encoder, sample_ids):
    """
    全面评估模型性能
    """
    logger = logging.getLogger(__name__)
    
    # 预测概率
    y_pred_proba = model.predict(X_test)
    
    # 对于二分类，确保形状正确
    if y_pred_proba.shape[1] == 1:
        y_pred = (y_pred_proba > 0.5).astype(int).flatten()
    else:
        y_pred = np.argmax(y_pred_proba, axis=1)
    
    # 分类报告
    logger.info("\n分类报告:")
    report = classification_report(y_test, y_pred, 
                                   target_names=label_encoder.classes_,
                                   output_dict=True)
    logger.info(classification_report(y_test, y_pred, 
                                      target_names=label_encoder.classes_))
    
    # 混淆矩阵
    logger.info("\n混淆矩阵:")
    cm = confusion_matrix(y_test, y_pred)
    logger.info(f"\n{cm}")
    
    # AUC计算（仅适用于二分类）
    if len(label_encoder.classes_) == 2:
        try:
            auc = roc_auc_score(y_test, y_pred_proba)
            logger.info(f"AUC分数: {auc:.4f}")
        except:
            auc = None
            logger.warning("无法计算AUC分数")
    else:
        auc = None
    
    # 保存预测结果
    results_df = pd.DataFrame({
        'SampleID': sample_ids,
        'True_Label': label_encoder.inverse_transform(y_test),
        'Predicted_Label': label_encoder.inverse_transform(y_pred),
        'Prediction_Score': y_pred_proba.max(axis=1) if y_pred_proba.shape[1] > 1 else y_pred_proba.flatten()
    })
    
    return report, cm, auc, results_df

# 保存模型和结果
def save_results(model, best_model, report, cm, auc, test_accuracy, 
                test_loss, results_df, feature_names, model_name, label_encoder,
                train_auc=None):
    """
    保存模型和所有结果
    """
    logger = logging.getLogger(__name__)
    
    # 创建输出目录
    output_dir = f"results_{model_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"保存结果到目录: {output_dir}")
    
    # 1. 保存模型
    model_path = os.path.join(output_dir, f"{model_name}.keras")
    try:
        # 尝试保存为Keras格式
        best_model.save(model_path)
        logger.info(f"模型保存为Keras格式: {model_path}")
    except:
        # 回退到h5格式
        model_path = os.path.join(output_dir, f"{model_name}.h5")
        best_model.save(model_path)
        logger.info(f"模型保存为h5格式: {model_path}")
    
    # 2. 保存AutoKeras完整模型（包含搜索空间）
    ak_model_path = os.path.join(output_dir, f"{model_name}_autokeras")
    model.export_model().save(ak_model_path, save_format='tf')
    logger.info(f"AutoKeras完整模型保存到: {ak_model_path}")
    
    # 3. 保存性能指标
    metrics = {
        'model_name': model_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'test_accuracy': float(test_accuracy),
        'test_loss': float(test_loss),
        'confusion_matrix': cm.tolist(),
        'classification_report': report,
        'feature_count': len(feature_names),
        'label_mapping': dict(zip(label_encoder.classes_.tolist(), 
                                 label_encoder.transform(label_encoder.classes_).tolist()))
    }
    
    # 验证集 AUC
    if auc is not None:
        metrics['valid_auc_score'] = float(auc)
        # 为兼容旧字段，保留 auc_score 作为验证集 AUC
        metrics['auc_score'] = float(auc)
    
    # 训练集 AUC（用于监控过拟合）
    if train_auc is not None:
        metrics['train_auc_score'] = float(train_auc)
    
    metrics_path = os.path.join(output_dir, 'performance_metrics.json')
    with open(metrics_path, 'w', encoding='utf-8') as f:
        json.dump(metrics, f, indent=2, ensure_ascii=False)
    
    # 4. 保存预测结果
    results_path = os.path.join(output_dir, 'predictions.csv')
    results_df.to_csv(results_path, index=False, encoding='utf-8')
    
    # 5. 保存特征列表
    features_path = os.path.join(output_dir, 'feature_names.txt')
    with open(features_path, 'w', encoding='utf-8') as f:
        for feature in feature_names:
            f.write(f"{feature}\n")
    
    # 6. 保存训练摘要
    summary_path = os.path.join(output_dir, 'training_summary.txt')
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(f"模型训练摘要\n")
        f.write("="*50 + "\n")
        f.write(f"模型名称: {model_name}\n")
        f.write(f"训练时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"验证集准确率: {test_accuracy:.4f}\n")
        f.write(f"验证集损失: {test_loss:.4f}\n")
        if auc is not None:
            f.write(f"验证集AUC分数: {auc:.4f}\n")
        if train_auc is not None:
            f.write(f"训练集AUC分数: {train_auc:.4f}\n")
        f.write(f"特征数量: {len(feature_names)}\n")
        f.write(f"标签类别: {', '.join(label_encoder.classes_.tolist())}\n")
    
    logger.info(f"所有结果已保存到: {output_dir}")
    return output_dir

# 主函数
def main():
    parser = argparse.ArgumentParser(description='AutoKeras模型训练脚本')
    parser.add_argument('--data', type=str, required=True, help='特征数据文件路径(CSV格式)')
    parser.add_argument('--train', type=str, required=True, help='训练集分组文件路径(TSV格式)')
    parser.add_argument('--valid', type=str, required=True, help='验证集分组文件路径(TSV格式)')
    parser.add_argument('--model', type=str, required=True, help='模型名称')
    parser.add_argument('--max_trials', type=int, default=20, help='AutoKeras最大尝试次数(默认: 20)')
    parser.add_argument('--epochs', type=int, default=100, help='最大训练轮数(默认: 100)')
    parser.add_argument('--scaler', type=str, default='minmax', 
                       choices=['minmax', 'standard'], help='特征标准化方法(默认: minmax)')
    parser.add_argument('--gpu', action='store_true', help='启用GPU训练(如果可用)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1, 2], 
                       help='日志详细程度: 0=静默, 1=进度条, 2=每个epoch详情(默认: 1)')
    
    args = parser.parse_args()
    
    # 设置GPU（如果可用）
    if args.gpu:
        gpus = tf.config.list_physical_devices('GPU')
        if gpus:
            try:
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
                print(f"使用GPU: {[gpu.name for gpu in gpus]}")
            except RuntimeError as e:
                print(f"GPU设置错误: {e}")
        else:
            print("警告: 未检测到GPU，将使用CPU训练")
    
    # 设置日志
    logger = setup_logging(args.model)
    logger.info("="*60)
    logger.info("开始AutoKeras模型训练")
    logger.info(f"参数: {vars(args)}")
    logger.info("="*60)
    
    try:
        # 1. 加载和预处理数据
        logger.info("步骤1: 加载和预处理数据")
        (X_train, y_train, X_valid, y_valid, 
         train_ids, valid_ids, label_encoder, scaler, feature_names) = load_and_preprocess(
            args.data, args.train, args.valid, args.scaler
        )
        
        logger.info(f"训练集形状: {X_train.shape}")
        logger.info(f"验证集形状: {X_valid.shape}")
        logger.info(f"特征数量: {len(feature_names)}")
        logger.info(f"类别数量: {len(label_encoder.classes_)}")
        
        # 2. 训练模型
        logger.info("\n步骤2: 训练AutoKeras模型")
        model, best_model, history, test_loss, test_accuracy = train_autokeras_model(
            X_train, y_train, X_valid, y_valid,
            args.model, 
            num_classes=len(label_encoder.classes_),
            max_trials=args.max_trials,
            epochs=args.epochs
        )
        
        # 3. 评估模型
        logger.info("\n步骤3: 评估模型性能")
        # 3.1 验证集评估
        report, cm, auc, results_df = evaluate_model(
            best_model, X_valid, y_valid, label_encoder, valid_ids
        )
        
        # 3.2 训练集 AUC（用于监控过拟合，仅二分类时有效）
        train_auc = None
        if len(label_encoder.classes_) == 2:
            logger.info("额外计算训练集 AUC 以监控过拟合")
            try:
                _, _, train_auc, _ = evaluate_model(
                    best_model, X_train, y_train, label_encoder, train_ids
                )
            except Exception as e:
                logger.warning(f"计算训练集 AUC 时出错: {e}")
        
        # 4. 保存结果
        logger.info("\n步骤4: 保存模型和结果")
        output_dir = save_results(
            model, best_model, report, cm, auc, test_accuracy, test_loss,
            results_df, feature_names, args.model, label_encoder,
            train_auc=train_auc
        )
        
        logger.info(f"\n{'='*60}")
        logger.info("模型训练完成!")
        logger.info(f"验证集准确率: {test_accuracy:.4f}")
        logger.info(f"验证集损失: {test_loss:.4f}")
        if auc is not None:
            logger.info(f"验证集 AUC 分数: {auc:.4f}")
        if 'train_auc' in locals() and train_auc is not None:
            logger.info(f"训练集 AUC 分数: {train_auc:.4f}")
        logger.info(f"所有结果保存在: {output_dir}")
        logger.info(f"{'='*60}")
        
        # 返回准确率到标准输出（便于脚本调用）
        print(f"ACCURACY:{test_accuracy:.4f}")
        
    except Exception as e:
        logger.error(f"训练过程中发生错误: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()

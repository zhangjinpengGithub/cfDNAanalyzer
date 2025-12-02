import os
import glob
import csv
import sys
import argparse
import logging
import json
import pandas as pd

def convert_tsv_to_txt(root_dir):
    tsv_files = glob.glob(os.path.join(root_dir, '**', '*.tsv'), recursive=True)
    for tsv_file in tsv_files:
        txt_file = os.path.splitext(tsv_file)[0] + '.txt'
        try:
            with open(tsv_file, 'r', encoding='utf-8') as f_in, \
                 open(txt_file, 'w', encoding='utf-8', newline='') as f_out:
                for line_number, line in enumerate(f_in, 1):
                    stripped_line = line.strip()
                    cleaned_line = stripped_line.replace(' ', '\t')
                    f_out.write(cleaned_line + '\n')
        except Exception as e:
            print(f"Failed to convert {tsv_file} to TXT: {e}")

def load_config(config_path):
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        return config
    except Exception as e:
        print(f"Error loading configuration file: {e}")
        sys.exit(1)

def load_labels(label_path):
    try:
        labels = {}
        with open(label_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                labels[row['sample']] = row['label']
        return labels
    except Exception as e:
        print(f"Error loading label file: {e}")
        sys.exit(1)

def collect_small_features_per_sample(feature, details, sample_name, ROOT_DIR):
    small_features = set()
    sample_dir = os.path.join(ROOT_DIR, sample_name)
    feature_file_pattern = os.path.join(sample_dir, details['path'])
    feature_files = glob.glob(feature_file_pattern)
    if not feature_files:
        print(f"No files found for feature {feature} in sample {sample_name}, pattern: {feature_file_pattern}")
        return small_features

    for file in feature_files:
        try:
            sep = details.get('sep', '\t')
            header = details.get('has_header', True)
            remove_tabs = details.get('remove_tabs', False)

            id_columns = []
            if 'site_name' in details['columns']:
                id_columns = ['site_name']
            elif details['type'] == 'single':
                if header:
                    with open(file, 'r') as f:
                        headers = f.readline().strip().split(sep)
                        for col in details['columns'].keys():
                            if col in headers:
                                small_features.add(details['columns'][col])
                            else:
                                print(f"    Column {col} not found in file {file}")
                else:
                    for col in details['columns'].values():
                        small_features.add(col)
                continue
            elif details['type'] == 'region':
                if 'id_column' in details:
                    id_columns = [details['id_column']]
                else:
                    chrom_col = details.get('chrom_column', 'chr')
                    id_columns = [chrom_col, 'start', 'end']
            elif details['type'] == 'custom':
                id_columns = details.get('id_columns', [])
            elif details['type'] == 'non_region':
                id_columns = details.get('id_columns', [])
                if not id_columns:
                    id_columns = [0]
            else:
                continue

            if remove_tabs:
                temp_file = file + '.tmp'
                with open(file, 'r') as f_in, open(temp_file, 'w') as f_out:
                    for line in f_in:
                        f_out.write(line.replace('\\t', ''))
                file_to_read = temp_file
            else:
                file_to_read = file

            with open(file_to_read, 'r') as f:
                if header:
                    headers = f.readline().strip().split(sep)
                    col_indices = []
                    for col in id_columns:
                        if col in headers:
                            col_indices.append(headers.index(col))
                        else:
                            print(f"    Column {col} not found in file {file}")
                    if not col_indices:
                        print(f"    Could not find valid column indices, skipping file {file}")
                        continue
                else:
                    col_indices = []
                    for col in id_columns:
                        if str(col).isdigit():
                            col_indices.append(int(col))
                        else:
                            print(f"    Non-header file {file} requires numeric column indices, found non-numeric {col}")
                    if not col_indices:
                        print(f"    Could not find valid column indices, skipping file {file}")
                        continue

                for line in f:
                    parts = line.strip().split(sep)
                    if col_indices and len(parts) >= max(col_indices) + 1:
                        feature_id = '_'.join([parts[i] for i in col_indices])
                        small_features.add(feature_id)
            if remove_tabs:
                os.remove(file_to_read)
        except Exception as e:
            print(f"Error processing file {file}: {e}")
    return small_features

def write_headers(feature, small_features, details, OUTPUT_DIR):
    csv_path = os.path.join(OUTPUT_DIR, f"{feature}.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = ['sample', 'label']
        if small_features:
            for sf in sorted(small_features):
                header.append(sf)
        else:
            for original_col, new_col in details['columns'].items():
                header.append(new_col)
        writer.writerow(header)

def format_number(value):
    """格式化数字，保留最多6位小数，不足时不补零"""
    try:
        num = float(value)
        if num.is_integer():
            return str(int(num))
        # 使用rstrip('0')去掉多余的0，然后确保最多6位小数
        formatted = "{0:.6f}".format(num).rstrip('0').rstrip('.') if '.' in "{0:.6f}".format(num) else str(num)
        return formatted
    except ValueError:
        return value

def process_sample_feature(sample_name, label, feature, details, small_features, ROOT_DIR, OUTPUT_DIR):
    csv_path = os.path.join(OUTPUT_DIR, f"{feature}.csv")
    sample_dir = os.path.join(ROOT_DIR, sample_name)
    feature_file_pattern = os.path.join(sample_dir, details['path'])
    feature_files = glob.glob(feature_file_pattern)
    coverage_dict = {}

    if small_features:
        for sf in sorted(small_features):
            coverage_dict[sf] = "NA"
    else:
        for original_col, new_col in details['columns'].items():
            coverage_dict[new_col] = "NA"

    if not feature_files:
        print(f"  No files found for feature {feature} in sample {sample_name}, pattern: {details['path']}")
    else:
        for file in feature_files:
            try:
                sep = details.get('sep', '\t')
                header = details.get('has_header', True)
                remove_tabs = details.get('remove_tabs', False)
                id_columns = []

                if 'site_name' in details['columns']:
                    id_columns = ['site_name']
                elif details['type'] == 'single':
                    with open(file, 'r') as f:
                        if header:
                            headers = f.readline().strip().split(sep)
                            data_line = f.readline().strip().split(sep)
                            for original_col, new_col in details['columns'].items():
                                if original_col in headers:
                                    idx = headers.index(original_col)
                                    value = data_line[idx] if idx < len(data_line) else "NA"
                                    coverage_dict[new_col] = format_number(value)
                                else:
                                    print(f"    Column {original_col} not found in file {file}")
                        else:
                            data_line = f.readline().strip().split(sep)
                            for idx, (original_col, new_col) in enumerate(details['columns'].items()):
                                value = data_line[idx] if idx < len(data_line) else "NA"
                                coverage_dict[new_col] = format_number(value)
                    continue
                elif details['type'] == 'region':
                    if 'id_column' in details:
                        id_columns = [details['id_column']]
                    else:
                        chrom_col = details.get('chrom_column', 'chr')
                        id_columns = [chrom_col, 'start', 'end']
                elif details['type'] == 'custom':
                    id_columns = details.get('id_columns', [])
                elif details['type'] == 'non_region':
                    id_columns = details.get('id_columns', [])
                    if not id_columns:
                        id_columns = [0]
                else:
                    pass

                if remove_tabs:
                    temp_file = file + '.tmp'
                    with open(file, 'r') as f_in, open(temp_file, 'w') as f_out:
                        for line in f_in:
                            f_out.write(line.replace('\\t', ''))
                    file_to_read = temp_file
                else:
                    file_to_read = file

                with open(file_to_read, 'r') as f:
                    if header:
                        headers = f.readline().strip().split(sep)
                        col_indices = []
                        for col in id_columns:
                            if col in headers:
                                col_indices.append(headers.index(col))
                            else:
                                print(f"    Column {col} not found in file {file}")
                        if not col_indices:
                            print(f"    Could not find valid column indices, skipping file {file}")
                            continue
                        data_indices = {}
                        for original_col, new_col in details['columns'].items():
                            if original_col in headers:
                                data_indices[original_col] = headers.index(original_col)
                            else:
                                print(f"    Column {original_col} not found in file {file}")
                    else:
                        col_indices = []
                        for col in id_columns:
                            if str(col).isdigit():
                                col_indices.append(int(col))
                            else:
                                print(f"    Non-header file requires numeric column indices, found non-numeric {col}")
                        if not col_indices:
                            print(f"    Could not find valid column indices, skipping file {file}")
                            continue
                        data_indices = {}
                        for original_col, new_col in details['columns'].items():
                            if str(original_col).isdigit():
                                data_indices[original_col] = int(original_col)
                            else:
                                print(f"    Non-header file requires numeric column indices, found non-numeric {original_col}")

                    for line in f:
                        parts = line.strip().split(sep)
                        if col_indices and len(parts) >= max(col_indices) + 1:
                            feature_id = '_'.join([parts[i] for i in col_indices])
                            if feature_id in small_features:
                                for original_col, new_col in details['columns'].items():
                                    if original_col not in id_columns and original_col not in ['site_name']:
                                        idx = data_indices.get(original_col, None)
                                        if idx is not None and idx < len(parts):
                                            coverage = parts[idx]
                                            coverage_dict[feature_id] = format_number(coverage)
                if remove_tabs:
                    os.remove(file_to_read)
            except Exception as e:
                print(f"    Error processing file {file}: {e}")

    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        row = [sample_name, label]
        if coverage_dict:
            for key in sorted(coverage_dict.keys()):
                row.append(coverage_dict.get(key, "NA"))
        else:
            for original_col, new_col in details['columns'].items():
                row.append("NA")
        writer.writerow(row)

def remove_all_na_samples(csv_file):
    """删除所有特征值都是NA的样本行"""
    try:
        # 读取CSV文件
        df = pd.read_csv(csv_file)
        
        if df.empty:
            return df
        
        # 找出特征列（排除sample和label列）
        feature_columns = [col for col in df.columns if col not in ['sample', 'label']]
        
        if not feature_columns:
            return df
        
        # 标记哪些行所有特征都是NA
        all_na_mask = df[feature_columns].isna().all(axis=1) | (df[feature_columns] == 'NA').all(axis=1)
        
        # 删除所有特征都是NA的行
        rows_before = len(df)
        df_cleaned = df[~all_na_mask]
        rows_after = len(df_cleaned)
        
        if rows_before > rows_after:
            print(f"  Removed {rows_before - rows_after} all-NA samples from {csv_file}")
            
            # 保存清理后的文件
            df_cleaned.to_csv(csv_file, index=False)
            return df_cleaned
        else:
            return df
            
    except Exception as e:
        print(f"Error processing {csv_file}: {e}")
        return None

def main():
    args = parse_arguments()

    CONFIG_FILE = args.config
    ROOT_DIR = args.root
    LABEL_FILE = os.path.join(ROOT_DIR, 'label.txt')
    OUTPUT_DIR = os.path.join(ROOT_DIR, args.output)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    config = load_config(CONFIG_FILE)
    labels_dict = load_labels(LABEL_FILE)

    convert_tsv_to_txt(ROOT_DIR)

    for feature, details in config.items():
        print(f"Processing feature: {feature}")
        sample_to_small_features = {}
        valid_samples = []  # 记录有该feature的样本
        
        for sample_name in labels_dict.keys():
            small_features = collect_small_features_per_sample(feature, details, sample_name, ROOT_DIR)
            sample_to_small_features[sample_name] = small_features
            # 如果样本有该feature，添加到valid_samples
            if small_features or details['type'] == 'single':
                valid_samples.append(sample_name)

        # 修改这里的逻辑：使用所有样本的交集，或者至少使用有该feature的样本的交集
        if valid_samples:
            # 只使用有该feature的样本的交集
            valid_features_sets = [sample_to_small_features[sample] for sample in valid_samples]
            if valid_features_sets:
                all_small_features = set.intersection(*valid_features_sets)
            else:
                all_small_features = set()
        else:
            all_small_features = set()

        # 即使没有small features，对于single类型的feature仍然继续
        if not all_small_features and details['type'] != 'single':
            print(f"Feature {feature} has no common small features, but will continue with NA values.")
            all_small_features = set()  # 设为空集但继续处理

        write_headers(feature, all_small_features, details, OUTPUT_DIR)

        # 处理所有样本，包括那些没有该feature的样本
        for sample_name, label in labels_dict.items():
            print(f"  Processing sample: {sample_name}")
            process_sample_feature(sample_name, label, feature, details, all_small_features, ROOT_DIR, OUTPUT_DIR)
        
        # 在每个feature处理完成后，删除所有特征值都是NA的样本
        csv_file = os.path.join(OUTPUT_DIR, f"{feature}.csv")
        if os.path.exists(csv_file):
            print(f"  Cleaning all-NA samples from {feature}.csv")
            remove_all_na_samples(csv_file)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Integrate feature coverage into individual CSV files.")
    parser.add_argument('--config', type=str, default='config.json', help='Path to the configuration JSON file.')
    parser.add_argument('--root', type=str, required=True, help='Input directory containing sample folders.')
    parser.add_argument('--output', type=str, default='features_csv', help='Output directory for the CSV files.')
    return parser.parse_args()

if __name__ == "__main__":
    main()
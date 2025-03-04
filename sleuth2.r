#!/usr/bin/env Rscript

##############################
# 加載所需的庫
library(sleuth)
library(dplyr)

# 創建樣本表
sample_table <- data.frame(
  sample = c("L_1", "L_2", "L_3", "S_1", "S_2", "S_3"),
  condition = c("L", "L", "L", "S", "S", "S"),
  path = c(
    "/work/u1432917/output_L1_root/L1abundance.h5",
    "/work/u1432917/output_L2_root/L2abundance.h5",
    "/work/u1432917/output_L3_root/L3abundance.h5",
    "/work/u1432917/output_S1_root/S1abundance.h5",
    "/work/u1432917/output_S2_root/S2abundance.h5",
    "/work/u1432917/output_S3_root/S3abundance.h5"
  )
)

# 創建sleuth對象
so <- sleuth_prep(sample_table, extra_bootstrap_summary = TRUE)

# 擬合模型
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')

# 執行似然比檢驗
so <- sleuth_lrt(so, 'reduced', 'full')

# 獲取結果
results <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# 載入所需套件
library(dplyr)

# 文件路徑
L1_txt <- "/work/u1432917/output_L1_root/L1abundance.tsv"
L2_txt <- "/work/u1432917/output_L2_root/L2abundance.tsv"
L3_txt <- "/work/u1432917/output_L3_root/L3abundance.tsv"
S1_txt <- "/work/u1432917/output_S1_root/S1abundance.tsv"
S2_txt <- "/work/u1432917/output_S2_root/S2abundance.tsv"
S3_txt <- "/work/u1432917/output_S3_root/S3abundance.tsv"
N1_txt <- "/work/u1432917/output_N1_root/N1abundance.tsv"
N2_txt <- "/work/u1432917/output_N2_root/N2abundance.tsv"
N3_txt <- "/work/u1432917/output_N3_root/N3abundance.tsv"
I1_txt <- "/work/u1432917/output_I1_root/I1abundance.tsv"
I2_txt <- "/work/u1432917/output_I2_root/I2abundance.tsv"
I3_txt <- "/work/u1432917/output_I3_root/I3abundance.tsv"


# 定義讀取 Kallisto .txt 檔案並提取 TPM 數據的函數
# 定義讀取 Kallisto .txt 檔案並提取 target_id 和 TPM 數據的函數
read_tpm_with_id <- function(file_path) {
  data <- read.table(file_path, header = TRUE)  # 讀取 .txt 文件
  return(data[, c("target_id", "tpm")])  # 提取 'target_id' 和 'tpm' 欄位
}

# 讀取每個樣本的 TPM 數據
L1_tpm <- read_tpm_with_id(L1_txt)
L2_tpm <- read_tpm_with_id(L2_txt)
L3_tpm <- read_tpm_with_id(L3_txt)
S1_tpm <- read_tpm_with_id(S1_txt)
S2_tpm <- read_tpm_with_id(S2_txt)
S3_tpm <- read_tpm_with_id(S3_txt)
N1_tpm <- read_tpm_with_id(N1_txt)
N2_tpm <- read_tpm_with_id(N2_txt)
N3_tpm <- read_tpm_with_id(N3_txt)
I1_tpm <- read_tpm_with_id(I1_txt)
I2_tpm <- read_tpm_with_id(I2_txt)
I3_tpm <- read_tpm_with_id(I3_txt)

# 合併所有樣本的 TPM 數據
all_tpm <- full_join(L1_tpm, L2_tpm, by = "target_id", suffix = c("_L1", "_L2")) %>%
  full_join(L3_tpm, by = "target_id") %>%
  rename(L3_tpm = tpm) %>%
  full_join(S1_tpm, by = "target_id") %>%
  rename(S1_tpm = tpm) %>%
  full_join(S2_tpm, by = "target_id") %>%
  rename(S2_tpm = tpm) %>%
  full_join(S3_tpm, by = "target_id") %>%
  rename(S3_tpm = tpm)  %>%
    full_join(N1_tpm, by = "target_id") %>%
  rename(N1_tpm = tpm) %>%
  full_join(N2_tpm, by = "target_id") %>%
  rename(N2_tpm = tpm) %>%
  full_join(N3_tpm, by = "target_id") %>%
  rename(N3_tpm = tpm) %>%
  full_join(I1_tpm, by = "target_id") %>%
  rename(I1_tpm = tpm)   %>%
    full_join(I2_tpm, by = "target_id") %>%
  rename(I2_tpm = tpm)  %>%
    full_join(I3_tpm, by = "target_id") %>%
  rename(I3_tpm = tpm)











# 將 TPM 數據加入 sleuth_table
results <- results %>%
  left_join(all_tpm, by = "target_id")


# 過濾顯著的結果 (例如，q值 < 0.05)
significant_results <- results %>% filter(pval < 0.05)
significant05_results <- results %>% filter(pval < 0.05)
significant01_results <- results %>% filter(pval < 0.01)
significant001_results <- results %>% filter(pval < 0.001)

# 查看結果
print(head(significant_results))

# 保存結果

write.csv(significant_results, "differential_expression_results.csv", row.names = FALSE)
write.csv(significant05_results, "differential_expression_results05.csv", row.names = FALSE)
write.csv(significant01_results, "differential_expression_results01.csv", row.names = FALSE)
write.csv(significant001_results, "differential_expression_results001.csv", row.names = FALSE)


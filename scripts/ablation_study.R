## ==================== TCGA内部验证集的消融研究 ====================



cat("\n开始TCGA内部验证集消融研究分析...\n")

# 读取训练好的模型系数
model_coef <- read.csv("train_model_coefficients.csv", header = TRUE)
model_genes <- as.character(model_coef$gene)

cat("模型包含基因:", paste(model_genes, collapse = ", "), "\n")

# 根据Lasso系数绝对值确定最重要的基因
model_coef$abs_coef <- abs(model_coef$coefficient)
model_coef <- model_coef[order(-model_coef$abs_coef), ]

cat("\n基因重要性排名 (按系数绝对值):\n")
for(i in 1:nrow(model_coef)) {
  cat(sprintf("%d. %s: %.4f\n", i, model_coef$gene[i], model_coef$coefficient[i]))
}

# 选择前两个最重要的基因进行消融研究
genes_to_remove <- model_coef$gene[1:2]
cat("\n选择进行消融研究的基因:", paste(genes_to_remove, collapse = " 和 "), "\n")

# 读取TCGA测试集数据
test_data <- read.csv('test_muticoxdata.csv', row.names = 1)

# 检查数据基本信息
cat("\n=== TCGA测试集基本信息 ===\n")
cat("样本数:", nrow(test_data), "\n")
cat("事件数:", sum(test_data$status), "(", round(mean(test_data$status)*100, 1), "%)\n")
cat("时间范围:", round(range(test_data$time, na.rm = TRUE), 2), "年\n")
cat("可用模型基因:", sum(model_genes %in% colnames(test_data)), "/", length(model_genes), "\n")

# 风险评分计算函数
calculate_risk_score <- function(data, coefficients) {
  risk_scores <- rep(0, nrow(data))
  for(i in 1:nrow(coefficients)) {
    gene <- coefficients$gene[i]
    coef_value <- coefficients$coefficient[i]
    if(gene %in% colnames(data)) {
      risk_scores <- risk_scores + coef_value * data[[gene]]
    }
  }
  return(risk_scores)
}

## ==================== 消融分析执行 ====================

# 1. 完整模型（基准）
cat("\n1. 计算完整模型性能...\n")
full_risk_scores <- calculate_risk_score(test_data, model_coef)
full_roc_1yr <- timeROC(T = test_data$time, 
                        delta = test_data$status,
                        marker = full_risk_scores, 
                        times = 1, cause = 1, iid = TRUE)
full_roc_2yr <- timeROC(T = test_data$time, 
                        delta = test_data$status,
                        marker = full_risk_scores, 
                        times = 2, cause = 1, iid = TRUE)
full_roc_3yr <- timeROC(T = test_data$time, 
                        delta = test_data$status,
                        marker = full_risk_scores, 
                        times = 3, cause = 1, iid = TRUE)

# 2. 移除最重要的基因 (Gene_top1)
cat("2. 计算移除", genes_to_remove[1], "后的性能...\n")
coef_no_gene1 <- model_coef[model_coef$gene != genes_to_remove[1], ]
risk_scores_no_gene1 <- calculate_risk_score(test_data, coef_no_gene1)
roc_no_gene1_1yr <- timeROC(T = test_data$time, 
                            delta = test_data$status,
                            marker = risk_scores_no_gene1, 
                            times = 1, cause = 1, iid = TRUE)
roc_no_gene1_2yr <- timeROC(T = test_data$time, 
                            delta = test_data$status,
                            marker = risk_scores_no_gene1, 
                            times = 2, cause = 1, iid = TRUE)
roc_no_gene1_3yr <- timeROC(T = test_data$time, 
                            delta = test_data$status,
                            marker = risk_scores_no_gene1, 
                            times = 3, cause = 1, iid = TRUE)

# 3. 移除前两个最重要的基因 (Gene_top1 + Gene_top2)
cat("3. 计算移除", paste(genes_to_remove, collapse = "和"), "后的性能...\n")
coef_no_gene12 <- model_coef[!model_coef$gene %in% genes_to_remove, ]
risk_scores_no_gene12 <- calculate_risk_score(test_data, coef_no_gene12)
roc_no_gene12_1yr <- timeROC(T = test_data$time, 
                             delta = test_data$status,
                             marker = risk_scores_no_gene12, 
                             times = 1, cause = 1, iid = TRUE)
roc_no_gene12_2yr <- timeROC(T = test_data$time, 
                             delta = test_data$status,
                             marker = risk_scores_no_gene12, 
                             times = 2, cause = 1, iid = TRUE)
roc_no_gene12_3yr <- timeROC(T = test_data$time, 
                             delta = test_data$status,
                             marker = risk_scores_no_gene12, 
                             times = 3, cause = 1, iid = TRUE)

## ==================== 结果汇总 ====================

# 创建结果汇总表
ablation_results <- data.frame(
  Model = c(paste0("完整模型 (", length(model_genes), "基因)"),
            paste0("移除 ", genes_to_remove[1], " (", length(model_genes)-1, "基因)"),
            paste0("移除 ", paste(genes_to_remove, collapse = "+"), " (", length(model_genes)-2, "基因)")),
  AUC_1yr = c(full_roc_1yr$AUC[2], roc_no_gene1_1yr$AUC[2], roc_no_gene12_1yr$AUC[2]),
  AUC_2yr = c(full_roc_2yr$AUC[2], roc_no_gene1_2yr$AUC[2], roc_no_gene12_2yr$AUC[2]),
  AUC_3yr = c(full_roc_3yr$AUC[2], roc_no_gene1_3yr$AUC[2], roc_no_gene12_3yr$AUC[2])
)

cat("\n" + "="*60 + "\n")
cat("            TCGA内部验证集消融研究结果汇总\n")
cat("="*60 + "\n")
print(ablation_results)

# 保存结果
write.csv(ablation_results, "TCGA_ablation_study_results.csv", row.names = FALSE)

## ==================== 性能变化分析 ====================

cat("\n" + "="*60 + "\n")
cat("                    性能变化分析\n")
cat("="*60 + "\n")

# 计算AUC变化
calculate_performance_change <- function(full_auc, reduced_auc) {
  change <- reduced_auc - full_auc
  percent_change <- (change / full_auc) * 100
  return(list(change = change, percent_change = percent_change))
}

cat("\n移除", genes_to_remove[1], "后的性能变化:\n")
change_1yr_gene1 <- calculate_performance_change(full_roc_1yr$AUC[2], roc_no_gene1_1yr$AUC[2])
change_2yr_gene1 <- calculate_performance_change(full_roc_2yr$AUC[2], roc_no_gene1_2yr$AUC[2])
change_3yr_gene1 <- calculate_performance_change(full_roc_3yr$AUC[2], roc_no_gene1_3yr$AUC[2])

cat(sprintf("  1年AUC: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
            full_roc_1yr$AUC[2], roc_no_gene1_1yr$AUC[2],
            change_1yr_gene1$change, change_1yr_gene1$percent_change))
cat(sprintf("  2年AUC: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
            full_roc_2yr$AUC[2], roc_no_gene1_2yr$AUC[2],
            change_2yr_gene1$change, change_2yr_gene1$percent_change))
cat(sprintf("  3年AUC: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
            full_roc_3yr$AUC[2], roc_no_gene1_3yr$AUC[2],
            change_3yr_gene1$change, change_3yr_gene1$percent_change))

cat("\n移除", paste(genes_to_remove, collapse = "和"), "后的性能变化:\n")
change_1yr_gene12 <- calculate_performance_change(full_roc_1yr$AUC[2], roc_no_gene12_1yr$AUC[2])
change_2yr_gene12 <- calculate_performance_change(full_roc_2yr$AUC[2], roc_no_gene12_2yr$AUC[2])
change_3yr_gene12 <- calculate_performance_change(full_roc_3yr$AUC[2], roc_no_gene12_3yr$AUC[2])

cat(sprintf("  1年AUC: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
            full_roc_1yr$AUC[2], roc_no_gene12_1yr$AUC[2],
            change_1yr_gene12$change, change_1yr_gene12$percent_change))
cat(sprintf("  2年AUC: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
            full_roc_2yr$AUC[2], roc_no_gene12_2yr$AUC[2],
            change_2yr_gene12$change, change_2yr_gene12$percent_change))
cat(sprintf("  3年AUC: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
            full_roc_3yr$AUC[2], roc_no_gene12_3yr$AUC[2],
            change_3yr_gene12$change, change_3yr_gene12$percent_change))

## ==================== 可视化结果 ====================

cat("\n生成可视化结果...\n")

library(ggplot2)
library(reshape2)

# 准备绘图数据
ablation_melted <- melt(ablation_results, id.vars = "Model", 
                        variable.name = "Time", value.name = "AUC")
ablation_melted$Time <- factor(ablation_melted$Time, 
                               levels = c("AUC_1yr", "AUC_2yr", "AUC_3yr"),
                               labels = c("1年", "2年", "3年"))

# 绘制消融研究结果图
ablation_plot <- ggplot(ablation_melted, aes(x = Time, y = AUC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = round(AUC, 3)), 
            position = position_dodge(0.8), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "TCGA内部验证集 - 消融研究分析",
       subtitle = paste("移除基因:", paste(genes_to_remove, collapse = ", ")),
       x = "随访时间",
       y = "AUC值",
       fill = "模型版本",
       caption = paste("TCGA测试集: n =", nrow(test_data), "个样本")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  ylim(0, 1)

print(ablation_plot)
ggsave("TCGA_ablation_study_plot.png", width = 10, height = 6, dpi = 300)

## ==================== 模型稳健性评估 ====================

cat("\n" + "="*60 + "\n")
cat("                    模型稳健性评估\n")
cat("="*60 + "\n")

# 计算性能保留百分比
performance_retention <- function(original_auc, new_auc) {
  return((new_auc / original_auc) * 100)
}

cat("\n性能保留分析:\n")
retention_1yr <- performance_retention(full_roc_1yr$AUC[2], roc_no_gene12_1yr$AUC[2])
retention_2yr <- performance_retention(full_roc_2yr$AUC[2], roc_no_gene12_2yr$AUC[2])
retention_3yr <- performance_retention(full_roc_3yr$AUC[2], roc_no_gene12_3yr$AUC[2])

cat(sprintf("移除前2个基因后的性能保留:\n"))
cat(sprintf("  1年: %.1f%%\n", retention_1yr))
cat(sprintf("  2年: %.1f%%\n", retention_2yr))
cat(sprintf("  3年: %.1f%%\n", retention_3yr))

# 稳健性判断
cat("\n稳健性结论:\n")
avg_retention <- mean(c(retention_1yr, retention_2yr, retention_3yr))
cat(sprintf("平均性能保留: %.1f%%\n", avg_retention))

if(avg_retention > 90) {
  cat("✅ 模型非常稳健：移除两个最重要基因后，平均性能保留超过90%\n")
} else if(avg_retention > 80) {
  cat("✅ 模型较为稳健：移除两个最重要基因后，平均性能保留超过80%\n")
} else if(avg_retention > 70) {
  cat("⚠️  模型稳健性一般：移除两个最重要基因后，平均性能保留超过70%\n")
} else {
  cat("❌ 模型稳健性较差：移除两个最重要基因后，性能下降较多\n")
}

## ==================== 生成完整报告 ====================

cat("\n" + "="*60 + "\n")
cat("                    TCGA消融研究完整报告\n")
cat("="*60 + "\n")

cat("分析数据集: TCGA内部验证集 (严格分割的测试集)\n")
cat("样本数量:", nrow(test_data), "\n")
cat("分析基因:", length(model_genes), "个\n")
cat("目标基因:", paste(genes_to_remove, collapse = ", "), "\n\n")

cat("主要发现:\n")
cat("1. 完整模型性能: 1年=", round(full_roc_1yr$AUC[2], 3), 
    ", 2年=", round(full_roc_2yr$AUC[2], 3), 
    ", 3年=", round(full_roc_3yr$AUC[2], 3), "\n")
cat("2. 移除前2基因后平均性能保留:", round(avg_retention, 1), "%\n")
cat("3. 最关键的时间点(2年)性能保留:", round(retention_2yr, 1), "%\n\n")

cat("结论:\n")
if(avg_retention > 85) {
  cat("该多基因特征在TCGA内部验证集上表现出优异的稳健性，\n")
  cat("不依赖于单一或少数基因，证明其生物学合理性和可靠性。\n")
} else {
  cat("该特征对关键基因存在一定依赖性，建议进一步优化或\n")
  cat("在更大样本中验证其稳健性。\n")
}

# 保存详细结果
detailed_results <- data.frame(
  Time = rep(c("1年", "2年", "3年"), each = 3),
  Model = rep(ablation_results$Model, 3),
  AUC = c(ablation_results$AUC_1yr, ablation_results$AUC_2yr, ablation_results$AUC_3yr)
)

write.csv(detailed_results, "TCGA_ablation_detailed_results.csv", row.names = FALSE)

cat("\n" + "="*60 + "\n")
cat("TCGA内部验证集消融研究完成！\n")
cat("结果文件已保存:\n")
cat("- TCGA_ablation_study_results.csv\n")
cat("- TCGA_ablation_detailed_results.csv\n")
cat("- TCGA_ablation_study_plot.png\n")
cat("="*60 + "\n")






##重新
## ==================== 消融研究数据诊断 ====================

cat("\n=== 消融研究数据一致性诊断 ===\n")

# 1. 检查当前使用的数据
cat("1. 当前消融研究使用的数据:\n")
cat("数据文件: test_data\n")
cat("样本数:", nrow(test_data), "\n")
cat("时间范围:", range(test_data$time), "\n")
cat("事件率:", round(mean(test_data$status)*100, 1), "%\n")

# 2. 重新计算完整模型的AUC进行验证
cat("\n2. 重新验证完整模型AUC:\n")

# 使用与之前训练集相同的计算方法
full_risk_scores_verify <- predict(step_cox, type = "risk", newdata = test_data)

roc_verify_1yr <- timeROC(T = test_data$time, 
                          delta = test_data$status,
                          marker = full_risk_scores_verify, 
                          times = 1, cause = 1, iid = TRUE)

roc_verify_2yr <- timeROC(T = test_data$time, 
                          delta = test_data$status,
                          marker = full_risk_scores_verify, 
                          times = 2, cause = 1, iid = TRUE)

roc_verify_3yr <- timeROC(T = test_data$time, 
                          delta = test_data$status,
                          marker = full_risk_scores_verify, 
                          times = 3, cause = 1, iid = TRUE)

cat("验证计算 - 完整模型AUC:\n")
cat(sprintf("  1年: %.3f\n", roc_verify_1yr$AUC[2]))
cat(sprintf("  2年: %.3f\n", roc_verify_2yr$AUC[2]))
cat(sprintf("  3年: %.3f\n", roc_verify_3yr$AUC[2]))

# 3. 对比您之前的结果
cat("\n3. 数据对比:\n")
cat("您之前报告的测试集AUC: 0.704, 0.674, 0.656\n")
cat("当前计算的测试集AUC:", paste(round(c(roc_verify_1yr$AUC[2], roc_verify_2yr$AUC[2], roc_verify_3yr$AUC[2]), 3), collapse = ", "), "\n")

# 4. 检查风险评分分布
cat("\n4. 风险评分分布:\n")
cat("完整模型风险评分范围:", round(range(full_risk_scores_verify), 3), "\n")
cat("风险评分中位数:", round(median(full_risk_scores_verify), 3), "\n")

# 5. 检查关键基因的表达
cat("\n5. 关键基因表达情况:\n")
if("SPC25" %in% colnames(test_data)) {
  cat("SPC25表达范围:", round(range(test_data$SPC25), 3), "\n")
}
if("KITLG" %in% colnames(test_data)) {
  cat("KITLG表达范围:", round(range(test_data$KITLG), 3), "\n")
}

# 6. 重新运行可靠的消融研究
cat("\n6. 重新运行消融研究...\n")

reliable_ablation <- function(test_data, step_cox) {
  # 完整模型
  full_scores <- predict(step_cox, type = "risk", newdata = test_data)
  
  # 移除SPC25
  test_data_no_spc25 <- test_data
  test_data_no_spc25$SPC25 <- mean(test_data_no_spc25$SPC25)  # 设置为均值而非移除
  scores_no_spc25 <- predict(step_cox, type = "risk", newdata = test_data_no_spc25)
  
  # 移除SPC25和KITLG
  test_data_no_both <- test_data
  test_data_no_both$SPC25 <- mean(test_data_no_both$SPC25)
  test_data_no_both$KITLG <- mean(test_data_no_both$KITLG)
  scores_no_both <- predict(step_cox, type = "risk", newdata = test_data_no_both)
  
  # 计算AUC
  calculate_auc <- function(scores, time, status) {
    c(
      timeROC(T = time, delta = status, marker = scores, times = 1, cause = 1)$AUC[2],
      timeROC(T = time, delta = status, marker = scores, times = 2, cause = 1)$AUC[2],
      timeROC(T = time, delta = status, marker = scores, times = 3, cause = 1)$AUC[2]
    )
  }
  
  results <- list()
  results$full <- calculate_auc(full_scores, test_data$time, test_data$status)
  results$no_spc25 <- calculate_auc(scores_no_spc25, test_data$time, test_data$status)
  results$no_both <- calculate_auc(scores_no_both, test_data$time, test_data$status)
  
  return(results)
}

# 运行可靠的消融分析
reliable_results <- reliable_ablation(test_data, step_cox)

cat("\n=== 可靠的消融研究结果 ===\n")
cat("完整模型AUC (1,2,3年):", paste(round(reliable_results$full, 3), collapse = ", "), "\n")
cat("排除SPC25后AUC (1,2,3年):", paste(round(reliable_results$no_spc25, 3), collapse = ", "), "\n")
cat("排除SPC25+KITLG后AUC (1,2,3年):", paste(round(reliable_results$no_both, 3), collapse = ", "), "\n")

# 计算性能变化
cat("\n=== 性能变化分析 ===\n")
calculate_change <- function(full, reduced) {
  change <- reduced - full
  percent <- (change / full) * 100
  return(list(change = change, percent = percent))
}

cat("排除SPC25后的变化:\n")
for(i in 1:3) {
  change <- calculate_change(reliable_results$full[i], reliable_results$no_spc25[i])
  cat(sprintf("  %d年: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
              i, reliable_results$full[i], reliable_results$no_spc25[i],
              change$change, change$percent))
}

cat("\n排除SPC25+KITLG后的变化:\n")
for(i in 1:3) {
  change <- calculate_change(reliable_results$full[i], reliable_results$no_both[i])
  cat(sprintf("  %d年: %.3f -> %.3f (变化: %+.3f, %+.1f%%)\n", 
              i, reliable_results$full[i], reliable_results$no_both[i],
              change$change, change$percent))
}



library(tidyverse)
library(survival)
install.packages('survminer')
library(survminer)
install.packages('timeROC')
library(timeROC)
library(glmnet)
install.packages('forestplot')
library(forestplot)
BiocManager::install('ComplexHeatmap')
library(ComplexHeatmap)

setwd('C:\\SYH\\Study\\ccRCC\\KIRC\\2.cox')
unicox1 <- read.csv('coxinput.csv', header = T, row.names = 1) %>% 
  na.omit()
unicox1$time <- unicox1$time/365

##Cox;X9i
set.seed(20240612)
train.Idx <- sample(1:dim(unicox1)[1], floor(1/2*dim(unicox1)[1]))
test.Idx <- setdiff(1:dim(unicox1)[1], train.Idx)
x.train <- unicox1[train.Idx ,]
x.test <- unicox1[test.Idx ,]

##5%RrKXCox;X9i
library(survival)
mergeEXP = log2(x.train[,3:ncol(x.train)]+1)
merge = cbind(x.train[,1:2],mergeEXP)

coxR = data.frame()
coxf <- function(x){
  fmla1 <- as.formula(Surv(time,status) ~ merge[, x]) 
  mycox <- coxph(fmla1,data=merge)} 

for(a in colnames(merge[,3:ncol(merge)])){
  mycox = coxf(a)
  coxResult = summary(mycox)
  coxR = rbind(coxR,cbind(mergename=a,
                        HR = as.numeric(coxResult$coefficients[,"exp(coef)"])[1],
                        z = as.numeric(coxResult$coefficients[,"z"])[1],
                        pvalue = as.numeric(coxResult$coefficients[,"Pr(>|z|)"])[1],
                        lower = as.numeric(coxResult$conf.int[,3][1]),
                        upper = as.numeric(coxResult$conf.int[,4][1])))  
}

##I8Q!P<0.05
coxRR <- coxR %>% 
  arrange(., coxR$pvalue) %>% 
  filter(., pvalue < 0.05)
coxRR1 = sapply(coxRR[2:6], as.numeric)
coxRR = cbind(coxRR[1], coxRR1) %>% 
  head(., n = 20)

##I-AVM<
library(forestplot)
hz <- paste(round(coxRR$HR,3),
            "(",round(coxRR$lower,3),
            "-",round(coxRR$upper,3),")",sep = "")

tabletext <- cbind(c(NA,NA,coxRR$mergename),
                   c(NA,"pvalue",ifelse(coxRR$pvalue < 0.001,"P < 0.001", round(coxRR$pvalue,3))),
                   c(NA,"Hazard Ratio",hz))

forestplot(labeltext = tabletext, 
           graph.pos = 4,  #PvalueOdO_KyTZN;VC
           col = fpColors(box = "red", lines = "darkblue", zero = "gray50"),
           mean = c(NA,NA,coxRR$HR),
           lower = c(NA,NA,coxRR$lower), #95%VCPEGx<dOBO^
           upper = c(NA,NA,coxRR$upper), #95%VCPEGx<dIOO^
           boxsize = 0.3,
           lwd.ci = 0.7,   #OdWS5D4sP!#,O_5D?m6H
           ci.vertices = TRUE, #VCPEGx<dSCO_!"8_!"PM
           zero = 1,
           colgap = unit(3,"mm"),    #AP<dO6
           xticks = c(0, 0.5, 1,1.5, 2.0), #:aWx1j?L6H
           lineheight = unit(0.8,"cm"), #9L6(PP8_
           graphwidth = unit(.2,"npc"), #M<TZ1mVP5D?m6H1H@}
           cex = 1.5, 
           mar = unit(rep(0.5, times = 4), "cm"),#M<PNR31_>`
           txt_gp = fpTxtGp(label=gpar(cex = 0.5),
                            ticks=gpar(cex = 0.5),
                            xlab=gpar(cex = 0.5)),
           xlab="Hazard Ratio") 


##:sPxlasso
coxR1 <- coxR %>% 
  filter(pvalue < 0.05)

trainlassodat <- merge[, colnames(merge) %in% coxR1$mergename]
trainlassodat1 <- cbind(merge[,1:2], trainlassodat) %>% 
  na.omit()
write.csv(trainlassodat1, 'train_lasso_data.csv')

##lasso;X9iGzO_
v1 <- as.matrix(trainlassodat1[,c(3:ncol(trainlassodat1))])
v2 <- as.matrix(Surv(trainlassodat1$time,trainlassodat1$status))

###lamdaGzO_;fVF
library(glmnet)
myfit <- glmnet(v1, v2, family = "cox")
plot(myfit, xvar = "lambda", label = TRUE, cex = 1.5) 

###=;2fQiV$:s;yRrJ}A?
myfit2 <- cv.glmnet(v1, v2, family = "cox")
plot(myfit2, cex = 1.0)

###LaH!;yRr
coe <- coef(myfit, s = myfit2$lambda.min)
act_index <- which(coe != 0)
vec <- as.vector(row.names(coe)[act_index])

###W<186`RrKXcox;X9iJ}>]
muticoxdat  <- trainlassodat1[, colnames(trainlassodat1) %in% c('status', 'time', vec)] %>% 
  na.omit()
write.csv(muticoxdat, 'muti_coxdata.csv')

##Wv6`RrKXCox;X9i
fmla1 <- as.formula(Surv(time, status) ~ .)
mycox <- coxph(fmla1,data = muticoxdat)
summary(mycox)

###Vp2=;X9i
step_cox <- step(mycox ,direction = "both") 
summary(step_cox)

vec2 <- as.vector(row.names(summary(step_cox)$coefficients))

step_muticoxdat <-  trainlassodat1[, colnames(trainlassodat1) %in% c('status', 'time', vec2)] 
write.csv(step_muticoxdat, 'step_muticoxdata.csv')

##I-AVM<
library(survminer)
ggforest(step_cox, data = step_muticoxdat, fontsize = 1.5)
summary(step_cox)$coefficients

###Q5A7</<FKcN#OU5C7V2"IhVC8_5M7gOUWi
risk_score <- predict(step_cox, type = "risk", newdata = step_muticoxdat)
risk_level <- as.factor(ifelse(risk_score > median(risk_score), "High", "Low"))
mg = cbind(id = rownames(cbind(step_muticoxdat[,1:2], 
                               risk_score , risk_level)), 
           cbind(step_muticoxdat[,1:2],risk_score,risk_level))

write.table(mg, "traindata_riskscore.txt", sep = "\t" , quote = F, row.names = F)

##;fVFQ5A7</Iz4fGzO_
inputRisk <- read.table('traindata_riskscore.txt', header = T , sep = "\t" , check.names = F,row.names = 1)
kms2 <- survfit(Surv(time,status) ~ risk_level,data = inputRisk)
ggsurvplot(kms2, 
           data = inputRisk,
           pval = TRUE, 
           surv.median.line = "hv",
           legend = c(0.9,0.75),
           legend.labs = c('High Risk', 'Low Risk'))

##2bJT;z<FKcN#OU5C7V2"IhVC8_5M7gOUWi
tdrs <- select(x.test,  vec2)
tdrs1 = log2(tdrs+1)
tdrs2 <- cbind(x.test[,c(1:2)], tdrs1)
write.csv(tdrs2, 'test_muticoxdata.csv')

fmla12 <- as.formula(Surv(time,status)~.)
mycox2 <- coxph(fmla12,data = tdrs2)
risk_score2 <- predict(mycox2, type = "risk", newdata = tdrs2)
risk_level2 <- as.factor(ifelse(risk_score2 > median(risk_score2),"High","Low"))
mg2 = cbind(id = rownames(cbind(tdrs2[,1:2], 
                                risk_score2, risk_level2)), 
            cbind(tdrs2[,1:2],risk_score2,risk_level2))
write.table(mg2, "testdata_riskscore.txt",sep = "\t",quote = F,row.names = F)

##;fVF2bJT</Iz4fGzO_
inputRisk2 <- read.table('testdata_riskscore.txt', header = T , sep = "\t" , check.names = F,row.names = 1)
kms22 <- survfit(Surv(time,status) ~ risk_level2, data = inputRisk2)

##ggsurvplot WwM<
ggsurvplot(kms22, 
           data = inputRisk2,
           pval = TRUE, 
           surv.median.line = "hv",
           legend = c(0.9,0.75),
           legend.labs = c('High Risk', 'Low Risk'))


##2bJT</ROCGzO_
ROC1 <-timeROC(T = inputRisk2$time,
               delta = inputRisk2$status,
               marker = inputRisk2$risk_score, cause = 1,
               weighting = "marginal",
               times = c(1, 2, 3),
               iid = TRUE, 
               ROC = TRUE)

plot(ROC1, 
     time=1, col="red", lwd=1.5, title = "")  
plot(ROC1,
     time=2, col="blue", add=TRUE, lwd=1.5)    
plot(ROC1,
     time=3, col="orange", add=TRUE, lwd=1.5)

legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC1[["AUC"]][1],3)), 
         paste0("AUC at 2 year: ",round(ROC1[["AUC"]][2],3)), 
         paste0("AUC at 3 year: ",round(ROC1[["AUC"]][3],3))),
       col=c("red", "blue", "orange"),
       lty = 1, lwd = 1.5 ,bty = "n") 
print(confint(ROC1))

##Q5A7</ROCGzO_
ROC <-timeROC(T = inputRisk$time,
              delta = inputRisk$status,
              marker = inputRisk$risk_score, cause = 1,
              weighting = "marginal",
              times = c(1, 2, 3),
              iid = TRUE,  # <- 就是加上这一小段
              ROC = TRUE)

plot(ROC, 
     time=1, col="red", lwd=1.5, title = "")  
plot(ROC,
     time=2, col="blue", add=TRUE, lwd=1.5)    
plot(ROC,
     time=3, col="orange", add=TRUE, lwd=1.5)

legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
         paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],3)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "orange"),
       lty = 1, lwd = 2 ,bty = "n") 
print(confint(ROC))

##Q5A7</7gOURrWS
mg = arrange(mg, risk_score)
mg$patientid = c(1:nrow(mg))
ggplot(mg,aes(x = patientid,y = risk_score)) + 
  geom_point(aes(color = risk_level))+
  labs(x = "Patient (increasing risk score)", y = "Risk score")+
  geom_hline(yintercept = median(mg$risk_score),colour="black", linetype = "dotted",size = 1)+
  geom_vline(xintercept = nrow(mg)/2,colour="black", linetype = "dotted",size = 1)+
  theme_bw()


##Q5A7</Iz4fGx<d
mg$event = ifelse(mg$status==0,'Alive','Dead')
ggplot(mg,aes(x = patientid,y = time))+geom_point(aes(col = event))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(nrow(mg)/2),colour="black", linetype="dotted",size = 1)+
  theme_bw()

##2bJT</7gOURrWS
mg2 = arrange(mg2, risk_score2)
mg2$patientid = c(1:nrow(mg2))
names(mg2)[5] <- 'risk_level2'
ggplot(mg2,aes(x = patientid,y = risk_score2)) + 
  geom_point(aes(color = risk_level2))+
  labs(x = "Patient (increasing risk score)", y = "Risk score")+
  geom_vline(xintercept = nrow(mg2)/2,colour="black", linetype = "dotted",size = 1)+
  theme_bw()
 

##2bJT</Iz4fGx<d
mg2$event = ifelse(mg2$status==0,'Alive','Dead')
ggplot(mg2,aes(x = patientid,y = time))+geom_point(aes(col = event))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(nrow(mg2)/2),colour="black", linetype="dotted",size = 1)+
  theme_bw()


##Q5A7</;yRr1m4oGi?v
library(ggpubr)
merge <- step_muticoxdat %>% 
  mutate(group = inputRisk$risk_level) %>% 
  select(-c(1, 2))

library(reshape2)
merge_train <- melt(merge, id.vars = 'group')

ggplot(merge_train,aes(x = variable, 
             y = value,
             fill = group))+
  geom_boxplot(color = 'black',
               aes(fill = group),
               position = position_dodge(0.6),
               size = 0.8,
               width = 0.6)+
  theme_classic()+
  ylim(0,20)+
  ylab('Gene expression')+
  xlab('')+
  stat_compare_means(label = 'p.signif',
                     size = 5,
                     colour = 'Black')

##2bJT</;yRr1m4oGi?v
merge.1 <- tdrs2 %>% 
  mutate(group = inputRisk2$risk_level) %>% 
  select(-c(1, 2))

library(reshape2)
merge_test <- melt(merge.1, id.vars = 'group')

ggplot(merge_test,aes(x = variable, 
                       y = value,
                       fill = group))+
  geom_boxplot(color = 'black',
               aes(fill = group),
               position = position_dodge(0.6),
               size = 0.8,
               width = 0.6)+
  theme_classic()+
  ylim(0,20)+
  ylab('Gene expression')+
  xlab('')+
  stat_compare_means(label = 'p.signif',
                     size = 5,
                     colour = 'Black')


##Q5A7</84TSHHM<
library(ComplexHeatmap)
merge <- step_muticoxdat %>% 
  mutate(group = inputRisk$risk_level) %>% 
  select(-c(1, 2)) %>% 
  arrange(desc(group))
  
annol.train <- merge %>% 
  select(group)

mat.train <- merge %>% 
  select(-group) %>% 
  t()

annocol <- list(group = c('Low' = "darkblue", 'High' = "pink"))

pheatmap(mat.train,
         color = colorRampPalette(c("blue", "white", "red"))(1000),
         annotation_col = annol.train,
         annotation_colors = annocol,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F, 
         legend = T,
         scale = 'row') 


##2bJT</84TSHHM<
merge.1 <- tdrs2 %>% 
  mutate(group = inputRisk2$risk_level2) %>% 
  select(-c(1, 2)) %>% 
  arrange(desc(group))

annol.test <- merge.1 %>% 
  select(group)

mat.test <- merge.1 %>% 
  select(-group) %>% 
  t()

annocol <- list(group = c('Low' = "darkblue", 'High' = "pink"))

pheatmap(mat.test,
         color = colorRampPalette(c("blue", "white", "darkred"))(1000),
         annotation_col = annol.test,
         annotation_colors = annocol,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         legend = T,
         scale = 'row') 

## ==================== 保存训练模型 ==================== 
## 保存训练好的最终模型
saveRDS(step_cox, "final_train_model.rds")

## 保存模型系数
model_coef <- coef(step_cox)
model_info <- data.frame(
  gene = names(model_coef),
  coefficient = as.numeric(model_coef)
)
write.csv(model_info, "train_model_coefficients.csv", row.names = FALSE)

## 保存模型使用的基因列表
model_genes <- names(coef(step_cox))
write.csv(data.frame(gene = model_genes), "model_genes.csv", row.names = FALSE)

cat("训练模型已保存！\n")
cat("模型包含基因:", paste(model_genes, collapse = ", "), "\n")



##第二步
# ==================== 完整的外部验证代码 - 包含置信区间 ====================

## 设置工作目录
setwd('C:\\SYH\\Study\\ccRCC\\KIRC\\2.cox')

## 加载包
library(GEOquery)
library(tidyverse)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)

## 设置选项
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_ALL", "C")

# ==================== 第一步：检查文件 ====================

cat("检查必要文件...\n")
required_files <- c("model_genes.csv", "train_model_coefficients.csv", 
                    "GSE167573_Processed_normalized_count_matrix.txt")

for(file in required_files) {
  if(file.exists(file)) {
    cat("???", file, "存在\n")
  } else {
    cat("???", file, "缺失\n")
    stop("缺少必要文件！")
  }
}

# ==================== 第二步：读取模型信息 ====================

cat("\n读取训练集模型信息...\n")
gene_info <- read.csv('model_genes.csv', header = TRUE)
model_genes <- gene_info$gene
cat("模型包含基因:", paste(model_genes, collapse = ", "), "\n")

model_coef <- read.csv('train_model_coefficients.csv', header = TRUE)
cat("模型系数已加载\n")

# ==================== 第三步：处理临床数据 ====================

cat("\n处理验证集临床数据...\n")

gse <- getGEO('GSE167573',
              AnnotGPL = TRUE,
              getGPL = TRUE,
              GSEMatrix = TRUE,
              parseCharacteristics = FALSE)

eSet = gse[[1]] 

## 提取临床信息
phenoDat <- pData(eSet) %>% 
  select(title, geo_accession, characteristics_ch1, 'OS status:ch1', 'os time (months):ch1') %>% 
  filter(characteristics_ch1 == 'tissue: Renal cell carcinoma') %>% 
  rename(status = 'OS status:ch1') %>% 
  rename(time = 'os time (months):ch1') %>% 
  filter(status != '-') %>% 
  filter(time != '-')

cat("原始临床数据样本数:", nrow(phenoDat), "\n")

# ==================== 第四步：处理表达数据 ====================

cat("\n读取表达数据...\n")

a <- read.table('GSE167573_Processed_normalized_count_matrix.txt', 
                sep = '\t', 
                header = TRUE, 
                check.names = FALSE)

cat("表达矩阵原始维度:", dim(a), "\n")

## 处理基因名重复
a_processed <- a %>% 
  filter(!duplicated(GeneSymbol)) %>% 
  column_to_rownames(var = 'GeneSymbol')

cat("去重后表达矩阵维度:", dim(a_processed), "\n")

# ==================== 第五步：修正样本名匹配 ====================

cat("\n修正样本名匹配...\n")

# 将临床数据中的点转换为连字符
phenoDat$matched_title <- gsub("\\.", "-", phenoDat$title)

cat("临床样本名转换示例:\n")
cat("原始:", head(phenoDat$title, 3), "-> 转换:", head(phenoDat$matched_title, 3), "\n")

cat("表达数据样本名示例:", head(colnames(a_processed), 3), "\n")

## 提取模型基因表达数据
available_genes <- model_genes[model_genes %in% rownames(a_processed)]
missing_genes <- model_genes[!model_genes %in% rownames(a_processed)]

cat("可用的模型基因:", length(available_genes), "/", length(model_genes), "\n")
if(length(missing_genes) > 0) {
  cat("缺失的基因:", paste(missing_genes, collapse = ", "), "\n")
}

if(length(available_genes) == 0) {
  stop("没有模型基因在表达矩阵中找到！")
}

## 提取并转置表达数据
exprs_subset <- a_processed[available_genes, , drop = FALSE]
exprs_t <- t(exprs_subset) %>% as.data.frame()
exprs_t$sample_id <- rownames(exprs_t)

## 匹配临床数据
clinical_matched <- phenoDat %>%
  mutate(sample_id = matched_title) %>%
  inner_join(exprs_t, by = "sample_id")

cat("匹配后的样本数:", nrow(clinical_matched), "\n")

if(nrow(clinical_matched) == 0) {
  # 详细调试信息
  cat("\n=== 详细调试信息 ===\n")
  cat("临床数据样本名 (前5个):", head(phenoDat$matched_title, 5), "\n")
  cat("表达数据样本名 (前5个):", head(colnames(a_processed), 5), "\n")
  cat("共同样本:", intersect(phenoDat$matched_title, colnames(a_processed)), "\n")
  stop("样本匹配失败！")
}

# ==================== 第六步：数据预处理 ====================

cat("\n数据预处理...\n")

## 提取表达数据并进行log2转换
expression_data <- clinical_matched[, available_genes, drop = FALSE]
expression_data_log <- log2(expression_data + 1)

## 创建最终数据框
final_data <- data.frame(
  sample_id = clinical_matched$sample_id,
  time = as.numeric(clinical_matched$time),
  status = as.numeric(clinical_matched$status)
)

## 添加表达数据
final_data <- cbind(final_data, expression_data_log)
rownames(final_data) <- final_data$sample_id

## 转换为数值型
final_data$time <- as.numeric(final_data$time)
final_data$status <- as.numeric(final_data$status)

cat("最终数据维度:", dim(final_data), "\n")
cat("时间范围:", range(final_data$time, na.rm = TRUE), "个月\n")
cat("状态分布: 事件=", sum(final_data$status), ", 删失=", sum(final_data$status == 0), "\n")

# ==================== 第七步：计算风险分数 ====================

cat("\n计算风险分数...\n")

## 手动计算风险分数
risk_score <- rep(0, nrow(final_data))

for(i in 1:nrow(model_coef)) {
  gene <- model_coef$gene[i]
  coef_value <- model_coef$coefficient[i]
  if(gene %in% colnames(final_data)) {
    risk_score <- risk_score + coef_value * final_data[[gene]]
  }
}

## 计算风险等级
risk_level <- as.factor(ifelse(risk_score > median(risk_score), "High", "Low"))

## 创建结果数据框
results <- data.frame(
  id = final_data$sample_id,
  time = final_data$time / 12,  # 转换为年
  status = final_data$status,
  risk_score = risk_score,
  risk_level = risk_level
)

## 移除缺失值
results <- na.omit(results)

cat("最终验证样本数:", nrow(results), "\n")
cat("高风险组:", sum(results$risk_level == "High"), "\n")
cat("低风险组:", sum(results$risk_level == "Low"), "\n")

# ==================== 第八步：保存结果 ====================

write.table(results, "GEO_validation_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(final_data, 'phe_GEO_167573_validation.csv')

cat("\n数据预处理完成！开始生存分析...\n")

# ==================== 第九步：生存分析和置信区间 ====================

if(nrow(results) >= 10) {
  
  ## 生存分析并显示置信区间
  kms <- survfit(Surv(time, status) ~ risk_level, data = results)
  
  # 查看详细的置信区间
  cat("\n=== 生存率及其95%置信区间 ===\n")
  surv_summary <- summary(kms, times = c(1, 2, 3))  # 1年、2年、3年
  surv_ci <- data.frame(
    time = surv_summary$time,
    strata = surv_summary$strata,
    survival = surv_summary$surv,
    lower_ci = surv_summary$lower,
    upper_ci = surv_summary$upper
  )
  print(surv_ci)
  
  # 绘制带置信区间的生存曲线
  surv_plot <- ggsurvplot(kms, 
                          data = results,
                          pval = TRUE,
                          conf.int = TRUE,
                          conf.int.style = "ribbon",
                          conf.int.alpha = 0.2,
                          risk.table = TRUE,
                          surv.median.line = "hv",
                          legend = c(0.8, 0.8),
                          legend.labs = c('High Risk', 'Low Risk'),
                          title = "External Validation - Survival Analysis with CI")
  print(surv_plot)
  
  # ==================== 第十步：ROC分析和置信区间 ====================
  
  ## ROC分析
  ROC <- timeROC(T = results$time,
                 delta = results$status,
                 marker = results$risk_score,
                 cause = 1,
                 weighting = "marginal",
                 times = c(1, 2, 3),
                 iid = TRUE,
                 ROC = TRUE)
  
  # 计算AUC的置信区间
  roc_ci <- ci(ROC)
  cat("\n=== ROC AUC及其95%置信区间 ===\n")
  print(roc_ci)
  
  # 绘制ROC曲线
  plot(ROC, 
       time = 1, col = "red", lwd = 2, title = "ROC with Confidence Intervals")  
  plot(ROC,
       time = 2, col = "blue", add = TRUE, lwd = 2)    
  plot(ROC,
       time = 3, col = "green", add = TRUE, lwd = 2)
  
  # 添加置信区间到图例
  legend("bottomright",
         c(paste0("1-year AUC: ", round(ROC$AUC[1], 3), 
                  " (", round(roc_ci$CI_AUC[1, 1], 3), "-", round(roc_ci$CI_AUC[1, 2], 3), ")"),
           paste0("2-year AUC: ", round(ROC$AUC[2], 3),
                  " (", round(roc_ci$CI_AUC[2, 1], 3), "-", round(roc_ci$CI_AUC[2, 2], 3), ")"),
           paste0("3-year AUC: ", round(ROC$AUC[3], 3),
                  " (", round(roc_ci$CI_AUC[3, 1], 3), "-", round(roc_ci$CI_AUC[3, 2], 3), ")")),
         col = c("red", "blue", "green"),
         lty = 1, lwd = 2, bty = "n")
  
  # ==================== 第十一步：风险分数的置信区间 ====================
  
  ## 风险分数的描述性统计和置信区间
  cat("\n=== 风险分数的95%置信区间 ===\n")
  
  risk_stats <- results %>%
    group_by(risk_level) %>%
    summarise(
      n = n(),
      mean_risk = mean(risk_score),
      sd_risk = sd(risk_score),
      se_risk = sd_risk / sqrt(n),
      lower_ci = mean_risk - qt(0.975, n-1) * se_risk,
      upper_ci = mean_risk + qt(0.975, n-1) * se_risk
    )
  
  print(risk_stats)
  
  # 绘制带置信区间的箱线图
  library(ggplot2)
  risk_plot <- ggplot(results, aes(x = risk_level, y = risk_score, fill = risk_level)) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun.data = mean_cl_normal, geom = "pointrange", 
                 color = "red", size = 0.8) +
    labs(title = "Risk Score by Group with 95% CI",
         x = "Risk Group", y = "Risk Score") +
    theme_bw()
  print(risk_plot)
  
  # ==================== 第十二步：完整置信区间分析 ====================
  
  calculate_all_ci <- function(results) {
    cat("\n" + "="*50 + "\n")
    cat("            COMPLETE CONFIDENCE INTERVAL ANALYSIS\n")
    cat("="*50 + "\n")
    
    # 1. 生存分析CI
    kms <- survfit(Surv(time, status) ~ risk_level, data = results)
    surv_ci <- summary(kms, times = c(1,2,3))
    
    cat("\n1. SURVIVAL PROBABILITIES WITH 95% CI:\n")
    for(i in 1:length(surv_ci$time)) {
      cat(sprintf("   %s at %.0f year: %.3f (95%% CI: %.3f-%.3f)\n",
                  surv_ci$strata[i], surv_ci$time[i],
                  surv_ci$surv[i], surv_ci$lower[i], surv_ci$upper[i]))
    }
    
    # 2. ROC AUC CI
    ROC <- timeROC(T = results$time, delta = results$status,
                   marker = results$risk_score, cause = 1,
                   times = c(1,2,3), iid = TRUE)
    roc_ci <- ci(ROC)
    
    cat("\n2. ROC AUC WITH 95% CI:\n")
    for(i in 1:3) {
      cat(sprintf("   %d-year AUC: %.3f (95%% CI: %.3f-%.3f)\n",
                  i, ROC$AUC[i], roc_ci$CI_AUC[i,1], roc_ci$CI_AUC[i,2]))
    }
    
    # 3. 风险分数CI
    cat("\n3. RISK SCORE BY GROUP WITH 95% CI:\n")
    risk_summary <- results %>%
      group_by(risk_level) %>%
      summarise(
        mean = mean(risk_score),
        ci_lower = t.test(risk_score)$conf.int[1],
        ci_upper = t.test(risk_score)$conf.int[2]
      )
    print(risk_summary)
  }
  
  # 运行完整的置信区间分析
  calculate_all_ci(results)
  
  # ==================== 第十三步：输出总结 ====================
  
  cat("\n=== 外部验证结果总结 ===\n")
  cat("总样本数:", nrow(results), "\n")
  cat("事件数:", sum(results$status), "\n")
  cat("高风险组:", sum(results$risk_level == "High"), "\n")
  cat("低风险组:", sum(results$risk_level == "Low"), "\n")
  cat("1年AUC:", round(ROC$AUC[1], 3), "\n")
  cat("2年AUC:", round(ROC$AUC[2], 3), "\n")
  cat("3年AUC:", round(ROC$AUC[3], 3), "\n")
  
  # 生存差异检验
  surv_diff <- survdiff(Surv(time, status) ~ risk_level, data = results)
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  cat("Log-rank检验p值:", round(p_value, 4), "\n")
  
} else {
  cat("样本数不足，跳过生存分析\n")
}

cat("\n外部验证完成！所有置信区间分析已完成！\n")

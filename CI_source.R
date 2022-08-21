library(tidyverse)
library(ggplot2)
library(magrittr)
#install.packages('patchwork')
library('patchwork')
#install.packages("cobalt")
library(cobalt)
#install.packages('Matching')
#library(Matching)
#install.packages('rgenoud')
library(rgenoud)
#library(CBPS)
#library(gesttools)

#preprocessing df
preprocess = function(df){
  df <- df %>% mutate(X2_A = ifelse(X2 == "A", 1, 0)) %>% 
    mutate(X2_B = ifelse(X2 == "B", 1, 0)) %>% 
    mutate(X4_A = ifelse(X4 == "A", 1, 0)) %>%
    mutate(X4_B = ifelse(X4 == "B", 1, 0))
  df <- df[, (c(1, 11:15, 2, 23:24, 4, 25:26, 6:10, 16:22))]
  df["Tr"] = ifelse(df["Z"] == 1 & df["post"] == 1, 1, 0)
  y1 = df %>% filter(year == 1) %>% dplyr::select(Tr, Y, V1_avg, V2_avg, V3_avg, V4_avg, V5_A_avg, V5_B_avg, V5_C_avg)
  colnames(y1) = paste("lag", colnames(y1), sep = "")
  y2 = df %>% filter(year == 2) %>% dplyr::select(Tr, Y, V1_avg, V2_avg, V3_avg, V4_avg, V5_A_avg, V5_B_avg, V5_C_avg)
  colnames(y2) = paste("lag", colnames(y2), sep = "")
  y3 = df %>% filter(year == 3) %>% dplyr::select(Tr, Y, V1_avg, V2_avg, V3_avg, V4_avg, V5_A_avg, V5_B_avg, V5_C_avg)
  colnames(y3) = paste("lag", colnames(y3), sep = "")
  df2 = cbind(df %>% filter(year == 2), y1)
  df3 = cbind(df %>% filter(year == 3), y2)
  df4 = cbind(df %>% filter(year == 4), y3)
  df1 = cbind(df %>% filter(year == 1), matrix(0, nrow = 500, ncol = 9))
  colnames(df1) = c(colnames(df1)[1:25], colnames(y1))
  df_new = rbind(df1, df2, df3, df4) %>% arrange(id.practice)
  return(df_new)
}

#continuous to binary with median cutoff. use as df %>% summarise_all(med)
med = function(x){ifelse(x > median(x), 1, 0)}

#g_estimation
gest = function(df, S){
  df_new = preprocess(df)
  df1 = df_new %>% dplyr::select(Y, lagY, Tr, year, starts_with("X"))
  df1 = df_new %>% dplyr::select(starts_with("V"), starts_with("l")&!ends_with("Y")) %>% summarise_all(med) %>% cbind(df1)  
  df11 = df1 %>% dplyr::select(year, !starts_with("V") & !starts_with("Y"))
  #n = df_new["n.patients"] / sum(df_new["n.patients"])
  V1_glm = glm(df1$V1_avg ~., data = df11, family = "binomial")
  V2_glm = glm(df1$V2_avg ~., data = df11, family = "binomial")
  V3_glm = glm(df1$V3_avg ~., data = df11, family = "binomial")
  V4_glm = glm(df1$V4_avg ~., data = df11, family = "binomial")
  V5A_glm = glm(df1$V5_A_avg ~., data = df11, family = "binomial")
  V5B_glm = glm(df1$V5_B_avg ~., data = df11, family = "binomial")
  V5C_glm = glm(df1$V5_C_avg ~., data = df11, family = "binomial") 
  y1_glm = glm(Y ~ ., data = df1[, colnames(df1) != "year"])
  coef_v = cbind(V1_glm$coefficients, V2_glm$coefficients, V3_glm$coefficients, V4_glm$coefficients, V5A_glm$coefficients, V5B_glm$coefficients, V5C_glm$coefficients)
  Y3 = Y4 = numeric(500)
  for (s in 1:S){
    #year1
    df_1 = cbind(df["id.practice"], df11) %>% filter(year == 1)
    df_1["id.practice"] = rep(1, 500);colnames(df_1) = c("int", colnames(df_1)[2:23])
    z = df_1 %>% as.matrix() %*% coef_v
    z = exp(z) / (1 + exp(z))
    for (j in 1:ncol(z)) for (i in 1:nrow(z)) z[i, j] = rbernoulli(1, p = z[i, j])
    V_1 = z
    y_1 = cbind(df_1["int"], V_1, df_1["lagTr"], df_1[, 4:23]) %>% as.matrix() %*% y1_glm$coefficients
    #year2
    df_1 = cbind(df["id.practice"], df11) %>% filter(year == 2)
    df_1["id.practice"] = rep(1, 500); colnames(df_1) = c("int", colnames(df_1)[2:23])
    df_1[, 4:10] = V_1
    df_1["lagY"] = y_1
    V_2 = df_1 %>% as.matrix() %*% coef_v
    V_2 = exp(V_2) / (1 + exp(V_2))
    for (j in 1:ncol(V_2)) for (i in 1:nrow(V_2)) V_2[i, j] = rbernoulli(1, p = V_2[i, j])
    y_2 = cbind(df_1["int"], V_2, df_1["lagTr"], df_1[, 4:23]) %>% as.matrix() %*% y1_glm$coefficients
    #year3
    df_1 = cbind(df["id.practice"], df11) %>% filter(year == 3)
    df_1["Tr"] = 1 - df_1["Tr"]
    df_1["id.practice"] = rep(1, 500);colnames(df_1) = c("int", colnames(df_1)[2:23])
    df_1[, 4:10] = V_2
    df_1["lagY"] = y_2
    z = df_1 %>% as.matrix() %*% coef_v
    z = exp(z) / (1 + exp(z))
    for (j in 1:ncol(z)) for (i in 1:nrow(z)) z[i, j] = rbernoulli(1, p = z[i, j])
    V_3 = z
    y_3 = cbind(df_1["int"], V_3, df_1["lagTr"], df_1[, 4:23]) %>% as.matrix() %*% y1_glm$coefficients
    Y3 = Y3 + y_3
    #year4
    df_1 = cbind(df["id.practice"], df11) %>% filter(year == 4)
    df_1["Tr"] = 1 - df_1["Tr"]
    df_1["lagTr"] = 1 - df_1["lagTr"]
    df_1["id.practice"] = rep(1, 500);colnames(df_1) = c("int", colnames(df_1)[2:23])
    df_1[, 4:10] = V_3
    df_1["lagY"] = y_3
    z = df_1 %>% as.matrix() %*% coef_v
    z = exp(z) / (1 + exp(z))
    for (j in 1:ncol(z)) for (i in 1:nrow(z)) z[i, j] = rbernoulli(1, p = z[i, j])
    V_4 = z
    y_4 = cbind(df_1["int"], V_4, df_1["lagTr"], df_1[, 4:23]) %>% as.matrix() %*% y1_glm$coefficients
    Y4 = Y4 + y_4
  }
  Y3 = Y3 / S; Y4 = Y4 / S
  diff3 = df_new %>% filter(year == 3) %>% dplyr::select(n.patients, Z, Y) %>% cbind(Y3) %>% 
    filter(Z == 1) %>% mutate(diff = (Y - Y3) * n.patients) %>% 
    dplyr::select(n.patients, diff) %>% apply(2, sum)
  diff4 = df_new %>% filter(year == 4) %>% dplyr::select(n.patients, Z, Y) %>% cbind(Y4) %>%
    filter(Z == 1) %>% mutate(diff = (Y - Y4) * n.patients) %>%
    dplyr::select(n.patients, diff) %>% apply(2, sum)
  SATT = (diff3[2] + diff4[2]) / (diff3[1] + diff4[1])
  SATT_t = c(diff3[2] / diff3[1], diff4[2] / diff4[1])
  out = list(SATT, SATT_t)
  return(out)
}

ps_g = function(df){
  df["ps"] = rep(0, dim(df)[1])
  for (i in 1:4){
    X = df %>% filter(year == i) %>% dplyr::select(Z, starts_with("X"), starts_with("V"))
    glm1 = glm(Z ~ ., family = "binomial", data = X)
    df[df["year"] == i, "ps"] = glm1$fitted.values
  }
  return(df)
}

satt = function(df, t, g, yname = "Y", gname = "Z", psname = "ps"){
  Z = df[df["year"] == g, "Z"]
  ps_g = df[df["year"] == g, "ps"]
  Yt = df[df["year"] == t, "Y"]
  Yg_1 = df[df["year"] == g-1, "Y"]
  ds = ps_g * (1 - Z) / (1 - ps_g)
  SATT = (Z / mean(Z) - ds / mean(ds)) * (Yt - Yg_1)
  return(SATT)
}

att = function(df, t = 3:4, g = 3, M = 10000){
  ATT = numeric(M)
  ATT_t = matrix(0, nrow = M, ncol = length(t)) %>% data.frame()
  colnames(ATT_t) = paste("ATT(", g, ", ", t, ")", sep = "")
  for (i in 1:M){
    df_b = data.frame()
    for (ii in 1:4){
      df_b = df_b %>% rbind(df[df["year"] == ii & df["Z"] == 0, ][sample(1:303, 303, replace = T),])
      df_b = df_b %>% rbind(df[df["year"] == ii & df["Z"] == 1, ][sample(1:197, 197, replace = T),])
    }
    df_b = ps_g(df_b)
    N = 0
    for (j in 1:length(t)){
      ATT[i] = ATT[i] + sum(satt(df_b, t[j], g) * df_b[df_b["year"] == t[j], "n.patients"])
      N = N + sum(df_b[df_b["year"] == t[j], "n.patients"])
      ATT_t[i, j] = sum(satt(df_b, t[j], g) * df_b[df_b["year"] == t[j], "n.patients"]) / sum(df_b[df_b["year"] == t[j], "n.patients"])
    }
    ATT[i] = ATT[i] / N
    #print(ATT[i,])
  }
  CI_L = function(x, alpha = 0.05) sort(x)[round(length(x)*alpha / 2)]
  CI_U = function(x, alpha = 0.05) sort(x)[round(length(x)*(1 - alpha / 2))]
  
  out = rbind(apply(ATT_t, 2, mean), apply(ATT_t, 2, sd),  apply(ATT_t, 2, CI_L), apply(ATT_t, 2, CI_U))
  rownames(out) = c("ATT", "SE", "LB 95% CI", "UB 95% CI")
  r = list(ATT = data.frame(method = "IPW w DID", att = mean(ATT), sd = sd(ATT)), ATT_t = out)
  return(r)
}

did = function(df, t = 3:4, g = 3, M = 10000, yname = "Y", gname = "Z"){
  ATT = numeric(M)
  ATT_t = matrix(0, nrow = M, ncol = length(t)) %>% data.frame()
  colnames(ATT_t) = paste("ATT(", g, ", ", t, ")", sep = "")
  for (i in 1:M){
    df_b = data.frame()
    for (ii in 1:4){
      df_b = df_b %>% rbind(df[df["year"] == ii & df["Z"] == 0, ][sample(1:303, 303, replace = T),])
      df_b = df_b %>% rbind(df[df["year"] == ii & df["Z"] == 1, ][sample(1:197, 197, replace = T),])
    }
    Yg_1 = df_b[df_b["year"] == g-1 & df_b["Z"] == 1, "Y"]
    Yg_0 = df_b[df_b["year"] == g-1 & df_b["Z"] == 0, "Y"]
    for (j in 1:length(t)){
      Yt_1 = df_b[df_b["year"] == t[j] & df_b["Z"] == 1, "Y"]
      N_1 = df_b[df_b["year"] == t[j] & df_b["Z"] == 1, "n.patients"]
      Yt_0 = df_b[df_b["year"] == t[j] & df_b["Z"] == 0, "Y"]
      N_0 = df_b[df_b["year"] == t[j] & df_b["Z"] == 0, "n.patients"]
      ATT_t[i, j] = sum((Yt_1 - Yg_1) * N_1) / sum(N_1) - sum((Yt_0 - Yg_0) * N_0) / sum(N_0)
    }
    ATT[i] = sum(ATT_t[i, ] * df_b[df_b["year"] == t,] %>% group_by(year) %>% summarise(n = sum(n.patients)) %>% dplyr::select(n) %>% as.matrix()) / sum(df_b[df_b["year"] == t, "n.patients"])
  }
  CI_L = function(x, alpha = 0.05) sort(x)[round(length(x)*alpha / 2)]
  CI_U = function(x, alpha = 0.05) sort(x)[round(length(x)*(1 - alpha / 2))]
  
  out = rbind(apply(ATT_t, 2, mean), apply(ATT_t, 2, sd),  apply(ATT_t, 2, CI_L), apply(ATT_t, 2, CI_U))
  rownames(out) = c("ATT", "SE", "LB 95% CI", "UB 95% CI")
  r = list(ATT = data.frame(method = "DID", att = mean(ATT), sd = sd(ATT)), ATT_t = out)
  return(r)
}
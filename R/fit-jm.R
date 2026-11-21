# 결합모형(JMbayes2)을 적합하고 작은 요약만 outputs/jm-01.rds에 저장한다.
# 적합에 수십 분이 걸리고 적합 객체가 커서, 분석 문서 안에서 돌리지 않는다.
#
#   Rscript R/fit-jm.R              기본 반복 수
#   Rscript R/fit-jm.R 20000 4000   반복 수, 버림 수를 바꿔서

suppressPackageStartupMessages({
  library(data.table); library(survival); library(nlme); library(JMbayes2)
})
args <- commandArgs(trailingOnly = TRUE)
# 기본값은 모의 데이터에서 수렴한 반복 수다(6000회에서는 value·slope의 Rhat이 1.1을 넘었다)
N_ITER   <- if (length(args) >= 1) as.integer(args[1]) else 15000L
N_BURNIN <- if (length(args) >= 2) as.integer(args[2]) else 3000L

ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(ROOT, "R", "oai.R"))

src <- pick_data(ROOT)
d <- build_analysis(src$dir)
j <- jm_data(d$knees, d$long)

t0 <- Sys.time()
# 종단: KL 등급마다 수준과 기울기가 다르고, 무릎마다 무작위 절편·기울기
lme_fit <- lme(womac_func ~ year * kl, random = ~ year | knee, data = j$long,
               control = lmeControl(opt = "optim"))
# 생존: 기준 공변량 + 그 시점의 참 기능 점수(value)와 변화율(slope)
cox_fit <- coxph(Surv(time_y, event) ~ age_c + female + bmi_c + kl, data = j$surv)
jm_fit <- jm(cox_fit, lme_fit, time_var = "year",
             functional_forms = ~ value(womac_func) + slope(womac_func),
             n_iter = N_ITER, n_burnin = N_BURNIN, n_chains = 3L, cores = 3L, seed = 20261118L)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

s <- summary(jm_fit)
tidy <- function(tab, expo) {
  f <- if (expo) exp else identity
  data.table(term = rownames(tab), est = f(tab$Mean), lo = f(tab$`2.5%`), hi = f(tab$`97.5%`), Rhat = tab$Rhat)
}
res <- list(
  survival     = tidy(s$Survival, expo = TRUE),
  longitudinal = tidy(s$Outcome1, expo = FALSE),
  n_knees = nrow(j$surv), n_events = sum(j$surv$event), n_obs = nrow(j$long),
  n_iter = N_ITER, n_burnin = N_BURNIN, minutes = minutes,
  real_data = src$real, fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M"))

dir.create(file.path(ROOT, "outputs"), showWarnings = FALSE)
saveRDS(res, file.path(ROOT, "outputs", "jm-01.rds"))
print(res$survival, digits = 3)
cat(sprintf("%.1f분, 최대 Rhat %.3f\n", minutes, max(c(res$survival$Rhat, res$longitudinal$Rhat), na.rm = TRUE)))

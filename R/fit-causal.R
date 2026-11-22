# 근력운동과 걷기의 시작(과 중단)을 처치로 한 모사 시험의 평균 효과를 추정하고
# 요약만 outputs/causal-07.rds에 저장한다(docs/causal-inference-plan.md 6장. 장에서는 쓰지 않는 쌓은 시험의 추정이다).
#
#   Rscript R/fit-causal.R          나무 2,000그루(10분 안팎)
#   Rscript R/fit-causal.R 500      나무 수를 줄여 빨리 돌린다
#
# 추정량
#   raw   보정하지 않은 평균 차이
#   ols   공변량을 넣은 선형 회귀(결측은 중앙값으로 채우고 결측 표시 변수를 더한다)
#   ate   인과 숲의 성향점수와 결과 회귀로 만든 이중 강건(AIPW) 평균 처치 효과
#   att   같은 방법으로 구한 처치군에서의 효과
#   ato   겹침 가중 효과(성향점수가 0이나 1에 가까운 관측의 비중을 줄인다)
# 표준오차는 모두 사람 단위 군집이다. 결과 방문의 탈락은 관측 확률의 역수로 가중한다.

suppressPackageStartupMessages({ library(data.table); library(grf) })
args <- commandArgs(trailingOnly = TRUE)
TREES <- if (length(args) >= 1) as.integer(args[1]) else 2000L
SEED <- 20261122L

ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(ROOT, "R", "oai.R"))

src <- pick_data(ROOT)
panel <- tryCatch(exercise_panel(src$dir), error = function(e) {
  if (src$real) stop(e)
  make_synthetic(src$dir)                       # PASE 항목이 없는 옛 시뮬레이션 데이터
  exercise_panel(src$dir)
})

FIXED <- c("age_i", "female", "race", "cohort", "site", "edu", "income", "kl", "align_r", "align_l",
           "injury", "surgery", "trial", "gap")
# 지표 방문의 PASE 총점은 처치 항목을 포함하므로 넣지 않는다
covariate_sets <- function(treatment) list(
  prev = c(FIXED, paste0(TRIAL_COVARIATES, "_p"), "pain_pp"),
  both = c(FIXED, paste0(TRIAL_COVARIATES, "_p"), "pain_pp",
           paste0(setdiff(TRIAL_COVARIATES, c(treatment, "pase")), "_i")))

OUTCOMES <- c(pain = "y_pain", func = "y_func", pace = "y_pace", chair = "y_chair", hand = "y_hand")

# 값이 하나뿐이거나 모두 결측인 열을 뺀 숫자 행렬
design <- function(d, vars) {
  X <- as.matrix(d[, vars, with = FALSE]); storage.mode(X) <- "double"
  X[, apply(X, 2, function(x) { u <- unique(x[!is.na(x)]); length(u) > 1 }), drop = FALSE]
}

# 군집 강건 표준오차를 붙인 가중 최소제곱의 처치 계수
cluster_wls <- function(y, a, Z, id, w = rep(1, length(y))) {
  X <- cbind(1, a, Z)
  X <- X[, c(TRUE, TRUE, if (is.null(Z)) NULL else apply(Z, 2, sd) > 0), drop = FALSE]
  f <- lm.wfit(X, y, w)
  keep <- !is.na(f$coefficients); X <- X[, keep, drop = FALSE]
  B <- solve(crossprod(X * sqrt(w)))
  S <- rowsum(X * (w * f$residuals), id)
  V <- B %*% crossprod(S) %*% B * (nrow(S) / (nrow(S) - 1))
  c(est = unname(f$coefficients[2]), se = sqrt(V[2, 2]))
}

impute <- function(X) {
  miss <- is.na(X)
  ind <- miss[, colMeans(miss) > 0.01, drop = FALSE]
  if (ncol(ind)) colnames(ind) <- paste0(colnames(ind), "_na")
  for (j in seq_len(ncol(X))) X[miss[, j], j] <- median(X[, j], na.rm = TRUE)
  cbind(X, ind * 1)
}

smd <- function(x, a, w = rep(1, length(x))) {
  ok <- !is.na(x)
  m <- function(g) weighted.mean(x[ok & a == g], w[ok & a == g])
  s <- sqrt((var(x[ok & a == 1]) + var(x[ok & a == 0])) / 2)
  if (is.na(s) || s == 0) NA_real_ else (m(1) - m(0)) / s
}

# 한 설계(처치, 방향, 공변량 집합, 결과)의 추정. y가 관측된 관측만 쓰고 탈락은 역확률로 가중한다.
fit_one <- function(d, vars, ycol, censor = TRUE) {
  X <- design(d, vars)
  y <- d[[ycol]]; a <- d$a; id <- d$ID
  obs <- !is.na(y)
  w <- rep(1, nrow(d))
  if (censor && any(!obs)) {
    cf <- regression_forest(cbind(X, a = a), as.numeric(obs), clusters = id, num.trees = 500, seed = SEED)
    w <- 1 / pmax(predict(cf)$predictions, 0.05)
  }
  X <- X[obs, , drop = FALSE]; y <- y[obs]; a <- a[obs]; id <- id[obs]; w <- w[obs]
  f <- causal_forest(X, y, a, clusters = id, sample.weights = w, num.trees = TREES, seed = SEED)
  e <- f$W.hat
  est <- rbind(
    raw = cluster_wls(y, a, NULL, id),
    ols = cluster_wls(y, a, impute(X), id, w),
    ate = average_treatment_effect(f, target.sample = "all"),
    att = average_treatment_effect(f, target.sample = "treated"),
    ato = average_treatment_effect(f, target.sample = "overlap"))
  est <- data.table(estimator = rownames(est), est = est[, 1], se = est[, 2])
  ipw <- w * ifelse(a == 1, 1 / e, 1 / (1 - e))
  bal <- data.table(variable = colnames(X),
                    smd_raw = apply(X, 2, smd, a = a),
                    smd_ipw = apply(X, 2, smd, a = a, w = ipw))
  overlap <- data.table(arm = c("처치", "비교"),
                        rbind(quantile(e[a == 1], c(0, .01, .05, .5, .95, .99, 1)),
                              quantile(e[a == 0], c(0, .01, .05, .5, .95, .99, 1))),
                        outside = c(mean(e[a == 1] < 0.05 | e[a == 1] > 0.95), mean(e[a == 0] < 0.05 | e[a == 0] > 0.95)))
  list(est = est, bal = bal, overlap = overlap,
       info = data.table(n = length(y), treated = sum(a), persons = uniqueN(id), dropped = sum(!obs),
                         y_mean0 = weighted.mean(y[a == 0], w[a == 0]), y_sd = sd(y), max_w = max(w)))
}

t0 <- Sys.time()
est <- list(); bal <- list(); overlap <- list(); info <- list(); counts <- list()
for (treatment in c("strength", "walk")) for (direction in c("start", "stop")) {
  d <- trial_table(src$dir, treatment, direction, panel = panel)
  key0 <- list(treatment = treatment, direction = direction)
  counts[[length(counts) + 1]] <- cbind(as.data.table(key0),
    d[, .(n = .N, treated = sum(a), tkr_out = sum(tkr_out), y_obs = sum(!is.na(y_pain))), by = trial])
  sets <- covariate_sets(treatment)
  for (cs in names(sets)) {
    jobs <- lapply(names(OUTCOMES), function(o) list(kind = if (o == "hand") "negative-control" else "outcome",
                                                     outcome = o, ycol = OUTCOMES[[o]], vars = sets[[cs]], data = d))
    if (cs == "prev" && direction == "start") {
      # 위약: 직전 방문의 공변량이 같을 때, 시작한 사람과 하지 않은 사람의 그 전 방문 통증이 다른가(과거의 추세)
      jobs[[length(jobs) + 1]] <- list(kind = "placebo", outcome = "pain_pp", ycol = "pain_pp",
                                       vars = setdiff(sets$prev, "pain_pp"), data = d[!is.na(pain_pp)], censor = FALSE)
      # 동시점: 지표 방문의 통증 차이. 시작 전의 교란과 시작 뒤의 효과가 섞인 값이다
      jobs[[length(jobs) + 1]] <- list(kind = "concurrent", outcome = "pain_i", ycol = "pain_i",
                                       vars = sets$prev, data = d, censor = FALSE)
    }
    for (j in jobs) {
      if (all(is.na(j$data[[j$ycol]])) || uniqueN(j$data[[j$ycol]], na.rm = TRUE) < 2) next
      r <- fit_one(j$data, j$vars, j$ycol, censor = !identical(j$censor, FALSE))
      key <- c(key0, list(covariates = cs, kind = j$kind, outcome = j$outcome))
      est[[length(est) + 1]] <- cbind(as.data.table(key), r$est)
      bal[[length(bal) + 1]] <- cbind(as.data.table(key), r$bal)
      overlap[[length(overlap) + 1]] <- cbind(as.data.table(key), r$overlap)
      info[[length(info) + 1]] <- cbind(as.data.table(key), r$info)
      cat(sprintf("%s %s %s %s: ATE %.3f (%.3f), n %d, %.1f분\n", treatment, direction, cs, j$outcome,
                  r$est[estimator == "ate", est], r$est[estimator == "ate", se], r$info$n,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
}

out <- list(est = rbindlist(est), balance = rbindlist(bal), overlap = rbindlist(overlap), info = rbindlist(info),
            counts = rbindlist(counts), real = src$real, trees = TREES, seed = SEED,
            minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), date = Sys.Date())
out$est[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
dir.create(file.path(ROOT, "outputs"), showWarnings = FALSE)
saveRDS(out, file.path(ROOT, "outputs", "causal-07.rds"))
cat(sprintf("저장: outputs/causal-07.rds (%.1f분, %s 데이터)\n", out$minutes, if (src$real) "OAI" else "시뮬레이션"))

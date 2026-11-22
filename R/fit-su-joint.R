# SU정규 결합모형 장의 다변량 정규 결합모형을 적합하고 요약만 outputs/su-joint-06.rds 에 저장한다.
# 종단 하위모형의 무작위 효과와 변환한 치환술 시간이 함께 다변량 정규분포를 따른다.
# 치환술 시간의 주변분포는 로그정규 또는 SU정규(Johnson S_U) 가속수명 모형이다.
#
#   Rscript R/fit-su-joint.R        전체: 적합, 분포의 비교, 동적 예측, 모의 검증 20회(약 25분)
#   Rscript R/fit-su-joint.R 0      모의 검증 없이(약 5분)
#
# 표본은 결합모형과 동적 예측 장 I장의 결합모형과 같은 정의(jm_data)다. 동적 예측은 그 장의 JMbayes2와 같은 무릎으로 적합하고
# 같은 무릎에서 평가한다(outputs/jm-01.rds 가 있을 때).

suppressPackageStartupMessages({library(data.table); library(survival)})
args <- commandArgs(trailingOnly = TRUE)
N_REP <- if (length(args) >= 1) as.integer(args[1]) else 20L
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }

# ---------------------------------------------------------------- 모형
# long: knee, year, womac_func 와 고정효과 설계행렬 X.  surv: knee, time_y, event 와 설계행렬 W(상수항 포함)
prep <- function(surv, long, form_l = ~ year * kl, form_s = ~ age_c + female + bmi_c + kl) {
  surv <- as.data.frame(surv); long <- as.data.frame(long)
  g <- match(long$knee, surv$knee)
  X <- model.matrix(form_l, long); W <- model.matrix(form_s, surv)
  n_i <- tabulate(g, nrow(surv))
  list(X = X, W = W, y = long$womac_func, t = long$year, g = g, n = n_i,
       st = as.vector(rowsum(long$year, factor(g, levels = seq_len(nrow(surv))))),
       stt = as.vector(rowsum(long$year^2, factor(g, levels = seq_len(nrow(surv))))),
       logT = log(surv$time_y), d = surv$event, N = nrow(surv), pl = ncol(X), ps = ncol(W),
       xn = colnames(X), wn = colnames(W), knee = as.character(surv$knee))
}

# 파라미터 벡터를 이름 붙은 목록으로. su = TRUE 이면 형상 파라미터 lambda, log theta 가 끝에 붙는다
unpack <- function(p, D, su) {
  i <- 0; take <- function(k) { v <- p[i + seq_len(k)]; i <<- i + k; v }
  o <- list(beta = take(D$pl), ls = take(1), ch = take(3), gamma = take(D$ps), lt = take(1), a = take(2))
  if (su) { o$lam <- take(1); o$lth <- take(1) }
  L11 <- exp(o$ch[1]); L21 <- o$ch[2]; L22 <- exp(o$ch[3])
  o$d11 <- L11^2; o$d12 <- L11 * L21; o$d22 <- L21^2 + L22^2
  o
}

# 무릎별 조건부 분포의 재료: 무작위 효과의 사후 평균(mb1, mb2), M 행렬, 종단 로그우도
long_parts <- function(o, D, rows = NULL) {
  r <- as.vector(D$y - D$X %*% o$beta)
  if (!is.null(rows)) {                       # 일부 측정만 쓴다(동적 예측)
    S <- rowsum(cbind(1, D$t, D$t^2, r, D$t * r, r^2)[rows, , drop = FALSE], factor(D$g[rows], levels = seq_len(D$N)))
    n <- S[, 1]; st <- S[, 2]; stt <- S[, 3]; sr <- S[, 4]; str <- S[, 5]; srr <- S[, 6]
  } else {
    S <- rowsum(cbind(r, D$t * r, r^2), factor(D$g, levels = seq_len(D$N)))
    n <- D$n; st <- D$st; stt <- D$stt; sr <- S[, 1]; str <- S[, 2]; srr <- S[, 3]
  }
  s2 <- exp(2 * o$ls); detD <- o$d11 * o$d22 - o$d12^2
  M11 <- s2 * o$d22 / detD + n; M12 <- -s2 * o$d12 / detD + st; M22 <- s2 * o$d11 / detD + stt
  detM <- M11 * M22 - M12^2
  mb1 <- (M22 * sr - M12 * str) / detM; mb2 <- (-M12 * sr + M11 * str) / detM
  quad <- (srr - (sr * mb1 + str * mb2)) / s2
  ll <- -0.5 * (n * log(2 * pi) + (n - 2) * log(s2) + log(detM) + log(detD) + quad)
  list(ll = ll, mb1 = mb1, mb2 = mb2, M11 = M11, M12 = M12, M22 = M22, detM = detM, s2 = s2)
}

# 표준화한 오차 g(Z): 로그정규이면 Z, SU정규이면 (sinh(lambda + theta Z) - sinh lambda) / (theta cosh lambda)
g_of <- function(o, z, su) if (su) (sinh(o$lam + exp(o$lth) * z) - sinh(o$lam)) / (exp(o$lth) * cosh(o$lam)) else z

# 변환한 사건 시간 z = g^{-1}((ln t - eta) / tau) 와 야코비안의 로그
z_of <- function(o, logT, eta, su) {
  tau <- exp(o$lt); e <- (logT - eta) / tau
  if (!su) return(list(z = e, lj = -o$lt - logT))
  th <- exp(o$lth); eo <- sinh(o$lam) + th * cosh(o$lam) * e
  list(z = (asinh(eo) - o$lam) / th, lj = log(cosh(o$lam)) - 0.5 * log1p(eo^2) - o$lt - logT)
}

# 측정이 주어졌을 때 Z 의 조건부 평균과 분산
z_cond <- function(o, L) {
  a1 <- o$a[1]; a2 <- o$a[2]
  k2 <- 1 + a1^2 * o$d11 + 2 * a1 * a2 * o$d12 + a2^2 * o$d22
  aSa <- L$s2 * (a1^2 * L$M22 - 2 * a1 * a2 * L$M12 + a2^2 * L$M11) / L$detM
  list(mu = (a1 * L$mb1 + a2 * L$mb2) / sqrt(k2), v = (aSa + 1) / k2, k2 = k2)
}

negll <- function(p, D, su) {
  o <- unpack(p, D, su)
  if (abs(o$lt) > 8 || (su && (abs(o$lam) > 8 || abs(o$lth) > 8))) return(1e12)
  L <- long_parts(o, D); zc <- z_cond(o, L)
  zz <- z_of(o, D$logT, as.vector(D$W %*% o$gamma), su)
  q <- (zz$z - zc$mu) / sqrt(zc$v)
  ls <- ifelse(D$d == 1, dnorm(q, log = TRUE) - 0.5 * log(zc$v) + zz$lj, pnorm(q, lower.tail = FALSE, log.p = TRUE))
  v <- -sum(L$ll) - sum(ls)
  if (is.finite(v)) v else 1e12
}

start_values <- function(D, su) {
  b <- qr.solve(D$X, D$y); r <- D$y - D$X %*% b
  mk <- as.vector(rowsum(r, D$g)) / D$n
  ls <- 0.5 * log(mean((r - mk[D$g])^2))
  f <- survreg(Surv(exp(D$logT), D$d) ~ D$W - 1, dist = "lognormal")
  c(b, ls, log(sd(mk)), 0, log(1), coef(f), log(f$scale), 0, 0, if (su) c(0, log(0.2)))
}

par_names <- function(D, su) c(paste0("long:", D$xn), "log_sigma", "chol11", "chol21", "chol22", paste0("surv:", D$wn),
                                "log_tau", "a0", "a1", if (su) c("lambda", "log_theta"))

fit <- function(D, su, start = NULL, hessian = TRUE) {
  p0 <- if (is.null(start)) start_values(D, su) else start
  t0 <- Sys.time()
  r <- optim(p0, negll, D = D, su = su, method = "BFGS", control = list(maxit = 2000, reltol = 1e-11))
  r <- optim(r$par, negll, D = D, su = su, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  out <- list(par = r$par, ll = -r$value, conv = r$convergence, su = su, k = length(r$par))
  if (hessian) {
    H <- optimHess(r$par, negll, D = D, su = su)
    out$pd <- all(eigen((H + t(H)) / 2, only.values = TRUE)$values > 1e-8)
    out$vcov <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
    out$se <- sqrt(abs(diag(out$vcov)))
  }
  out$sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  nm <- par_names(D, su); names(out$par) <- nm
  if (hessian) { names(out$se) <- nm; dimnames(out$vcov) <- list(nm, nm) }
  out
}

# 동적 예측: 시각 s 까지의 측정이 주어졌을 때 (s, s + h] 안에 사건이 일어날 확률
predict_risk <- function(f, D, s, h) {
  o <- unpack(f$par, D, f$su)
  L <- long_parts(o, D, rows = which(D$t <= s)); zc <- z_cond(o, L)
  eta <- as.vector(D$W %*% o$gamma)
  q <- function(tt) (z_of(o, rep(log(tt), D$N), eta, f$su)$z - zc$mu) / sqrt(zc$v)
  1 - exp(pnorm(q(s + h), lower.tail = FALSE, log.p = TRUE) - pnorm(q(s), lower.tail = FALSE, log.p = TRUE))
}

# 주변 생존함수: Z 의 주변분포가 N(0, 1) 이므로 S(t | w) = Phi(-z(t))
marg_surv <- function(f, D, times, rows = seq_len(D$N)) {
  o <- unpack(f$par, D, f$su); eta <- as.vector(D$W[rows, , drop = FALSE] %*% o$gamma)
  sapply(times, function(tt) mean(pnorm(-z_of(o, rep(log(tt), length(rows)), eta, f$su)$z)))
}

# 무작위 효과가 delta 만큼 다를 때 조건부 중위값의 비
time_ratio <- function(f, D, delta) {
  o <- unpack(f$par, D, f$su)
  k <- sqrt(z_cond(o, list(s2 = 1, M11 = 1, M12 = 0, M22 = 1, detM = 1, mb1 = 0, mb2 = 0))$k2)
  exp(exp(o$lt) * (g_of(o, sum(o$a * delta) / k, f$su) - g_of(o, 0, f$su)))
}

# 공변량 w 에서 치환술 시간의 주변 분위수: T_p = exp(w'gamma + tau g(Phi^{-1}(p)))
marg_quantile <- function(f, D, w, p) {
  o <- unpack(f$par, D, f$su); exp(sum(w * o$gamma) + exp(o$lt) * g_of(o, qnorm(p), f$su))
}

# 모형에서 데이터를 생성한다. 방문 시각과 공변량은 D 의 것을 쓰고, 중도절단 시각은 cens 에서 받는다
simulate <- function(p, D, su, cens, seed = 1) {
  set.seed(seed); o <- unpack(p, D, su)
  L11 <- exp(o$ch[1]); L21 <- o$ch[2]; L22 <- exp(o$ch[3])
  u1 <- rnorm(D$N); u2 <- rnorm(D$N); b0 <- L11 * u1; b1 <- L21 * u1 + L22 * u2
  k <- sqrt(1 + o$a[1]^2 * o$d11 + 2 * o$a[1] * o$a[2] * o$d12 + o$a[2]^2 * o$d22)
  Z <- (o$a[1] * b0 + o$a[2] * b1 + rnorm(D$N)) / k
  T <- exp(as.vector(D$W %*% o$gamma) + exp(o$lt) * g_of(o, Z, su))
  y <- as.vector(D$X %*% o$beta) + b0[D$g] + b1[D$g] * D$t + rnorm(length(D$t), sd = exp(o$ls))
  time <- pmin(T, cens); keep <- D$t < time[D$g] | D$t == 0
  list(time = time, event = as.integer(T <= cens), y = y, keep = keep)
}

subset_jm <- function(j, keep_knees) {
  s <- j$surv[as.character(j$surv$knee) %in% keep_knees, ]; s$knee <- droplevels(s$knee)
  l <- j$long[as.character(j$long$knee) %in% keep_knees, ]; l$knee <- factor(as.character(l$knee), levels = levels(s$knee))
  list(surv = s, long = l)
}

# ---------------------------------------------------------------- 실행
if (sys.nframe() == 0) {
  ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
  source(file.path(ROOT, "R", "oai.R"))
  src <- pick_data(ROOT); d <- build_analysis(src$dir); j <- jm_data(d$knees, d$long)
  D <- prep(j$surv, j$long)
  res <- list(real_data = src$real, fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M"),
              n_knees = D$N, n_events = sum(D$d), n_obs = length(D$y))
  tab <- function(f) data.table(term = names(f$par), est = unname(f$par), se = unname(f$se))
  say("표본: 무릎", D$N, "치환술", sum(D$d), "관측", length(D$y))

  # --- 1. 무릎 전체로 적합: 독립(a = 0), 결합 로그정규, 결합 SU정규
  f_ln <- fit(D, su = FALSE)
  f_su <- fit(D, su = TRUE, start = c(f_ln$par, 0, log(0.2)))
  s0 <- f_ln$par; s0[c("a0", "a1")] <- 0; free <- !names(s0) %in% c("a0", "a1")
  f_ind <- optim(s0[free], function(q) { p <- s0; p[free] <- q; negll(p, D, FALSE) },
                 method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  say(sprintf("적합 끝: 로그정규 %.0f초, SU정규 %.0f초", f_ln$sec, f_su$sec))
  res$fits <- list(
    independent = list(ll = -f_ind$value, k = sum(free)),
    lognormal = c(f_ln[c("ll", "k", "sec", "conv", "pd")], list(coef = tab(f_ln))),
    su = c(f_su[c("ll", "k", "sec", "conv", "pd")], list(coef = tab(f_su))))
  res$time_ratio <- sapply(list(lognormal = f_ln, su = f_su), function(f)
    c(intercept10 = time_ratio(f, D, c(10, 0)), slope1 = time_ratio(f, D, c(0, 1))))
  # Z 와 무작위 효과의 상관, 무작위 효과가 설명하는 Z 분산의 몫
  res$assoc <- rbindlist(lapply(list(lognormal = f_ln, su = f_su), function(f) {
    o <- unpack(f$par, D, f$su); k2 <- z_cond(o, list(s2 = 1, M11 = 1, M12 = 0, M22 = 1, detM = 1, mb1 = 0, mb2 = 0))$k2
    data.table(sd_b0 = sqrt(o$d11), sd_b1 = sqrt(o$d22), cor_b = o$d12 / sqrt(o$d11 * o$d22),
               cor_z_b0 = (o$a[1] * o$d11 + o$a[2] * o$d12) / sqrt(k2 * o$d11),
               cor_z_b1 = (o$a[1] * o$d12 + o$a[2] * o$d22) / sqrt(k2 * o$d22), share = 1 - 1 / k2)
  }), idcol = "model")

  # --- 2. 주변분포: KL 등급 묶음별 누적 치환술(모형)과 Kaplan-Meier, 기준 프로필의 분위수
  grid <- seq(0.25, 16, by = 0.25)
  klg <- cut(as.integer(as.character(j$surv$kl)), c(-1, 1, 2, 4), labels = c("KL 0-1", "KL 2", "KL 3-4"))
  res$marginal <- rbindlist(lapply(levels(klg), function(g) {
    rows <- which(klg == g)
    km <- survfit(Surv(time_y, event) ~ 1, data = j$surv[rows, ])
    data.table(group = g, t = grid,
               km = 1 - summary(km, times = grid, extend = TRUE)$surv,
               lognormal = 1 - marg_surv(f_ln, D, grid, rows), su = 1 - marg_surv(f_su, D, grid, rows))
  }))
  qs <- c(0.05, 0.10, 0.25, 0.50)
  res$quantiles <- rbindlist(lapply(c("0", "2", "3", "4"), function(k) {
    w <- setNames(rep(0, D$ps), D$wn); w["(Intercept)"] <- 1; if (k != "0") w[paste0("kl", k)] <- 1
    rbindlist(lapply(list(lognormal = f_ln, su = f_su), function(f)
      data.table(p = qs, years = sapply(qs, function(p) marg_quantile(f, D, w, p)))), idcol = "model")[, kl := k]
  }))
  # 오차 분포의 밀도: 표준화한 오차 epsilon = g(Z) 의 밀도
  # z = (asinh(u) - lambda) / theta, u = sinh lambda + theta cosh lambda epsilon, dz/d epsilon = cosh lambda / sqrt(1 + u^2)
  o_su <- unpack(f_su$par, D, TRUE); ee <- seq(-5, 3, by = 0.02)
  u <- sinh(o_su$lam) + exp(o_su$lth) * cosh(o_su$lam) * ee
  res$density <- data.table(eps = ee, normal = dnorm(ee),
                            su = dnorm((asinh(u) - o_su$lam) / exp(o_su$lth)) * cosh(o_su$lam) / sqrt(1 + u^2))
  say("주변분포 계산 끝")

  # --- 3. 동적 예측: 결합모형과 동적 예측 장의 JMbayes2와 같은 무릎으로 적합하고, 같은 무릎에서 4년 시점의 2년 위험을 예측한다
  jm_path <- file.path(ROOT, "outputs", "jm-01.rds")
  jmr <- if (file.exists(jm_path)) readRDS(jm_path) else NULL
  if (!is.null(jmr) && isTRUE(jmr$real_data == src$real)) {
    fit_knees <- jmr$fit_knees; test_knees <- jmr$landmark$pred$knee
  } else {
    set.seed(20261121L); all_k <- as.character(j$surv$knee)
    fit_knees <- sample(all_k, min(4000L, length(all_k)))
    rest <- setdiff(as.character(j$surv$knee[j$surv$time_y > 4]), fit_knees)
    test_knees <- sample(rest, min(2500L, length(rest))); jmr <- NULL
  }
  tr <- subset_jm(j, fit_knees); te <- subset_jm(j, test_knees)
  Dtr <- prep(tr$surv, tr$long); Dte <- prep(te$surv, te$long)
  g_ln <- fit(Dtr, FALSE, hessian = FALSE)
  g_su <- fit(Dtr, TRUE, start = c(g_ln$par, 0, log(0.2)), hessian = FALSE)
  g_base <- g_ln; g_base$par[c("a0", "a1")] <- 0
  LM <- 4; H <- 2
  last_f <- function(sv, lg) { x <- as.data.table(lg)[year <= LM, .(f_lm = womac_func[which.max(year)]), by = knee]
                               x$f_lm[match(as.character(sv$knee), as.character(x$knee))] }
  ltr <- as.data.table(tr$surv)[time_y > LM]; ltr$f_lm <- last_f(ltr, tr$long)
  ltr[, `:=`(tt = pmin(time_y, LM + H) - LM, ev = as.integer(event == 1 & time_y <= LM + H))]
  cox <- coxph(Surv(tt, ev) ~ age_c + female + bmi_c + kl + f_lm, data = ltr)
  lte <- as.data.table(te$surv); lte$f_lm <- last_f(lte, te$long)
  hz <- copy(lte); hz[, `:=`(tt = H, ev = 0)]
  pred <- data.table(knee = as.character(te$surv$knee), time_y = te$surv$time_y, event = te$surv$event,
                     su = predict_risk(g_su, Dte, LM, H), lognormal = predict_risk(g_ln, Dte, LM, H),
                     baseline = predict_risk(g_base, Dte, LM, H),
                     landmark_cox = 1 - exp(-predict(cox, hz, type = "expected")))
  if (!is.null(jmr)) pred <- merge(pred, jmr$landmark$pred[, .(knee, jmbayes2 = risk)], by = "knee", all.x = TRUE)
  res$prediction <- list(pred = pred, fit_sec = c(lognormal = g_ln$sec, su = g_su$sec), n_fit = Dtr$N, n_fit_events = sum(Dtr$d),
                         jm_minutes = if (!is.null(jmr)) jmr$minutes else NA_real_, t0 = LM, horizon = H,
                         ll = c(lognormal = g_ln$ll, su = g_su$ll))
  say("동적 예측 끝:", nrow(pred), "무릎")

  # --- 4. 모의 검증: SU정규 적합값을 참값으로 두고 데이터를 생성해 다시 추정한다
  if (N_REP > 0) {
    rec <- rbindlist(lapply(seq_len(N_REP), function(sd) {
      set.seed(1000 + sd)
      cens <- ifelse(D$d == 1, sample(exp(D$logT[D$d == 0]), D$N, replace = TRUE), exp(D$logT))
      sm <- simulate(f_su$par, D, TRUE, cens, seed = sd)
      Ds <- D; Ds$logT <- log(sm$time); Ds$d <- sm$event
      kp <- sm$keep; Ds$X <- D$X[kp, , drop = FALSE]; Ds$y <- sm$y[kp]; Ds$t <- D$t[kp]; Ds$g <- D$g[kp]
      Ds$n <- tabulate(Ds$g, D$N)
      Ds$st <- as.vector(rowsum(Ds$t, factor(Ds$g, levels = 1:D$N))); Ds$stt <- as.vector(rowsum(Ds$t^2, factor(Ds$g, levels = 1:D$N)))
      fs <- fit(Ds, TRUE, start = f_su$par)
      say("모의 검증", sd, "/", N_REP)
      data.table(rep = sd, term = names(fs$par), est = unname(fs$par), se = unname(fs$se), truth = unname(f_su$par),
                 events = sum(sm$event), pd = fs$pd)
    }))
    res$recovery <- rec[, .(truth = truth[1], mean = mean(est), sd = sd(est), se = median(se),
                            cover = mean(abs(est - truth) <= 1.96 * se)), by = term]
    res$recovery_info <- rec[, .(events = events[1], pd = pd[1]), by = rep]
  }

  dir.create(file.path(ROOT, "outputs"), showWarnings = FALSE)
  saveRDS(res, file.path(ROOT, "outputs", "su-joint-06.rds"))
  say("저장: outputs/su-joint-06.rds")
}

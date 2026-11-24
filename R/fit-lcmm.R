# WOMAC 통증 궤적의 잠재계층 혼합모형(lcmm::hlme)을 계층 수 1~K로 피팅하고
# 작은 요약만 outputs/lcmm-02.rds에 저장한다. 계층 수가 늘수록 오래 걸려 분석 문서 안에서 돌리지 않는다.
#
#   Rscript R/fit-lcmm.R        계층 1~5
#   Rscript R/fit-lcmm.R 3      계층 1~3

suppressPackageStartupMessages({ library(data.table); library(lcmm) })
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.integer(args[1]) else 5L

ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(ROOT, "R", "oai.R"))

src <- pick_data(ROOT)
d <- build_analysis(src$dir)
pl <- pain_long(d$knees, d$long)

# 고정효과와 계층별 효과: 연수의 2차식. 무작위 효과: 무릎마다 절편과 기울기
fit_k <- function(k, init = NULL) {
  t0 <- Sys.time()
  m <- if (k == 1) {
    hlme(womac_pain ~ year + I(year^2), random = ~ year, subject = "kid", ng = 1, data = pl, verbose = FALSE)
  } else {
    # 국소 최댓값을 피하려고 무작위 초기값 10개로 짧게 돌린 뒤 가장 나은 것에서 끝까지 돌린다
    gridsearch(rep = 10, maxiter = 15, minit = init,
               hlme(womac_pain ~ year + I(year^2), mixture = ~ year + I(year^2), random = ~ year,
                    subject = "kid", ng = k, data = pl, verbose = FALSE))
  }
  cat(sprintf("계층 %d: %.1f분, 수렴 %d, BIC %.1f\n", k,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), m$conv, m$BIC))
  m
}

t0 <- Sys.time()
fits <- list(fit_k(1))
for (k in seq_len(K)[-1]) fits[[k]] <- fit_k(k, init = fits[[1]])
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# 계층 분리의 정도: 상대 엔트로피(1에 가까울수록 잘 갈린다)
entropy <- function(m) {
  if (m$ng == 1) return(NA_real_)
  p <- as.matrix(m$pprob[, grep("^prob", names(m$pprob))])
  p <- pmax(p, 1e-12)
  1 + sum(p * log(p)) / (nrow(p) * log(ncol(p)))
}

grid <- data.frame(year = seq(0, 16, by = 0.25))
summ <- rbindlist(lapply(fits, function(m) {
  cls <- if (m$ng == 1) rep(1L, nrow(m$pprob)) else m$pprob$class
  data.table(ng = m$ng, conv = m$conv, loglik = m$loglik, npar = length(m$best), BIC = m$BIC,
             entropy = entropy(m),
             smallest = min(tabulate(cls, m$ng)) / length(cls))
}))
curves <- rbindlist(lapply(fits, function(m) {
  pr <- predictY(m, newdata = grid, var.time = "year")$pred
  data.table(ng = m$ng, year = grid$year, pr)[
    , melt(.SD, id.vars = c("ng", "year"), variable.name = "class", value.name = "pain")][
    , class := as.integer(gsub("\\D", "", class))]
}))
posterior <- rbindlist(lapply(fits[-1], function(m) {
  p <- as.data.table(m$pprob)
  p[, .(ng = m$ng, kid, class, pmax = do.call(pmax, .SD)), .SDcols = patterns("^prob")]
}))
posterior <- merge(posterior, unique(as.data.table(pl)[, .(kid, knee)]), by = "kid")[, kid := NULL]

res <- list(summary = summ, curves = curves, posterior = posterior,
            n_knees = uniqueN(pl$knee), n_obs = nrow(pl), K = K, minutes = minutes,
            real_data = src$real, fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M"))
dir.create(file.path(ROOT, "outputs"), showWarnings = FALSE)
saveRDS(res, file.path(ROOT, "outputs", "lcmm-02.rds"))
print(summ, digits = 4)
cat(sprintf("전체 %.1f분\n", minutes))

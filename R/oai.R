# OAI(Osteoarthritis Initiative) 공개 임상 데이터를 읽어 분석용 표로 만든다.
#
# 변수 이름의 형식과 확인 상태는 docs/oai-variables.md에 있다. 이 파일의
# VISIT_MONTHS와 VARIABLES는 그 문서의 3장·4장과 같은 값이어야 한다.
#
# 공개 함수
#   read_oai_txt      세로막대 구분 ASCII 파일 하나를 읽는다
#   resolve_columns   변수 형식을 실제 열 이름과 대조한다
#   outcome_columns   Outcomes99.txt의 V99E 열
#   read_oai_cols     파일에서 필요한 열만 읽는다
#   build_analysis    원데이터 폴더에서 분석용 표 세 개를 만든다
#   cohort_data       데이터 1장용 표: 사람, 무릎, 무릎·방문(치환술 뒤 측정 표시)
#   event_table       생존 분석용 무릎의 사건 시간 표(치환술, 사망, 중도절단)
#   exercise_panel    운동 항목과 시간가변 공변량의 사람·방문 표
#   trial_table       운동의 시작(또는 중단)을 처치로 한 모사 시험의 사람·시험 표
#   make_synthetic    OAI와 같은 형식의 시뮬레이션 데이터를 만든다

suppressPackageStartupMessages(library(data.table))

DAYS_PER_MONTH <- 30.4375

# 방문 번호 -> 개월. Outcomes99.txt의 V99RNTCNT 부호 이름으로 확인했다(docs/oai-variables.md 3장).
VISIT_MONTHS <- c("00" = 0, "01" = 12, "02" = 18, "03" = 24, "04" = 30, "05" = 36,
                  "06" = 48, "07" = 60, "08" = 72, "09" = 84, "10" = 96, "11" = 108,
                  "12" = 120, "13" = 168, "14" = 192)

# 치환술 판정 부호. 3 판정 확정만 사건이다(0 판정 불확인, 1 자기 보고).
TKR_CONFIRMED <- 3

# 개념 -> 열 이름의 정규식. {v}는 방문, {s}는 무릎(R/L).
VARIABLES <- c(
  sex         = "P02SEX",
  age         = "V00AGE",
  bmi         = "P01BMI",
  kl          = "V00XRKL",
  womac_func  = "V{v}WOMADL?{s}",
  womac_pain  = "V{v}WOMKP{s}",
  womac_stiff = "V{v}WOMSTF?{s}",
  pase        = "V{v}PASE",
  bmi_visit   = "V{v}BMI",
  tkr_event   = "V99E{s}KRPCF",
  tkr_days    = "V99E{s}KDAYS"
)

SIDES <- c(R = 1L, L = 2L)


# ---------------------------------------------------------------- 읽기

# '1: Male' -> 1, '.: Missing…' -> NA, '27.3' -> 27.3
decode <- function(x) {
  if (!is.character(x)) return(x)
  s <- trimws(x)
  code <- sub("^\\s*(-?[0-9]+(\\.[0-9]+)?)\\s*:.*$", "\\1", s)
  coded <- grepl("^\\s*-?[0-9]+(\\.[0-9]+)?\\s*:", s)
  val <- suppressWarnings(as.numeric(ifelse(coded, code, s)))
  # 숫자도 부호도 아닌 문자열 열(예: 날짜)은 그대로 둔다
  if (all(is.na(val)) && any(!is.na(s)) && !all(startsWith(s[!is.na(s)], "."))) return(x)
  if (!anyNA(val) && all(val == round(val))) val <- as.integer(val)
  val
}

read_oai_txt <- function(path) {
  dt <- fread(path, sep = "|", colClasses = "character", na.strings = "",
              encoding = "Latin-1", showProgress = FALSE)
  setnames(dt, toupper(trimws(names(dt))))
  dt[, names(dt) := lapply(.SD, decode)]
  dt
}

# 대소문자를 가리지 않고 파일을 찾는다(kXR_SQ_BU00.txt와 KXR_SQ_BU00.txt).
find_file <- function(folder, name) {
  f <- list.files(folder, full.names = TRUE)
  hit <- f[tolower(basename(f)) == tolower(name)]
  if (length(hit)) hit[1] else NA_character_
}

header <- function(path) toupper(trimws(names(fread(path, sep = "|", nrows = 0))))


# ---------------------------------------------------------------- 대조

fill <- function(concept, v = "00", s = "R") {
  gsub("{s}", s, gsub("{v}", v, VARIABLES[[concept]], fixed = TRUE), fixed = TRUE)
}

match_column <- function(columns, concept, v = "00", s = "R") {
  hit <- grep(paste0("^", fill(concept, v, s), "$"), columns, value = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

# 맞는 열이 없을 때 고칠 실마리로 쓸 비슷한 열 이름
suggest <- function(columns, concept, k = 8) {
  core <- gsub("[^A-Z]", "", gsub("\\{[vs]\\}", "", VARIABLES[[concept]]))
  head(grep(substr(core, 1, 4), columns, value = TRUE, fixed = TRUE), k)
}

resolve_columns <- function(raw_dir) {
  plan <- data.table(
    concept = c("sex", "age", "bmi", "kl", "womac_func", "womac_func", "womac_pain", "tkr_event", "tkr_days"),
    file    = c("Enrollees.txt", "AllClinical00.txt", "AllClinical00.txt", "kXR_SQ_BU00.txt",
                rep("AllClinical00.txt", 3), "Outcomes99.txt", "Outcomes99.txt"),
    side    = c("R", "R", "R", "R", "R", "L", "R", "R", "R"))
  plan[, c("expected", "found", "hints") := {
    out <- lapply(seq_len(.N), function(i) {
      exp <- fill(concept[i], "00", side[i])
      p <- find_file(raw_dir, file[i])
      if (is.na(p)) return(list(exp, NA_character_, "(파일 없음)"))
      cols <- header(p)
      f <- match_column(cols, concept[i], "00", side[i])
      list(exp, f, if (is.na(f)) paste(suggest(cols, concept[i]), collapse = ", ") else "")
    })
    list(sapply(out, `[[`, 1), sapply(out, `[[`, 2), sapply(out, `[[`, 3))
  }]
  plan[, side := NULL][]
}

outcome_columns <- function(raw_dir) {
  p <- find_file(raw_dir, "Outcomes99.txt")
  if (is.na(p)) return(character())
  grep("^V99E", header(p), value = TRUE)
}


# ---------------------------------------------------------------- 분석용 표

require_column <- function(columns, concept, v = "00", s = "R") {
  col <- match_column(columns, concept, v, s)
  if (is.na(col)) stop(sprintf("%s: '%s'에 맞는 열이 없다. 비슷한 열: %s",
                               concept, fill(concept, v, s), paste(suggest(columns, concept), collapse = ", ")))
  col
}

# 넓은 형식(오른쪽·왼쪽 열)을 무릎 하나에 한 줄로 바꾼다
wide_to_knee <- function(dt, concepts, v) {
  parts <- lapply(names(SIDES), function(s) {
    cols <- vapply(concepts, function(cn) match_column(names(dt), cn, v, s), "")
    cols <- cols[!is.na(cols)]
    if (!length(cols)) return(NULL)
    part <- dt[, c("ID", cols), with = FALSE]
    setnames(part, cols, names(cols))
    part[, SIDE := SIDES[[s]]]
  })
  rbindlist(parts, fill = TRUE)
}

# 원데이터 폴더에서 분석용 표 세 개를 만든다.
#   knees  무릎 하나에 한 줄. 기준 공변량, 사건 여부(event), 사건·중도절단 개월(time)
#   long   무릎·방문 하나에 한 줄. 사건 뒤 측정은 뺀다
#   tv     시간가변 Cox용 구간(start, stop, event, 직전 WOMAC)
build_analysis <- function(raw_dir) {
  # 사람 단위
  enr <- read_oai_txt(find_file(raw_dir, "Enrollees.txt"))
  c00 <- read_oai_txt(find_file(raw_dir, "AllClinical00.txt"))
  person <- enr[, .(ID, sex = get(require_column(names(enr), "sex")))]
  person <- merge(person,
                  c00[, .(ID, age = get(require_column(names(c00), "age")),
                          bmi = get(require_column(names(c00), "bmi")))],
                  by = "ID", all.x = TRUE)
  person[, female := as.numeric(sex == 2)]

  # 무릎 단위 기준 KL. 판독 사업이 여럿이면 15를 먼저 쓴다
  xr <- read_oai_txt(find_file(raw_dir, "kXR_SQ_BU00.txt"))
  klcol <- require_column(names(xr), "kl")
  xr <- xr[!is.na(get(klcol))]
  if ("READPRJ" %in% names(xr)) xr <- xr[order(ID, SIDE, READPRJ != 15)]
  kl <- unique(xr, by = c("ID", "SIDE"))[, .(ID, SIDE, kl = get(klcol))]

  # 방문마다 WOMAC
  long <- rbindlist(lapply(names(VISIT_MONTHS), function(v) {
    p <- find_file(raw_dir, sprintf("AllClinical%s.txt", v))
    if (is.na(p)) return(NULL)
    dt <- if (v == "00") c00 else read_oai_txt(p)
    part <- wide_to_knee(dt, c("womac_func", "womac_pain"), v)
    if (is.null(part) || !nrow(part)) return(NULL)
    part[, `:=`(visit = v, month = VISIT_MONTHS[[v]])]
  }), fill = TRUE)
  long <- long[!is.na(womac_func)]

  # 사건
  out <- read_oai_txt(find_file(raw_dir, "Outcomes99.txt"))
  ev <- wide_to_knee(out, c("tkr_event", "tkr_days"), "99")
  if (is.null(ev) || !all(c("tkr_event", "tkr_days") %in% names(ev)))
    stop("사건 변수를 찾지 못했다. outcome_columns()로 열을 보고 VARIABLES를 고친다.")
  ev[, `:=`(tkr = as.integer(tkr_event %in% TKR_CONFIRMED), tkr_month = tkr_days / DAYS_PER_MONTH)]

  # 무릎 표
  base <- long[visit == "00", .(ID, SIDE, func0 = womac_func, pain0 = womac_pain)]
  knees <- merge(base, kl, by = c("ID", "SIDE"))
  knees <- merge(knees, person, by = "ID", all.x = TRUE)
  knees <- merge(knees, ev[, .(ID, SIDE, tkr, tkr_month)], by = c("ID", "SIDE"), all.x = TRUE)
  knees[is.na(tkr), tkr := 0L]
  knees <- knees[!(tkr == 1 & tkr_month <= 0)]                       # 기준 전 치환술

  # 사건 뒤 측정을 빼고, 중도절단 시각을 정한다
  long <- merge(long, knees[, .(ID, SIDE, tkr, tkr_month)], by = c("ID", "SIDE"))
  long <- long[!(tkr == 1 & month >= tkr_month)]
  knees <- merge(knees, long[, .(last_month = max(month)), by = .(ID, SIDE)], by = c("ID", "SIDE"), all.x = TRUE)
  knees[, `:=`(event = tkr, time = fifelse(tkr == 1, tkr_month, as.numeric(last_month)))]
  knees <- knees[time > 0]
  long <- merge(long, knees[, .(ID, SIDE)], by = c("ID", "SIDE"))[, c("tkr", "tkr_month") := NULL]
  setorder(long, ID, SIDE, month)

  knees[, knee := paste(ID, SIDE, sep = "-")]
  long[, knee := paste(ID, SIDE, sep = "-")]
  list(knees = knees[], long = long[], tv = counting_process(knees, long))
}

# 직전 방문의 WOMAC 기능을 다음 방문까지 유지하는 (start, stop] 구간
counting_process <- function(knees, long) {
  tv <- merge(long[, .(ID, SIDE, knee, start = as.numeric(month), func = womac_func)],
              knees[, .(knee, t_end = time, ev = event)], by = "knee")
  setorder(tv, knee, start)
  tv[, stop := pmin(shift(start, -1, fill = Inf), t_end), by = knee]
  tv <- tv[stop > start]
  tv[, event := as.integer(stop == t_end & ev == 1)]
  tv[, .(ID, SIDE, knee, start, stop, event, func)]
}


# 활동량과 통증의 사람·방문 표. PASE가 있는 방문만 담는다.
#   pase      PASE 총점(지난 7일의 신체활동, 클수록 활동적)
#   pain      두 무릎 가운데 더 아픈 쪽의 WOMAC 통증(0~20)
#   bmi_t     그 방문의 BMI(기준 방문은 P01BMI)
#   기준 공변량: age, female, bmi0, kl(두 무릎 가운데 높은 등급)
# 어느 한쪽이라도 확정된 치환술을 받은 뒤의 방문은 뺀다.
activity_table <- function(raw_dir) {
  enr <- read_oai_txt(find_file(raw_dir, "Enrollees.txt"))
  person <- enr[, .(ID, female = as.numeric(get(require_column(names(enr), "sex")) == 2))]

  xr <- read_oai_txt(find_file(raw_dir, "kXR_SQ_BU00.txt"))
  klcol <- require_column(names(xr), "kl")
  person <- merge(person, xr[!is.na(get(klcol)), .(kl = max(get(klcol))), by = ID], by = "ID", all.x = TRUE)

  out <- read_oai_txt(find_file(raw_dir, "Outcomes99.txt"))
  ev <- wide_to_knee(out, c("tkr_event", "tkr_days"), "99")
  tkr <- ev[tkr_event %in% TKR_CONFIRMED, .(tkr_month = min(tkr_days / DAYS_PER_MONTH)), by = ID]
  person <- merge(person, tkr, by = "ID", all.x = TRUE)

  pv <- rbindlist(lapply(names(VISIT_MONTHS), function(v) {
    p <- find_file(raw_dir, sprintf("AllClinical%s.txt", v))
    if (is.na(p)) return(NULL)
    h <- header(p)
    pc <- match_column(h, "pase", v)
    if (is.na(pc)) return(NULL)
    dt <- read_oai_txt(p)
    pr <- match_column(h, "womac_pain", v, "R"); pl <- match_column(h, "womac_pain", v, "L")
    bc <- if (v == "00") match_column(h, "bmi") else match_column(h, "bmi_visit", v)
    x <- data.table(ID = dt$ID, visit = v, month = VISIT_MONTHS[[v]], pase = as.numeric(dt[[pc]]),
                    pain_r = as.numeric(dt[[pr]]), pain_l = as.numeric(dt[[pl]]),
                    bmi_t = if (is.na(bc)) NA_real_ else as.numeric(dt[[bc]]))
    if (v == "00") x[, age := as.numeric(dt[[require_column(h, "age")]])]
    x
  }), fill = TRUE)
  if (!nrow(pv)) stop("PASE 열이 없다. 시뮬레이션 데이터라면 make_synthetic()으로 다시 만든다.")
  pv[, pain := pmax(pain_r, pain_l, na.rm = TRUE)]
  pv <- merge(pv, pv[visit == "00", .(ID, age0 = age, bmi0 = bmi_t)], by = "ID")
  pv <- merge(pv[, age := NULL], person, by = "ID")
  setnames(pv, "age0", "age")
  pv <- pv[is.na(tkr_month) | month < tkr_month]
  pv <- pv[!is.na(pase) & !is.na(pain)]
  setorder(pv, ID, month)
  pv[, .(ID, visit, month, year = month / 12, pase, pain, pain_r, pain_l, bmi_t, age, female, bmi0, kl)]
}


# 필요한 열만 골라 읽는다. 방문 파일은 열이 천 개가 넘어 전부 읽으면 느리다.
read_oai_cols <- function(path, pattern) {
  h <- header(path)
  idx <- which(h == "ID" | grepl(pattern, h))
  dt <- fread(path, sep = "|", select = idx, colClasses = "character", na.strings = "",
              encoding = "Latin-1", showProgress = FALSE)
  setnames(dt, h[idx])
  dt[, names(dt) := lapply(.SD, decode)]
  dt
}

# 데이터 1장(데이터와 측정)용 표. 치환술 뒤의 측정도 지우지 않고 표시만 한다.
#   person  사람 하나에 한 줄. 하위 코호트, 기관, 기준 특성, 마지막 접촉 방문(개월), 사망 판정과 시각(개월)
#   knees   무릎 하나에 한 줄. 기준 KL, 치환술 여부와 시각(개월)
#   kv      무릎·방문 하나에 한 줄. WOMAC 통증·기능·뻣뻣함, 치환술 뒤 측정 표시(after_tkr)
# 하위 코호트, 기관, 마지막 접촉, 사망 열이 없는 데이터(옛 시뮬레이션 데이터)에서는 그 값이 NA다.
cohort_data <- function(raw_dir) {
  enr <- read_oai_txt(find_file(raw_dir, "Enrollees.txt"))
  opt <- function(dt, col) if (col %in% names(dt)) dt[[col]] else NA
  person <- data.table(ID = enr$ID, female = as.numeric(enr[[require_column(names(enr), "sex")]] == 2),
                       cohort = opt(enr, "V00COHORT"), site = opt(enr, "V00SITE"))

  c00 <- read_oai_cols(find_file(raw_dir, "AllClinical00.txt"), "^(V00AGE|P01BMI|V00PASE|V00EVDATE)$")
  person <- merge(person, data.table(ID = c00$ID, age = as.numeric(opt(c00, "V00AGE")), bmi = as.numeric(opt(c00, "P01BMI")),
                                     pase = as.numeric(opt(c00, "V00PASE")),
                                     enroll_date = suppressWarnings(as.Date(as.character(opt(c00, "V00EVDATE")), "%m/%d/%Y"))),
                  by = "ID", all.x = TRUE)

  out <- read_oai_txt(find_file(raw_dir, "Outcomes99.txt"))
  lc <- suppressWarnings(as.integer(opt(out, "V99RNTCNT")))
  last_month <- unname(VISIT_MONTHS[sprintf("%02d", lc)])
  person <- merge(person, data.table(ID = out$ID, last_month = last_month,
                                     died = as.integer(opt(out, "V99EDDCF") %in% 2),
                                     death_date = suppressWarnings(as.Date(as.character(opt(out, "V99EDDDATE")), "%m/%d/%y"))),
                  by = "ID", all.x = TRUE)
  if (!"V99RNTCNT" %in% names(out)) person[, died := NA_integer_]
  # 사망 시각(등록부터의 개월). 사건 파일의 치환술 일수와 같게 날짜 차이에 1일을 더한다
  person[, death_month := fifelse(died %in% 1, as.numeric(death_date - enroll_date + 1) / DAYS_PER_MONTH, NA_real_)]
  person[, c("enroll_date", "death_date") := NULL]

  xr <- read_oai_txt(find_file(raw_dir, "kXR_SQ_BU00.txt"))
  klcol <- require_column(names(xr), "kl")
  xr <- xr[!is.na(get(klcol))]
  if ("READPRJ" %in% names(xr)) xr <- xr[order(ID, SIDE, READPRJ != 15)]
  kl <- unique(xr, by = c("ID", "SIDE"))[, .(ID, SIDE, kl = get(klcol))]

  ev <- wide_to_knee(out, c("tkr_event", "tkr_days"), "99")
  ev[, `:=`(tkr = as.integer(tkr_event %in% TKR_CONFIRMED), tkr_month = tkr_days / DAYS_PER_MONTH)]
  knees <- CJ(ID = person$ID, SIDE = unname(SIDES))
  knees <- merge(knees, kl, by = c("ID", "SIDE"), all.x = TRUE)
  knees <- merge(knees, ev[, .(ID, SIDE, tkr, tkr_month)], by = c("ID", "SIDE"), all.x = TRUE)
  knees[is.na(tkr), tkr := 0L]
  knees[tkr == 0, tkr_month := NA_real_]

  kv <- rbindlist(lapply(names(VISIT_MONTHS), function(v) {
    p <- find_file(raw_dir, sprintf("AllClinical%s.txt", v))
    if (is.na(p)) return(NULL)
    dt <- read_oai_cols(p, sprintf("^V%sWOM(KP|ADL|STF)[RL]$", v))
    part <- wide_to_knee(dt, c("womac_pain", "womac_func", "womac_stiff"), v)
    if (is.null(part) || !nrow(part)) return(NULL)
    part[, `:=`(visit = v, month = VISIT_MONTHS[[v]])]
  }), fill = TRUE)
  setnames(kv, c("womac_pain", "womac_func", "womac_stiff"), c("pain", "func", "stiff"), skip_absent = TRUE)
  if (!"stiff" %in% names(kv)) kv[, stiff := NA_real_]
  kv <- kv[!(is.na(pain) & is.na(func) & is.na(stiff))]
  kv <- merge(kv, knees[, .(ID, SIDE, tkr_month)], by = c("ID", "SIDE"), all.x = TRUE)
  kv[, after_tkr := !is.na(tkr_month) & month >= tkr_month][, tkr_month := NULL]
  setorder(kv, ID, SIDE, month)
  list(person = person[], knees = knees[], kv = kv[])
}


# 무릎의 사건 시간 표. 생존 분석의 장들이 함께 쓴다. cohort_data()의 결과를 받는다.
#   time     사건 또는 중도절단까지의 개월
#   status   0 중도절단, 1 치환술, 2 사망
#   year, klf, func10, age10, bmi5, knee   분석에 쓰는 변환
# 대상: 기준 KL 판독, 기준 WOMAC 기능, 나이, 성별, BMI가 있고 등록 전에 치환술을 받지 않은 무릎.
# 중도절단: 그 사람의 마지막 접촉 방문과 그 무릎의 마지막 WOMAC 방문 가운데 늦은 쪽. 추적은 마지막 방문 개월에서 끝낸다.
# 사망은 마지막 접촉 뒤 다음 예정 방문 전에 일어났을 때만 센다. 그보다 늦은 사망은 마지막 접촉에서 중도절단한다.
event_table <- function(d) {
  P <- d$person; K <- d$knees; KV <- d$kv
  CAP <- max(VISIT_MONTHS)
  base <- KV[month == 0 & !after_tkr, .(ID, SIDE, pain0 = pain, func0 = func)]
  S <- merge(K[!is.na(kl) & (is.na(tkr_month) | tkr_month > 0)], base, by = c("ID", "SIDE"))
  S <- merge(S, P[, .(ID, age, female, bmi, last_month, died, death_month)], by = "ID")
  S <- S[!is.na(func0) & !is.na(age) & !is.na(bmi) & !is.na(female)]
  lastw <- KV[!after_tkr & !is.na(func), .(lastw = max(month)), by = .(ID, SIDE)]
  S <- merge(S, lastw, by = c("ID", "SIDE"), all.x = TRUE)
  S[, cens := pmax(fifelse(is.na(last_month), 0, last_month), fifelse(is.na(lastw), 0, lastw))]
  next_visit <- function(m) vapply(m, function(x) { v <- VISIT_MONTHS[VISIT_MONTHS > x]; if (length(v)) min(v) else CAP }, 1)
  S[, nextv := next_visit(cens)]
  S[, `:=`(time = cens, status = 0L)]
  S[!is.na(death_month) & death_month <= nextv, `:=`(time = death_month, status = 2L)]
  S[tkr == 1 & (status != 2 | tkr_month <= death_month), `:=`(time = tkr_month, status = 1L)]
  S[time > CAP, `:=`(time = CAP, status = 0L)]
  S <- S[time > 0]
  S[, `:=`(year = time / 12, klf = factor(kl, 0:4), func10 = func0 / 10, age10 = (age - 61) / 10, bmi5 = (bmi - 28.5) / 5,
           knee = paste(ID, SIDE, sep = "-"))]
  setorder(S, ID, SIDE)
  S[]
}

# ---------------------------------------------------------------- 모사 시험 (docs/causal-inference-plan.md 1장, 6장)

# PASE를 잰 0~96개월 방문. 이웃한 세 방문이 직전·지표·결과 방문이 된다.
EXERCISE_VISITS <- c("00", "01", "03", "05", "06", "08", "10")

# 방문마다 읽는 사람 단위 변수: 이름 -> V{v} 뒤의 열 이름. 그 방문에 열이 없으면 NA다.
# 기준 방문의 BMI와 손 통증은 접두어가 P01이다.
EXERCISE_VARS <- c(
  pase = "PASE", sit = "PASE1", walk = "PASE2", sport_l = "PASE3", sport_m = "PASE4", sport_s = "PASE5",
  strength = "PASE6", strength_hr = "PASE6HR",
  pain_r = "WOMKPR", pain_l = "WOMKPL", func_r = "WOMADLR", func_l = "WOMADLL",
  stiff_r = "WOMSTFR", stiff_l = "WOMSTFL", koosqol = "KOOSQOL",
  cesd = "CESD", comorb = "COMORB", bmi = "BMI", pcs = "HSPSS", mcs = "HSMSS",
  nsaid = "NSAIDS", pnmed = "PNMEDT", narc = "RXNARC", inj = "KNINJ",
  pace = "20MPACE", chair = "CSPACE", quad_r = "REMAXF", quad_l = "LEMAXF",
  employed = "CUREMP", hand_r = "OJPNRH", hand_l = "OJPNLH")

# 시험 표에 직전 방문(_p)과 지표 방문(_i)의 값으로 넣는 시간가변 공변량
TRIAL_COVARIATES <- c("pase", "sit", "walk", "sport_l", "sport_m", "sport_s", "strength",
                      "pain", "func", "stiff", "koosqol", "cesd", "comorb", "bmi", "pcs", "mcs",
                      "nsaid", "pnmed", "narc", "inj", "pace", "chair", "quad", "employed", "hand")

# 운동 항목과 시간가변 공변량의 사람·방문 표.
#   person  사람 하나에 한 줄. 고정 특성과 가장 이른 확정 치환술의 개월(tkr_month)
#   pv      사람·방문 하나에 한 줄. pain, func, stiff는 두 무릎 가운데 나쁜 쪽, quad는 두 다리의 평균,
#           hand는 어느 한 손이라도 통증이 있으면 1, 없으면 0. 동반질환 지수는 격년으로 재므로 직전 값을 유지한다
exercise_panel <- function(raw_dir) {
  opt <- function(dt, col) if (col %in% names(dt)) dt[[col]] else NA
  enr <- read_oai_txt(find_file(raw_dir, "Enrollees.txt"))
  person <- data.table(ID = enr$ID, female = as.numeric(enr[[require_column(names(enr), "sex")]] == 2),
                       race = as.numeric(opt(enr, "P02RACE")), cohort = as.numeric(opt(enr, "V00COHORT")),
                       site = as.numeric(factor(as.character(opt(enr, "V00SITE")))))
  c00 <- read_oai_cols(find_file(raw_dir, "AllClinical00.txt"),
                       "^(V00AGE|V00EDCV|V00INCOME|V00[RL]KALNMT|P01INJ[RL]|P01KSURG[RL])$")
  num <- function(col) as.numeric(opt(c00, col))
  worse <- function(a, b) { x <- pmax(a, b, na.rm = TRUE); x[is.infinite(x)] <- NA; x }
  person <- merge(person, data.table(ID = c00$ID, age = num("V00AGE"), edu = num("V00EDCV"), income = num("V00INCOME"),
                                     align_r = num("V00RKALNMT"), align_l = num("V00LKALNMT"),
                                     injury = worse(num("P01INJR"), num("P01INJL")),
                                     surgery = worse(num("P01KSURGR"), num("P01KSURGL"))),
                  by = "ID", all.x = TRUE)

  xr <- read_oai_txt(find_file(raw_dir, "kXR_SQ_BU00.txt"))
  klcol <- require_column(names(xr), "kl")
  person <- merge(person, xr[!is.na(get(klcol)), .(kl = as.numeric(max(get(klcol)))), by = ID], by = "ID", all.x = TRUE)

  out <- read_oai_txt(find_file(raw_dir, "Outcomes99.txt"))
  ev <- wide_to_knee(out, c("tkr_event", "tkr_days"), "99")
  person <- merge(person, ev[tkr_event %in% TKR_CONFIRMED, .(tkr_month = min(tkr_days / DAYS_PER_MONTH)), by = ID],
                  by = "ID", all.x = TRUE)

  pv <- rbindlist(lapply(EXERCISE_VISITS, function(v) {
    p <- find_file(raw_dir, sprintf("AllClinical%s.txt", v))
    if (is.na(p)) return(NULL)
    pre <- if (v == "00") "(V00|P01)" else paste0("V", v)
    dt <- read_oai_cols(p, sprintf("^%s(%s)$", pre, paste(EXERCISE_VARS, collapse = "|")))
    x <- data.table(ID = dt$ID, visit = v, month = VISIT_MONTHS[[v]])
    for (nm in names(EXERCISE_VARS)) {
      col <- intersect(paste0(c(paste0("V", v), "P01"), EXERCISE_VARS[[nm]]), names(dt))
      x[, (nm) := if (length(col)) as.numeric(dt[[col[1]]]) else NA_real_]
    }
    x
  }))
  if (all(is.na(pv$strength)) || all(is.na(pv$walk)))
    stop("PASE 항목(PASE2, PASE6) 열이 없다. 시뮬레이션 데이터라면 make_synthetic()으로 다시 만든다.")
  pv[, `:=`(pain = worse(pain_r, pain_l), func = worse(func_r, func_l), stiff = worse(stiff_r, stiff_l),
            hand = worse(hand_r, hand_l), quad = rowMeans(cbind(quad_r, quad_l), na.rm = TRUE))]
  pv[is.nan(quad), quad := NA_real_]
  # 다른 관절의 통증은 표시한 칸만 1이고 표시하지 않은 칸은 결측 부호다. 그 방문에 응답한 사람의 빈칸은 0으로 본다
  pv[, hand := fifelse(!is.na(hand), 1, fifelse(!is.na(pase) | !is.na(pain), 0, NA_real_))]
  setorder(pv, ID, month)
  pv[, comorb := nafill(comorb, type = "locf"), by = ID]
  list(person = person[], pv = pv[])
}

# 모사 시험의 사람·시험 표. 이웃한 세 방문(직전, 지표, 결과)마다 시험 하나를 만들어 쌓는다.
#   treatment  "strength": 근력운동 항목이 0이다가 1 이상(주 1일 이상)이 됨
#              "walk": 운동 삼아 걷기 항목이 1 이하(주 2일 이하)이다가 2 이상(주 3일 이상)이 됨
#   direction  "start": 직전 방문에 하지 않던 사람 가운데 지표 방문에 하는 사람이 처치군(a = 1)
#              "stop": 직전 방문에 하던 사람 가운데 지표 방문에 하지 않는 사람이 처치군
# 대상: 직전 방문의 상태가 조건에 맞고, 지표 방문의 그 항목과 통증이 관측되고, 지표 방문까지 치환술이 없다.
# 결과(y_*)는 결과 방문의 값이며, 결과 방문까지 치환술을 받았으면 NA이고 tkr_out이 1이다.
# 공변량은 직전 방문(_p)과 지표 방문(_i)의 값이다. pain_pp는 직전 방문보다 한 방문 앞의 통증이다(첫 시험은 NA).
trial_table <- function(raw_dir, treatment = c("strength", "walk"), direction = c("start", "stop"), panel = NULL) {
  treatment <- match.arg(treatment); direction <- match.arg(direction)
  on <- c(strength = 1, walk = 2)[[treatment]]
  if (is.null(panel)) panel <- exercise_panel(raw_dir)
  pv <- panel$pv
  at <- function(v, vars, suffix) {
    x <- pv[visit == v, c("ID", vars), with = FALSE]
    setnames(x, vars, paste0(vars, suffix))
    x
  }
  V <- EXERCISE_VISITS
  tr <- rbindlist(lapply(2:(length(V) - 1), function(k) {
    x <- merge(at(V[k - 1], TRIAL_COVARIATES, "_p"), at(V[k], TRIAL_COVARIATES, "_i"), by = "ID")
    x <- merge(x, at(V[k + 1], c("pain", "func", "pace", "chair", "hand"), "_y"), by = "ID", all.x = TRUE)
    x <- merge(x, if (k > 2) at(V[k - 2], "pain", "_pp") else data.table(ID = x$ID, pain_pp = NA_real_), by = "ID", all.x = TRUE)
    x[, `:=`(trial = VISIT_MONTHS[[V[k]]], month_p = VISIT_MONTHS[[V[k - 1]]], month_y = VISIT_MONTHS[[V[k + 1]]])]
  }))
  tr <- merge(tr, panel$person, by = "ID")
  prev <- tr[[paste0(treatment, "_p")]]; idx <- tr[[paste0(treatment, "_i")]]
  tr[, a := if (direction == "start") as.integer(idx >= on) else as.integer(idx < on)]
  tr <- tr[!is.na(prev) & (if (direction == "start") prev < on else prev >= on) &
             !is.na(a) & !is.na(pain_i) & (is.na(tkr_month) | tkr_month > trial)]
  tr[, tkr_out := as.integer(!is.na(tkr_month) & tkr_month <= month_y)]
  for (y in c("pain", "func", "pace", "chair", "hand")) {
    tr[, (paste0("y_", y)) := fifelse(tkr_out == 1, NA_real_, get(paste0(y, "_y")))]
    tr[, (paste0(y, "_y")) := NULL]
  }
  tr[, `:=`(age_i = age + trial / 12, gap = (month_y - trial) / 12)]
  setorder(tr, trial, ID)
  tr[]
}

# ---------------------------------------------------------------- 공용 준비

# data/raw에 OAI 파일이 있으면 그것을, 없으면 시뮬레이션 데이터를 쓴다
pick_data <- function(root) {
  raw <- file.path(root, "data", "raw")
  if (dir.exists(raw) && !is.na(find_file(raw, "AllClinical00.txt")))
    return(list(dir = raw, real = TRUE))
  syn <- file.path(root, "data", "synthetic")
  if (is.na(find_file(syn, "AllClinical00.txt"))) make_synthetic(syn)
  list(dir = syn, real = FALSE)
}

# 결합모형용 표. 시간은 연, 나이·BMI는 중심화, 무릎은 요인
jm_data <- function(knees, long) {
  k <- knees[complete.cases(knees[, .(age, bmi, female)])]
  kj <- k[, .(knee = factor(knee), time_y = time / 12, event,
              age_c = age - 61, female, bmi_c = bmi - 28.5, kl = factor(kl))]
  lj <- merge(long[knee %in% k$knee, .(knee, year = month / 12, womac_func)],
              k[, .(knee, kl = factor(kl))], by = "knee")
  lj[, knee := factor(knee, levels = levels(kj$knee))]
  setorder(lj, knee, year)
  list(surv = as.data.frame(kj), long = as.data.frame(lj))
}

# 통증 궤적용 긴 표. lcmm은 정수 대상자 번호(kid)를 쓴다
pain_long <- function(knees, long) {
  pl <- long[!is.na(womac_pain), .(knee, year = month / 12, womac_pain)]
  pl <- pl[knee %in% knees$knee]
  pl[, kid := as.integer(factor(knee))]
  setorder(pl, kid, year)
  as.data.frame(pl)
}


# ---------------------------------------------------------------- 시뮬레이션 데이터

# OAI와 같은 파일·열·부호 형식의 시뮬레이션 데이터를 만든다.
# 참값: 기능 점수는 KL 등급에 따라 수준과 기울기가 다른 선형 궤적이고,
# 치환술 위험은 그 시점의 참 기능 점수와 KL 등급에 따라 오른다.
# 수치는 그럴듯하게 고른 것이며 OAI에서 추정한 값이 아니다.
make_synthetic <- function(out_dir, n = 1500, seed = 20261118) {
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  miss <- ".: Missing Form/Incomplete Workbook"

  ids <- 9000000L + seq_len(n)
  sex <- sample(1:2, n, replace = TRUE)
  age <- round(pmin(pmax(rnorm(n, 61, 9), 45), 79))
  bmi <- round(pmin(pmax(rnorm(n, 28.5, 4.8), 17), 48), 1)

  knees <- lapply(names(SIDES), function(s) {
    kl <- sample(0:4, n, replace = TRUE, prob = c(0.38, 0.18, 0.22, 0.15, 0.07))
    list(side = s, kl = kl,
         b0 = -4 + 5.5 * kl + 0.25 * (bmi - 28.5) + rnorm(n, 0, 7),   # 기준 참 기능
         b1 = 0.2 + 0.45 * kl + rnorm(n, 0, 1.2))                     # 연간 변화
  })
  true_func <- function(k, months) pmin(pmax(outer(k$b0, rep(1, length(months))) +
                                              outer(k$b1, months / 12), 0), 68)

  # 치환술: 월 단위 이산 위험
  grid <- seq(0.5, 119.5, by = 1)
  dropout <- ifelse(runif(n) < 0.25, runif(n, 12, 110), Inf)
  for (i in seq_along(knees)) {
    k <- knees[[i]]
    h <- exp(-9.6 + 0.07 * true_func(k, grid) + 0.35 * k$kl + 0.02 * (age - 61))
    hit <- matrix(runif(length(h)) < h, nrow = n)
    first <- ifelse(rowSums(hit) > 0, grid[max.col(hit, ties.method = "first")], Inf)
    knees[[i]]$tkr_month <- ifelse(first < dropout, first, Inf)
  }

  # 방문 파일
  sub <- runif(n) < 0.3                                              # 18·30개월 일부 대상자
  for (v in names(VISIT_MONTHS)[VISIT_MONTHS <= 96]) {                # 시뮬레이션 데이터는 96개월까지
    m <- VISIT_MONTHS[[v]]
    dt <- data.table(ID = ids)
    if (v == "00") dt[, `:=`(V00AGE = age, P01BMI = bmi)]
    attend <- (m < dropout) & (runif(n) > if (v == "00") 0 else 0.08)
    if (v %in% c("02", "04")) attend <- attend & sub
    for (k in knees) {
      obs <- pmin(pmax(round(true_func(k, m)[, 1] + rnorm(n, 0, 4)), 0), 68)
      pain <- pmin(pmax(round(obs / 68 * 20 + rnorm(n, 0, 2)), 0), 20)
      ok <- attend & (m < k$tkr_month)
      dt[, (sprintf("V%sWOMADL%s", v, k$side)) := ifelse(ok, as.character(obs), miss)]
      dt[, (sprintf("V%sWOMKP%s", v, k$side)) := ifelse(ok, as.character(pain), miss)]
    }
    # PASE와 방문별 BMI. 주 난수 흐름을 건드리지 않도록 따로 씨앗을 준다(기존 시뮬레이션 값이 바뀌지 않는다).
    # 참값: 활동량은 사람마다 수준이 다르고 해마다 2점 줄며, 그 방문의 통증 1점에 3점 줄어든다.
    # 활동량이 통증에 주는 효과는 없다.
    if (!v %in% c("02", "04")) {
      keep <- .Random.seed
      set.seed(seed + 1000L + as.integer(v))
      worst <- pmax(suppressWarnings(as.numeric(dt[[sprintf("V%sWOMKPR", v)]])),
                    suppressWarnings(as.numeric(dt[[sprintf("V%sWOMKPL", v)]])), na.rm = TRUE)
      set.seed(seed + 999L); level <- rnorm(n, 170, 60)
      set.seed(seed + 1000L + as.integer(v))
      pase <- pmax(round(level - 2 * m / 12 - 3 * ifelse(is.finite(worst), worst, 0) + rnorm(n, 0, 45)), 0)
      dt[, (sprintf("V%sPASE", v)) := ifelse(attend, as.character(pase), miss)]
      if (v != "00") dt[, (sprintf("V%sBMI", v)) := ifelse(attend, as.character(round(bmi + rnorm(n, 0, 0.8), 1)), miss)]
      # PASE 항목(걷기 PASE2, 근력운동 PASE6). 사람마다 성향이 있고 그 방문의 통증이 높으면 덜 한다.
      # 운동이 통증에 주는 효과는 없다.
      lab <- c("0: Never", "1: Seldom (1-2 days)", "2: Sometimes (3-4 days)", "3: Often (5-7 days)")
      for (item in c(2L, 6L)) {
        set.seed(seed + 3000L + item); liking <- rnorm(n)
        set.seed(seed + 3100L + 20L * item + as.integer(v))
        z <- liking + rnorm(n, 0, 0.8) - 0.05 * ifelse(is.finite(worst), worst, 0) + if (item == 2L) 0.7 else -0.3
        dt[, (sprintf("V%sPASE%d", v, item)) := ifelse(attend, lab[findInterval(z, c(-0.2, 0.5, 1.2)) + 1L], miss)]
      }
      .Random.seed <<- keep
    }
    fwrite(dt, file.path(out_dir, sprintf("AllClinical%s.txt", v)), sep = "|")
  }

  # 하위 코호트, 기관, 마지막 접촉, 사망. 주 난수 흐름을 건드리지 않도록 따로 씨앗을 준다.
  keep <- .Random.seed
  set.seed(seed + 2000L)
  cohort <- sample(c("1: Progression", "2: Incidence", "3: Non-exposed control group"), n, replace = TRUE,
                   prob = c(0.29, 0.68, 0.03))
  site <- sample(LETTERS[1:5], n, replace = TRUE)
  vm <- VISIT_MONTHS[VISIT_MONTHS <= 96]
  last_idx <- vapply(dropout, function(d) max(as.integer(names(vm))[vm < d]), 1L)
  last_code <- sprintf("%d: %d-month", last_idx, vm[sprintf("%02d", last_idx)])
  died <- is.finite(dropout) & runif(n) < 0.3
  .Random.seed <<- keep

  fwrite(data.table(ID = ids, P02SEX = ifelse(sex == 1, "1: Male", "2: Female"), V00COHORT = cohort, V00SITE = site),
         file.path(out_dir, "Enrollees.txt"), sep = "|")

  xr <- rbindlist(lapply(knees, function(k)
    data.table(ID = ids, SIDE = if (k$side == "R") "1: Right" else "2: Left",
               READPRJ = 15L, V00XRKL = sprintf("%d: %d", k$kl, k$kl))))
  fwrite(xr, file.path(out_dir, "kXR_SQ_BU00.txt"), sep = "|")

  o <- data.table(ID = ids)
  for (k in knees) {
    hit <- is.finite(k$tkr_month)
    o[, (sprintf("V99E%sKRPCF", k$side)) := ifelse(hit, "3: Replacement adjudicated, confirmed", miss)]
    o[, (sprintf("V99E%sKDAYS", k$side)) := ifelse(hit, as.character(round(k$tkr_month * DAYS_PER_MONTH)),
                                                   ".A: Not Applicable")]
  }
  o[, `:=`(V99RNTCNT = last_code, V99EDDCF = ifelse(died, "2: Death adjudicated, confirmed", miss))]
  fwrite(o, file.path(out_dir, "Outcomes99.txt"), sep = "|")
  invisible(out_dir)
}

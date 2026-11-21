# 렌더된 GitHub용 MD를 손본다. chapters/_quarto.yml의 post-render가 부른다.
#
# pandoc은 원문의 \~를 맨 ~로 내보낸다. GitHub은 한 줄에 맨 ~가 둘 이상이면
# 그 사이를 취소선으로 그리므로(research-writing-principles.md §II.3),
# 코드와 수식 밖의 ~를 다시 \~로 바꾼다.

escape_tilde <- function(line) {
  # 인라인 코드와 수식 구간은 건드리지 않는다
  parts <- regmatches(line, gregexpr("`[^`]*`|\\$\\$[^$]*\\$\\$|\\$[^$]+\\$", line), invert = NA)[[1]]
  for (i in seq(1, length(parts), by = 2)) parts[i] <- gsub("(?<!\\\\)~", "\\\\~", parts[i], perl = TRUE)
  paste(parts, collapse = "")
}

fix_file <- function(path) {
  x <- readLines(path, encoding = "UTF-8", warn = FALSE)
  fence <- cumsum(grepl("^```", x)) %% 2 == 1 | grepl("^```", x)
  img <- grepl("^!\\[", x)
  x[!fence & !img] <- vapply(x[!fence & !img], escape_tilde, "")
  writeLines(enc2utf8(x), path, useBytes = TRUE)
}

files <- strsplit(Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES"), "\n")[[1]]
files <- files[grepl("\\.md$", files)]
for (f in files) fix_file(f)

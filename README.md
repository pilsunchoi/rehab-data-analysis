# 재활 데이터 분석 (집필 중)

재활과 관련된 공개 자료를 내용 중심으로 분석하는 저서의 저장소다. 장마다 재활의 질문 하나를 세우고, 답하는 데 필요한 방법을 수식과 가정까지 설명한 뒤, 자료로 추정하고 해석하는 데이터 분석 교과서로 쓴다. 첫 편은 무릎 골관절염이고 Osteoarthritis Initiative(OAI) 공개 코호트를 쓴다. 분석은 모두 R로 한다.

소개 대시보드: <https://pilsunchoi.github.io/rehab-data-analysis/>

장 구성과 서술 형식은 `docs/book-plan.md`에 있다.

## 구조

```
rehab-data-analysis/
├── README.md
├── chapters/                                 무릎 골관절염 편(5개 장)
│   ├── _quarto.yml                           렌더 설정(렌더 뒤 MD 손질)
│   ├── 01-data-and-measures.qmd              1장 원본(R 코드 포함)
│   ├── 01-data-and-measures.md               1장 렌더 결과(GitHub에서 읽는 본문)
│   ├── 01-data-and-measures_files/           1장 그림
│   ├── 02-structure-and-symptoms.qmd         2장 원본
│   ├── 02-structure-and-symptoms.md          2장 렌더 결과
│   └── 02-structure-and-symptoms_files/      2장 그림
├── docs/
│   ├── index.html                            소개 대시보드(GitHub Pages)
│   ├── img/                                  대시보드의 그림(장의 그림을 복사)
│   ├── book-plan.md                          장 구성 계획
│   ├── oai-variables.md                      받을 파일, 변수 이름의 꼴, 사건 정의
│   └── stroke-data-sources.md                뇌졸중 재활 데이터 후보와 접근 조건
├── R/
│   ├── oai.R                                 OAI 파일 읽기, 분석용 표, 모의 데이터
│   ├── postrender-md.R                       렌더된 MD를 GitHub용으로 손질(물결 이스케이프)
│   ├── fit-jm.R                              결합모형 적합(수십 분), 요약을 outputs/에 저장
│   └── fit-lcmm.R                            통증 궤적 계층 적합(몇 시간), 요약을 outputs/에 저장
├── archive/                                  2026-11-18~21의 첫 분석 세 편(고치지 않는다)
├── outputs/                                  오래 걸리는 적합의 요약(공개하지 않는다)
└── data/                                     공개하지 않는다
    ├── raw/                                  OAI에서 받은 파일(이름을 바꾸지 않는다)
    └── synthetic/                            모의 데이터
```

`archive/`의 세 편은 장 구성을 정하기 전에 쓴 첫 분석이다. 그 내용은 2\~5장으로 옮겨 다시 쓴다(`docs/book-plan.md` 2장).

## 렌더하기

장의 `.md`는 `.qmd` 원본을 렌더해 만든 것이므로 손으로 고치지 않는다. `chapters/`에서 렌더한다. R이 PATH에 없으면 `QUARTO_R`로 R의 `bin` 폴더를 알려 준다. 1장은 1분쯤 걸린다.

```bash
cd chapters
quarto render 01-data-and-measures.qmd
```

OAI 데이터가 없어도 장은 끝까지 돈다. `data/raw/`가 비어 있으면 같은 파일·열·부호 꼴의 모의 데이터를 만들어 읽고, 문서 첫머리에 그 사실을 적는다. 모의 데이터의 수치는 해석하지 않는다.

필요한 패키지는 다음과 같다.

```r
install.packages(c("data.table", "survival", "lme4", "nlme", "JMbayes2", "lcmm", "ggplot2", "gridExtra", "scales", "knitr", "rmarkdown"))
```

## 원데이터로 바꾸기

1. NIMH Data Archive(NDA)에 Research Auth Service로 로그인하고 OAI 이용 약관에 동의한다.
2. `Enrollees.txt`, `AllClinical00.txt`\~`AllClinical14.txt`, `kXR_SQ_BU00.txt`, `Outcomes99.txt`를 받아 `data/raw/`에 둔다. ASCII(세로막대 구분) 판을 받는다.
3. 변수 이름이 맞지 않으면 `R/oai.R`의 `VARIABLES`를 고친다. 맞는 열이 없을 때 함수가 비슷한 열 이름을 보여 준다.

## 주의

- 변수 이름, 사건 부호, 방문과 개월의 대응은 2026-11-20에 받은 OAI 데이터로 확인했다. 확인 상태와 남은 결정은 `docs/oai-variables.md`에 있다.
- 사건은 판정 확정(부호 3)된 치환술만 센다.
- OAI 데이터는 이용 약관에 따라 다시 배포하지 않는다. `data/`와 `outputs/`는 git에 올리지 않는다.
- `analysis/Choi/`(생존분포 방법 연구)는 공개 전이라 git에 올리지 않는다.

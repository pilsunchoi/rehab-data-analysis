# 재활 데이터 분석 I (집필 중)

재활과 관련된 공개 데이터로 분석 방법을 가르치는 교재의 저장소다. 데이터에 기본적인 분석 방법을 하나씩 적용하고, 방법마다 수식과 가정, 추정 절차, 결과를 읽는 법을 자세히 설명한다. 장은 방법의 계열(회귀, 종단 분석, 생존 분석, 인과 추론, 그래프 모형)로 묶는다. 그 첫 번째인 I권은 무릎 골관절염 공개 코호트인 Osteoarthritis Initiative(OAI)를 다룬다. 분석은 모두 R로 한다.

소개 대시보드: <https://pilsunchoi.github.io/rehab-data-analysis/>

## 구조

```
rehab-data-analysis/
├── README.md
├── chapters/                                 무릎 골관절염 편(12개 장)
│   ├── _quarto.yml                           렌더 설정(렌더 뒤 MD 손질)
│   ├── 01-data-and-measures.qmd              데이터 1장 원본(R 코드 포함)
│   ├── 01-data-and-measures.md               데이터 1장 렌더 결과(GitHub에서 읽는 본문)
│   ├── 01-data-and-measures_files/           데이터 1장 그림
│   ├── 02-structure-and-symptoms.qmd         회귀 1장 원본
│   ├── 02-structure-and-symptoms.md          회귀 1장 렌더 결과
│   ├── 02-structure-and-symptoms_files/      회귀 1장 그림
│   ├── 03-pain-trajectories.qmd              종단 분석 1장 원본(잠재계층 적합은 outputs/lcmm-02.rds를 읽는다)
│   ├── 03-pain-trajectories.md               종단 분석 1장 렌더 결과
│   ├── 03-pain-trajectories_files/           종단 분석 1장 그림
│   ├── 04-activity-and-pain.qmd              인과 추론 1장 원본
│   ├── 04-activity-and-pain.md               인과 추론 1장 렌더 결과(그림 없음)
│   ├── 05-time-to-replacement.qmd            생존 분석 1장 원본
│   ├── 05-time-to-replacement.md             생존 분석 1장 렌더 결과
│   ├── 05-time-to-replacement_files/         생존 분석 1장 그림
│   ├── 06-su-normal-joint-model.qmd          SU정규 결합모형 장 원본(적합은 outputs/su-joint-06.rds를 읽는다)
│   ├── 06-su-normal-joint-model.md           SU정규 결합모형 장 렌더 결과
│   ├── 06-su-normal-joint-model_files/       SU정규 결합모형 장 그림
│   ├── 07-propensity-score.qmd               인과 추론 2장 원본(성향점수와 이중 강건 추정)
│   ├── 07-propensity-score.md                인과 추론 2장 렌더 결과
│   ├── 07-propensity-score_files/            인과 추론 2장 그림
│   ├── 08-causal-ml.qmd                      인과 추론 3장 원본(교차 적합과 인과 숲)
│   ├── 08-causal-ml.md                       인과 추론 3장 렌더 결과
│   ├── 08-causal-ml_files/                   인과 추론 3장 그림
│   ├── 09-competing-risks.qmd                생존 분석 2장 원본(경쟁 위험)
│   ├── 09-competing-risks.md                 생존 분석 2장 렌더 결과(그림 없음)
│   ├── 10-joint-model.qmd                    결합모형과 동적 예측 장 원본(outputs/jm-01.rds를 읽는다)
│   ├── 10-joint-model.md                     결합모형과 동적 예측 장 렌더 결과(그림 없음)
│   ├── 11-graphical-model.qmd                그래프 모형 1장 원본(부분상관 네트워크)
│   ├── 11-graphical-model.md                 그래프 모형 1장 렌더 결과
│   ├── 11-graphical-model_files/             그래프 모형 1장 그림
│   ├── 12-logistic-regression.qmd            회귀 2장 원본(로지스틱 회귀)
│   ├── 12-logistic-regression.md             회귀 2장 렌더 결과
│   └── 12-logistic-regression_files/         회귀 2장 그림
├── docs/
│   ├── index.html                            소개 대시보드(GitHub Pages)
│   ├── img/                                  대시보드의 그림(장의 그림을 복사)
│   ├── oai-variables.md                      받을 파일, 변수 이름의 형식, 사건 정의
│   └── stroke-data-sources.md                뇌졸중 재활 데이터 후보와 접근 조건(I권에서는 쓰지 않는다)
├── R/
│   ├── oai.R                                 OAI 파일 읽기, 분석용 표, 시뮬레이션 데이터
│   ├── postrender-md.R                       렌더된 MD를 GitHub용으로 손질(물결 이스케이프)
│   ├── fit-jm.R                              결합모형 적합(수십 분), 요약을 outputs/에 저장
│   ├── fit-causal.R                          운동의 시작·중단을 처치로 쌓은 시험의 평균 효과(13분). 장에서는 쓰지 않는다
│   ├── fit-su-joint.R                        SU정규 결합모형 장의 적합과 모의 검증(25분 안팎)
│   └── fit-lcmm.R                            통증 궤적 계층 적합(몇 시간), 요약을 outputs/에 저장
├── archive/                                  2026-11-18~21의 첫 분석 세 편(고치지 않는다)
├── outputs/                                  오래 걸리는 적합의 요약(공개하지 않는다)
└── data/                                     공개하지 않는다
    ├── raw/                                  OAI에서 받은 파일(이름을 바꾸지 않는다)
    └── synthetic/                            시뮬레이션 데이터
```

`archive/`의 세 편은 장 구성을 정하기 전에 쓴 첫 분석이다. 그 내용은 회귀, 종단 분석, 인과 추론, 생존 분석의 장으로 옮겨 다시 썼다.

## 렌더하기

장의 `.md`는 `.qmd` 원본을 렌더해 만든 것이므로 손으로 고치지 않는다. `chapters/`에서 렌더한다. R이 PATH에 없으면 `QUARTO_R`로 R의 `bin` 폴더를 알려 준다. 데이터 1장은 1분쯤 걸린다.

```bash
cd chapters
quarto render 01-data-and-measures.qmd
```

종단 분석 1장의 잠재계층 궤적(몇 시간)과 결합모형과 동적 예측 장의 결합모형(무릎 4,000개로 20분 안팎)은 적합이 오래 걸려 장 안에서 하지 않는다. 먼저 `Rscript R/fit-lcmm.R`과 `Rscript R/fit-jm.R`로 적합하면 요약이 `outputs/lcmm-02.rds`와 `outputs/jm-01.rds`에 저장되고, 장은 그 파일을 읽는다. SU정규 결합모형 장은 `Rscript R/fit-su-joint.R`의 요약(`outputs/su-joint-06.rds`)을 읽는다. 적합 자체는 1분 남짓이지만 모의 검증 20회가 함께 돈다. 생존 분석 1장과 같은 무릎으로 예측을 비교하므로 `R/fit-jm.R`을 먼저 돌린다. 파일이 없으면 그 절은 안내만 남기고 나머지는 그대로 돈다.

OAI 데이터가 없어도 장은 끝까지 돈다. `data/raw/`가 비어 있으면 같은 파일·열·부호 형식의 시뮬레이션 데이터를 만들어 읽고, 문서 첫머리에 그 사실을 적는다. 시뮬레이션 데이터의 수치는 해석하지 않는다.

필요한 패키지는 다음과 같다.

```r
install.packages(c("data.table", "survival", "lme4", "nlme", "JMbayes2", "lcmm", "grf", "MatchIt", "glasso", "ggplot2", "gridExtra", "scales", "knitr", "rmarkdown"))
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

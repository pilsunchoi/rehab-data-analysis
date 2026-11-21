# 재활 데이터 분석 (집필 중)

재활과 관련된 공개 자료를 내용 중심으로 분석하는 저서의 저장소다. 첫 주제는 무릎 골관절염이고, Osteoarthritis Initiative(OAI) 공개 코호트로 기능과 통증의 궤적, 신체활동, 인공관절 치환술을 본다. 분석은 모두 R로 한다.

소개 대시보드: <https://pilsunchoi.github.io/rehab-data-analysis/>

`analysis/Choi/`(생존분포 방법 연구, Python)와 `outputs/`(적합 결과)는 공개 전이라 git에 올리지 않는다(`.gitignore`). 아래 구조도에는 로컬 작업본 기준으로 함께 적었다.

## 구조

```
rehab-data-analysis/
├── README.md
├── docs/
│   ├── index.html                            소개 대시보드(GitHub Pages)
│   ├── img/                                  대시보드의 그림(분석 문서의 그림을 복사)
│   ├── oai-variables.md                      받을 파일, 변수 이름의 꼴, 사건 정의
│   └── stroke-data-sources.md                뇌졸중 재활 데이터 후보와 접근 조건
├── R/
│   ├── oai.R                                 OAI 파일 읽기, 분석용 표(무릎 단위, 사람·방문 단위), 모의 데이터
│   ├── fit-jm.R                              결합모형 적합(수십 분), 요약을 outputs/에 저장
│   └── fit-lcmm.R                            통증 궤적 계층 적합(몇 시간), 요약을 outputs/에 저장
├── analysis/
│   ├── 01-function-trajectory-tkr-method.md  방법 해설(수식, 인용)
│   ├── 01-function-trajectory-tkr.qmd        첫 분석의 원본(R 코드 포함)
│   ├── 01-function-trajectory-tkr.md         렌더 결과(글·표·그림)
│   ├── 01-function-trajectory-tkr_files/     렌더 결과의 그림
│   ├── 02-knee-pain-trajectory.qmd           통증 분석의 원본: 구조와 통증, 궤적 계층
│   ├── 02-knee-pain-trajectory.md            렌더 결과
│   ├── 02-knee-pain-trajectory_files/        렌더 결과의 그림
│   ├── 03-activity-pain.qmd                  활동량과 통증 분석의 원본: PASE와 WOMAC 통증
│   ├── 03-activity-pain.md                   렌더 결과
│   ├── 03-activity-pain_files/               렌더 결과의 그림
│   └── Choi/                                 TSU 생존모형(Python)
│       ├── *.docx                            2024년 원고와 메모(원본, 고치지 않는다)
│       ├── TSU Survival (vs KM).ipynb        2024년 노트북(원본)
│       ├── tsu-distribution.md               모형: 공변량 없는 TSU 분포
│       ├── tsu-regression.md                 모형: 기준 분포와 TSU 회귀
│       ├── tsu-simulation.md                 시뮬레이션의 설계와 결과
│       ├── tsu-reg-simulation.md             회귀모형의 소규모 시뮬레이션
│       ├── tsu-aft-simulation.md             수정한 회귀모형(TSU-AFT, SAS-AFT)의 비교
│       ├── su-multivariate-ideas.md          메모: 다변량 SU를 쓸 수 있는 분야
│       ├── su-joint-model.md                 다변량 정규 결합모형의 시험(OAI)
│       ├── su-joint.R                        그 코드(R, 약 25분)
│       ├── su-joint-literature.md            메모: 결합모형의 선행 문헌
│       ├── tsu.py                            모형 함수
│       ├── tsu-check.py                      수식 검증
│       ├── tsu-sim.py                        시뮬레이션(약 1시간 반)
│       ├── tsu-reg-sim.py                    회귀모형 시뮬레이션(약 7분)
│       ├── tsu-reg-ci.py                     회귀모형의 구간추정 시뮬레이션(약 45분)
│       └── tsu-aft-sim.py                    수정한 회귀모형의 비교(약 16분)
├── outputs/
│   ├── jm-01.rds                             결합모형 요약
│   ├── lcmm-02.rds                           통증 궤적 계층 요약
│   ├── tsu-sim-summary.csv                   TSU 시뮬레이션 요약(문서의 표)
│   ├── tsu-sim-raw.csv.gz                    TSU 시뮬레이션의 반복별 결과
│   ├── tsu-reg-sim-summary.csv               TSU 회귀 시뮬레이션 요약
│   ├── tsu-reg-sim-raw.csv.gz                TSU 회귀 시뮬레이션의 반복별 결과
│   ├── tsu-reg-ci-summary.csv                TSU 회귀 구간추정 요약
│   ├── tsu-reg-ci-raw.csv.gz                 TSU 회귀 구간추정의 반복별 결과
│   ├── tsu-aft-sim-summary.csv               수정한 회귀모형 비교 요약
│   ├── tsu-aft-sim-raw.csv.gz                수정한 회귀모형 비교의 반복별 결과
│   └── su-joint-01.rds                       다변량 정규 결합모형 시험의 요약
└── data/                                     공개하지 않는다
    ├── raw/                                  OAI에서 받은 파일(이름을 바꾸지 않는다)
    └── synthetic/                            분석 문서가 만드는 모의 데이터
```

## 시작하기

OAI 데이터 없이도 분석 문서는 끝까지 돈다. `data/raw/`가 비어 있으면 같은 파일·열·부호 꼴의 모의 데이터를 만들어 읽는다. 모의 데이터의 수치는 해석하지 않는다.

모형과 정의의 설명은 `analysis/01-function-trajectory-tkr-method.md`에 있다. 이 문서는 손으로 쓰고 고치는 문서이며 추정값을 옮겨 적지 않는다. 결과는 `analysis/01-function-trajectory-tkr.md`에서 바로 읽는다. 이 MD는 `.qmd` 원본을 렌더해 만든 것이므로 손으로 고치지 않고, 원본을 고쳐 다시 렌더한다.

결합모형은 적합에 수십 분이 걸리므로 먼저 따로 돌린다. 그다음 RStudio에서 `analysis/01-function-trajectory-tkr.qmd`를 열고 Render를 누른다. 명령줄에서는 다음과 같다(R이 PATH에 없으면 `QUARTO_R`로 R의 `bin` 폴더를 알려 준다).

```bash
Rscript R/fit-jm.R
```

두 번째 분석(통증 궤적)은 계층 적합이 따로 있다. 계층 1~5에 약 3시간 반이 걸렸다(2026-11-20).

```bash
Rscript R/fit-lcmm.R
```

```bash
quarto render analysis/01-function-trajectory-tkr.qmd
```

세 번째 분석(활동량과 통증)은 따로 돌릴 적합이 없고 렌더에 1분쯤 걸린다.

```bash
quarto render analysis/03-activity-pain.qmd
```

필요한 패키지는 다음과 같다.

```r
install.packages(c("data.table", "survival", "lme4", "nlme", "JMbayes2", "lcmm", "ggplot2", "gridExtra", "knitr", "rmarkdown"))
```

## 원데이터로 바꾸기

1. NIMH Data Archive(NDA)에 Research Auth Service로 로그인하고 OAI 이용 약관에 동의한다.
2. `Enrollees.txt`, `AllClinical00.txt`\~`AllClinical10.txt`, `kXR_SQ_BU00.txt`, `Outcomes99.txt`를 받아 `data/raw/`에 둔다. ASCII(세로막대 구분) 판을 받는다.
3. 분석 문서 1절의 변수 대조가 통과할 때까지 `R/oai.R`의 `VARIABLES`를 고친다.

## TSU 생존모형

`analysis/Choi/`는 0에서 절단한 Johnson $S_U$ 분포(TSU)로 생존함수를 추정하는 방법 연구다. 2024년의 원고와 메모(DOCX), 노트북이 원본이고, 2026-11-20에 모형을 MD로 옮겨 완성하고 시뮬레이션을 새로 설계했다. OAI 데이터를 쓰지 않는다.

모형은 `tsu-distribution.md`(분포, 우도, 식별)와 `tsu-regression.md`(공변량이 있는 모형)에, 시뮬레이션의 설계와 결과는 `tsu-simulation.md`에, 회귀모형의 소규모 시뮬레이션은 `tsu-reg-simulation.md`에, 그 결과로 모형을 고쳐 일반화 감마와 비교한 것은 `tsu-aft-simulation.md`에 있다. 지금까지의 판정은 마지막 문서의 IV장에 있다. 원본에서 고친 곳은 두 모형 문서의 마지막 장에 적었다. 표의 수치는 `outputs/`의 요약 파일 넷(`tsu-sim-summary.csv`, `tsu-reg-sim-summary.csv`, `tsu-reg-ci-summary.csv`, `tsu-aft-sim-summary.csv`)에서 옮긴 것이므로, 시뮬레이션을 다시 돌리면 문서의 표와 본문 수치도 함께 고친다.

수식 검증은 몇 초, 분포 시뮬레이션은 프로세스 7개로 약 1시간 반, 회귀 시뮬레이션은 약 7분, 구간추정 시뮬레이션은 약 45분, 수정한 모형의 비교는 약 16분이 걸린다. `analysis/Choi/`에서 실행한다.

```bash
python tsu-check.py
```

```bash
python tsu-sim.py
```

```bash
python tsu-reg-sim.py
```

```bash
python tsu-reg-ci.py
```

```bash
python tsu-aft-sim.py
```

필요한 패키지는 `numpy`, `scipy`, `pandas`이고 `tsu-check.py`는 `lifelines`도 쓴다.

## 주의

- 변수 이름, 사건 부호, 방문과 개월의 대응은 2026-11-20에 받은 OAI 데이터로 확인했다. 확인 상태와 남은 결정은 `docs/oai-variables.md`에 있다.
- 사건은 판정 확정(부호 3)된 치환술만 센다.
- OAI 데이터는 이용 약관에 따라 다시 배포하지 않는다. `data/`는 git에 올리지 않는다.

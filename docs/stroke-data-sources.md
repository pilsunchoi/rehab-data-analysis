# 뇌졸중 재활 데이터 후보 — 공개·신청형 데이터의 접근 조건

작성일: 2026-11-20
대상: 뇌졸중 재활 연구에 쓸 수 있는 공개 데이터와 신청형 데이터
목적: 두 질문(아래 0장)에 맞는 데이터를 접근 난이도별로 고른다. 접근 조건은 2026-11-20 웹 검색으로 확인한 것이며, 각 기관 페이지를 직접 열어 확인하지 않은 것은 미확인으로 적는다.

---

## 0. 목적과 범위

이 목록은 두 질문을 기준으로 데이터를 고른다.

| 질문 | 필요한 것 |
|---|---|
| A. 회복 궤적과 비례 회복 규칙 | 같은 사람의 운동 기능 점수(Fugl-Meyer, 보행 속도 등)를 여러 시점에 잰 것 |
| B. 재활 시점과 강도의 인과 효과 | 재활을 언제, 얼마나 받았는지와 그 뒤의 결과. 교란을 보정할 기준 상태 |

A는 무작위 시험 데이터만으로도 할 수 있다. B는 처치가 무작위가 아닌 관찰 데이터가 필요하고, 치료 시간이나 수가 기록이 있어야 한다.

판단이 갈리면 다음 순서를 따른다. 첫째, 질문에 맞는 측정이 있는가. 둘째, 지금 받을 수 있는가. 셋째, 표본이 큰가.

---

## 1. 바로 받을 수 있는 것 (2026-11-20 확인)

가입이나 약정서만 내면 받는다. 심사 기간이 거의 없다.

| 데이터 | 내용 | 접근 | 질문 |
|---|---|---|---|
| IST-1 (International Stroke Trial) | 급성 뇌졸중 19,435명, 아스피린·헤파린 시험. 6개월 의존 여부, 추적 99% | Edinburgh DataShare에서 공개(DOI 10.7488/ds/104) | 재활 측정 없음. 순서형 결과와 결측 방법론 연습용 |
| RAND HRS Longitudinal File | 미국 50세 이상 패널, 2년 간격. 의사 진단 뇌졸중(자기·대리 보고, 발생 연월)과 ADL 다섯 항목 | HRS 계정 가입 뒤 내려받기 | A의 변형. 뇌졸중 전후 ADL 궤적(사건 연구 설계) |
| 고령화연구패널(KLoSA) | 한국 45세 이상 패널, 2006년부터 2년 간격. 건강 상태, 의료 이용, ADL | KOSSDA, 한국고용정보원 KEIS, Gateway to Global Aging(가입) | HRS와 같은 설계의 한국판. 뇌졸중 변수의 정확한 이름은 미확인 |
| ATLAS v2.0 | 뇌졸중 병변 T1 MRI와 수작업 병변 마스크. 논문은 1,271건, ICPSR은 공개 655건으로 적는다 | ICPSR 36684. 데이터 이용 약정서를 내야 한다 | 임상 결과가 거의 없어 단독으로는 맞지 않는다. 병변으로 회복을 예측할 때 짝으로 쓴다 |

HRS와 KLoSA는 재활 처치가 없고 기능도 자기 보고 ADL이다. 대신 뇌졸중이 일어나기 전의 궤적까지 있다는 점이 다른 데이터에 없는 강점이다. 뇌졸중 몇 해 전부터 ADL 저하 속도가 달라진다는 HRS 연구가 이미 있다.

---

## 2. 신청해 승인받는 것

심사가 있다. 기간과 비용은 대부분 미확인이다.

| 데이터 | 내용 | 접근 | 질문 |
|---|---|---|---|
| ICARE (NINDS 보관) | 상지 재활 3상 무작위 시험. 1차 결과 Wolf Motor Function Test(1년), 상지 Fugl-Meyer를 기준·중재 후·6·12개월에 측정 | NINDS Archived Clinical Research Datasets. 이용 조건에 동의하고 신청. 절차 세부와 파일 구성은 미확인 | A에 가장 잘 맞는다 |
| LEAPS (NINDS 보관) | 보행 재활 3상 무작위 시험(Locomotor Experience Applied Post-Stroke) | ICARE와 같은 경로. 세부 미확인 | A. 하지·보행 회복 |
| VISTA-Rehab | 재활 시험을 모은 보관소. 한 개관은 10,194명으로 적는다 | 운영위원회에 분석 계획을 내 승인받는다. 데이터는 글래스고대 Robertson Centre가 보관. 근거 문서가 2009\~2013년이라 현재 절차는 미확인 | A. 여러 시험을 합친 큰 표본 |
| IST-3 | 혈전용해 시험 | Edinburgh DataShare의 통제 접근 신청 | B와는 거리가 있다 |
| MIMIC-IV | 미국 한 병원의 중환자실·입원 전자의무기록 | PhysioNet 자격 심사. CITI "Data or Specimens Only Research" 교육 수료 보고서와 이용 약정, 보통 1\~2주 | 입원 중 재활 협진 시점. 뇌졸중 환자와 재활 기록이 얼마나 있는지는 미확인 |
| SSNAP (영국 국가 뇌졸중 감사) | 영국 뇌졸중 입원 전수. 물리·작업·언어치료를 실제로 받은 시간을 28일 단위로 기록 | 질문지를 내고, 1차 심사를 통과하면 심사비를 낸다. 이어 HQIP의 월례 심사(DARG)를 받는다 | B에 가장 잘 맞는다. 치료 용량-반응 |
| 국민건강보험 표본코호트 DB | 자격, 진료, 건강검진, 요양기관 DB. 2014\~2023년 100만 명 저위험 표본데이터를 시범 제공한다는 보도가 있다 | 건강보험 빅데이터 플랫폼 신청. 비용과 재활 수가 코드의 포함 여부는 미확인 | B의 한국판. 재활 시작 시점과 횟수를 처치로, 재입원·요양시설·사망을 결과로 |

---

## 3. 협업해야 하는 것

| 데이터 | 내용 | 접근 |
|---|---|---|
| KOSCO (Korean Stroke Cohort for Functioning and Rehabilitation) | 첫 뇌졸중 환자 7,858명(9개 지역 병원, 2012\~2015년 등록), 10년 추적. K-MBI, Fugl-Meyer, FAC, K-MMSE 등을 반복 측정 | 원데이터 비공개. 질병관리청 국립보건연구원 규정에 따라 연구진과 협업해야 한다 |

A와 B 모두에 가장 이상적인 데이터다. 기능 점수를 반복 측정했고, 입원 재활 강도와 장애의 관계를 다룬 논문도 이미 나왔다. 다른 데이터로 방법을 먼저 세우고 그 결과를 들고 협업을 제안하는 순서가 현실적이다.

---

## 4. 권하는 순서

1. HRS와 KLoSA를 지금 받는다. 뇌졸중 전후 ADL 궤적으로 사건 연구 설계와 궤적 모형을 세운다. 무릎 연구의 궤적·바닥 효과 처리를 그대로 쓴다.
2. 같은 때에 NINDS에 ICARE와 LEAPS를 신청한다. 승인되면 질문 A(비례 회복 규칙)를 시험 데이터로 다룬다.
3. 질문 B는 건강보험 표본코호트의 재활 수가 코드를 확인한 뒤 정한다. 국외라면 SSNAP이 치료 시간 기록으로 가장 낫지만 심사비와 영국 기관 요건이 걸린다.
4. 1\~2의 결과로 KOSCO 연구진에 협업을 제안한다.

---

## 5. 미확인 사항

1. NINDS 보관 데이터의 신청 절차, 그리고 ICARE·LEAPS 파일에 실제로 들어 있는 변수.
2. VISTA-Rehab의 현재 신청 절차와 연락처.
3. KLoSA와 HRS의 뇌졸중 변수 이름. KLoSA는 진단 문항의 범위(뇌혈관질환 전체인지)도 확인한다.
4. 건강보험 표본코호트의 비용과 재활 수가 코드 포함 여부.
5. MIMIC-IV에 있는 뇌졸중 환자 수와 재활 협진 기록.

---

## 6. 출처 (2026-11-20 검색)

- IST-1: [The International Stroke Trial database (PMC3320519)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3320519/), [에든버러대 공개 소식](https://medicine-vet-medicine.ed.ac.uk/news-events/latest-news/data-from-international-stroke-trial-released)
- IST-3: [IST-3 데이터 공개 소식](https://clinical-brain-sciences.ed.ac.uk/news/latest-news/ist3-trial-data-available)
- HRS: [RAND HRS Longitudinal File](https://www.rand.org/health/surveys/hrs/hrs-data.html), [뇌졸중 전 ADL 저하 연구](https://experts.umn.edu/en/publications/stroke-associated-differences-in-rates-of-activity-of-daily-livin/)
- KLoSA: [KOSSDA](https://kossda.snu.ac.kr/handle/20.500.12236/14660), [Gateway to Global Aging 안내](https://g2aging.org/hrd/get-data/0e7d72286d2840a7a7ffe9672a9f24aa/instruction)
- ATLAS: [ICPSR 36684](https://www.icpsr.umich.edu/web/ICPSR/studies/36684/versions/V4), [데이터 논문 (PMC9203460)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9203460)
- NINDS: [Archived Clinical Research Datasets](https://www.ninds.nih.gov/current-research/research-funded-ninds/clinical-research/archived-clinical-research-datasets), [ICARE Fugl-Meyer 재척도화 논문 (PMC6517555)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6517555)
- VISTA-Rehab: [VISTA-Rehab: a resource for stroke rehabilitation trials](https://ueaeprints.uea.ac.uk/id/eprint/43719/), [2009년 포스터](https://www.gcu.ac.uk/__data/assets/pdf_file/0024/22389/vista.pdf)
- MIMIC-IV: [MIMIC 접근 안내](https://mimic.mit.edu/docs/faq/how-to-get-access.html)
- SSNAP: [데이터 요청](https://strokeaudit.org/Resources/Using-SSNAP-data-in-research/Data-requests.aspx), [치료 빈도 안내](https://ssnap.zendesk.com/hc/en-us/sections/23534755876893-Domain-6-Therapy-frequency)
- 건강보험: [저위험 표본데이터 설명회 보도](https://medifonews.com/news/article.html?no=214824)
- KOSCO: [장기 기능 회복 (PMC9508656)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9508656), [설계 (PMC5698661)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5698661)

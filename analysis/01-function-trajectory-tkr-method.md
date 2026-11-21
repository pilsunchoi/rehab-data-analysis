# 무릎 기능 궤적과 인공관절 치환술: 방법 해설

이 문서는 첫 분석(`analysis/01-function-trajectory-tkr.qmd`와 그 렌더 결과 `01-function-trajectory-tkr.md`)이 쓰는 데이터, 정의, 모형을 수식과 함께 설명한다. 추정값은 렌더 결과에 있고 이 문서에는 옮겨 적지 않는다. 데이터를 바꾸어 다시 렌더하면 수치는 바뀌지만 이 문서의 설명은 그대로 유지된다. 파일과 변수 이름은 `docs/oai-variables.md`가, 코드는 `R/oai.R`과 `R/fit-jm.R`이 정한다.

---

## I. 질문과 데이터

### 1. 질문

무릎 골관절염에서 인공관절 치환술(total knee replacement, TKR)은 증상과 기능이 보존적 치료로 더 이상 관리되지 않을 때 고려되는 최종 단계의 치료다. 따라서 치환술의 시점은 영상에서 보이는 구조적 손상만이 아니라 환자가 겪는 통증과 기능 저하의 경과에 의해서도 결정된다. 이 연구의 질문은 하나다. 반복 측정한 무릎 기능 점수의 궤적, 곧 어느 시점의 수준과 그 변화 속도가 치환술의 위험을 얼마나 예측하는가. 이 질문은 종단 데이터(기능 점수)와 사건 시간 데이터(치환술까지의 시간)를 함께 다루어야 답할 수 있고, IV장의 결합모형이 그 표준적인 도구다.

### 2. Osteoarthritis Initiative

데이터는 미국 국립보건원이 지원한 Osteoarthritis Initiative(OAI)다. OAI는 무릎 골관절염이 있거나 생길 위험이 높은 45\~79세 성인을 등록해 매년 추적한 다기관 전향적 코호트이며, 임상 평가, 자기 보고 설문, 방사선 영상과 MRI를 방문마다 수집했다(Nevitt et al., 2006). 데이터는 NIMH Data Archive에서 이용 약관에 동의한 연구자에게 공개된다. 방문은 등록 시점(V00)부터 번호가 붙고, 치환술 같은 사건은 방문과 별도로 판정해 `Outcomes99` 파일에 모은다. 방문 번호와 개월의 대응, 변수 이름 가운데 확인한 것과 확인하지 못한 것은 `docs/oai-variables.md` 3장과 4장에 있다.

### 3. 측정

기능 점수는 Western Ontario and McMaster Universities Osteoarthritis Index(WOMAC)의 신체 기능 하위척도다(Bellamy et al., 1988). 계단 오르내리기, 의자에서 일어서기, 양말 신기 같은 일상 활동 17개 문항에 대해 지난 48시간 동안의 곤란 정도를 0(없음)\~4(극심) 다섯 단계로 답하고, 그 합을 점수로 쓴다. 따라서 점수는 0\~68의 정수이고 클수록 기능이 나쁘다. 같은 설문의 통증 하위척도는 5개 문항으로 0\~20, 뻣뻣함 하위척도는 2개 문항으로 0\~8이다. OAI는 이 설문을 무릎마다 따로 받았으므로 기능 점수는 무릎 단위의 반복 측정이 된다.

구조적 중증도는 방사선 영상의 Kellgren-Lawrence(KL) 등급으로 측정한다(Kellgren and Lawrence, 1957). 0은 골관절염 소견이 없음, 1은 의심스러운 골극, 2는 확실한 골극, 3은 중등도의 관절 간격 협소, 4는 심한 협소와 골경화를 뜻하는 순서 척도다. 첫 분석은 등록 시점의 KL 등급을 기준 공변량으로 쓰고, 등급 사이의 간격이 같다고 가정하지 않도록 범주형으로 넣는다.

### 4. 분석 단위와 표기

분석 단위는 무릎이다. 무릎을 $i = 1, \dots, n$으로, 무릎 $i$가 속한 사람을 $p(i)$로 적는다. 한 사람에게 무릎이 둘이므로 같은 $p$를 공유하는 두 무릎은 독립이 아니며, 이 의존은 III.2절과 V.2절에서 다룬다. 시간 $t$는 등록부터 측정한 연수다.

무릎 $i$의 기능 점수는 방문 시각 $t_{i1} < t_{i2} < \dots < t_{in_i}$에서 관측되며, 관측값을 $y_{ij} = y_i(t_{ij})$로 적는다. 치환술까지의 참 시간을 $T_i^*$, 중도절단 시간을 $C_i$라 하면 관측되는 것은 $T_i = \min(T_i^*, C_i)$와 사건 지시자 $\delta_i = I(T_i^* \le C_i)$다. 중도절단 시간은 그 무릎의 기능 점수가 마지막으로 관측된 방문의 시각으로 정했다(`docs/oai-variables.md` 5장). 치환술 뒤의 기능 점수는 인공관절을 측정한 값이므로 쓰지 않는다. 따라서 모든 $j$에 대해 $t_{ij} < T_i$다. 기준 공변량(나이, 성별, BMI, KL 등급)은 벡터 $w_i$로 적는다.

---

## II. 기술 분석

### 1. 바닥 효과

분석 문서 3절의 그림은 기준 시점 WOMAC 기능 점수의 분포를 보여 준다. 왼쪽은 전체 무릎의 히스토그램이고, 오른쪽은 KL 등급별 경험 누적분포함수 $\hat F_k(y) = n_k^{-1} \sum_{i: \mathrm{KL}_i = k} I(y_i \le y)$다. 골관절염이 없거나 가벼운 무릎은 일상 활동에 곤란이 없다고 답하는 경우가 많아 0에 확률 질량이 몰린다. 이것이 척도의 바닥 효과다.

바닥 효과는 IV장에서 정규 오차를 가정하는 종단 모형과 충돌한다. 점수가 0에 막혀 있으면 관측값의 조건부 분포는 0에서 잘린 꼴이 되고, 0 근처 무릎의 평균 궤적은 실제보다 완만하게 추정된다. 한 가지 처리는 관측되지 않는 잠재 점수 $y^*$를 두고 $y = \max(0, y^*)$로 관측된다고 보는 검열 정규 모형(Tobin, 1958)이다. 이 모형에서 $y_{ij} = 0$인 관측의 우도 기여는 $\Phi(-m_{ij}/\sigma)$, $y_{ij} > 0$인 관측의 기여는 $\sigma^{-1}\phi\{(y_{ij} - m_{ij})/\sigma\}$다. 여기서 $m_{ij}$는 잠재 점수의 평균, $\Phi$와 $\phi$는 표준정규 분포함수와 밀도다. 첫 분석은 이 처리를 하지 않았고, V.1절에서 다음 단계로 둔다.

### 2. 궤적

분석 문서 4절의 왼쪽 그림은 KL 등급 $k$와 방문 시각 $t$마다 관측된 기능 점수의 평균 $\bar y_k(t)$를 이은 것이다. 이 평균은 그 방문까지 치환술을 받지 않고 추적이 유지된 무릎만으로 계산되므로, 위험이 큰 무릎이 차례로 빠져나가는 선택이 섞여 있다. 그림의 기울기를 개별 무릎의 변화 속도로 읽으면 안 되는 이유가 이것이다.

오른쪽 그림은 치환술을 받은 무릎만 골라 시간축을 수술 시점에 맞춘 것이다. 무릎 $i$의 관측을 $u_{ij} = t_{ij} - T_i$에 놓으면 $u_{ij} < 0$이고, 회색 선은 무릎별 궤적, 빨간 선은 $u$를 1년 구간으로 나눈 평균이다. 이 그림은 미래의 사건으로 표본을 고른 것이므로 예측 모형의 근거가 되지는 않는다. 수술에 가까워질수록 기능이 빠르게 나빠지는지라는, 결합모형이 잡아야 할 신호의 모양을 보여 주는 기술 통계로만 쓴다.

### 3. Kaplan-Meier 추정

분석 문서 5절의 그림은 기준 KL 등급별로 치환술 없이 지낸 비율 $S(t) = P(T^* > t)$의 Kaplan-Meier 추정량(Kaplan and Meier, 1958)이다. 서로 다른 사건 시각을 $\tau_1 < \tau_2 < \dots$, 각 시각 직전에 위험에 놓인 무릎 수를 $r_m$, 그 시각의 사건 수를 $d_m$이라 하면 다음과 같다.

$$\hat S(t) = \prod_{m: \tau_m \le t} \left(1 - \frac{d_m}{r_m}\right)$$

이 추정량은 중도절단이 사건과 독립이라는 가정 아래 일치한다. 첫 분석은 사망을 중도절단으로 다루었는데, 사망한 사람은 치환술을 받을 수 없으므로 사망은 독립 중도절단이 아니라 경쟁 위험이다. 이때 $1 - \hat S(t)$는 실제 치환술 누적 발생률보다 크게 나온다(Putter et al., 2007). V.3절에서 다시 다룬다.

---

## III. Cox 모형

### 1. 비례위험모형

III장의 세 모형은 모두 Cox (1972)의 비례위험모형이다. 무릎 $i$의 위험함수를 다음과 같이 둔다.

$$h_i(t) = h_0(t) \exp(x_i(t)^\top \beta)$$

$h_0(t)$는 모양을 정하지 않은 기준 위험이고 $x_i(t)$는 공변량이다. $\beta$는 부분우도를 최대화해 추정한다. 사건 시각에서 위험집합을 $R(\tau_m) = \{i : T_i \ge \tau_m\}$이라 하면 동시 사건이 없을 때 부분우도는 다음과 같다.

$$L(\beta) = \prod_{i: \delta_i = 1} \frac{\exp(x_i(T_i)^\top \beta)}{\sum_{l \in R(T_i)} \exp(x_l(T_i)^\top \beta)}$$

부분우도는 각 사건 시각에 사건이 난 무릎이 위험집합 안의 다른 무릎보다 먼저 사건을 겪을 조건부 확률의 곱이므로 $h_0(t)$가 약분되어 사라지고, 기준 위험의 모양을 정하지 않고도 $\beta$를 추정할 수 있다. 사건 시각을 개월이나 일로 반올림하면 동시 사건이 생기며, 분석은 이를 Efron (1977)의 근사로 처리한다. $\exp(\beta_k)$는 다른 공변량이 같을 때 $x_k$가 한 단위 커지는 데 따른 위험비다.

### 2. 군집 강건 분산

한 사람의 두 무릎은 같은 유전, 체중, 활동량, 통증 민감성을 공유하므로 사건 시간이 양의 상관을 가질 가능성이 크다. 첫 분석은 무릎을 독립으로 놓은 작업 독립 부분우도로 $\beta$를 추정하고, 분산만 사람 단위로 보정한다(Lin and Wei, 1989). 무릎 $i$의 점수 잔차를 $U_i(\hat\beta)$, 관측 정보행렬을 $I(\hat\beta)$라 하면 군집 강건 분산은 다음과 같다.

$$\widehat{\mathrm{Var}}(\hat\beta) = I(\hat\beta)^{-1} \left\{ \sum_{p} \Big( \sum_{i: p(i) = p} U_i(\hat\beta) \Big) \Big( \sum_{i: p(i) = p} U_i(\hat\beta) \Big)^\top \right\} I(\hat\beta)^{-1}$$

이 보정은 표준오차만 고치고 계수의 해석은 바꾸지 않는다. 계수는 사람 안의 상관을 모형으로 다룬 조건부 효과가 아니라 무릎 집단 수준의 효과다(Therneau and Grambsch, 2000). R에서는 `coxph(..., cluster = ID)`가 이 추정량을 준다.

### 3. 모형 A: 기준 공변량

모형 A의 공변량은 기준 시점의 값뿐이다.

$$h_i(t) = h_0(t) \exp\big(w_i^\top \gamma + \beta_y \, y_{i1}/10\big)$$

$y_{i1}$은 등록 시점의 기능 점수이고, 10으로 나눈 것은 계수를 10점당 위험비로 읽기 위해서다. $w_i$의 KL 등급은 0등급을 기준으로 한 지시변수 4개로 들어간다. 이 모형은 등록 시점에 이미 알 수 있는 정보로 이후의 위험을 설명하며, 기능이 그 뒤에 어떻게 변했는지는 쓰지 않는다. 기준 점수 하나에도 측정오차가 있으므로 $\beta_y$는 III.4절과 같은 이유로 0 쪽으로 편향될 수 있다.

### 4. 모형 B: 시간가변 공변량

모형 B는 기능 점수를 시간에 따라 바뀌는 공변량으로 넣는다. 시각 $t$ 직전의 마지막 방문을 $j(t) = \max\{j : t_{ij} \le t\}$라 하고, 그 방문의 값을 다음 방문까지 유지한다.

$$h_i(t) = h_0(t) \exp\big(w_i^\top \gamma + \beta_y \, y_{i, j(t)}/10\big)$$

추정은 계수 과정 표현(Andersen and Gill, 1982)으로 한다. 무릎마다 방문 사이 구간 $(t_{ij}, t_{i,j+1}]$을 한 줄로 만들고, 마지막 구간의 끝은 $T_i$로 자르며, 사건은 마지막 구간에만 표시한다. 분석 문서의 `tv` 표가 이 형식이고, R에서는 `Surv(start, stop, event)`로 적는다.

이 모형에는 문제가 셋 있다. 첫째, $y_{ij}$는 참 기능 수준 $m_i(t_{ij})$에 측정오차 $\varepsilon_{ij}$가 더해진 값이다. 오차가 섞인 공변량을 그대로 넣으면 계수가 0 쪽으로 편향된다(Prentice, 1982). 둘째, 값을 다음 방문까지 유지하는 것은 방문 사이의 변화를 무시하는 것이다. OAI의 방문 간격은 1\~2년이어서 이 근사가 거칠다. 셋째, 기능 점수는 무릎이 살아 있어야, 곧 치환술 전이어야 측정되는 내부 공변량이다(Kalbfleisch and Prentice, 2002). 내부 공변량에서는 위험함수와 생존함수 사이의 관계 $S_i(t) = \exp\{-\int_0^t h_i(s)\, ds\}$가 공변량 경로를 조건으로 해서는 성립하지 않으므로, 모형 B로는 개인별 생존 확률을 예측할 수 없다.

### 5. 모형 C: 2단계 추정

모형 C는 먼저 기능 점수에 선형 혼합모형(Laird and Ware, 1982)을 적합한다.

$$y_{ij} = (\beta_0 + b_{0i}) + (\beta_1 + b_{1i}) t_{ij} + \varepsilon_{ij}, \qquad b_i = (b_{0i}, b_{1i})^\top \sim N(0, D), \quad \varepsilon_{ij} \sim N(0, \sigma^2)$$

무작위 효과의 예측값은 최량 선형 불편 예측량(BLUP)이다. 무릎 $i$의 관측 벡터를 $y_i$, 고정효과와 무작위 효과의 설계행렬을 $X_i$, $Z_i$, 주변 공분산을 $V_i = Z_i D Z_i^\top + \sigma^2 I$라 하면 다음과 같다.

$$\hat b_i = \hat D Z_i^\top \hat V_i^{-1} (y_i - X_i \hat\beta)$$

관측이 적거나 잡음이 큰 무릎일수록 $\hat b_i$는 0 쪽으로 수축된다. 무릎별 연간 기울기 $\hat\beta_1 + \hat b_{1i}$을 얻으면, 이것을 모형 A에 공변량으로 더해 Cox 모형을 적합한다. 분석은 이 혼합모형을 `lme4`로 적합했다(Bates et al., 2015).

이 2단계 추정은 모형 B의 측정오차를 줄이지만 다른 문제를 만든다. 기울기를 추적 전체의 관측으로 추정하므로 시각 $t$의 위험을 설명하는 데 $t$ 이후의 정보가 섞인다. 첫 단계의 추정 불확실성이 둘째 단계로 넘어가지 않아 표준오차가 작게 나온다. 혼합모형은 추적 중단이 관측된 값에만 의존한다고 가정하는데, 치환술로 인한 중단은 관측되지 않은 기능 악화와 관련되므로 이 가정이 깨진다. 2단계 접근의 편향은 Ye et al. (2008)과 Sweeting and Thompson (2011)이 모의실험으로 보였다.

---

## IV. 결합모형

### 1. 구조

결합모형은 종단 하위모형과 생존 하위모형을 공통의 무작위 효과로 묶는다(Wulfsohn and Tsiatis, 1997; Henderson et al., 2000; Tsiatis and Davidian, 2004; Rizopoulos, 2012). 종단 하위모형은 관측 점수를 참 궤적 $m_i(t)$와 측정오차로 분해한다.

$$y_i(t) = m_i(t) + \varepsilon_i(t), \qquad m_i(t) = x_i(t)^\top \beta + z_i(t)^\top b_i, \qquad b_i \sim N(0, D), \quad \varepsilon_i(t) \sim N(0, \sigma^2)$$

첫 분석에서 고정효과는 시간, KL 등급, 그 상호작용이고 무작위 효과는 무릎별 절편과 기울기다. KL 등급 $k$의 무릎에 대해 적으면 다음과 같다.

$$m_i(t) = \beta_0 + \beta_{0k} + (\beta_1 + \beta_{1k}) t + b_{0i} + b_{1i} t, \qquad \beta_{00} = \beta_{10} = 0$$

생존 하위모형은 시각 $t$의 위험을 그 시각의 참 궤적 값과 변화율에 연결한다.

$$h_i(t \mid \mathcal M_i(t), w_i) = h_0(t) \exp\big\{ w_i^\top \gamma + \alpha_1 m_i(t) + \alpha_2 m_i'(t) \big\}$$

$\mathcal M_i(t) = \{m_i(s) : 0 \le s < t\}$는 시각 $t$까지의 참 궤적이다. $\exp(\alpha_1)$은 같은 시각에 참 기능 점수가 1점 나쁜 데 따른 위험비이고, $\exp(\alpha_2)$는 연간 악화 속도가 1점 빠른 데 따른 위험비다. JMbayes2에서는 이 두 항을 `functional_forms = ~ value(womac_func) + slope(womac_func)`로 지정한다. 궤적이 시간의 1차식이므로 변화율은 $m_i'(t) = \beta_1 + \beta_{1k} + b_{1i}$로, 무릎마다 시간에 따라 변하지 않는 값이다.

모형 B와 비교하면 차이는 공변량 자리에 관측값 $y_{i,j(t)}$ 대신 참 값 $m_i(t)$가 들어간다는 데 있다. $m_i(t)$는 측정오차가 없고 방문 사이에서도 연속적으로 정의되므로 III.4절의 첫째와 둘째 문제가 해결된다. 셋째 문제도 해결된다. 생존 하위모형이 무작위 효과를 조건으로 정의되므로 $S_i(t \mid b_i) = \exp\{-\int_0^t h_i(s \mid b_i)\, ds\}$가 성립하고, 이 관계가 V.4절의 동적 예측을 가능하게 한다.

### 2. 우도

두 하위모형은 무작위 효과가 주어지면 조건부 독립이라고 가정한다. 무릎 $i$의 우도 기여는 다음과 같다.

$$p(y_i, T_i, \delta_i) = \int p(y_i \mid b_i)\, h_i(T_i \mid b_i)^{\delta_i} \exp\Big\{ -\int_0^{T_i} h_i(s \mid b_i)\, ds \Big\}\, p(b_i)\, db_i$$

$p(y_i \mid b_i)$는 독립 정규 밀도의 곱이다. 치환술로 인한 추적 중단은 공유된 $b_i$를 통해 종단 하위모형과 연결되므로, 중단이 관측되지 않은 궤적에 의존하는 경우도 이 우도 안에서 다루어진다. III.5절의 2단계 추정이 이 점을 다루지 못한 것과 다르다.

누적 위험의 적분은 닫힌 꼴이 없어 수치적분으로 계산한다. JMbayes2(Rizopoulos et al., 2026)는 기본값으로 15점 Gauss-Kronrod 구적을 쓰고, 로그 기준 위험을 B-스플라인으로 전개한다. 설치된 판(0.6-5)의 기본 설정을 확인한 결과, 스플라인은 2차이고 구간은 9개이며, 계수에 2차 차분 벌점에 해당하는 평활 사전분포를 둔다. 따라서 Cox 모형과 마찬가지로 기준 위험의 모양을 미리 정하지 않는다.

### 3. KL 등급을 종단 하위모형에 넣는 이유

KL 등급을 종단 하위모형에서 빼면 $m_i(t) = \beta_0 + \beta_1 t + b_{0i} + b_{1i} t$가 되고, 등급에 따른 수준 차이는 무작위 절편이 떠안는다. 등급 $k$의 평균 수준 차이를 $\mu_k$라 하면 $b_{0i} \approx \mu_k + \tilde b_{0i}$꼴이 되어, 생존 하위모형의 선형 예측식에서 등급의 효과가 $\gamma_k$와 $\alpha_1 \mu_k$ 두 자리에 동시에 들어간다. 사후분포는 $\gamma_k + \alpha_1 \mu_k$가 일정한 방향으로 길게 늘어난 능선이 되고, 마르코프 연쇄는 이 능선을 따라 느리게 움직인다. 첫 적합에서 KL 계수의 $\hat R$이 가장 나빴던 것이 이 때문이다. 등급을 종단 하위모형의 고정효과로 넣으면 무작위 효과는 같은 등급 안의 편차가 되어 $\gamma_k$와 분리되고, 같은 반복 수에서 KL 계수의 $\hat R$은 기준 1.05 아래로 내려갔다.

### 4. 베이즈 추정과 수렴 진단

JMbayes2는 사후분포를 마르코프 연쇄 몬테카를로(MCMC)로 추출한다. 분석은 연쇄 3개를 돌리고 각 연쇄의 앞부분을 버린 뒤 나머지로 사후 평균과 95% 신용구간을 낸다. 위험비와 그 구간은 생존 하위모형 계수의 사후 표본에 지수를 취해 얻는다.

수렴은 잠재 척도 축소 인자 $\hat R$로 진단한다(Gelman and Rubin, 1992). 연쇄 $M$개가 각각 표본 $N$개를 가질 때, 연쇄 안 분산의 평균을 $W$, 연쇄 평균들 사이의 분산에 $N$을 곱한 값을 $B$라 하면 다음과 같다.

$$\hat R = \sqrt{\frac{\frac{N-1}{N} W + \frac{1}{N} B}{W}}$$

연쇄들이 같은 분포에 도달했으면 $\hat R$은 1에 가깝다. 분석 문서는 1.05를 기준으로 삼았다. Vehtari et al. (2021)은 순위 정규화를 적용한 개선된 $\hat R$과 더 엄격한 기준 1.01을 권고한다. 첫 분석에서 가장 늦게 섞인 것은 $\alpha_1$과 $\alpha_2$였다. 두 항은 모두 $b_{1i}$를 통해 같은 무작위 효과에 의존하므로 사후 상관이 높다. 반복 수를 늘려 기준을 만족시켰고, 그 반복 수를 `R/fit-jm.R`의 기본값으로 두었다.

### 5. 모의 데이터에서의 점검

원데이터를 받기 전까지 분석은 `R/oai.R`의 `make_synthetic`이 만든 모의 데이터로 돈다. 모의 데이터의 참 궤적은 KL 등급 $k$, BMI, 무릎별 무작위 효과로 정해지는 선형식을 $[0, 68]$로 자른 것이다.

$$m_i^*(t) = \min\Big\{68, \max\big(0,\; -4 + 5.5k + 0.25(\mathrm{BMI}_i - 28.5) + u_{0i} + (0.2 + 0.45k + u_{1i})\, t\big)\Big\}$$

$u_{0i} \sim N(0, 7^2)$, $u_{1i} \sim N(0, 1.2^2)$이고, 관측값은 여기에 $N(0, 4^2)$ 오차를 더해 반올림하고 $[0, 68]$로 자른 것이다. 치환술은 월 단위 이산 위험으로 생성했다.

$$h_i^*(t) = \exp\{-9.6 + 0.07\, m_i^*(t) + 0.35k + 0.02(\mathrm{age}_i - 61)\}$$

따라서 결합모형의 참값은 $\alpha_1 = 0.07$(위험비 $e^{0.07} \approx 1.073$)이고 $\alpha_2 = 0$이다. 참 궤적이 0과 68에서 잘려 있고 관측값도 잘려 있으므로, 정규 오차를 가정한 종단 하위모형은 이 데이터를 조금 잘못 지정한다. 점검의 기준은 둘이다. 결합모형의 $\alpha_1$과 $\alpha_2$ 신용구간이 참값을 담는지, 그리고 모형 B의 시간가변 계수가 III.4절의 측정오차 때문에 참값보다 0 쪽에 있는지다. 수치는 분석 문서 6절과 7절의 표에 있다. 모의 데이터의 결과는 모형과 코드의 점검일 뿐이며 OAI에 대한 결론이 아니다.

---

## V. 한계와 다음 단계

### 1. 바닥 효과

II.1절의 검열 정규 모형을 종단 하위모형에 넣으면, 0으로 관측된 점수는 잠재 점수가 0 이하라는 정보로 우도에 들어간다. 이렇게 하면 기능이 좋은 무릎의 궤적과 그 변화 속도를 바닥 효과로 인한 편향 없이 추정할 수 있다. JMbayes2가 지원하는 분포족에 이 모형이 있는지 먼저 확인하고, 없거나 맞지 않으면 Stan으로 직접 작성한다. 0에서의 질량이 크고 0이 아닌 값의 분포가 정규에서 멀면, 0 여부와 양의 값을 따로 모형화하는 두 부분 모형도 대안이다.

### 2. 한 사람의 두 무릎

III.2절의 군집 강건 분산은 Cox 모형의 표준오차만 보정하며, 결합모형은 두 무릎을 독립으로 다룬다. 사람 단위의 무작위 효과 $v_{p}$를 두 하위모형에 함께 넣으면 의존을 모형 안에서 다룰 수 있다. 종단 쪽은 $m_i(t)$에 $v_{p(i)}$를 더하고, 생존 쪽은 공유 취약성 $\exp(\omega v_{p(i)})$를 곱하는 꼴이다. 더 간단한 민감도 분석은 한 사람에서 무릎 하나만 고르는 것이다. 예를 들어 기준 KL 등급이 높은 쪽을 고르고, 결과가 바뀌는지 본다.

### 3. 사망과 경쟁 위험

사망을 중도절단으로 두면 사망한 사람이 계속 살았다면 치환술을 받았을 것처럼 다루게 된다. 원인별 위험 $h_1(t)$(치환술)과 $h_2(t)$(사망)를 따로 두면 치환술의 누적 발생률은 다음과 같다.

$$F_1(t) = \int_0^t S(s^-)\, h_1(s)\, ds, \qquad S(t) = \exp\Big\{ -\int_0^t \big(h_1(s) + h_2(s)\big) ds \Big\}$$

이 값은 II.3절의 $1 - \hat S(t)$보다 작다. 공변량 효과를 누적 발생률에 직접 연결하려면 Fine and Gray (1999)의 부분분포 위험모형을 쓴다. 결합모형에서는 생존 하위모형을 원인별 위험 두 개로 확장한다. 먼저 `Outcomes99`에 사망 일자 변수가 있는지 확인해야 한다.

### 4. 동적 예측

결합모형의 쓸모가 가장 크게 드러나는 것은 예측이다. 시각 $t$까지 치환술 없이 지냈고 그때까지의 점수 $\mathcal Y_i(t) = \{y_{ij} : t_{ij} \le t\}$가 관측된 무릎에 대해, $u > t$까지 치환술 없이 지낼 확률은 다음과 같다(Rizopoulos, 2011).

$$\pi_i(u \mid t) = P\big(T_i^* > u \mid T_i^* > t, \mathcal Y_i(t)\big) = \int \frac{S_i(u \mid b_i)}{S_i(t \mid b_i)}\, p\big(b_i \mid T_i^* > t, \mathcal Y_i(t)\big)\, db_i$$

새 방문에서 점수가 하나 더 관측되면 $b_i$의 사후분포가 갱신되고 예측도 함께 갱신된다. 예를 들어 앞으로 2년 안의 치환술 확률을 방문마다 다시 계산하면 물리치료사가 운동치료의 강도나 의뢰 시점을 정할 때 참고할 수 있는 숫자가 된다. 예측의 판별력은 시간의존 ROC 곡선과 AUC로 평가하고(Heagerty et al., 2000), 보정은 예측 확률과 관측 발생률을 구간별로 비교해 평가한다.

---

## 참고문헌

Andersen, P. K., and Gill, R. D. (1982). "Cox's regression model for counting processes: A large sample study." 『The Annals of Statistics』 10(4), 1100–1120.

Bates, D., Mächler, M., Bolker, B., and Walker, S. (2015). "Fitting linear mixed-effects models using lme4." 『Journal of Statistical Software』 67(1), 1–48.

Bellamy, N., Buchanan, W. W., Goldsmith, C. H., Campbell, J., and Stitt, L. W. (1988). "Validation study of WOMAC: A health status instrument for measuring clinically important patient relevant outcomes to antirheumatic drug therapy in patients with osteoarthritis of the hip or knee." 『The Journal of Rheumatology』 15(12), 1833–1840.

Cox, D. R. (1972). "Regression models and life-tables." 『Journal of the Royal Statistical Society: Series B』 34(2), 187–220.

Efron, B. (1977). "The efficiency of Cox's likelihood function for censored data." 『Journal of the American Statistical Association』 72(359), 557–565.

Fine, J. P., and Gray, R. J. (1999). "A proportional hazards model for the subdistribution of a competing risk." 『Journal of the American Statistical Association』 94(446), 496–509.

Gelman, A., and Rubin, D. B. (1992). "Inference from iterative simulation using multiple sequences." 『Statistical Science』 7(4), 457–472.

Heagerty, P. J., Lumley, T., and Pepe, M. S. (2000). "Time-dependent ROC curves for censored survival data and a diagnostic marker." 『Biometrics』 56(2), 337–344.

Henderson, R., Diggle, P., and Dobson, A. (2000). "Joint modelling of longitudinal measurements and event time data." 『Biostatistics』 1(4), 465–480.

Kalbfleisch, J. D., and Prentice, R. L. (2002). 『The Statistical Analysis of Failure Time Data』 (2nd ed.). Wiley.

Kaplan, E. L., and Meier, P. (1958). "Nonparametric estimation from incomplete observations." 『Journal of the American Statistical Association』 53(282), 457–481.

Kellgren, J. H., and Lawrence, J. S. (1957). "Radiological assessment of osteo-arthrosis." 『Annals of the Rheumatic Diseases』 16(4), 494–502.

Laird, N. M., and Ware, J. H. (1982). "Random-effects models for longitudinal data." 『Biometrics』 38(4), 963–974.

Lin, D. Y., and Wei, L. J. (1989). "The robust inference for the Cox proportional hazards model." 『Journal of the American Statistical Association』 84(408), 1074–1078.

Nevitt, M. C., Felson, D. T., and Lester, G. (2006). 『The Osteoarthritis Initiative: Protocol for the Cohort Study』 (Version 1.1). National Institute of Arthritis and Musculoskeletal and Skin Diseases.

Prentice, R. L. (1982). "Covariate measurement errors and parameter estimation in a failure time regression model." 『Biometrika』 69(2), 331–342.

Putter, H., Fiocco, M., and Geskus, R. B. (2007). "Tutorial in biostatistics: Competing risks and multi-state models." 『Statistics in Medicine』 26(11), 2389–2430.

Rizopoulos, D. (2011). "Dynamic predictions and prospective accuracy in joint models for longitudinal and time-to-event data." 『Biometrics』 67(3), 819–829.

Rizopoulos, D. (2012). 『Joint Models for Longitudinal and Time-to-Event Data: With Applications in R』. Chapman and Hall/CRC.

Rizopoulos, D. (2016). "The R package JMbayes for fitting joint models for longitudinal and time-to-event data using MCMC." 『Journal of Statistical Software』 72(7), 1–45.

Rizopoulos, D., Miranda-Afonso, P., and Papageorgiou, G. (2026). 『JMbayes2: Extended Joint Models for Longitudinal and Time-to-Event Data』. R package version 0.6-5. https://doi.org/10.32614/CRAN.package.JMbayes2

Sweeting, M. J., and Thompson, S. G. (2011). "Joint modelling of longitudinal and time-to-event data with application to predicting abdominal aortic aneurysm growth and rupture." 『Biometrical Journal』 53(5), 750–763.

Therneau, T. M., and Grambsch, P. M. (2000). 『Modeling Survival Data: Extending the Cox Model』. Springer.

Tobin, J. (1958). "Estimation of relationships for limited dependent variables." 『Econometrica』 26(1), 24–36.

Tsiatis, A. A., and Davidian, M. (2004). "Joint modeling of longitudinal and time-to-event data: An overview." 『Statistica Sinica』 14(3), 809–834.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Bürkner, P.-C. (2021). "Rank-normalization, folding, and localization: An improved $\widehat{R}$ for assessing convergence of MCMC." 『Bayesian Analysis』 16(2), 667–718.

Wulfsohn, M. S., and Tsiatis, A. A. (1997). "A joint model for survival and longitudinal data measured with error." 『Biometrics』 53(1), 330–339.

Ye, W., Lin, X., and Taylor, J. M. G. (2008). "Semiparametric modeling of longitudinal measurements and time-to-event data: A two-stage regression calibration approach." 『Biometrics』 64(4), 1238–1246.

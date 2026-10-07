# Ordering 메커니즘 실험 결과

원본 캠페인: `/Users/albertlee0125/Repos/qhat/codex_mechanism_experiment/runs/mechanism-t06qcgac`

## 검증 상태

- 완료: 129개 설정; baseline 108개, 통제 21개.
- 정확한 시간별 오차 분해: 63개.
- 기존 물리적 오차 anchor 일치: 9/9.
- 최대 상태벡터 재구성 잔차: 9.857e-14.
- 최대 복소 overlap 재구성 잔차: 9.842e-14.
- 최대 정규화 drift: 3.530e-12.
- exact 상태의 최대 spin-sector 누출: 1.461e-34.
- 분석 바닥 E<=1e-12: 0개. 바닥 이하 값은 성능비/수렴률 주장에 쓰지 않음.

E=1-|정규화 overlap|, I=1-|정규화 overlap|². 두 지표를 혼용하지 않음.
exact와 Trotter는 동일한 identity-free Pauli Hamiltonian을 사용한다.
도달 가능 공간은 각 Pauli 항에 대해 닫혀 있고 입자수 누출 상태를 제거하지 않는다.

## 기준 결과: T=1, r=100

|분자|순서|E|I 중 구역 밖 비율|국소 오차 총량 D|누적 계수 C=I/D|
|---|---|---:|---:|---:|---:|
|f2|fermionic_signed|1.058153e-07|0.000%|2.492994e-09|84.8901|
|f2|jw_signed|4.345745e-08|99.829%|4.263807e-08|2.03843|
|f2|jw_magnitude|1.059283e-07|0.169%|2.651991e-09|79.8859|
|nh3|fermionic_signed|9.631484e-10|0.000%|2.111384e-11|91.2338|
|nh3|jw_signed|1.039576e-09|0.970%|3.210019e-11|64.7707|
|nh3|jw_magnitude|9.948193e-10|0.000%|2.115362e-11|94.0566|
|li2|fermionic_signed|2.855298e-10|0.000%|4.372809e-11|13.0593|
|li2|jw_signed|1.142612e-08|96.494%|1.856264e-08|1.23109|
|li2|jw_magnitude|3.727978e-10|0.000%|3.917423e-11|19.0328|

## F2 순위 역전의 정량적 분해

- JW signed / Fermionic signed의 초기 BCH norm 비: 4.13579.
- JW signed / Fermionic signed의 국소 projected power D 비: 17.1032.
- JW signed / Fermionic signed의 누적 계수 C 비: 0.0240126.
- 두 비의 곱 = 최종 infidelity 비: 0.410691.

D와 C의 곱 분해는 정확한 항등식이지 새로운 예측 이론의 증명은 아니다.
C가 작아도 C>1이면 교차항의 순효과는 보강이다. 이를 음의 순상쇄라고 부르면 안 된다.
전파된 국소 오차의 방향/위상 누적을 무시한 초기-HF norm만으로 최종 순위를 단정할 수 없다.

### F2 시간 변화 (r=100)

|T|순서|I|D|C|교차항 I-D|
|---:|---|---:|---:|---:|---:|
|0.25|fermionic_signed|9.63979e-10|9.73757e-12|98.996|9.54242e-10|
|0.25|jw_signed|1.21012e-08|1.66557e-10|72.655|1.19346e-08|
|0.25|jw_magnitude|1.00987e-09|1.03586e-11|97.491|9.99511e-10|
|0.5|fermionic_signed|1.49623e-08|1.55804e-10|96.033|1.48065e-08|
|0.5|jw_signed|6.98853e-08|2.66491e-09|26.224|6.72204e-08|
|0.5|jw_magnitude|1.52217e-08|1.65741e-10|91.84|1.50559e-08|
|1|fermionic_signed|2.11631e-07|2.49299e-09|84.89|2.09138e-07|
|1|jw_signed|8.69149e-08|4.26381e-08|2.0384|4.42768e-08|
|1|jw_magnitude|2.11857e-07|2.65199e-09|79.886|2.09205e-07|

## 모든 grid 결과와 통제 실험

아래 winner 횟수는 동일 분자의 여러 설정을 세는 기술 통계이며 독립 표본 수가 아니다.
- f2: {'fermionic_signed': 8, 'jw_signed': 4}
- nh3: {'fermionic_signed': 12}
- li2: {'jw_magnitude': 8, 'fermionic_signed': 4}

|분자|통제 순서|Fermionic 대비 E 비|누적 계수 C|
|---|---|---:|---:|
|f2|round_robin|0.332749|5.97727|
|f2|within_blocks_20260917|1|84.8901|
|f2|blocks_random_1101|1.70987|84.6816|
|f2|blocks_random_1102|13.9375|84.9217|
|f2|blocks_random_1103|13.7338|84.9067|
|f2|blocks_random_1104|4.07646|84.9039|
|f2|blocks_random_1105|20.1619|84.9259|
|nh3|round_robin|6.63387|94.0533|
|nh3|within_blocks_20260917|1|91.2338|
|nh3|blocks_random_1101|185.282|97.0771|
|nh3|blocks_random_1102|151.907|97.1227|
|nh3|blocks_random_1103|219.412|97.3433|
|nh3|blocks_random_1104|206.269|97.1422|
|nh3|blocks_random_1105|397.447|97.1273|
|li2|round_robin|12.4313|1.94291|
|li2|within_blocks_20260917|1|13.0593|
|li2|blocks_random_1101|693.562|62.8894|
|li2|blocks_random_1102|201.801|63.1783|
|li2|blocks_random_1103|44.3737|87.4559|
|li2|blocks_random_1104|528.472|62.3271|
|li2|blocks_random_1105|22.9948|64.1743|

## 수렴률

분석 바닥 이상인 가장 큰 두 step 수로 구한 log(E)/log(r) 기울기. 일반적인 first-order 상태 오차의 제곱 지표는 -2에 접근할 수 있으나 이를 강제하지 않았다.

|분자|T|순서|r 구간|기울기|
|---|---:|---|---|---:|
|f2|0.25|fermionic_signed|100–200|-2.0000|
|f2|0.25|jw_signed|100–200|-2.0000|
|f2|0.25|jw_magnitude|100–200|-2.0000|
|f2|0.5|fermionic_signed|100–200|-2.0000|
|f2|0.5|jw_signed|100–200|-2.0001|
|f2|0.5|jw_magnitude|100–200|-2.0000|
|f2|1|fermionic_signed|100–200|-2.0000|
|f2|1|jw_signed|100–200|-2.0006|
|f2|1|jw_magnitude|100–200|-2.0000|
|nh3|0.25|fermionic_signed|100–200|-2.0000|
|nh3|0.25|jw_signed|100–200|-2.0000|
|nh3|0.25|jw_magnitude|100–200|-2.0000|
|nh3|0.5|fermionic_signed|100–200|-2.0000|
|nh3|0.5|jw_signed|100–200|-2.0000|
|nh3|0.5|jw_magnitude|100–200|-2.0000|
|nh3|1|fermionic_signed|100–200|-2.0000|
|nh3|1|jw_signed|100–200|-2.0001|
|nh3|1|jw_magnitude|100–200|-1.9999|
|li2|0.25|fermionic_signed|100–200|-2.0000|
|li2|0.25|jw_signed|100–200|-2.0001|
|li2|0.25|jw_magnitude|100–200|-2.0000|
|li2|0.5|fermionic_signed|100–200|-2.0000|
|li2|0.5|jw_signed|100–200|-2.0004|
|li2|0.5|jw_magnitude|100–200|-1.9999|
|li2|1|fermionic_signed|100–200|-2.0000|
|li2|1|jw_signed|100–200|-2.0015|
|li2|1|jw_magnitude|100–200|-1.9999|

## 이 실험으로 해결되지 않은 것

- 과거 tensor/코드 불일치의 원인은 이 실험으로 판별하지 않았다.
- 세 사례는 이미 결과를 본 사례다. 외부·전향 검증이 아니며 일반적 우월성을 입증하지 않는다.
- 시간별 분해는 설명 도구다. exact trajectory가 필요한 이 진단을 효율적인 ordering 선택 알고리즘으로 주장할 수 없다.
- 5개 block random seed와 round robin은 구조적 통제이나 서로 같은 크기의 perturbation은 아니다.
- 1-body/2-body의 독립 효과, orbital-gauge 견고성, 동일 정확도 회로 비용은 검증하지 않았다.
- 작은 누출 값의 소수점 수준 차이보다 전체 오차와 수치 잔차를 우선 해석한다.
- conditional sector metric은 사후 선택의 비용을 무시한 성능 주장에 사용하지 않는다.

## 산출물

- `baseline_grid.png`: 시간·step 변화와 세 ordering 비교.
- `f2_decomposition.png`: F2의 국소 오차 크기와 누적 계수.
- `sector_controls.png`: 구조적 통제에서 구역 안/밖 오차.
- `analysis.json`: 표·검증 통계의 원본.
- 상위 캠페인의 `inputs/`, `manifest.json`, `*.setup.json`, `*.trace.json`, `results.json`: 재현 입력과 수치 기록.

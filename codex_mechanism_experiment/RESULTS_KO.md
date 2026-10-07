# 실험 완료: 초기 BCH norm과 최종 정확도의 불일치를 설명하기

## 결론

F2의 ordering 순위 역전은 scalar offset 불일치나 반복 불안정성으로
사라지는 현상이 아니었다. 동일한 Hamiltonian에서 국소 오차의 크기와
시간에 따른 누적을 분리하면 역전을 정량적으로 설명할 수 있었고,
시간전파한 leading-BCH 모델도 맞춤 계수 없이 실제 오차를 재현했다.

이는 F2에 대한 설명·검증의 연결을 확보한 결과다. 세 분자의 일반적
우월성, 새로운 오차 간섭 이론, 논문 신규성, 실용적 회로 비용 절감까지
입증한 결과는 아니다.

## 수행 범위와 보존

- 코드 기준: QHAT L-sweep, 7cd519d520662260f13e16e213eacc8f74c91c7f.
- 고정 입력: 이전 재현성 캠페인의 F2, NH3, Li2 tensor. 각 캠페인에 실제
  입력을 복사하고 Hamiltonian/parent/물리적 순서의 해시를 대조했다.
- 본실험: 108개 baseline 설정 + 21개 구조적 통제 = 129개 완료.
- 기본 grid: 세 분자 × 세 ordering × T={0.25,0.5,1} × r={25,50,100,200}.
- 정확한 시간별 오차 분해: 63개. 나머지는 최종 오차와 수렴률을 측정.
- 별도 탐색적 실험: F2의 T=0.75,r=100 세 ordering 추가. 예측을 먼저
  파일에 저장한 뒤 실제 product-formula 상태를 계산했다.
- 수치 검증용 4-qubit smoke 6개와 별도 프로세스 F2 anchor 3개는 위의
  새로운 물리 설정 수에 포함하지 않는다.
- HF, first-order, coefficient cutoff 1e-12. 기존 연구 코드·과거 결과는
  수정하지 않았다. commit/push는 하지 않았다.

## 1. 시간에 따라 실제로 우열이 바뀐다

아래 비는 r=100에서 E(JW signed)/E(Fermionic signed)이다.
E=1-|정규화 overlap|이며, 표준 infidelity I=1-|정규화 overlap|²와 다르다.

|F2의 T|오차 비|더 정확한 순서|
|---:|---:|---|
|0.25|12.553|Fermionic signed|
|0.5|4.671|Fermionic signed|
|0.75 (추가 검증)|약 0.774|JW signed|
|1|0.411|JW signed|

원래 grid의 시간별 우열은 r=25,50,100,200에서도 유지됐다. Li2에서도
T=0.25,0.5에서는 JW magnitude, T=1에서는 Fermionic signed가 우세했다.
NH3에서는 조사한 grid 모두 Fermionic signed가 우세했다. 이 36개의
분자·시간·step 조합은 세 분자에서 얻은 반복 조건이지 독립 분자 36개가 아니다.

분석 바닥 이상인 가장 큰 두 r로 측정한 E의 log-log 수렴 기울기는
-2.00154에서 -1.99986이었다. 따라서 조사 구간의 순위 역전은 first-order
제곱 오차의 수렴 형태가 무너지는 현상으로 설명되지 않는다.

## 2. 국소 오차는 크지만 덜 누적된다

정규화한 exact 최종 상태를 phi, Q=I-|phi><phi|라 하자. 각 step에서 생긴
defect를 최종 시간으로 전파한 벡터 v_j에 대해

    I = ||sum_j Q v_j||² / ||psi_Trotter||²
      = D_local + X = D_local * C

여기서 D_local은 각 ||Q v_j||²의 합, X는 교차항, C는 누적 계수다.
이 항등식 자체는 새로운 이론이나 사전 예측 법칙이 아니다.

F2, r=100에서 JW signed/Fermionic signed의 D_local 비는 세 시간 모두
약 17.10이었다. 반면 C의 비는 다음처럼 달라졌다.

|T|D_local 비|C 비|최종 infidelity 비|
|---:|---:|---:|---:|
|0.25|17.105|0.7339|12.553|
|0.5|17.104|0.2731|4.671|
|1|17.103|0.02401|0.4107|

따라서 초기 BCH norm은 국소 오차의 큰 차이를 반영하지만, 유한 시간에
전파된 오차의 방향·위상 누적을 반영하지 못해 최종 순위를 놓친다.

중요한 표현상의 제한: T=1,r=100의 C는 Fermionic 약 84.89, JW signed
약 2.038로 둘 다 1보다 크다. 이 설정은 **음의 순상쇄**보다 **훨씬 약한
보강 누적**으로 기술하는 것이 정확하다. 반면 T=1,r=25의 JW signed는
C≈0.514로 음의 순교차항을 보였다. C는 step 분할에 의존하므로 다른 r를
비교할 때 C/r 또는 같은 r에서의 ordering 비를 같이 보아야 한다.

## 3. 설명을 넘어 leading-order 예측도 확인했다

F2에서 실제 순서의 한-step 오차 계수 B_pi를

    S_pi(dt) - exp(-i H dt) = dt² B_pi + O(dt³)

로 정의하고, 다음 선형 비균질 방정식을 계산했다.

    eta'(t) = -i H eta(t) + B_pi psi_exact(t),  eta(0)=0
    I_pred(T,r) = (T/r)² ||Q_T eta(T)||²

Hamiltonian과 실제 Pauli 순서에서 계산했으며, 측정된 Trotter 오차에
맞추는 회귀 계수나 보정 계수는 사용하지 않았다. 이는 표준적인 leading
error 전파 구성의 적용이며 새 정리를 주장하는 것이 아니다.

- 기존 grid와 새 시간의 총 39개 F2 예측 비교: 상대 차이 최대 0.954%.
- r=100의 모든 비교: 상대 차이 최대 0.0598%.
- r=200: 상대 차이 최대 0.0150%.
- 새 시간 T=0.75,r=100에서 실제 E는 Fermionic 3.59899e-8,
  JW signed 2.78728e-8, JW magnitude 3.60272e-8이었다.
- 새 시간에서도 세 순서의 순위를 맞췄고, infidelity 상대 차이는 최대
  0.0219%였다. 예측 파일은 이 세 실제 상태를 계산하기 전에 저장됐다.

이 보조 실험은 F2 결과를 본 뒤 설계한 **탐색적 검증**이다. 새 분자에
대한 전향 검증으로 부르면 안 되며, exact 진화가 필요한 이 진단을
저비용 ordering 선택 알고리즘으로 주장할 수도 없다.

## 4. Parent 구조의 역할과 한계

세 입력의 실제 parent-owned block 모두 내부 Pauli들이 서로 교환했다.
block Hamiltonian과 N, Nalpha, Nbeta의 commutator 계수 l1 잔차 최대값은
F2 1.03e-16, NH3 5.42e-19, Li2 1.44e-17이었다. 작은 계수를 자동 삭제하지
않는 방식으로 검증했다. fallback Pauli는 없었다.

이는 이번 입력에서 다음을 설명하는 충분조건이다.

- block 내부 순열: 교환하는 지수들의 재배치이므로 결과가 변하지 않음.
- 완전한 block 단위 재배치: 입자수·스핀 보존은 유지되지만 block 간
  순서 오차는 달라질 수 있음.

실제로 내부 순열은 세 사례 모두 기준 결과와 같았다. 5개 고정 random
seed의 block 재배치는 모두 기준보다 나빴지만, 이는 작은 표본의 기술적
통제 결과이지 일반 최적성의 증명은 아니다.

반면 block을 서로 끼워 넣는 round robin은 Fermionic 기준 대비:

|분자|E 비|효과|
|---|---:|---|
|F2|0.333|약 3.0배 더 작은 오차|
|NH3|6.634|오차 증가|
|Li2|12.431|오차 증가|

따라서 **parent 구조는 대칭성 보존을 설명하지만, 그 보존만으로 최종
정확도의 우월성을 보장하지 않는다.** F2의 JW signed는 전체 infidelity의
약 99.83%가 구역 밖 오차인데도, 구역 내부 오차가 훨씬 작아 전체에서는
Fermionic signed보다 정확했다. 조건부 구역 오차는 진단용이며 공짜
postselection 성능으로 해석하지 않았다.

## 검증과 남은 과제

- 단위·행렬·위상·예측 검증 테스트 13개 통과.
- 기존 9개 물리적 오차 anchor 모두 일치.
- 별도 프로세스 F2 anchor의 비교 18개 통과.
- 기존 전체 Hilbert 공간의 한-step 결과와 축약 계산의 잔차: 0.
- 시간별 defect 합의 최대 상태벡터 재구성 잔차: 9.86e-14.
- 재구성한 projected infidelity와 직접 계산의 최대 차이: 3.98e-19.
- 최대 상태 norm drift: 3.53e-12. 물리 지표는 정규화한 직교 잔차로
  계산해 1에 가까운 overlap의 뺄셈 및 norm drift 영향을 줄였다.
- 입력 tensor는 새 캠페인에 포함됐고, 새 폴더의 ignore 규칙은 이 입력을
  제외하지 않는다. 계산 캐시는 제외한다.

논문을 향한 다음 핵심 과제는 (1) 시간전파 오차 진단을 아직 분석하지 않은
분자에서 검증, (2) 기존 error-interference 및 symmetry-preserving
partitioning 연구와 구체적 기여를 대조하는 것이다. 원한다면 실용성 주장에
앞서 동일 정확도 비용도 측정해야 한다. 과거 자료의 계보 불일치는 여전히
별도 문제이며 이번 결과로 해결됐다고 하지 않는다.

## 실제 파일

- 본실험: `runs/mechanism-t06qcgac/`
- 최종 표·그림: `runs/mechanism-t06qcgac/analysis-dl3ym9tx/`
- 구조·예측 검증: `runs/mechanism-t06qcgac/exploratory-kdaizs_h/`
- 추가 전 원래 계획: `PROTOCOL.md`
- 탐색적 추가 계획: `EXPLORATORY_ADDENDUM.md`
- 테스트 기록: `VALIDATION_TESTS.txt`

앞선 `analysis-jr4926ap/`는 같은 수치의 첫 렌더링이다. 최종 그림은 위의
`analysis-dl3ym9tx/`를 사용한다. 과거 결과와 중간 산출물은 삭제하지 않았다.

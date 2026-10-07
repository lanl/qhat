# 선행연구 대비 신규성 감사

검토일: 2026-09-17. 아래 판단은 공개 원문의 관련 절을 확인한 범위의 연구 판단이며,
포괄적인 우선권 증명이나 게재 보장이 아니다. 새로운 수치 결과와 새로운 이론을 구별한다.
검색은 논문 제목, molecular Trotter ordering, fermionic/signed-coefficient ordering,
first occurrence, state-dependent error, error interference, orbital transformations로 수행했다.
검색에서 동일 문구가 없다는 사실은 신규성의 증거로 사용하지 않았다.

## 결론

**새 분자에 대한 예측 검증은 확보했다. 그러나 “fermionic parent ordering 자체”,
“BCH norm과 실제 오차가 다르다”, “시간적 오차 간섭이 중요하다”를 새 발견으로 내세울 수는 없다.**

기여 후보는 좁혀야 한다: 고정된 merged Pauli 목록에서 parent 기반 스케줄을 명확히 정의하고,
대칭 섹터 내부 오차·누출·시간적 누적을 분리하여 순위의 성공과 실패를 재현 가능한 방식으로 설명하는
진단/벤치마크 연구. 이 결합의 구체적 가치도 아래 문헌 대비 직접 비교로 입증해야 한다.

## 주장별 판정

| 주장 | 선행연구와 겹치는 부분 | 현재 판정 |
|---|---|---|
| 항 ordering으로 Trotter 오차를 줄인다 | Hastings 2014 §4; Tranter 2019 §§III–V | 알려진 문제·접근 |
| fermionic 구조에 맞춰 묶거나 순서를 정한다 | Hastings의 interleaved 항 순서; Martínez-Martínez의 fermionic fragment; Kronenberger 2026의 Hermitian fermionic 항 | 일반 원리는 신규 아님 |
| commuting block 내부 재배열은 오차를 바꾸지 않는다 | Tomesh의 그룹 내부 gate-cancellation 재배열 | 알려진 대수적 성질 |
| 대칭 보존을 이용한다 | Martínez-Martínez §2.5의 입자수·스핀 제한 | 알려진 분석 도구 |
| norm-bound/단순 descriptor가 실제 오차를 잘 못 설명할 수 있다 | Babbush 2015; Kronenberger 2026 | 일반적 관찰은 신규 아님 |
| 시간적으로 운반된 오차가 간섭한다 | Tran 2020; Layden 2022; Chen 2024 preprint/2026 PRL | 일반 현상·이론은 명백히 기존 연구 |
| global phase를 제거하고 상태 fidelity를 본다 | Yi–Crosson의 spectral analysis | 신규 metric 아님 |
| leading BCH를 시간 전파해 finite-time 오차를 예측한다 | effective Hamiltonian·섭동 분석의 표준 계열 | 이번 구현은 응용/검증, 새로운 일반이론 아님 |
| merged Pauli를 처음 등장하는 signed parent에 배정한다 | 조사한 관련 절에서 이 정확한 규칙의 동일 구현은 확인하지 못함 | 우선권 미확정; 규칙의 사소한 차이만으로 논문 기여가 충분하지 않음 |
| 대칭 보존과 fidelity 우월성이 다름을 새 분자로 정량 검증한다 | 기존 이론으로 예상 가능한 구분 | 응용·벤치마크 기여 후보; 보편적 신현상으로 표현 금지 |

## 반드시 대조해야 할 원문

### 1. Hastings et al., Improving Quantum Algorithms for Quantum Chemistry (2014)

§4는 HF 조건을 이용해 one-body hopping과 관련 interaction을 함께 처리하는 interleaved
순서를 제시한다. §3에는 fermionic 항에서 생긴 Pauli subterm 순서 조정도 있다.
따라서 “fermionic parent를 활용하는 ordering을 처음 생각했다”는 표현은 부적절하다.
우리 first-occurrence merged-coefficient 스케줄과 동일한 알고리즘이라는 뜻은 아니다.
[원문 §§3–4](https://arxiv.org/html/1403.1539v2)

### 2. Babbush et al., Chemical Basis of Trotter-Suzuki Errors in Quantum Chemistry Simulation (2015)

오차 연산자 norm과 실제 ground-state energy error의 큰 차이, 전자구조와 orbital filling의
영향을 이미 분석했다. 우리의 HF 동역학 infidelity와 목적함수는 다르지만,
“norm이 실제 물리 오차를 그대로 뜻하지 않는다”는 일반 명제는 기존 배경이다.
[원문](https://arxiv.org/html/1410.8159v2)

### 3. Tranter et al., Ordering of Trotterization: Impact on Errors in Quantum Simulation of Electronic Structure (2019)

44 Hamiltonian benchmark, magnitude ordering, graph-coloring 및 error-operator 기반 순서가
포함된다. 우리 세 가지 baseline 비교만으로 ordering benchmark 분야에 새 방법을 제시했다고
말하기는 어렵다. Published algorithm과 동일 목적함수·같은 입력에서의 직접 비교가 필요하다.
[원문 §§IV–V](https://arxiv.org/html/1912.07555v1)

### 4. Tran et al., Destructive Error Interference in Product-Formula Lattice Simulation (2020)

서로 다른 Trotter step의 destructive interference를 분석하고 기존 오차 상계보다 나은
scaling을 제시했다. 주요 설정은 lattice model이다. 분자에의 응용 차이는 남지만
“step 간 상쇄” 자체는 새 이론이 아니다.
[원문](https://arxiv.org/html/1912.11047v1)

### 5. Layden, First-Order Trotter Error from a Second-Order Perspective (2022; preprint 2021)

두 항 분해에서 first-order와 second-order 회로의 관계로 간섭과 시간 scaling을 설명한다.
우리 많은 Pauli 항의 상태별 실험과 적용 범위가 같지는 않다. 그러나 간섭 설명이나 시간별
오차 변화 그 자체를 신규성으로 삼아서는 안 된다.
[원문, 특히 Eqs. 6–14](https://arxiv.org/html/2107.08032v1)

### 6. Yi and Crosson, Spectral Analysis of Product Formulas for Quantum Simulation (2022; preprint 2021)

effective Hamiltonian의 eigenvalue/eigenvector 섭동 및 global phase와 fidelity error의
분리를 사용한다. 우리의 위상 불변 지표·leading perturbation 계산에 가까운 개념적 선행연구다.
초기상태 가정과 주요 응용(QPE/DAS)은 구별해야 한다.
[원문 §I, Eqs. 2–3](https://arxiv.org/html/2102.12655v1)

### 7. Tomesh et al., Optimized Quantum Program Execution Ordering to Mitigate Errors in Simulations of Quantum Systems (ICRC 2021; preprint 2022)

commuting groups를 만든 뒤 그룹 내부를 TSP로 재배열해 gate cancellation을 개선한다.
그룹 내 순서 불변성과 오류·회로 비용 공동 최적화를 이미 다룬다.
우리 within-bucket shuffle 결과는 이 원리를 확인하는 control이지 새로운 발견이 아니다.
[원문 §IV](https://arxiv.org/html/2203.12713v1)

### 8. Martínez-Martínez, Yen and Izmaylov, Assessment of various Hamiltonian partitionings… (Quantum 2023)

fermionic/qubit fragment 분해, 대칭 제한 norm, T-gate 비용을 비교한다.
특히 §2.5는 particle number 및 spin symmetry를 이용한다.
우리는 완전히 같은 merged Pauli rotations를 재배열하므로 fragment 분해를 바꾸는 비교와 다르다.
이 차이를 Methods에서 명시하고 두 방법을 같은 것으로 혼용하지 않아야 한다.
[원문 §§2.3–2.5](https://arxiv.org/html/2210.10189v2)

### 9. Chen et al., General Framework for Error Interference in Quantum Simulation (PRL 136, 200601, 2026-05-18)

2024 preprint 제목은 Error Interference in Quantum Simulation이다.
effective Hamiltonian, 간섭의 필요충분 조건, approximate interference 및 상태 의존 조건을 다룬다.
“처음으로 오차 간섭의 일반 틀을 제시했다”는 주장은 특히 이 연구와 충돌한다.
우리 finite-T에서 C=I/D<1이라는 진단은 그 논문의 asymptotic spectral-norm 정의와 동일하지 않다.
[출판본 정보](https://journals.aps.org/prl/abstract/10.1103/g2n5-qdxh),
[preprint 원문 §II, Theorem 1, Corollary 3](https://arxiv.org/html/2411.03255v2)

### 10. Kronenberger, Erakovic and Reiher, Trotter Error and Orbital Transformations in Quantum Phase Estimation (2026)

가장 가까운 최근 molecular 선행연구 중 하나다. Hermitian fermionic 항과 qubit 표현,
orbital 변화, magnitude/index ordering, perturbative energy estimate와 descriptor의 한계를
분석한다. §4.2.2는 orbital 효과와 ordering 효과의 혼입을 직접 통제한다.
우리 고정 Pauli 목록·HF finite-time leakage 분석은 다르지만, “fermionic vs Pauli + 설명 지표”만으로는
충분한 차별점이 아니다. orbital gauge/활성공간 민감성도 제출 전 확인해야 한다.
[원문 §§2.1, 2.3, 3, 4.2.2](https://arxiv.org/html/2602.18913v1),
[출판 DOI](https://doi.org/10.1080/00268976.2026.2681062)

### 11. Tate, Aktar and Eidenbenz, An Analysis of Commutation-Based Trotter Ordering Strategies on Heisenberg-Style Hamiltonians (2026 preprint)

저자 소속이 LANL이다. commutation-based ordering, random controls, group permutations,
first/second-order 및 fidelity 비교를 수행한다. 이 연구와 현재 작업의 내부적 관계는
확인하지 않았으므로 별개 연구라고 단정하지 않는다. 투고 전 지도·공동연구진과 경계를 정해야 한다.
우리 분자 parent 구조와 symmetry-channel 설명은 Heisenberg benchmark와 구분 가능한 축이다.
[원문 §§III–V](https://arxiv.org/html/2604.23138v1)

### 12. Bay-Smidt et al., Quantum simulation of nanographenes and Trotter error cancellation (2026 preprint)

분자계의 application-specific error와 에너지 차이의 상쇄를 정량화한다.
우리 HF 상태 infidelity 및 time-step 간섭과 관측량이 다르므로 직접 같은 현상이라고 하면 안 된다.
“분자 Trotter error cancellation 연구가 없다”는 주장도 피해야 한다.
[원문 §§IV–V](https://arxiv.org/html/2605.00745v1)

## 우리 실험에 맞는 수학적 위치

S_pi(dt)=exp(-iH dt)+dt^2 B_pi+O(dt^3)라 쓰면, fixed-T에서

eta(T) = integral_0^T exp[-iH(T-s)] B_pi exp(-iHs)|HF> ds,

I_pred = (T/r)^2 ||Q_T eta(T)||^2.

이것은 leading error를 운반한 1차 섭동식이며 이번 코드의 augmented ODE와 동치이다.
정확한 고전 시간발전을 필요로 하므로 큰 분자에서 저비용으로 ordering을 정하는 새 알고리즘이라는
주장은 하지 않는다. Q_T는 정확한 최종 상태 방향을 제거하는 projector다.

P가 정확히 보존되는 particle/spin sector이고 exact state가 그 안에 있으면
I = ||P Q_T psi_Trotter||^2 + ||(1-P)psi_Trotter||^2 (정규화된 상태).
따라서 leakage=0은 둘째 항만 없애며, 첫째 항을 최소화한다는 보장은 없다.
이 직교 분해도 표준 대수다. 기여는 이 구분이 ordering 판단에서 얼마나 필요한지 보여주는 증거다.

## 현재 데이터가 지지하는 논문 방향과 남은 관문

가능한 제목: *Symmetry preservation does not determine finite-time Trotter ordering accuracy
in molecular Hamiltonians*.

이는 작업용 주제이며 새 현상/일반 정리를 보증하는 제목이 아니다.
현재는 “일반적으로 좋은 fermionic ordering” 논문보다 “ordering 선택 지표의 한계와 channel-resolved
진단” 논문이 정직하다. 새 panel은 dynamic model의 오차값 정확도를 지지하지만, static BCH도
72/72 순위를 맞혔으므로 더 나은 ordering selector라는 신규성은 아직 입증되지 않았다.

제출 수준의 기여를 강화할 다음 관문은 두 가지다.

1. **동일 조건 선행 알고리즘 비교:** Tranter의 graph-group 전략과 fermionic magnitude 등
   강한 comparator를 원문대로 구현하고, state error와 실제 회로 비용을 구분해 비교한다.
   현재 random-block control은 이 comparator를 대체하지 않는다.
2. **진단의 실제 효용:** orbital gauge/degenerate subspace/활성공간 크기를 통제한 후,
   static 지표가 틀리는 조건을 사전 규칙으로 예측하는지 평가한다. 더 긴 시간 구간을 이미 본
   4분자에 추가하면 exploratory로 표시하고, 새 확인용 panel에서 다시 검증한다.
   단순히 좋은 사례를 늘리거나 결과를 본 뒤 시간 범위를 골라 “blind validation”이라 부르지 않는다.

이 관문 없이 현재 모델 정확도만을 “새 이론”으로 제출하는 것은 추천하지 않는다.
소규모지만 정직한 재현성·진단 benchmark로 원고를 구성할 가능성은 있으며, 적절한 학술적
기여 판단은 위 비교와 공동연구진 검토를 거쳐야 한다.

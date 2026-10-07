# 새 분자 검증 결과

결론: **새 분자로의 수치 예측 전이는 확인했다. Fermionic signed의 일반적 정확도 우월성은
확인되지 않았고, 오히려 다수의 반례를 얻었다. 신규성 검토 결과도 일반 ordering/간섭 이론은 기존 연구다.**

## 무엇을 실행했나

- HCN, CH2O, H2S, SiH4: 기존 저장소 tensor 목록의 14분자군에 없던 네 분자군.
- 각 이상화 기하구조와 전체 좌표를 1.25배 늘린 구조: 총 8 Hamiltonian.
- 6-31g/RHF, CAS(4e,5 spatial orbitals)=10 qubits. 아래는 이 잘린 활성공간의 결과다.
- HF, first-order, T={0.25,0.75,1}, r={50,100,200}, 세 ordering: 216조건.
- T=1,r=100에서 7개 구조 control×8입력: 56조건.
- **총 272조건을 두 독립 Python 프로세스에서 실행**했다. 544개 서로 다른 조건은 아니다.
- 모든 예측을 216개 target 계산 전에 저장했다. protocol은 분자 생성 전 작성했다.
  공개 preregistration은 아니고 local code/protocol/hash 동결이다.

기준 campaign: `runs/transfer-vtrgy3ds`.
독립 반복: `runs/transfer-hux1k985`.
최종 분석: `runs/transfer-vtrgy3ds/analysis-x8cr3me0`.

## 예측 검증

오차 지표 I=1-|overlap|^2, E=1-|overlap|. 아래 “상대 차이”는 I 예측의 상대 차이다.

| 항목 | 결과 |
|---|---:|
| r=100 최대 상대 차이 | 0.10150% |
| r=200 최대 상대 차이 | 0.05057% |
| 사전 기준 ≤5% 통과 | 143/143 eligible rows |
| 수치 floor로 제외한 baseline row | 1/216 |
| r=100 순위 일치 (actual gap>1%) | dynamic 72/72; static HF-BCH도 72/72 |
| r=100,200 상대 차이 중앙값 | dynamic 0.002302%; static 4.33748% |
| E의 r-scaling 기울기 범위 | -2.000630 ~ -1.998892 |

제외된 row는 SiH4 s=1.00 / JW signed / T=.25 / r=200, E=3.9124e-13이다.
사전에 정한 E<=1e-12 기준에 따라 상대 오차·순위 판정에서 제외했고 원자료는 보존했다.
72개 pairwise 비교는 4분자×2기하×3시간의 각 세 쌍 비교다. 독립 분자 72개로 해석하지 않는다.

**중요한 부정 결과:** 새 panel에서 dynamic predictor는 오차값을 더 정확히 예측했지만,
static HF-BCH보다 ordering을 더 잘 고른다는 근거는 얻지 못했다.
각 기하구조 안에서 세 시간 동안 최우수 ordering이 바뀌는 현상도 없었다.
앞선 F2의 시간별 순위 역전을 새 네 분자에서도 재현했다고 주장할 수 없다.

## 정확도 비교: T=1,r=100의 모든 사례

E=1-|overlap|, 낮을수록 좋다.

| 분자·좌표 scale | Fermionic signed | JW signed | JW magnitude | 최우수 |
|---|---:|---:|---:|---|
| HCN 1.00 | 3.35373e-7 | 8.03305e-8 | 2.94998e-7 | JW signed |
| HCN 1.25 | 1.00638e-7 | 4.08194e-7 | 1.10576e-7 | Fermionic signed |
| CH2O 1.00 | 1.87183e-7 | 1.57621e-7 | 4.81440e-7 | JW signed |
| CH2O 1.25 | 1.27784e-7 | 2.41721e-7 | 3.61926e-7 | Fermionic signed |
| H2S 1.00 | 4.31809e-8 | 4.67899e-9 | 4.26298e-8 | JW signed |
| H2S 1.25 | 1.24158e-7 | 1.94141e-8 | 8.41206e-8 | JW signed |
| SiH4 1.00 | 4.40692e-8 | 3.45201e-10 | 4.92689e-8 | JW signed |
| SiH4 1.25 | 6.66818e-8 | 5.21954e-9 | 6.80857e-8 | JW signed |

JW signed 6/8, fermionic signed 2/8. 동일 family의 두 기하구조는 독립 분자가 아니다.
이 비율은 사전 선정한 작은 panel의 기술 통계이며 모든 분자의 승률 추정이 아니다.

## 구조·간섭·control

- 8입력 모두 first-occurrence bucket 내부 noncommuting pair=0, fallback=0.
- bucket의 N_alpha,N_beta,N commutator coefficient-l1 상계 최대 6.94e-18.
- fermionic baseline spin-sector leakage 최대 5.04e-37: 수치적으로 0.
  JW signed baseline 최대 1.17e-7. 그럼에도 JW의 total I가 더 작을 수 있다.
- baseline의 r=100에서 C=I/D 범위는 83.20~99.77, C<1은 0개다.
  즉 이번 panel을 “순 destructive cancellation을 새 분자에서 확인했다”고 쓰면 안 된다.
- within-bucket shuffle의 E 비율은 1에서 최대 약 2e-12 상대 차이: commuting 구조와 일치.
- whole-bucket random 40개 중 13개는 기준 fermionic signed보다 작았다.
  이전 세 분자에서 random block이 모두 나빴던 현상은 일반화되지 않는다.
  예: H2S s=1.00 seed1102의 E는 기준의 0.1866배.
- round-robin도 일부에서 더 작고 일부에서 더 크다. bucket 보존 자체가 accuracy 최적성은 아니다.

## 수치·재현성 점검

17개 unit test 통과 (기존13+신규4). 모든 입력 SCF 수렴.
PySCF HF energy와 active tensor+core/nuclear constant의 HF expectation 차이 최대 1.42e-12 Hartree.
전체 Hilbert 공간과 Pauli-flip 폐공간의 one-step 결과 및 Hamiltonian action을 비교했다.
Trotter 상태를 고정 입자수 sector에 투영해 누출을 지우지 않았다.

| 점검 | 최대값 |
|---|---:|
| state norm drift | 1.32e-12 |
| 전체 error-vector defect 재구성 잔차 | 5.15e-14 |
| overlap 재구성 잔차 | 5.02e-14 |
| projected I 재구성 차이 | 1.11e-17 |

두 프로세스의 입력·ordering hash 및 예측 파일은 완전히 일치한다.
물리 지표 2,288개는 rtol=1e-8, atol=1e-20 이내 일치하며,
|값|>1e-20인 지표들의 실제 최대 상대 차이는 3.23e-13이다.
**결과 JSON 전체가 bitwise 동일한 것은 아니다.** 극소 부동소수점 차이가 있어
최초 분석의 exact-equality assert가 중단됐고, 비교 방식을 위 수치 허용오차로 바꿨다.
실험·입력·예측·5% 판정 기준은 수정하지 않았다. 첫 분석의 부분 산출물도 보존했다.
`expm_multiply` 내부 norm estimation에 따른 계산 경로 변화가 가능한 설명이나 원인을 별도로 확정하지는 않았다.

## 해석 한계와 논문 방향

네 분자 모두 처음 생성한 입력이지만, 한 basis/한 active-space 크기/HF/first-order에 국한된다.
기하구조는 최적화되지 않았고 화학적 활성공간 수렴성도 검증하지 않았다.
특히 축퇴 orbital의 회전 자유도와 active/frozen 경계를 가로지르는 축퇴에 민감할 수 있다.
orbital sign을 고정하고 실제 orbital 행렬을 저장했지만 이 물리적 선택의 강건성까지 해결한 것은 아니다.

현재 권할 중심 문장:

> 같은 Pauli 항 목록을 재배열할 때 대칭 보존, 섹터 내부 오차, 시간적 오차 누적은
> 서로 다른 평가 축이며, 대칭을 보존하는 parent ordering이 항상 가장 정확하지는 않다.
> 고정 입력의 오차 채널을 분해하고 unfitted leading-error propagation을 이용하여
> 이 차이를 정량적으로 진단한다.

이 문장은 현재 연구의 초점이지 우선권 확정 선언이 아니다.
신규성 감사와 제출 전 필요한 차별화는 [NOVELTY_KO.md](NOVELTY_KO.md)를 참고한다.

기존 tracked 코드와 Overleaf는 수정하지 않았다. commit/push도 하지 않았다.

# 첫 실험: 동일 입력에서 기준 결과 재현

목적은 F₂·NH₃의 서로 다른 과거 결과 중 무엇이 **현재 입력과 코드로 재현되는지** 확인하는 것입니다. Li₂는 대조군입니다. 새로운 ordering의 우월성, 과거 tensor의 동일성, 데이터 차이의 원인을 이 실험 하나로 증명하지는 않습니다.

## 기본 규모

- F₂/HGBS-5/14+2, NH₃/HGBS-5/6+10, Li₂/HGBS-5/6+6.
- Fermionic signed, JW signed, JW magnitude의 세 순서.
- 서로 다른 Python 프로세스로 각 사례를 두 번 실행: 총 18개 결과.
- HF 초기 상태, first-order Trotter, T=1, r=100, coefficient tolerance=1e-12.
- 작업은 순차 실행하며 기본 한 스레드만 사용합니다. 각 사례·반복 작업의 제한 시간은 30분입니다. 이는 예상 시간이 아니며 전체 실행 시간 상한은 기본 여섯 작업 기준 약 3시간과 준비 시간입니다.
- Codex 남은 사용량 13%는 이 준비 작업의 제약입니다. 아래 로컬 실행의 계산 자원 13%를 의미하지 않습니다.

## 실행

QHAT 환경을 활성화하고 스크립트가 있는 폴더에서 실행합니다. QHAT 원본에 복사할 필요가 없습니다.

```bash
conda activate qhat
cd /Users/albertlee0125/Repos/qhat/codex_first_experiment
python first_reproducibility_experiment.py
```

이 명령은 입력 파일·해시·실험 구성을 확인만 하며 파일을 생성하거나 실험을 시작하지 않습니다. 실행하려면:

```bash
python first_reproducibility_experiment.py --execute
```

4-qubit B₂ 입력으로 가장 작은 기능 검사만 할 때:

```bash
python first_reproducibility_experiment.py --cases smoke --execute
```

F₂만 먼저 검사할 때는 `--cases f2 --execute`를 사용합니다. 다른 컴퓨터에서는 `--repo /path/to/qhat`을 지정합니다. 현재 QHAT의 `analysis/benchmark_b2_signed_coefficient_baseline.py`와 해당 입력·과거 CSV가 필요합니다.

## 보존과 실패 처리

실행마다 `repro_runs/first-repro-.../`라는 새 폴더를 만듭니다. tensor와 과거 비교 CSV를 복사해 보존하고, 입력 파일·Hamiltonian·parent·실제 Pauli 순서와 계수의 해시, QHAT 커밋·소스 해시·수치 라이브러리 버전을 기록합니다. 기존 짧은 순서 인덱스 해시끼리는 비교하지 않습니다.

QHAT 소스, 원 tensor, 과거 CSV는 수정하지 않습니다. 입력이 없으면 재생성하지 않고 중단합니다. 코드나 복사한 입력이 실행 중 변경되어도 실패합니다. 로그와 계산 캐시는 새 결과 폴더에 저장합니다. 실행 중 QHAT 코드를 변경하거나 환경을 업데이트하지 마세요.

중단되면 이전 결과에 이어 쓰지 않습니다. 동일 명령을 다시 실행하면 새 캠페인이 만들어집니다. 이는 잘못된 입력으로 이전 결과를 재사용하는 문제를 피하기 위한 의도적인 제한입니다.

## 결과 읽기

먼저 `SUMMARY.md`를 읽으세요. `manifest.json`에는 입력과 설정, `*.repeat*.json`에는 수치 결과, `*.operators.json`에는 실제 연산자·순서, `comparisons.json`에는 비교 상세가 있습니다. 실패 시 해당 `.log`를 확인하세요.

- 반복 결과 PASS: 현재 입력·코드에 대한 반복성이 확인됨. 과거 결과의 원인까지 해명된 것은 아님.
- 반복 PASS + 과거 불일치: tensor/orbital/code 버전 추적이 다음 작업. 과거 결과를 새 값으로 덮어쓰지 않기.
- 반복 FAIL 또는 실행 실패: 새 sweep 전에 환경·입력·수치 허용오차를 점검.

반복 비교는 오차에 절대 허용오차 1e-13, BCH norm에 1e-10, 상대 허용오차 1e-5를 사용합니다. 과거 비교의 상대 허용오차는 1e-3입니다. 항 개수와 해시는 정확히 일치해야 합니다. 오차 1e-12 이하의 비교는 `below_floor_*`로 별도 표시하며 성능 개선 비율로 해석하지 않습니다. 이 허용오차는 감사용 설정이지 통계적 유의성 기준이 아닙니다.

기존 QHAT 계산을 그대로 사용하므로 scalar/global-phase 처리는 이번 실험에서 고치지 않습니다. **물리적 주 지표는 `1 - abs(overlap)`이며 raw state-vector error는 판정에 사용하지 않습니다.** 시간 누적 상쇄의 위상 문제는 이 재현성 검사 이후 별도로 수정·검증해야 합니다.

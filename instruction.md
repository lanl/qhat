# Overnight Fermionic-Ordering Case Studies

이 문서는 fermionic parent-ordering 메커니즘을 보강하기 위해 다음 네 개의 HGBS-5, 18-qubit case를 밤새 실행하는 절차입니다.

| Molecule | Active space | Expected role |
|---|---:|---|
| H2O | 8 occupied + 10 vacant | strong-success polyatomic case |
| NH3 | 8 occupied + 10 vacant | strong-success polyatomic case |
| BeH2 | 6 occupied + 12 vacant | neutral/mild-improvement control |
| O2 | 8 occupied + 10 vacant | strong-failure control |

실행 순서는 다음과 같습니다.

1. 네 Hamiltonian tensor 생성
2. `r=100`에서 full structure ablation 실행
3. fermionic-signed와 JW-signed의 `r=25,50,100,200,400` step scaling 실행

Benchmark CSV는 append-only이며 이미 완료된 항목을 건너뜁니다. 중간에 중단되더라도 같은 명령을 다시 실행하면 이어서 진행됩니다.

## 1. 실행 전 확인

새 터미널을 열고 프로젝트에서 사용하는 Python 환경을 활성화한 뒤 아래 블록을 복사해 실행합니다.

```bash
cd /Users/albertlee0125/Repos/qhat

python -c "import numpy, scipy, openfermion, openfermionpyscf, pyscf, basis_set_exchange, mendeleev, numba, qhat; print('Python dependencies: OK')"

test -f hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/H2O_s-1.00_hgbs-5_as-008-010_JW.config
test -f hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-008-010_JW.config
test -f hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/BeH2_s-1.00_hgbs-5_as-006-012_JW.config
test -f hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/O-O_1.26_hgbs-5_as-008-010_JW.config
```

마지막 네 명령이 아무 메시지 없이 끝나면 config 파일이 모두 존재한다는 뜻입니다. dependency import나 `test`가 실패하면 overnight 명령을 시작하지 말고 환경 또는 파일 경로를 먼저 고쳐야 합니다.

## 2. 권장: 전체 작업을 백그라운드에서 한 번에 실행

아래 블록 전체를 그대로 복사해 실행합니다. macOS의 `caffeinate`가 계산 중 system sleep을 막고, 모든 출력은 `analysis/night_logs/fermionic_overnight.log`에 기록됩니다.

```bash
cd /Users/albertlee0125/Repos/qhat
mkdir -p analysis/night_logs

nohup caffeinate -dimsu zsh -c '
set -euo pipefail
cd /Users/albertlee0125/Repos/qhat

echo "=== Overnight run started ==="
date

echo "=== Generating H2O tensor ==="
python hamiltonian_generator/hamgen.py \
  hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/H2O_s-1.00_hgbs-5_as-008-010_JW.config

echo "=== Generating NH3 tensor ==="
python hamiltonian_generator/hamgen.py \
  hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-008-010_JW.config

echo "=== Generating BeH2 tensor ==="
python hamiltonian_generator/hamgen.py \
  hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/BeH2_s-1.00_hgbs-5_as-006-012_JW.config

echo "=== Generating O2 tensor ==="
python hamiltonian_generator/hamgen.py \
  hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/O-O_1.26_hgbs-5_as-008-010_JW.config

echo "=== Verifying generated tensors ==="
for tensor_path in \
  hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/H2O_s-1.00_hgbs-5_as-008-010.tensors.npz \
  hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-008-010.tensors.npz \
  hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/BeH2_s-1.00_hgbs-5_as-006-012.tensors.npz \
  hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/O-O_1.26_hgbs-5_as-008-010.tensors.npz
do
  if [[ ! -f "$tensor_path" ]]; then
    echo "Missing generated tensor: $tensor_path"
    exit 1
  fi
done

echo "=== Running r=100 structure ablation, 10 random samples per family ==="
python analysis/benchmark_fermionic_structure_ablation.py \
  --tensor hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/H2O_s-1.00_hgbs-5_as-008-010.tensors.npz \
  --tensor hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-008-010.tensors.npz \
  --tensor hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/BeH2_s-1.00_hgbs-5_as-006-012.tensors.npz \
  --tensor hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/O-O_1.26_hgbs-5_as-008-010.tensors.npz \
  --steps 100 \
  --time 1.0 \
  --samples 10 \
  --output analysis/fermionic_structure_ablation_extension_hgbs5.csv

echo "=== Running signed fermionic versus signed JW step scaling ==="
python analysis/benchmark_fermionic_structure_ablation.py \
  --tensor hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/H2O_s-1.00_hgbs-5_as-008-010.tensors.npz \
  --tensor hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-008-010.tensors.npz \
  --tensor hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/BeH2_s-1.00_hgbs-5_as-006-012.tensors.npz \
  --tensor hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/O-O_1.26_hgbs-5_as-008-010.tensors.npz \
  --steps 25 50 100 200 400 \
  --time 1.0 \
  --samples 0 \
  --schedules fermionic_signed_reference jw_signed_baseline \
  --output analysis/final_validation_step_scaling_extension_hgbs5.csv

echo "=== Overnight run completed ==="
date
' > analysis/night_logs/fermionic_overnight.log 2>&1 &

OVERNIGHT_PID=$!
echo "Overnight PID: $OVERNIGHT_PID"
echo "Log: /Users/albertlee0125/Repos/qhat/analysis/night_logs/fermionic_overnight.log"
```

터미널에 표시된 `Overnight PID`를 기록해 둡니다.

## 3. 진행 상황 확인

실시간 로그를 확인합니다.

```bash
cd /Users/albertlee0125/Repos/qhat
tail -f analysis/night_logs/fermionic_overnight.log
```

`tail` 화면에서 빠져나오려면 `Ctrl-C`를 누릅니다. 계산 작업 자체는 계속 실행됩니다.

프로세스가 아직 실행 중인지 확인하려면 기록해 둔 PID를 사용합니다.

```bash
ps -p OVERNIGHT_PID
```

위 명령의 `OVERNIGHT_PID`를 실제 숫자로 바꿉니다.

## 4. 아침에 결과 확인

```bash
cd /Users/albertlee0125/Repos/qhat

tail -50 analysis/night_logs/fermionic_overnight.log

wc -l analysis/fermionic_structure_ablation_extension_hgbs5.csv
wc -l analysis/final_validation_step_scaling_extension_hgbs5.csv

awk -F, 'NR > 1 {count[$1]++} END {for (status in count) print status, count[status]}' \
  analysis/fermionic_structure_ablation_extension_hgbs5.csv

awk -F, 'NR > 1 {count[$1]++} END {for (status in count) print status, count[status]}' \
  analysis/final_validation_step_scaling_extension_hgbs5.csv
```

빈 CSV에서 시작했다면 예상되는 결과 행 수는 다음과 같습니다.

- structure ablation: 4 cases x `(4 deterministic + 2 x 10 random)` = 96 data rows, header 포함 97 lines
- step scaling: 4 cases x 2 schedules x 5 step counts = 40 data rows, header 포함 41 lines

모든 data row의 첫 번째 열이 `success`여야 합니다.

## 5. 중단되었을 때 재시작

먼저 기존 작업이 실제로 종료되었는지 확인합니다. 같은 benchmark를 동시에 두 번 실행하면 안 됩니다.

종료가 확인되면 **2번의 전체 백그라운드 블록을 그대로 다시 실행**합니다. Hamiltonian generator는 생성된 intermediate를 재사용하고, benchmark는 output CSV에서 완료된 schedule을 읽어 건너뜁니다.

## 6. Random sample을 10개에서 20개로 확장

첫 번째 overnight run이 성공한 후 random distribution을 더 안정적으로 만들고 싶다면 아래 명령을 실행합니다. 같은 output CSV를 사용하므로 기존 sample 0--9는 건너뛰고 sample 10--19만 추가됩니다.

```bash
cd /Users/albertlee0125/Repos/qhat

python analysis/benchmark_fermionic_structure_ablation.py \
  --tensor hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/H2O_s-1.00_hgbs-5_as-008-010.tensors.npz \
  --tensor hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-008-010.tensors.npz \
  --tensor hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/BeH2_s-1.00_hgbs-5_as-006-012.tensors.npz \
  --tensor hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/O-O_1.26_hgbs-5_as-008-010.tensors.npz \
  --steps 100 \
  --time 1.0 \
  --samples 20 \
  --output analysis/fermionic_structure_ablation_extension_hgbs5.csv
```

20 samples까지 완료되면 structure-ablation CSV는 4 cases x `(4 + 2 x 20)` = 176 data rows, header 포함 177 lines가 됩니다.

## Output files

- `analysis/fermionic_structure_ablation_extension_hgbs5.csv`
- `analysis/final_validation_step_scaling_extension_hgbs5.csv`
- `analysis/night_logs/fermionic_overnight.log`

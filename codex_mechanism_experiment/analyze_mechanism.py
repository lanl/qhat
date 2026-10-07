#!/usr/bin/env python3
"""Read completed mechanism JSON; create an exclusive analysis directory."""
import argparse
import collections
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import numpy as np

BASE=['fermionic_signed','jw_signed','jw_magnitude']
LABELS={'fermionic_signed':'F signed','jw_signed':'JW signed','jw_magnitude':'JW magnitude',
        'round_robin':'Round robin','within_blocks_20260917':'Within-block'}
COLORS={'fermionic_signed':'#176b87','jw_signed':'#cf573c','jw_magnitude':'#8462a0'}


def save(path,obj):
    with path.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False)


def build_summary(rows):
    out={'rows':len(rows),'baseline_rows':sum(r['ordering'] in BASE for r in rows),
         'decomposed_rows':sum('coherence_factor' in r for r in rows),
         'floor_rows':sum(r['below_analysis_floor'] for r in rows),
         'old_anchor_matches':sum(r.get('old_anchor_match',False) for r in rows),
         'max_state_reconstruction_residual':max(r.get('state_reconstruction_residual',0) for r in rows),
         'max_overlap_reconstruction_residual':max(r.get('overlap_reconstruction_residual',0) for r in rows),
         'max_projected_reconstruction_absolute_difference':max(abs(r['infidelity']-r.get('reconstructed_projected_infidelity',r['infidelity'])) for r in rows),
         'max_norm_drift':max(max(r['exact_norm_drift'],r['trotter_norm_drift']) for r in rows),
         'max_exact_spin_leakage':max(r['exact_spin_leakage'] for r in rows),
         'max_sector_error_sum_residual':max(abs(r['infidelity']-r['inside_sector_error']-r['outside_sector_error']) for r in rows),
         'anchors':[],'f2_time_decomposition':[],'grid_winners':{},'controls':{},'slopes':[]}
    for case in ('f2','nh3','li2'):
        case_rows=[r for r in rows if r['case']==case]
        anchors=[r for r in case_rows if r['T']==1 and r['steps']==100]
        for r in anchors:
            if r['ordering'] in BASE:
                out['anchors'].append({**r,'outside_error_fraction':r['outside_sector_error']/r['infidelity']})
        wins=collections.Counter()
        for T in sorted({r['T'] for r in case_rows}):
            for steps in sorted({r['steps'] for r in case_rows}):
                group=[r for r in case_rows if r['T']==T and r['steps']==steps and r['ordering'] in BASE]
                if len(group)==3:
                    winner=min(group,key=lambda r:r['one_minus_overlap'])
                    wins[winner['ordering']]+=1
            for name in BASE:
                group=sorted([r for r in case_rows if r['T']==T and r['ordering']==name and not r['below_analysis_floor']],key=lambda r:r['steps'])
                if len(group)>=2:
                    x,y=group[-2:]
                    slope=math.log(y['one_minus_overlap']/x['one_minus_overlap'])/math.log(y['steps']/x['steps'])
                    out['slopes'].append({'case':case,'ordering':name,'T':T,'r_low':x['steps'],'r_high':y['steps'],'slope':slope})
        out['grid_winners'][case]=dict(wins)
        if anchors:
            ref=next(r for r in anchors if r['ordering']=='fermionic_signed')
            out['controls'][case]=[{**r,'error_ratio_to_fermionic':r['one_minus_overlap']/ref['one_minus_overlap']} for r in anchors if r['ordering'] not in BASE]
    for r in rows:
        if r['case']=='f2' and r['steps']==100 and r['ordering'] in BASE:
            out['f2_time_decomposition'].append({k:r[k] for k in ('ordering','T','infidelity','projected_local_power','coherence_factor','projected_interference','spin_sector_leakage','conditional_sector_infidelity')})
    return out


def plot(rows,dest):
    os.environ['MPLCONFIGDIR']=str(dest/'mpl_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':160,'savefig.facecolor':'white'})
    fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True)
    for i,case in enumerate(('f2','nh3','li2')):
        for j,T in enumerate((.25,.5,1.)):
            ax=axes[i,j]
            for name in BASE:
                subset=sorted([r for r in rows if r['case']==case and r['T']==T and r['ordering']==name],key=lambda r:r['steps'])
                ax.loglog([r['steps'] for r in subset],[r['one_minus_overlap'] for r in subset],'-o',ms=4,color=COLORS[name],label=LABELS[name])
                for r in subset:
                    if r['below_analysis_floor']:
                        ax.plot(r['steps'],r['one_minus_overlap'],'o',mfc='white',mec=COLORS[name],ms=7)
            ax.axhline(1e-12,color='#999999',ls=':',lw=.8)
            ax.set_xticks([25,50,100,200],['25','50','100','200'])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_title(f'{case.upper()} | T={T:g}')
            ax.grid(alpha=.15,which='both')
            if i==2:ax.set_xlabel('Trotter steps r')
            if j==0:ax.set_ylabel('1 - |normalized overlap|')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Fixed-input first-order ordering comparison\nDotted line: predeclared numerical-analysis floor',y=1.01)
    fig.tight_layout();fig.savefig(dest/'baseline_grid.png',bbox_inches='tight');plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    fields=[('infidelity','Final infidelity I'),('projected_local_power','Projected local power D'),('coherence_factor','Accumulation factor C = I / D')]
    for ax,(field,label) in zip(axes,fields):
        for name in BASE:
            subset=sorted([r for r in rows if r['case']=='f2' and r['steps']==100 and r['ordering']==name],key=lambda r:r['T'])
            ax.plot([r['T'] for r in subset],[r[field] for r in subset],'-o',color=COLORS[name],label=LABELS[name])
        ax.set_yscale('log');ax.set_xlabel('Evolution time T');ax.set_ylabel(label);ax.grid(alpha=.15)
        if field=='coherence_factor':ax.axhline(1,color='#555555',ls='--',lw=1)
    axes[0].legend(fontsize=8)
    fig.suptitle('F2, r=100: local error size and finite-time accumulation are distinct\nC < 1: net destructive interference; C > 1: net constructive interference',y=1.08)
    fig.tight_layout();fig.savefig(dest/'f2_decomposition.png',bbox_inches='tight');plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    for ax,case in zip(axes,('f2','nh3','li2')):
        subset=[r for r in rows if r['case']==case and r['T']==1 and r['steps']==100]
        labels=[LABELS.get(r['ordering'],r['ordering'].replace('blocks_random_','Seed ')) for r in subset]
        x=np.arange(len(subset))
        inside=np.array([r['inside_sector_error'] for r in subset]);outside=np.array([r['outside_sector_error'] for r in subset])
        ax.bar(x,inside,color='#176b87',label='Inside-sector error')
        ax.bar(x,outside,bottom=inside,color='#d58c3c',label='Outside-sector error')
        ax.set_yscale('log');ax.set_xticks(x,labels,rotation=65,ha='right',fontsize=8)
        ax.set_ylim(min(inside+outside)*.2,max(inside+outside)*2.5)
        ax.set_title(case.upper());ax.set_ylabel('Standard infidelity = inside + outside')
    axes[0].legend(fontsize=8)
    fig.suptitle('Structural controls at T=1, r=100\nFive fixed random seeds are descriptive controls, not independent molecules',y=1.01)
    fig.tight_layout();fig.savefig(dest/'sector_controls.png',bbox_inches='tight');plt.close(fig)


def report(summary,campaign,dest):
    s=summary
    lines=['# Ordering 메커니즘 실험 결과','',
           f'원본 캠페인: `{campaign}`','',
           '## 검증 상태','',
           f'- 완료: {s["rows"]}개 설정; baseline {s["baseline_rows"]}개, 통제 {s["rows"]-s["baseline_rows"]}개.',
           f'- 정확한 시간별 오차 분해: {s["decomposed_rows"]}개.',
           f'- 기존 물리적 오차 anchor 일치: {s["old_anchor_matches"]}/9.',
           f'- 최대 상태벡터 재구성 잔차: {s["max_state_reconstruction_residual"]:.3e}.',
           f'- 최대 복소 overlap 재구성 잔차: {s["max_overlap_reconstruction_residual"]:.3e}.',
           f'- 최대 정규화 drift: {s["max_norm_drift"]:.3e}.',
           f'- exact 상태의 최대 spin-sector 누출: {s["max_exact_spin_leakage"]:.3e}.',
           f'- 분석 바닥 E<=1e-12: {s["floor_rows"]}개. 바닥 이하 값은 성능비/수렴률 주장에 쓰지 않음.','',
           'E=1-|정규화 overlap|, I=1-|정규화 overlap|². 두 지표를 혼용하지 않음.',
           'exact와 Trotter는 동일한 identity-free Pauli Hamiltonian을 사용한다.',
           '도달 가능 공간은 각 Pauli 항에 대해 닫혀 있고 입자수 누출 상태를 제거하지 않는다.','',
           '## 기준 결과: T=1, r=100','',
           '|분자|순서|E|I 중 구역 밖 비율|국소 오차 총량 D|누적 계수 C=I/D|',
           '|---|---|---:|---:|---:|---:|']
    for r in s['anchors']:
        lines.append(f'|{r["case"]}|{r["ordering"]}|{r["one_minus_overlap"]:.6e}|{100*r["outside_error_fraction"]:.3f}%|{r["projected_local_power"]:.6e}|{r["coherence_factor"]:.6g}|')
    f=next(r for r in s['anchors'] if r['case']=='f2' and r['ordering']=='fermionic_signed')
    j=next(r for r in s['anchors'] if r['case']=='f2' and r['ordering']=='jw_signed')
    lines += ['', '## F2 순위 역전의 정량적 분해','',
              f'- JW signed / Fermionic signed의 초기 BCH norm 비: {j["bch_hf_norm"]/f["bch_hf_norm"]:.6g}.',
              f'- JW signed / Fermionic signed의 국소 projected power D 비: {j["projected_local_power"]/f["projected_local_power"]:.6g}.',
              f'- JW signed / Fermionic signed의 누적 계수 C 비: {j["coherence_factor"]/f["coherence_factor"]:.6g}.',
              f'- 두 비의 곱 = 최종 infidelity 비: {j["infidelity"]/f["infidelity"]:.6g}.','',
              'D와 C의 곱 분해는 정확한 항등식이지 새로운 예측 이론의 증명은 아니다.',
              'C가 작아도 C>1이면 교차항의 순효과는 보강이다. 이를 음의 순상쇄라고 부르면 안 된다.',
              '전파된 국소 오차의 방향/위상 누적을 무시한 초기-HF norm만으로 최종 순위를 단정할 수 없다.','',
              '### F2 시간 변화 (r=100)','',
              '|T|순서|I|D|C|교차항 I-D|','|---:|---|---:|---:|---:|---:|']
    for r in s['f2_time_decomposition']:
        lines.append(f'|{r["T"]:g}|{r["ordering"]}|{r["infidelity"]:.5e}|{r["projected_local_power"]:.5e}|{r["coherence_factor"]:.5g}|{r["projected_interference"]:.5e}|')
    lines += ['', '## 모든 grid 결과와 통제 실험','',
              '아래 winner 횟수는 동일 분자의 여러 설정을 세는 기술 통계이며 독립 표본 수가 아니다.']
    for case,wins in s['grid_winners'].items():
        lines += [f'- {case}: {wins}']
    lines += ['', '|분자|통제 순서|Fermionic 대비 E 비|누적 계수 C|', '|---|---|---:|---:|']
    for case,group in s['controls'].items():
        for r in group:
            lines.append(f'|{case}|{r["ordering"]}|{r["error_ratio_to_fermionic"]:.6g}|{r["coherence_factor"]:.6g}|')
    lines += ['', '## 수렴률','',
              '분석 바닥 이상인 가장 큰 두 step 수로 구한 log(E)/log(r) 기울기. 일반적인 first-order 상태 오차의 제곱 지표는 -2에 접근할 수 있으나 이를 강제하지 않았다.','',
              '|분자|T|순서|r 구간|기울기|','|---|---:|---|---|---:|']
    for r in s['slopes']:
        lines.append(f'|{r["case"]}|{r["T"]:g}|{r["ordering"]}|{r["r_low"]}–{r["r_high"]}|{r["slope"]:.4f}|')
    lines += ['', '## 이 실험으로 해결되지 않은 것','',
              '- 과거 tensor/코드 불일치의 원인은 이 실험으로 판별하지 않았다.',
              '- 세 사례는 이미 결과를 본 사례다. 외부·전향 검증이 아니며 일반적 우월성을 입증하지 않는다.',
              '- 시간별 분해는 설명 도구다. exact trajectory가 필요한 이 진단을 효율적인 ordering 선택 알고리즘으로 주장할 수 없다.',
              '- 5개 block random seed와 round robin은 구조적 통제이나 서로 같은 크기의 perturbation은 아니다.',
              '- 1-body/2-body의 독립 효과, orbital-gauge 견고성, 동일 정확도 회로 비용은 검증하지 않았다.',
              '- 작은 누출 값의 소수점 수준 차이보다 전체 오차와 수치 잔차를 우선 해석한다.',
              '- conditional sector metric은 사후 선택의 비용을 무시한 성능 주장에 사용하지 않는다.','',
              '## 산출물','',
              '- `baseline_grid.png`: 시간·step 변화와 세 ordering 비교.',
              '- `f2_decomposition.png`: F2의 국소 오차 크기와 누적 계수.',
              '- `sector_controls.png`: 구조적 통제에서 구역 안/밖 오차.',
              '- `analysis.json`: 표·검증 통계의 원본.',
              '- 상위 캠페인의 `inputs/`, `manifest.json`, `*.setup.json`, `*.trace.json`, `results.json`: 재현 입력과 수치 기록.','']
    with (dest/'REPORT.md').open('x') as stream:stream.write('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('campaign',type=Path)
    parser.add_argument('--repeat-campaign',type=Path)
    args=parser.parse_args()
    completion=json.loads((args.campaign/'completion.json').read_text())
    if completion['status']!='complete':raise RuntimeError('Campaign is incomplete')
    rows=json.loads((args.campaign/'results.json').read_text())
    destination=Path(tempfile.mkdtemp(prefix='analysis-',dir=args.campaign))
    s=build_summary(rows)
    s['analysis_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    setups=[json.loads(p.read_text()) for p in sorted(args.campaign.glob('*.setup.json'))]
    s['max_full_space_slice_residual']=max(v for d in setups for v in d['full_space_slice_residuals'].values())
    s['max_hamiltonian_action_residual']=max(d['hamiltonian_action_residual'] for d in setups)
    s['max_hamiltonian_sector_offblock_frobenius']=max(d['hamiltonian_sector_offblock_frobenius'] for d in setups)
    if args.repeat_campaign:
        previous=json.loads((args.repeat_campaign/'results.json').read_text())
        checks=[]
        for a in previous:
            b=next(r for r in rows if all(r[k]==a[k] for k in ('case','ordering','T','steps')))
            for key in ('one_minus_overlap','infidelity','projected_local_power','coherence_factor',
                        'inside_sector_error','outside_sector_error'):
                ok=math.isclose(a[key],b[key],rel_tol=1e-7,abs_tol=1e-15)
                checks.append({'case':a['case'],'ordering':a['ordering'],'metric':key,
                               'first':a[key],'fresh':b[key],'match':ok})
        s['independent_repeat_checks']=checks
        s['independent_repeat_pass']=all(x['match'] for x in checks)
        if not s['independent_repeat_pass']:raise RuntimeError('Independent repetition failed')
    save(destination/'analysis.json',s)
    report(s,args.campaign,destination);plot(rows,destination)
    print(destination)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Apply the fixed transfer criteria, retaining counterexamples and near ties."""
import argparse
import itertools
import json
from pathlib import Path
import tempfile
import external_validation as e


def rowkey(row):
    return tuple(row[k] for k in ('case','ordering','T','steps'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('campaign',type=Path)
    parser.add_argument('--repeat',type=Path);args=parser.parse_args()
    p=args.campaign
    destination=Path(tempfile.mkdtemp(prefix='analysis-',dir=p))
    e.m.init_numerics(destination/'cache');np=e.m.np
    rows=json.loads((p/'results.json').read_text())
    predicted=json.loads((p/'predictions_before_target_evolution.json').read_text())
    bykey={rowkey(r):r for r in rows}
    if len(rows)!=272 or len(bykey)!=272:raise RuntimeError('Missing or duplicate rows')
    comparisons=[]
    for pred in predicted:
        actual=bykey[rowkey(pred)]
        comparisons.append({**pred,'actual_infidelity':actual['infidelity'],
                            'eligible':not actual['below_analysis_floor'],
                            'dynamic_relative_error':pred['predicted_infidelity']/actual['infidelity']-1,
                            'static_relative_error':pred['static_bch_proxy_infidelity']/actual['infidelity']-1})
    e.m.save(destination/'prediction_comparisons.json',comparisons)
    matched={rowkey(r):r for r in comparisons}
    pairs=[];case_summary=[];slopes=[]
    for spec in e.specifications():
        case=spec['case'];wins=[]
        for T in e.TIMES:
            group=[bykey[(case,name,T,100)] for name in e.m.BASE]
            winners=sorted(group,key=lambda r:r['infidelity'])
            wins.append(winners[0]['ordering'])
            for a,b in itertools.combinations(group,2):
                gap=abs(a['infidelity']-b['infidelity'])/min(a['infidelity'],b['infidelity'])
                eligible=not(a['below_analysis_floor'] or b['below_analysis_floor']) and gap>.01
                pa,pb=matched[rowkey(a)],matched[rowkey(b)]
                sign=np.sign(a['infidelity']-b['infidelity'])
                pairs.append({'case':case,'T':T,'a':a['ordering'],'b':b['ordering'],
                              'relative_gap':gap,'eligible':eligible,
                              'dynamic_correct':bool(sign==np.sign(pa['predicted_infidelity']-pb['predicted_infidelity'])),
                              'static_correct':bool(sign==np.sign(pa['static_bch_proxy_infidelity']-pb['static_bch_proxy_infidelity']))})
            for name in e.m.BASE:
                group_r=[bykey[(case,name,T,r)] for r in e.STEPS]
                if all(not row['below_analysis_floor'] for row in group_r):
                    slopes.append(float(np.polyfit(np.log(e.STEPS),np.log([r['one_minus_overlap'] for r in group_r]),1)[0]))
        anchors={name:bykey[(case,name,1.,100)] for name in e.m.BASE}
        base=anchors['fermionic_signed']['one_minus_overlap']
        controls={name:bykey[(case,name,1.,100)]['one_minus_overlap']/base
                  for name in {r['ordering'] for r in rows}-set(e.m.BASE)}
        case_summary.append({'case':case,'family':spec['family'],'winners_by_time':wins,
                             'anchors':anchors,'control_error_ratio_to_fermionic':controls,
                             'max_relative_dynamic_error_r100':max(abs(r['dynamic_relative_error']) for r in comparisons
                                if r['case']==case and r['steps']==100 and r['eligible'])})
    eligible=[r for r in comparisons if r['eligible'] and r['steps'] in (100,200)]
    rankable=[r for r in pairs if r['eligible']]
    traces=[r for r in rows if 'state_reconstruction_residual' in r]
    summary={'campaign':str(p.resolve()),'unique_families':4,'geometries':8,
        'baseline_rows':216,'control_rows':56,'prediction_rows':len(comparisons),
        'excluded_below_floor_rows':sum(not r['eligible'] for r in comparisons),
        'primary_eligible_rows_r100_r200':len(eligible),
        'primary_prediction_failures':sum(abs(r['dynamic_relative_error'])>.05 for r in eligible),
        'r100_dynamic_max_relative_error':max(abs(r['dynamic_relative_error']) for r in eligible if r['steps']==100),
        'r200_dynamic_max_relative_error':max(abs(r['dynamic_relative_error']) for r in eligible if r['steps']==200),
        'dynamic_median_relative_error':float(np.median([abs(r['dynamic_relative_error']) for r in eligible])),
        'static_median_relative_error':float(np.median([abs(r['static_relative_error']) for r in eligible])),
        'static_max_relative_error':max(abs(r['static_relative_error']) for r in eligible),
        'rank_eligible_pairs':len(rankable),'rank_excluded_pairs':len(pairs)-len(rankable),
        'dynamic_rank_correct':sum(r['dynamic_correct'] for r in rankable),
        'static_rank_correct':sum(r['static_correct'] for r in rankable),
        'case_time_exact_winner_counts':{name:sum(c['winners_by_time'].count(name) for c in case_summary) for name in e.m.BASE},
        'r_scaling_slope_min_max':[min(slopes),max(slopes)],
        'max_norm_drift':max(max(r['exact_norm_drift'],r['trotter_norm_drift']) for r in rows),
        'max_defect_reconstruction':max(r['state_reconstruction_residual'] for r in traces),
        'max_overlap_reconstruction':max(r['overlap_reconstruction_residual'] for r in traces),
        'max_projected_reconstruction_difference':max(abs(r['infidelity']-r['reconstructed_projected_infidelity']) for r in traces),
        'fermionic_max_spin_leakage':max(r['spin_sector_leakage'] for r in rows if r['ordering']=='fermionic_signed'),
        'jw_signed_max_spin_leakage':max(r['spin_sector_leakage'] for r in rows if r['ordering']=='jw_signed'),
        'baseline_C_min_max':[min(r['coherence_factor'] for r in traces if r['ordering'] in e.m.BASE),
                               max(r['coherence_factor'] for r in traces if r['ordering'] in e.m.BASE)],
        'baseline_net_destructive_count':sum(r['coherence_factor']<1 for r in traces if r['ordering'] in e.m.BASE),
        'source_sha256':e.m.sha(Path(__file__))}
    if args.repeat:
        repeat=json.loads((args.repeat/'results.json').read_text())
        # expm_multiply's internal norm estimation can choose slightly different
        # floating-point paths. Compare physics, not bitwise residual diagnostics.
        physics=('one_minus_overlap','infidelity','inside_sector_error','outside_sector_error',
                 'particle_number_leakage','spin_sector_leakage','conditional_sector_infidelity',
                 'projected_local_power','projected_interference','coherence_factor')
        differences=[]
        for a,b in zip(rows,repeat):
            if rowkey(a)!=rowkey(b):raise RuntimeError('Repeat row index mismatch')
            for key in physics:
                if key not in a:continue
                np.testing.assert_allclose(a[key],b[key],rtol=1e-8,atol=1e-20,
                                           err_msg=f'{rowkey(a)}:{key}')
                differences.append({'metric':key,'absolute_difference':abs(a[key]-b[key]),
                                    'relative_difference':abs(a[key]-b[key])/abs(a[key]) if abs(a[key])>1e-20 else None})
        e.m.np.testing.assert_equal(predicted,json.loads((args.repeat/'predictions_before_target_evolution.json').read_text()))
        for spec in e.specifications():
            name=f"{spec['case']}.setup.json"
            s1=json.loads((p/name).read_text());s2=json.loads((args.repeat/name).read_text())
            if s1['fingerprints']!=s2['fingerprints'] or s1['tensor_sha256']!=s2['tensor_sha256']:
                raise RuntimeError('Repeat fingerprint mismatch')
        summary['repeat']={'path':str(args.repeat.resolve()),'all_272_result_rows_exactly_equal':rows==repeat,
                           'all_physical_metrics_within_rtol1e_8_atol1e_20':True,
                           'physical_metric_comparisons':len(differences),
                           'max_relative_difference_above_1e_20':max(d['relative_difference'] for d in differences if d['relative_difference'] is not None),
                           'all_216_predictions_exactly_equal':True,'all_input_and_order_hashes_equal':True}
    e.m.save(destination/'summary.json',summary)
    e.m.save(destination/'case_summary.json',case_summary)
    e.m.save(destination/'rank_pairs.json',pairs)
    lines=['# Prospective transfer results','',f"Campaign: `{p.name}`.",
           'Four new molecular families; two idealized geometries per family; 6-31g/CAS(4e,5o), HF.',
           'Rows/geometries/times are not independent molecules. No chemically converged molecular accuracy claim.','',
           f"Completed 216 baseline + 56 control rows. Prediction floor exclusions: {summary['excluded_below_floor_rows']}.",
           f"Primary 5% prediction failures: {summary['primary_prediction_failures']}/{len(eligible)}.",
           f"Max relative I error r=100: {summary['r100_dynamic_max_relative_error']:.6%}; r=200: {summary['r200_dynamic_max_relative_error']:.6%}.",
           f"Rank agreement (>1% actual margin, r=100): dynamic {summary['dynamic_rank_correct']}/{len(rankable)}, static HF-BCH {summary['static_rank_correct']}/{len(rankable)}; excluded pairs {summary['rank_excluded_pairs']}.",
           f"Median relative error (r=100,200): dynamic {summary['dynamic_median_relative_error']:.6%}, static {summary['static_median_relative_error']:.3%}.",
           '', '## All T=1,r=100 anchors (E=1-|overlap|)', '',
           '| Case | Fermionic signed | JW signed | JW magnitude | Best |', '|---|---:|---:|---:|---|']
    for c in case_summary:
        a=c['anchors'];lines.append('| '+c['case']+' | '+' | '.join(f"{a[n]['one_minus_overlap']:.9e}" for n in e.m.BASE)+' | '+c['winners_by_time'][-1]+' |')
    lines+=['','## Interpretation guards','',
            'Predictions are unfitted leading-error propagation, a standard perturbative construction, not a new general theory.',
            'The static comparator is a local HF BCH proxy, not an implementation of published best algorithms.',
            'Symmetry preservation is not a guarantee of lower total state error. Retain every contrary ranking.',
            f"Baseline C range at r=100: {summary['baseline_C_min_max']}; C<1 cases: {summary['baseline_net_destructive_count']}. A C reduction above 1 is NOT net cancellation.",
            '', '## Machine-readable checks','', '```json',json.dumps(summary,indent=2),'```','']
    with (destination/'REPORT.md').open('x') as f:f.write('\n'.join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for i,spec in enumerate(e.specifications()[::2]):
        selected=[r for r in comparisons if r['case'].startswith(spec['family']+'_') and r['steps']==100 and r['eligible']]
        axes[0].scatter([r['actual_infidelity'] for r in selected],[r['predicted_infidelity'] for r in selected],s=32,label=spec['family'],alpha=.8)
    lim=[min(r['actual_infidelity'] for r in comparisons if r['steps']==100)*.7,
         max(r['actual_infidelity'] for r in comparisons if r['steps']==100)*1.4]
    axes[0].plot(lim,lim,'k--',lw=1);axes[0].set(xscale='log',yscale='log',xlim=lim,ylim=lim,
        xlabel='Actual infidelity I',ylabel='Predicted infidelity I',title='Unfitted transfer: r = 100')
    axes[0].legend(fontsize=9)
    x=np.arange(8)
    for name,marker in [('jw_signed','o'),('jw_magnitude','s')]:
        axes[1].plot(x,[c['anchors'][name]['one_minus_overlap']/c['anchors']['fermionic_signed']['one_minus_overlap'] for c in case_summary],marker=marker,label=name)
    axes[1].axhline(1,color='k',ls='--',lw=1);axes[1].set(yscale='log',ylabel='E / E(fermionic signed)',title='All anchors: T = 1, r = 100')
    axes[1].set_xticks(x,[c['case'].replace('_s','\ns=') for c in case_summary],fontsize=8)
    axes[1].legend(fontsize=9);fig.savefig(destination/'transfer_validation.png',dpi=180);plt.close(fig)
    print('ANALYSIS',destination);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()

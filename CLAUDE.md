# Multi-Level Multi-Track Converter Optimization

MATLAB design-optimization code for a 6.25 kW, 1500 V → 48 V LLC-based converter
(GE Vernova project, MIT). Authors: Julia Estrin (magnetics), Qijia Li (switches).

## What the code does

Pipeline, in order:
1. **LLC tank** (`LLC Design/designLLC_v3.m`): turns ratio N, max Qe, Cr/Lr/Lm, worst-case `Ir_rms`.
2. **VIRT planar transformer** (`Transformer Design/optimize_VIRT_2LMT.m` → `calculateVIRT_round_2LMT.m`):
   `fmincon` minimizes core loss (Steinmetz + ferrite eddy loss) + copper loss (Dowell),
   subject to Bsat, footprint and temperature (cold-plate model).
3. **Resonant inductor** (`Lr Design/optimize_Lr_EQ.m` → `calculate_Lr_EQ.m`): loops over the number of turns N,
   with the air gap solved numerically including fringing.
4. **Primary switches** (`Switch Functions/analyzePriSwitches_v5.m`): GaN/SiC, 1–8 in parallel,
   with Tj < Tj,max − 30 °C, giving a loss-vs-area Pareto front.
5. **Secondary switches** (`Switch Functions/analyzeSecSwitches_v2.m`): conduction + gate loss.
6. **System efficiency** (`Switch Functions/calcEfficiency_v5.m`): primary × secondary Pareto combined with
   the magnetics loss/area, giving `effOut.overallParetoTable`.

Topology strings: `"2-level Multitrack"` (current focus), `"Multitrack"`, `"Multilevel Multitrack"`
(switch counts as [jj25, jj75] pairs), `"Fullbridge LLC"`.

Baseline specs: Vin 1500 V, Vo 48 V, nt = 2, Pmax 6.25 kW, fsw = f0 = 1 MHz, Ln = 7, ±5 % regulation,
±25 % frequency range, Fair-Rite 80 material. `designLLC_v3` gives N = 16, so np = N/nt = 8 turns per track.

## Current entry scripts (use these)

- `TwoLevelMutliTrack_design.m`: **main design script**. LLC → `k_area` sweep (0.1–2) of the round
  VIRT, keeps designs ≤ 4 in², best loss → subtracts transformer leakage (`calc_Llk`) from Lr →
  `optimize_Lr_EQ` → switches → `calcEfficiency_v5`. Has the newest transformer code. Slow run
  (191 transformer optimizations + 12 inductor optimizations).
- `main_only_switch_v2.m`: **switch-only reference** (Qijia's working version). Uses fixed
  magnetics placeholders (39.86 W, 4 in²).
- The switch section of `TwoLevelMutliTrack_design.m` is aligned with `main_only_switch_v2.m`:
  same functions, data files and arguments. Keep them aligned when either changes:
  - `analyzePriSwitches_v5(topology, f_per, Pmax, fsw, Ir_rms, 1, 8, [], [], 10000, 'GaN Data tf.xlsx', 'SiC Data tf.xlsx')`
  - `analyzeSecSwitches_v2(Pmax, f_per, f0, 1, 10, [4], [4], 'Sec Data Both.xlsx')`. **Always pass
    the file**: the default `'Sec Data v2.xlsx'` does not exist.
  - `calcEfficiency_v5(..., "pareto", "pareto", ...)`, then `extractBestPointSummary(effOut.overallParetoTable, ...)`

Broken or legacy scripts (do not use as references):
- `main.m`: calls `optimize_SPVIRT` / `SPVIRT3Dfigure_square`, which don't exist.
- `Ln_RegulationRange_Sweep_05_12_2026.m`: `design_params` is missing `Vin`, `k_area` and `sigma_1MHz`, so every
  iteration errors and the error is swallowed by `try/catch` (all-NaN plots); it also uses `effOut.points(1)`, which is
  not the min-loss point; and its switch calls are outdated.
- `TwoLevelMutliTrack_v2.m`: missing `Vin`.
- `compare_topology*.m`, `main_only_switch.m`, `design_summary_figure.m`: older APIs.

## Live function versions

| Purpose | Use | Superseded |
|---|---|---|
| LLC | `designLLC_v3` | `designLLC` (no 2-level branch); `designLLC_v4` only in Qijia's test (see below) |
| Transformer | `optimize_VIRT_2LMT`, `calculateVIRT_round_2LMT` | `optimize_VIRT`, `calculateVIRT_round`, `calculate_VIRT_design`, SPVIRT files |
| Copper loss (transformer) | `calculate_Pcond_2LMT` | `calculate_Pcond`, `calculate_Pcond_SPVIRT` |
| Leakage | `Transformer Design/calc_Llk.m` (4 args) | `old code/calc_Llk.m` (same name, different signature; must not be on the path) |
| Inductor | `optimize_Lr_EQ`, `calculate_Lr_EQ` | `optimize_Lr`, `calculate_Lr` |
| Primary switches | `analyzePriSwitches_v5`, `pickBestDevices_v3`, `packOneLocal_v3`, `resolveRankJJList_v3` | v1–v4 |
| Secondary switches | `analyzeSecSwitches_v2` (uses `pickBestDevices` / `packOneLocal` / `resolveRankJJList` v1) | `analyzeSecSwitches` |
| Efficiency | `calcEfficiency_v5` | v1–v4 |

Alternative transformer cores exist but are not in the active flow (commented blocks in the
design script): racetrack (`*_raceTrack_2LMT`), PQ (`*_PQ_2LMT`), round-EQ (`*_roundEQ_2LMT`).
`Lr Design/Fairrite80.m` is the 1 MHz Steinmetz fit used (hard-coded) in `calculate_Lr_EQ`.

`designLLC_v3` vs `v4`: identical results for these specs. v4 reports N per track (8) instead of
the total (16), but only for 2-level. If v4 is adopted, scripts must use `np = N` instead of `N/nt`, and
the other topologies should be made consistent.

## Known issues (open)

Transformer layer count is hard-wired in three places that must change together:
- `h_pcb = 6*h_prepreg + 5*h_core + 12*t_cu_pri` in `TwoLevelMutliTrack_design.m`
- `calculate_Pcond_2LMT.m`: fixed 8-primary / 2-secondary formula; ignores `stackup` (`'5Layer'` is a dead label)
- `calc_Llk.m`: hand-written 11-gap stack (S, shield, 8 P, shield, S)

Other modelling issues:
- Inductor counted as 2× for loss but 1× for area in the `calcEfficiency_v5` call of the design script.
- Flux density hard-coded as `Vin/(8*8*f*Ac)` in all `calculateVIRT_*_2LMT`; ignores `params.np`.
- Two different F80 loss models: the transformer uses `get_core_material('F80')` (with f^alpha), the inductor uses the hard-coded
  1 MHz fit. Center-post eddy-loss formula differs between files: `(π/4)·Ac` vs `(π/6)·b²`.
- No check that transformer leakage < Lr before `design_params.Lr = Llk - calc_Llk(...)`.
- `k_area` sweep: crashes if no design is ≤ 4 in²; doesn't filter on `is_valid`.
- `optimize_VIRT_2LMT` overwrites its `x0` argument; its square branch passes `(b,w,k_area)` to
  `calculateVIRT_square(a,b,w,...)`; the design script's square branch uses undefined `T_tx`.
- `calculate_Lr_EQ`: `fzero` bracket [1 µm, 10 mm] with the sign check commented out; can throw.
- `calculateVIRT_roundEQ_2LMT.m` line 5: typo `k_areSa`, so it doesn't run.
- `calculate_Lr` (rectangle shape): Rdc outer size uses `+w_winding` instead of `+2*w_winding`.
- Switch model assumptions to confirm: conduction `0.5·Ipk²·R` (half-bridge ≈ `Ipk²·R/4`);
  turn-off at peak resonant current rather than magnetizing current.
- `analyzePriSwitches_v5` calls a nonexistent `resolveRankJJList_v2` if `compare_list` is non-empty.
- `get_core_material('ML91S')` has no `Bsat` field.
- `one_device_check_epc.m`: `t_f` scaled by 1e-9 twice, so turn-off loss ≈ 0.
- Commented blocks in the design script are mislabelled (the "Round EQ" block calls the PQ optimizer, etc.).

## Repository / git workflow

- `main` = merged, current code (switch side from Qijia + magnetics from Julia, merged 2026-10-27).
- Branch naming: describe the work, lowercase-hyphenated, optional owner prefix,
  e.g. `jestrin/transformer-layer-reduction`. Avoid date-only names (old `jestrin_MM_DD_YYYY` branches).
- Active branches:
  - `jestrin/transformer-layer-reduction`: reducing transformer PCB layer count.
  - `cleanup-old-files`: empty, reserved for the file cleanup (merge `main` into it before resuming).
- Tag `pre-cleanup-2026-10-27` (pushed) = snapshot before any file removal. Restore a file with
  `git restore --source=pre-cleanup-2026-10-27 "<path>"`.
- Git notes describe the 16 commits merged on 2026-10-27; fetch them with
  `git fetch origin refs/notes/commits:refs/notes/commits`, push with `git push origin refs/notes/commits`.
- Repo lives inside Dropbox; paths contain spaces, so always quote them.

## Conventions

- **Edit files in place; do not create new `_vN` copies.** Use branches for experiments and commits for history.
- Qijia owns the switch-side files and device sheets (`Switch Functions/`, `Selected Switch Model/`,
  `one_device_check*.m`). Coordinate before changing or deleting them.
- Don't commit generated files: `.asv` autosaves, `sweep_results.mat`.
- When testing in MATLAB, run `restoredefaultpath` first so deleted or stale files can't mask errors.

## Planned file cleanup (paused)

Process: on `cleanup-old-files`: `git rm` files → add `*.asv` and `sweep_results.mat` to `.gitignore` →
test both entry scripts in MATLAB → commit → `git merge --ff-only` into `main` → push.
Agree on the list with Qijia first.
- Safe to delete: `old code/`; `main.m`; `main_only_switch.m`; `TwoLevelMutliTrack_v2.m`;
  `compare_topology*.m`; `analyzePriSwitches` v1–v4; `analyzeSecSwitches.m`; `calcEfficiency` v1–v4;
  `packOneLocal_v2`; `pickBestDevices_v2`; `try_with_functions*.m`; `compareParetoFreq.m`; `designLLC.m`;
  `test_designLLC.m`; old transformer model (`optimize_VIRT.m`, `calculateVIRT_round.m`,
  `calculate_VIRT_design.m`, `VIRT3Dfigure.m`, `test_*`, SPVIRT files, `calculate_Rdc.m`, `calc_Cps.m`);
  `optimize_Lr.m`, `calculate_Lr.m`, `Lr_AcWa.m` (unfinished); `.asv` files; `sweep_results.mat`.
- Keep: alternative core shapes, `Fairrite80.m`, the Ln sweep (needs fixing), Qijia's tools.
- Data sheets in use: `GaN Data tf.xlsx`, `SiC Data tf.xlsx`, `Sec Data Both.xlsx` (+ `GaN_05_12_26.xlsx`
  for the Ln sweep). Check with Qijia before removing the others.

## Working with Claude on this repo

- Claude has not been running MATLAB; the user runs scripts and reports errors and numbers.
  (Check whether `matlab -batch` is available to let Claude test changes.)
- Don't switch branches while Claude is reading or editing.
- For larger changes, ask for a plan first; commit after each working step.

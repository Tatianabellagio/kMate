# Savio HPC — submission reference for this project

Created 2026-05-28 during the Carnegie bsephc → Savio migration. Source of
truth on where jobs can land, what accounts/QoS to use, and where data
should live on Savio.

Upstream docs: https://docs-research-it.berkeley.edu/services/high-performance-computing/

---

## 1. Storage layout (use this everywhere)

| What | Path | Quota | Backup | Notes |
|---|---|---|---|---|
| `$HOME` | `/global/home/users/tbellg` | 50 GB | yes | Code, configs, binaries only. **Not for data.** |
| `$SCRATCH` (user) | `/global/scratch/users/tbellg` | soft 12 TB | no | Lustre. Active compute scratch. Project root lives here: `…/hapfire_sv/`. |
| Lab condo scratch | `/global/scratch/projects/co_moilab/` | 200 GB | no | Shared lab space (currently empty for us). |
| Lab FCA scratch | `/global/scratch/projects/fc_moilab/` | — | no | Shared lab space (in active use by lab; check before writing). |
| Convenience symlink | `/global/home/users/tbellg/scratch` → `/global/scratch/users/tbellg` | — | — | Pre-existing. Both paths work. |

**Rule of thumb:** project root is `/global/scratch/users/tbellg/hapfire_sv/`. All
runtime data, intermediate `.npz`, k-mer indexes, FASTQs, result TSVs go under
`$SCRATCH`. Only the code and the conda env live in `$HOME`.

Conda env on Savio: `/global/home/users/tbellg/miniforge3/envs/hapfm/` (already
exists — replaces the Carnegie path `/home/tbellagio/miniforge3/envs/hapfm/`).

---

## 2. Accounts the user has (from `sacctmgr show association`)

| Account | What it is | Partitions accessible | Notes |
|---|---|---|---|
| `co_moilab` | Moi-lab **condo** (purchased nodes) | All Savio partitions | Native QoS only on the condo's purchased hardware; everything else via `savio_lowprio` (preemptible). |
| `ac_moilab` | Moi-lab **faculty allocation** (paid SU pool) | All Savio partitions | Normal/debug QoS. Burns SUs from the lab's allocation. |
| `fc_moilab` | Moi-lab **FCA** (free annual allowance) | All Savio partitions | Normal/debug QoS. Burns the lab's free annual SUs. |
| `ucb` | UCB campus account | `ood-inter` only | Open OnDemand interactive nodes. |

The lab's condo purchased nodes (from `sacctmgr`):
- `savio4_htc` → QoS `moilab_htc4_normal` (NOT preempted; this is *our* hardware)
- `savio4_gpu` → QoS `a5k_gpu4_normal` (via `ac_moilab`; A5000 GPUs)
- `savio3_gpu` → QoS `a40_gpu3_normal`, `gtx2080_gpu3_normal` (via `ac_moilab`)

Everywhere else, the condo gets `savio_lowprio` only (preemptible — jobs can be
killed when the owning lab needs the node back).

---

## 3. Partitions on Savio (this user can reach all of them)

### Savio4 (newest, 2022+; scheduled **per-core** — request only what you need)
| Partition | Hardware | Cores/node | RAM/node | GPUs | Best for |
|---|---|---|---|---|---|
| `savio4_htc` | Xeon Gold 6330 | 56 | 256 or 512 GB | — | **Default for our CPU jobs.** Lab condo hardware here. |
| `savio4_gpu` | Xeon Gold 6326 / EPYC 9554 | 32–64 | 512–792 GB | 8× A5000 or 8× L40 | GPU jobs. Lab condo hardware here too. |

### Savio3 (2018+; mostly **per-node** unless noted)
| Partition | Hardware | Cores/node | RAM/node | GPUs | Best for |
|---|---|---|---|---|---|
| `savio3` | Skylake | 32–40 | 96 GB | — | General CPU, per-node. |
| `savio3_htc` | Skylake | 40 | 384 GB | — | Per-core CPU. |
| `savio3_bigmem` | Skylake | 32 | 384 GB | — | Memory-heavy per-node. |
| `savio3_xlmem` | Skylake | 32 | **1.5 TB** | — | Whole-genome assembly, huge in-memory ops. |
| `savio3_gpu` | Skylake | 8–64 | 65–515 GB | V100 / TITAN / 2080Ti / A40 | GPU; per-core. |

### Savio2 (2015; per-node; cheapest SU/hr but oldest)
| Partition | Hardware | Cores/node | RAM/node | GPUs | Best for |
|---|---|---|---|---|---|
| `savio2` | Haswell 12c×2 | 24 | 64 GB | — | Cheap CPU per-node. |
| `savio2_htc` | Haswell 12c (3.4 GHz) | 12 | 128 GB | — | Per-core, supports `savio_long` (10-day) QoS. |
| `savio2_bigmem` | Haswell | 24 | 128 GB | — | Per-node mem. |
| `savio2_1080ti` | — | 8 | 64 GB | 4× GTX 1080Ti | Cheap GPU. |
| `savio2_knl` | Xeon Phi (KNL) | 64 | 192 GB | — | Specialty — usually not what you want. |

### Other
| Partition | Notes |
|---|---|
| `ood-inter` | Open OnDemand interactive sessions (use `--account=ucb --qos=ood_interactive`). |
| `vector` / `bair` / `cloud4` | Other lab condos / specialty — not ours, ignore. |

---

## 4. QoS / submission recipes for *this* project

Pick by what you actually need. The recipes are ordered by what to try first.

### A. Default CPU job — use the lab's condo (no SU burn, no preemption)
```bash
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8       # savio4_htc is per-core — request only what you need
#SBATCH --mem=64G               # or omit to get default per-core RAM
#SBATCH --time=24:00:00         # condo QoS time limit — verify with `sacctmgr show qos moilab_htc4_normal`
```

### B. Large fan-out / opportunistic — condo lowprio on any partition
Use when you want to flood the cluster with hundreds of array tasks and you
don't mind some getting preempted (add `--requeue` so SLURM resubmits them).
```bash
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc     # or savio3, savio3_htc, savio2_htc, …
#SBATCH --qos=savio_lowprio
#SBATCH --requeue
#SBATCH --time=24:00:00
```
**Per-core partitions** (`savio4_htc`, `savio3_htc`, `savio2_htc`, `savio3_gpu`, `savio4_gpu`):
request exact `--cpus-per-task` + `--mem`. **Per-node partitions** (`savio2`,
`savio3`, `savio3_bigmem`, `savio3_xlmem`, `savio2_bigmem`, `savio2_knl`,
`savio2_1080ti`): you get the whole node, charged per-node.

### C. Time-sensitive runs — burn FCA or faculty-allocation SUs
Use when condo queue is backed up and you need it now.
```bash
#SBATCH --account=fc_moilab         # or ac_moilab
#SBATCH --partition=savio4_htc      # or savio2_htc / savio3_htc for cheaper SUs
#SBATCH --qos=savio_normal          # 72h max, up to 24 nodes
```

### D. Long jobs (>3 days, ≤10 days)
Only on `savio2_htc`, max 4 cores per job.
```bash
#SBATCH --account=fc_moilab         # or ac_moilab; condo lowprio also OK
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_long
#SBATCH --cpus-per-task=4
#SBATCH --time=10-00:00:00
```

### E. Debug / quick sanity checks (≤3h, ≤4 nodes)
```bash
#SBATCH --account=fc_moilab
#SBATCH --partition=savio4_htc      # or any partition
#SBATCH --qos=savio_debug
#SBATCH --time=3:00:00
```

### F. GPU on lab condo (A5000s on savio4_gpu)
```bash
#SBATCH --account=ac_moilab         # condo GPU access lives under ac_moilab
#SBATCH --partition=savio4_gpu
#SBATCH --qos=a5k_gpu4_normal
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:A5000:1
```

### Map of our existing job patterns → recommended Savio recipe

| Old Carnegie pattern | Savio recipe to use |
|---|---|
| Per-sample kMate (`run_site_array_perchrom.sh` — 8 CPU, 64 GB, hours) | **B** (lowprio, savio4_htc) for the 2,415-sample fan-out. Add `--requeue`. |
| `build_cn_full_*` / `build_cn_var_*` builders | **A** (condo, savio4_htc) — handful of jobs, fast turnaround matters. |
| SEEDMIX runs (`run_seedmix_*` — 8 CPU, 80 GB) | **A** (condo, savio4_htc). |
| PanGenie index / haploidize / merge VCFs | **A** or **C** on `savio4_htc`. |
| Anything needing >256 GB RAM | `savio3_bigmem` or `savio3_xlmem` (no per-node memory partition on savio4). |

---

## 5. Submission etiquette / gotchas

- **Per-core vs per-node:** on savio4_htc you pay per-core. On savio2/savio3 you
  pay per-node (so requesting `--cpus-per-task=4` on `savio3` still costs you all
  32 cores). Use savio4_htc and savio3_htc by default to avoid silent waste.
- **`savio_lowprio` jobs get preempted.** Always add `#SBATCH --requeue` so SLURM
  resubmits them. Driver scripts should be idempotent (skip already-done outputs).
- **No more `--partition=bse`.** That was Carnegie's bsephc cluster. Replace with
  `savio4_htc` (CPU default) or the right alternative from §4.
- **Log paths:** Carnegie scripts pointed `#SBATCH --output=` at absolute
  `/carnegie/...` paths. On Savio use relative paths (`logs/%x_%A_%a.out`) so
  the script works from any clone location.
- **Quotas:** `$HOME` is 50 GB and backed up — don't put data there. Use
  `$SCRATCH` (`/global/scratch/users/tbellg/...`). Lab condo space at
  `/global/scratch/projects/co_moilab/` is 200 GB shared.
- **Time limits:** check the QoS you're using with
  `sacctmgr show qos <qos_name> format=Name,MaxWall,MaxNodes,Priority`.

---

## 6. Useful inspection commands

```bash
# What can I submit to?
sacctmgr show association user=tbellg format=Account,Partition,QOS%50

# What's the time/node limit on a QoS?
sacctmgr show qos format=Name%24,MaxWall,MaxNodes,Priority | grep -E 'moilab|savio_'

# What's idle right now?
sinfo -o "%P %a %l %D %t" | grep -E 'idle|mix'

# My pending jobs and why they're queued
squeue -u tbellg --start
squeue -u tbellg -t PD -o "%i %P %j %r"
```

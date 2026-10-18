#!/usr/bin/env python
"""Submit every pool and fit that is not done yet, then the scoring.

Usage: python submit.py [FIGURE ...] [--dry-run]     (figures as named in runs.py; default all)

One Slurm job per missing pool (simulate.py), one per missing fit (fit.py, after its pool), and score.py after all of
them. A pool or fit is done when its directory has a DONE file.
"""
import os, sys, subprocess
import config as C
from runs import all_fits


def pool_done(pool_id):
    return os.path.exists(f"{C.POOLS}/{pool_id}/DONE")


def fit_done(fit_id):
    return os.path.exists(f"{C.FITS}/{fit_id}/DONE")


def sbatch(name, resources, command, after=(), dry=False, kind="afterok", where=C.SBATCH):
    dep = [f"--dependency={kind}:{':'.join(after)}"] if after else []
    cmd = ["sbatch", "--parsable", "--requeue", f"--job-name={name}", f"--output={C.LOGS}/{name}_%j.out",
           *where, *resources, *dep, f"--chdir={C.HERE}", f"--wrap={command}"]
    if dry:
        print(" ".join(cmd[:-1]), repr(command))
        return "DRY"
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip().split(";")[0]


def main(args):
    dry = "--dry-run" in args
    figures = [a for a in args if not a.startswith("--")] or None
    os.makedirs(C.LOGS, exist_ok=True)
    fits = all_fits(figures)
    pool_job, fit_jobs = {}, []
    for fit in fits.values():
        pid = fit.pool.id
        if pid not in pool_job:
            pool_job[pid] = None if pool_done(pid) else \
                sbatch(f"sim_{pid}", C.SIM_RES, f"{C.PYTHON} -u simulate.py {pid}", dry=dry)
        if not fit_done(fit.id):
            after = [pool_job[pid]] if pool_job[pid] else []
            where, res = (C.TIMED_SBATCH, C.TIMED_RES) if fit.timed else (C.SBATCH, C.FIT_RES)
            fit_jobs.append(sbatch(f"fit_{fit.id}", res, f"{C.PYTHON} -u fit.py {fit.id}", after, dry, where=where))
    n_sim = sum(j is not None for j in pool_job.values())
    print(f"{n_sim} of {len(pool_job)} pools and {len(fit_jobs)} of {len(fits)} fits submitted")
    if fit_jobs:
        print("score ->", sbatch("score", C.SCORE_RES, f"{C.PYTHON} -u score.py", fit_jobs, dry, kind="afterany"))


if __name__ == "__main__":
    main(sys.argv[1:])

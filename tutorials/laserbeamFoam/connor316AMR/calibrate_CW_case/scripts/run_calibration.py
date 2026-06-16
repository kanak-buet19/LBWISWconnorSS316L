#!/usr/bin/env python3
"""CW single-track 316L melt-pool calibration orchestrator (Bayesian).

Drives a Bayesian-optimization loop over the shared thermophysical parameters:
seed the search with a Latin-Hypercube design, then let Optuna's TPE surrogate
propose where the melt-pool error is likely lower. Candidates run as a parallel
swarm of laserbeamFoam sims, are scored against Hofmann (2026) experimental
depth/width, and the optimizer is told each result so the next proposals
improve. Proposals stream until the evaluation budget is spent.

Run via ./Allrun (local) or job_calibration.sh (HPC). All knobs live in
calibration_config.json; cores can be overridden by env CALIB_TOTAL_CORES /
CALIB_CORES_PER_SIM.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import caselib              # noqa: E402
import evaluate as ev       # noqa: E402
import plots                # noqa: E402
from optimizer import BOOptimizer  # noqa: E402

CONFIG = ROOT / "calibration_config.json"
TEMPLATE = ROOT / "template_case"
RUNS = ROOT / "runs"
RESULTS = ROOT / "results"
PENALTY = 9.99


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", flush=True)


def pspec_keys(pspec: dict) -> list[str]:
    return [k for k in pspec if not k.startswith("_")]


def run_command(case_dir: Path) -> str:
    """Per-sim shell command. Native OpenFOAM by default (local: OF sourced in the
    shell, exactly like cw_scantrack's Allrun -> Allrun_long). Set
    CALIB_USE_APPTAINER=1 (HPC, see job_calibration.sh) to run inside the image."""
    python = os.environ.get("PYTHON", sys.executable)
    if os.environ.get("CALIB_USE_APPTAINER", "0") == "1":
        of_image = os.environ.get("OF2506_IMAGE", str(Path.home() / "openfoam-dev_2506.sif"))
        of_user = os.environ.get("OF2506_USER", os.environ.get("USER", "user"))
        return (
            f"apptainer exec --cleanenv --env USER={of_user} {of_image} "
            f"bash -lc \"source /openfoam/bash.rc && export FOAM_SIGFPE=0 && "
            f"export PYTHON='{python}' && export MPLCONFIGDIR=/tmp && "
            f"cd '{case_dir}' && ./Allrun_long > log.run 2>&1\""
        )
    return (
        f"export FOAM_SIGFPE=0 && export PYTHON='{python}' && "
        f"export MPLCONFIGDIR=/tmp && cd '{case_dir}' && ./Allrun_long > log.run 2>&1"
    )


def sim_progress(case_dir: Path) -> tuple[float, float | None] | None:
    """Latest (Time, deltaT) from the solver log tail, or None if not started."""
    log_file = Path(case_dir) / "log.laserbeamFoam"
    if not log_file.exists():
        return None
    try:
        tail = log_file.read_bytes()[-40000:].decode("utf-8", "ignore")
    except OSError:
        return None
    t = dt = None
    for line in tail.splitlines():
        if line.startswith("Time = "):
            try:
                t = float(line[7:].strip())
            except ValueError:
                pass
        elif line.startswith("deltaT = "):
            try:
                dt = float(line[9:].strip())
            except ValueError:
                pass
    return (t, dt) if t is not None else None


def fmt_params(p: dict) -> str:
    parts = []
    for k, v in p.items():
        abbr = k[:5]
        parts.append(f"{abbr}={v:.3g}")
    return " ".join(parts)


def strip_heavy(case_dir: Path) -> int:
    """Delete the bulky run artifacts (decomposed + reconstructed time dirs and
    VTK) from a finished case, keeping the tiny post-processing CSV, logs and
    case_build.json. Returns number of dirs removed. Safe to call only after the
    sim has terminated."""
    cd = Path(case_dir)
    if not cd.exists():
        return 0
    removed = 0
    targets = list(cd.glob("processor*"))
    for name in ("VTK", "VTKs", "dynamicCode"):
        d = cd / name
        if d.is_dir():
            targets.append(d)
    # reconstructed root time dirs (numeric name, keep "0" = initial fields)
    for p in cd.iterdir():
        if (p.is_dir() and p.name != "0"
                and re.fullmatch(r"[0-9][0-9.eE+\-]*", p.name)):
            targets.append(p)
    for d in targets:
        shutil.rmtree(d, ignore_errors=True)
        removed += 1
    return removed


# --------------------------------------------------------------------------- #
# Job model
# --------------------------------------------------------------------------- #
class Job:
    def __init__(self, cand_id: int, params: dict, case_cfg: dict,
                 case_dir: Path, geom: dict, trial=None):
        self.cand_id = cand_id
        self.params = params
        self.case_cfg = case_cfg
        self.case_dir = case_dir
        self.trial = trial   # optuna Trial for sims proposed this run, else None
        self.surface_y_um = geom["surface_y_um"]
        self.proc: subprocess.Popen | None = None
        self.started = None
        self.status = "queued"   # queued|running|done|aborted|failed|cancelled
        self.note = ""
        self.result: dict | None = None

    @property
    def name(self) -> str:
        return self.case_cfg["name"]

    def start(self) -> None:
        self.started = time.time()
        self.status = "running"
        self.proc = subprocess.Popen(
            run_command(self.case_dir), shell=True, executable="/bin/bash",
            start_new_session=True,
        )

    def kill(self) -> None:
        if self.proc and self.proc.poll() is None:
            try:
                os.killpg(os.getpgid(self.proc.pid), signal.SIGTERM)
                time.sleep(2)
                if self.proc.poll() is None:
                    os.killpg(os.getpgid(self.proc.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass

    def evaluate(self, stability_tol: float) -> dict:
        return ev.evaluate(self.case_dir, self.case_cfg["exp_depth_um"],
                           self.case_cfg["exp_width_um"], stability_tol)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
class Calibrator:
    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.cases = cfg["cases"]
        self.geom = cfg["geometry"]
        self.control = cfg["control"]
        self.early = cfg["earlyAbort"]
        self.stability_tol = self.early.get("stabilityTol", 0.10)
        self.param_names = pspec_keys(cfg["parameters"])

        # disk saver: keep only the single best candidate's full data (VTK +
        # processor*). When a better candidate completes, the previous best is
        # stripped; every non-best candidate is stripped as soon as it finishes.
        clean = cfg.get("cleanup", {})
        self.cleanup_enabled = bool(clean.get("enabled", True))
        self.best_kept_cand: int | None = None
        self.best_kept_obj = float("inf")

        ex = cfg["execution"]
        self.total_cores = int(os.environ.get("CALIB_TOTAL_CORES", ex["totalCores"]))
        self.cores_per_sim = int(os.environ.get("CALIB_CORES_PER_SIM", ex["coresPerSim"]))
        self.max_parallel = max(1, self.total_cores // self.cores_per_sim)
        self.poll = ex.get("pollIntervalSec", 10)
        self.timeout = ex.get("perSimTimeoutSec", 14400)

        # Bayesian optimization: total evaluation budget + LHS init size
        opt = cfg["optimizer"]
        self.budget = int(opt["nSamples"])
        self.n_init = int(opt.get("initSamples", min(16, self.budget)))
        self.seed = int(opt.get("seed", 42))
        self.opt = BOOptimizer(cfg["parameters"], self.param_names,
                               self.n_init, self.seed)

        self.records: list[dict] = []   # grows as the optimizer proposes
        self.queue: list[Job] = []
        self.running: list[Job] = []
        self.resume_ids: list[int] = []  # incomplete prior candidates to rerun
        self.done_count = 0
        self.cancelled = False
        self.started_at = datetime.now().isoformat(timespec="seconds")

        self._load_prior()
        # fresh start: lay down the Latin-Hypercube initial design
        if not self.records:
            self.opt.seed_lhs(self.n_init, self.seed)

    # -- resume helpers ---------------------------------------------------- #
    def _cand_dir(self, cid: int) -> Path:
        return RUNS / f"cand_{cid:02d}"

    def _result_path(self, cid: int) -> Path:
        return self._cand_dir(cid) / "result.json"

    @staticmethod
    def _case_has_progress(case_dir: Path) -> bool:
        """True if a sim has written past t=0 (resumable, mirrors Allrun_long)."""
        p0 = case_dir / "processor0"
        if not p0.is_dir():
            return False
        for d in p0.iterdir():
            if d.is_dir() and d.name != "0" and re.fullmatch(r"[0-9][0-9.eE+\-]*", d.name):
                return True
        return False

    def _cand_has_vtk(self, cid: int) -> bool:
        base = self._cand_dir(cid)
        return base.exists() and any(
            (c / "VTK").is_dir() for c in base.iterdir() if c.is_dir())

    def _persist_candidate(self, cid: int) -> None:
        """Write a small result.json so a resubmitted job can skip/restore this
        candidate. Survives strip_heavy (which only deletes VTK/processor/time
        dirs)."""
        base = self._cand_dir(cid)
        if not base.exists():
            return
        rec = self.records[cid]
        payload = {
            "id": cid, "params": rec["params"], "objective": rec["objective"],
            "status": rec["status"], "disposed": bool(rec.get("_disposed")),
            "cases": {n: {k: v for k, v in c.items() if k != "series"}
                      for n, c in rec["cases"].items()},
        }
        try:
            self._result_path(cid).write_text(
                json.dumps(payload, indent=2, default=str))
        except OSError:
            pass

    def _new_record(self, params: dict) -> int:
        """Append a candidate record; its id is its position. Returns the id."""
        cid = len(self.records)
        self.records.append({
            "id": cid, "params": params, "params_pretty": fmt_params(params),
            "cases": {}, "objective": float("nan"), "status": "pending",
        })
        return cid

    def _load_prior(self) -> None:
        """Reload candidates from a previous job so we don't restart from
        scratch: completed ones are replayed into the surrogate and skipped;
        in-flight ones are queued for resume. A clean runs/ folder starts
        fresh. Candidate dirs are contiguous cand_00, cand_01, ..."""
        if not RUNS.exists():
            return
        loaded = completed = 0
        cid = 0
        while (RUNS / f"cand_{cid:02d}").exists():
            rp = self._result_path(cid)
            data = None
            if rp.exists():
                try:
                    data = json.loads(rp.read_text())
                except (OSError, ValueError):
                    data = None
            params = self._cand_params(cid, data)
            if params is None:
                break  # can't reconstruct -> stop (treat as the frontier)
            self._new_record(params)
            rec = self.records[cid]
            loaded += 1
            if data and data.get("status") == "complete":
                rec["objective"] = data.get("objective", float("nan"))
                rec["status"] = "complete"
                rec["cases"] = data.get("cases", {})
                rec["_told"] = True
                if data.get("disposed"):
                    rec["_disposed"] = True
                self.done_count += sum(
                    1 for c in rec["cases"].values()
                    if c.get("status") in ("done", "aborted", "failed", "skipped"))
                obj = rec["objective"]
                val = obj if (obj == obj and obj < PENALTY) else PENALTY
                self.opt.replay(params, val)       # teach surrogate
                completed += 1
                if obj == obj and obj < PENALTY and self._cand_has_vtk(cid):
                    if obj < self.best_kept_obj:
                        self.best_kept_obj, self.best_kept_cand = obj, cid
            else:
                self.resume_ids.append(cid)        # rerun / continue this one
            cid += 1
        if loaded:
            log(f"RESUME: restored {loaded} candidate(s) "
                f"({completed} complete, {len(self.resume_ids)} to resume); "
                f"best so far {self._best_objective()}")

    def _cand_params(self, cid: int, data: dict | None) -> dict | None:
        """Recover a candidate's params from result.json or case_build.json."""
        if data and data.get("params"):
            return data["params"]
        for cdir in (RUNS / f"cand_{cid:02d}").glob("*/case_build.json"):
            try:
                return json.loads(cdir.read_text())["params"]
            except (OSError, ValueError, KeyError):
                continue
        return None

    # -- scheduling -------------------------------------------------------- #
    def build_job(self, cand_id: int, case_cfg: dict, trial=None) -> Job:
        params = self.records[cand_id]["params"]
        case_dir = RUNS / f"cand_{cand_id:02d}" / case_cfg["name"]
        # resume an in-flight sim instead of wiping it; Allrun_long continues
        # from latestTime when processor dirs past t=0 already exist
        if self._case_has_progress(case_dir) and (case_dir / "case_build.json").exists():
            rec = json.loads((case_dir / "case_build.json").read_text())
            log(f"cand {cand_id:02d} {case_cfg['name']} RESUME (existing run)")
        else:
            rec = caselib.build_case(TEMPLATE, case_dir, params, case_cfg,
                                     self.geom, self.control, self.cores_per_sim)
        return Job(cand_id, params, case_cfg, case_dir,
                   {"surface_y_um": rec["surface_y_um"]}, trial=trial)

    def _enqueue_resumes(self) -> None:
        """Re-queue candidates that a prior job left in flight."""
        for cid in self.resume_ids:
            self.queue.append(self.build_job(cid, self.cases[0]))

    def _maybe_propose(self) -> None:
        """Keep the pipeline full: while a slot is free and budget remains, ask
        the optimizer for the next parameter set and queue it."""
        while (len(self.records) < self.budget
               and len(self.queue) + len(self.running) < self.max_parallel):
            params, trial = self.opt.ask()
            cid = self._new_record(params)
            phase = "INIT" if cid < self.n_init else "BO"
            log(f"cand {cid:02d} PROPOSE [{phase}] {fmt_params(params)}")
            self.queue.append(self.build_job(cid, self.cases[0], trial=trial))

    def maybe_enqueue_next(self, job: Job) -> None:
        try:
            current_idx = next(i for i, c in enumerate(self.cases) if c["name"] == job.name)
        except StopIteration:
            return
        next_idx = current_idx + 1
        if next_idx >= len(self.cases):
            return

        rec = self.records[job.cand_id]
        # Gate on case[0] error: if first case badly off, skip ALL remaining cases
        if self.early.get("skipSecondCaseIfFirstBad", True):
            first_err = rec["cases"].get(self.cases[0]["name"], {}).get("case_error")
            if (first_err is not None and first_err == first_err
                    and first_err > self.early["errorThreshold"]):
                for case_cfg in self.cases[next_idx:]:
                    log(f"cand {job.cand_id:02d} skip {case_cfg['name']} "
                        f"(case0 err {first_err*100:.0f}% > "
                        f"{self.early['errorThreshold']*100:.0f}%)")
                    rec["cases"][case_cfg["name"]] = {"status": "skipped"}
                self._update_objective(job.cand_id)
                return

        self.queue.append(self.build_job(job.cand_id, self.cases[next_idx]))

    # -- main loop --------------------------------------------------------- #
    def run(self) -> None:
        RESULTS.mkdir(parents=True, exist_ok=True)
        log(f"=== CW calibration (Bayesian): budget {self.budget} sims, "
            f"{self.n_init} LHS init, then TPE ===")
        log(f"cores: total={self.total_cores} perSim={self.cores_per_sim} "
            f"-> maxParallel={self.max_parallel}")
        self._enqueue_resumes()

        try:
            while (len(self.records) < self.budget or self.queue or self.running):
                self._maybe_propose()
                self._fill()
                self._poll_running()
                self._dispose_completed()
                self._print_progress()
                self._write_status()
                time.sleep(self.poll)
        except KeyboardInterrupt:
            # ignore further Ctrl-C so shutdown/finalize always completes
            signal.signal(signal.SIGINT, signal.SIG_IGN)
            self.cancelled = True
            log("!! interrupted - killing running sims and saving partial results")
            for job in self.running:
                job.kill()
                job.status = "cancelled"
                self._record_job(job, evaluate=False)

        self._finalize()

    def _print_progress(self) -> None:
        if not self.running:
            return
        end = self.control["endTime"]
        parts = []
        for job in self.running:
            tag = f"c{job.cand_id:02d}/{job.name.split('_')[-1]}"
            pr = sim_progress(job.case_dir)
            if pr is None:
                parts.append(f"{tag} init")
            else:
                t, dt = pr
                pct = 100.0 * t / end if end else 0.0
                dts = f"{dt:.1e}" if dt is not None else "?"
                parts.append(f"{tag} {pct:4.1f}% t={t:.2e} dt={dts}")
        log("   .. " + " | ".join(parts))

    def _fill(self) -> None:
        while self.queue and len(self.running) < self.max_parallel:
            job = self.queue.pop(0)
            job.start()
            self.running.append(job)
            best = self._best_objective()
            log(f"cand {job.cand_id:02d} {job.name} START [{fmt_params(job.params)}] "
                f"| running {len(self.running)}/{self.max_parallel} "
                f"| done {self.done_count} | best {best}")

    def _poll_running(self) -> None:
        still = []
        for job in self.running:
            rc = job.proc.poll()
            # early abort
            if rc is None and self.early.get("enabled", True):
                abort, why = ev.should_abort(
                    job.case_dir, job.case_cfg["exp_depth_um"],
                    job.case_cfg["exp_width_um"], self.early["errorThreshold"],
                    self.stability_tol, self.early.get("minPoints", 2))
                if abort:
                    job.kill()
                    job.status = "aborted"
                    job.note = why
                    log(f"cand {job.cand_id:02d} {job.name} ABORT - {why}")
                    self._record_job(job, evaluate=True)
                    self.maybe_enqueue_next(job)
                    self.done_count += 1
                    continue
            # timeout
            if rc is None and job.started and (time.time() - job.started) > self.timeout:
                job.kill()
                job.status = "failed"
                job.note = "timeout"
                log(f"cand {job.cand_id:02d} {job.name} TIMEOUT")
                self._record_job(job, evaluate=True)
                self.maybe_enqueue_next(job)
                self.done_count += 1
                continue
            if rc is None:
                still.append(job)
                continue
            # finished
            job.status = "done" if rc == 0 else "failed"
            if rc != 0:
                job.note = f"exit {rc}"
            self._record_job(job, evaluate=True)
            self._report_finish(job)
            self.maybe_enqueue_next(job)
            self.done_count += 1
        self.running = still

    def _record_job(self, job: Job, evaluate: bool) -> None:
        rec = self.records[job.cand_id]
        entry = {"status": job.status, "note": job.note}
        if evaluate:
            try:
                res = job.evaluate(self.stability_tol)
                entry.update({
                    "sim_depth_um": res["sim_depth_um"],
                    "sim_width_um": res["sim_width_um"],
                    "exp_depth_um": res["exp_depth_um"],
                    "exp_width_um": res["exp_width_um"],
                    "depth_err": res["depth_err"], "width_err": res["width_err"],
                    "case_error": res["case_error"],
                    "converged": res["converged"], "n_points": res["n_points"],
                    "series": res["series"],
                })
            except Exception as exc:
                entry["note"] = f"eval failed: {exc}"
        rec["cases"][job.name] = entry
        self._update_objective(job.cand_id)
        self._tell_if_complete(job)

    def _tell_if_complete(self, job: Job) -> None:
        """Feed a finished candidate's objective back to the surrogate so the
        next proposals improve. Tell the live trial if this sim was proposed
        this run, else replay (resumed sim)."""
        rec = self.records[job.cand_id]
        if rec["status"] != "complete" or rec.get("_told"):
            return
        rec["_told"] = True
        obj = rec["objective"]
        val = obj if (obj == obj and obj < PENALTY) else PENALTY
        try:
            if job.trial is not None:
                self.opt.tell(job.trial, val)
            else:
                self.opt.replay(rec["params"], val)
        except Exception as exc:   # never let optimizer bookkeeping kill the run
            log(f"cand {job.cand_id:02d} optimizer tell failed: {exc}")

    def _strip_candidate(self, cand_id: int) -> int:
        """Strip VTK + processor*/ + reconstructed time dirs from every case
        dir of a candidate. Tiny CSV/logs/build json are kept."""
        base = RUNS / f"cand_{cand_id:02d}"
        if not base.exists():
            return 0
        n = 0
        for cdir in sorted(base.iterdir()):
            if cdir.is_dir():
                n += strip_heavy(cdir)
        return n

    def _dispose_completed(self) -> None:
        """Keep only the best-so-far candidate's full data. When a candidate
        finishes: if it is the new best, strip the previous best; otherwise
        strip it now. Idempotent via the per-record '_disposed' flag."""
        if not self.cleanup_enabled:
            return
        for rec in self.records:
            if rec["status"] != "complete" or rec.get("_disposed"):
                continue
            cid = rec["id"]
            if cid == self.best_kept_cand:  # retained best, never re-examine
                continue
            obj = rec["objective"]
            valid = obj == obj and obj < PENALTY  # finite, real result
            if valid and obj < self.best_kept_obj:
                if self.best_kept_cand is not None:
                    m = self._strip_candidate(self.best_kept_cand)
                    log(f"cand {self.best_kept_cand:02d} STRIP (superseded by "
                        f"#{cid:02d} obj {obj:.3f}, {m} dirs)")
                    self.records[self.best_kept_cand]["_disposed"] = True
                    self._persist_candidate(self.best_kept_cand)
                self.best_kept_obj = obj
                self.best_kept_cand = cid
                log(f"cand {cid:02d} KEEP data (new best obj {obj:.3f})")
                # note: do NOT mark _disposed - this is the retained best
            else:
                m = self._strip_candidate(cid)
                rec["_disposed"] = True
                self._persist_candidate(cid)
                if m:
                    ostr = f"{obj:.3f}" if valid else "n/a"
                    bstr = (f"{self.best_kept_obj:.3f}"
                            if self.best_kept_obj < float("inf") else "n/a")
                    log(f"cand {cid:02d} STRIP ({m} dirs, obj {ostr} "
                        f"not best {bstr})")

    def _update_objective(self, cand_id: int) -> None:
        rec = self.records[cand_id]
        errs = [c.get("case_error") for c in rec["cases"].values()
                if c.get("case_error") == c.get("case_error")
                and c.get("case_error") is not None]
        n_expected = len(self.cases)
        n_have = len(rec["cases"])
        if errs:
            rec["objective"] = sum(errs) / len(errs)
        elif n_have >= n_expected:
            rec["objective"] = PENALTY
        # status summary
        was_complete = rec["status"] == "complete"
        if n_have >= n_expected or any(c.get("status") == "skipped"
                                       for c in rec["cases"].values()):
            rec["status"] = "complete"
        # persist on transition to complete so a resubmit can skip it
        if rec["status"] == "complete" and not was_complete:
            self._persist_candidate(cand_id)

    def _report_finish(self, job: Job) -> None:
        c = self.records[job.cand_id]["cases"].get(job.name, {})
        if "sim_depth_um" in c:
            d, w = c["sim_depth_um"], c["sim_width_um"]
            de = c["depth_err"] * 100 if c["depth_err"] == c["depth_err"] else float("nan")
            we = c["width_err"] * 100 if c["width_err"] == c["width_err"] else float("nan")
            conv = "conv" if c.get("converged") else "NOTconv"
            log(f"cand {job.cand_id:02d} {job.name} {job.status.upper()} "
                f"d={d:.1f}(exp{c['exp_depth_um']:.1f} {de:+.0f}%) "
                f"w={w:.1f}(exp{c['exp_width_um']:.1f} {we:+.0f}%) "
                f"caseErr {c['case_error']:.3f} [{conv}] | best {self._best_objective()}")
        else:
            log(f"cand {job.cand_id:02d} {job.name} {job.status.upper()} "
                f"({job.note}) | best {self._best_objective()}")

    def _best_record(self) -> dict | None:
        done = [r for r in self.records if r["objective"] == r["objective"]]
        return min(done, key=lambda r: r["objective"]) if done else None

    def _best_objective(self) -> str:
        b = self._best_record()
        return f"{b['objective']:.3f}(#{b['id']:02d})" if b else "n/a"

    # -- output ------------------------------------------------------------ #
    def _write_status(self) -> None:
        status = {
            "started_at": self.started_at,
            "updated_at": datetime.now().isoformat(timespec="seconds"),
            "cancelled": self.cancelled,
            "budget": self.budget,
            "candidates_proposed": len(self.records),
            "jobs_done": self.done_count,
            "jobs_running": len(self.running),
            "jobs_queued": len(self.queue),
            "running_now": [f"cand{j.cand_id:02d}/{j.name}" for j in self.running],
            "best": self._best_record(),
        }
        (RESULTS / "status.json").write_text(json.dumps(status, indent=2, default=str))

    def _finalize(self) -> None:
        (RESULTS / "plots").mkdir(parents=True, exist_ok=True)
        # strip heavy series from a compact copy for csv
        self._write_summary_csv()
        (RESULTS / "summary.json").write_text(
            json.dumps(self.records, indent=2, default=str))

        best = self._best_record()
        if best:
            bp = {
                "best_candidate_id": best["id"],
                "objective": best["objective"],
                "status": best["status"],
                "parameters": best["params"],
                "cases": {n: {k: v for k, v in c.items() if k != "series"}
                          for n, c in best["cases"].items()},
                "note": "Shared thermophysical params calibrated against both CW cases.",
            }
            (RESULTS / "best_params.json").write_text(json.dumps(bp, indent=2, default=str))
            log(f"BEST candidate #{best['id']} objective {best['objective']:.4f}")
            log(f"  params: {fmt_params(best['params'])}")

        plots.make_all(self.records, best, self.param_names, RESULTS / "plots")
        self._write_status()
        log(f"done. results in {RESULTS}")
        if self.cancelled:
            log("NOTE: run was cancelled - results are partial (see status.json).")

    def _write_summary_csv(self) -> None:
        cols = ["id", "status", "objective"] + self.param_names
        for case in self.cases:
            nm = case["name"]
            cols += [f"{nm}__status", f"{nm}__sim_depth_um", f"{nm}__sim_width_um",
                     f"{nm}__depth_err", f"{nm}__width_err", f"{nm}__case_error",
                     f"{nm}__converged"]
        with (RESULTS / "summary.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writeheader() if False else w.writerow(cols)
            for r in self.records:
                row = [r["id"], r["status"], r["objective"]]
                row += [r["params"][p] for p in self.param_names]
                for case in self.cases:
                    c = r["cases"].get(case["name"], {})
                    row += [c.get("status", ""), c.get("sim_depth_um", ""),
                            c.get("sim_width_um", ""), c.get("depth_err", ""),
                            c.get("width_err", ""), c.get("case_error", ""),
                            c.get("converged", "")]
                w.writerow(row)


def main() -> None:
    cfg = json.loads(CONFIG.read_text())
    Calibrator(cfg).run()


if __name__ == "__main__":
    main()

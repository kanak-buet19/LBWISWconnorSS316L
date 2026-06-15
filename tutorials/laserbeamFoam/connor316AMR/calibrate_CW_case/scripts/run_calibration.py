#!/usr/bin/env python3
"""CW single-track 316L melt-pool calibration orchestrator.

Samples thermophysical parameter sets (LHS + baseline), runs each candidate on
BOTH validation cases as a parallel swarm of laserbeamFoam sims, scores them
against Hofmann (2026) experimental melt-pool depth/width, and reports the best
shared parameter set.

Run via ./Allrun (local) or job_calibration.sh (HPC). All knobs live in
calibration_config.json; cores can be overridden by env CALIB_TOTAL_CORES /
CALIB_CORES_PER_SIM.
"""

from __future__ import annotations

import csv
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import caselib          # noqa: E402
import evaluate as ev   # noqa: E402
import plots            # noqa: E402

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


def sample_params(cfg: dict) -> list[dict]:
    """LHS samples in parameter ranges, optionally prepended with baseline."""
    pspec = cfg["parameters"]
    names = [k for k in pspec_keys(pspec)]
    lows = [pspec[k]["min"] for k in names]
    highs = [pspec[k]["max"] for k in names]

    opt = cfg["optimizer"]
    n_total = int(opt["nSamples"])
    include_baseline = bool(opt.get("includeBaseline", True))
    n_lhs = n_total - (1 if include_baseline else 0)

    from scipy.stats import qmc
    sampler = qmc.LatinHypercube(d=len(names), seed=opt.get("seed", 0))
    unit = sampler.random(n=n_lhs)
    scaled = qmc.scale(unit, lows, highs)

    candidates = []
    if include_baseline:
        candidates.append({k: pspec[k]["baseline"] for k in names})
    for row in scaled:
        candidates.append({k: float(v) for k, v in zip(names, row)})
    return candidates


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


# --------------------------------------------------------------------------- #
# Job model
# --------------------------------------------------------------------------- #
class Job:
    def __init__(self, cand_id: int, params: dict, case_cfg: dict,
                 case_dir: Path, geom: dict):
        self.cand_id = cand_id
        self.params = params
        self.case_cfg = case_cfg
        self.case_dir = case_dir
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

        ex = cfg["execution"]
        self.total_cores = int(os.environ.get("CALIB_TOTAL_CORES", ex["totalCores"]))
        self.cores_per_sim = int(os.environ.get("CALIB_CORES_PER_SIM", ex["coresPerSim"]))
        self.max_parallel = max(1, self.total_cores // self.cores_per_sim)
        self.poll = ex.get("pollIntervalSec", 10)
        self.timeout = ex.get("perSimTimeoutSec", 14400)

        self.candidates = sample_params(cfg)
        # per-candidate record
        self.records = [{
            "id": i, "params": p, "params_pretty": fmt_params(p),
            "cases": {}, "objective": float("nan"), "status": "pending",
        } for i, p in enumerate(self.candidates)]

        self.queue: list[Job] = []
        self.running: list[Job] = []
        self.done_count = 0
        self.cancelled = False
        self.started_at = datetime.now().isoformat(timespec="seconds")

    # -- scheduling -------------------------------------------------------- #
    def build_job(self, cand_id: int, case_cfg: dict) -> Job:
        case_dir = RUNS / f"cand_{cand_id:02d}" / case_cfg["name"]
        rec = caselib.build_case(TEMPLATE, case_dir, self.candidates[cand_id],
                                 case_cfg, self.geom, self.control,
                                 self.cores_per_sim)
        return Job(cand_id, self.candidates[cand_id], case_cfg, case_dir,
                   {"surface_y_um": rec["surface_y_um"]})

    def enqueue_first_cases(self) -> None:
        for cid in range(len(self.candidates)):
            self.queue.append(self.build_job(cid, self.cases[0]))

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
                return

        self.queue.append(self.build_job(job.cand_id, self.cases[next_idx]))

    # -- main loop --------------------------------------------------------- #
    def total_planned(self) -> int:
        # lower bound; second cases may be skipped
        return len(self.candidates) * len(self.cases)

    def run(self) -> None:
        RESULTS.mkdir(parents=True, exist_ok=True)
        log(f"=== CW calibration: {len(self.candidates)} candidates x "
            f"{len(self.cases)} cases ===")
        log(f"cores: total={self.total_cores} perSim={self.cores_per_sim} "
            f"-> maxParallel={self.max_parallel}")
        self.enqueue_first_cases()

        try:
            while self.queue or self.running:
                self._fill()
                self._poll_running()
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
        if n_have >= n_expected or any(c.get("status") == "skipped"
                                       for c in rec["cases"].values()):
            rec["status"] = "complete"

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
            "candidates_total": len(self.candidates),
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

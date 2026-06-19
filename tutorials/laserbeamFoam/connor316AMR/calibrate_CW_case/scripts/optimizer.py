"""Bayesian-optimization driver for CW melt-pool calibration.

Wraps Optuna's TPE sampler behind a small ask/tell API the orchestrator drives:

  1. Seed the search with `n_init` Latin-Hypercube points (broad exploration).
  2. After those, TPE fits a density-ratio surrogate over completed trials and
     proposes where the objective is likely lower (exploit + explore).
  3. `constant_liar` makes batch proposals (8-wide swarm) diversify instead of
     piling onto the same spot while siblings are still running.

The study is held in memory and rebuilt on resume by replaying the completed
trials the orchestrator reloads from disk (result.json), so no separate DB is
needed and the optimizer state always matches the sims that actually exist.
"""

from __future__ import annotations

import warnings

import optuna
from optuna.distributions import FloatDistribution
from optuna.exceptions import ExperimentalWarning
from optuna.trial import Trial, TrialState, create_trial

from scipy.stats import qmc

# TPE's multivariate/group/constant_liar knobs are flagged experimental but are
# stable in Optuna 4.x; silence the per-construction warnings to keep logs clean.
warnings.filterwarnings("ignore", category=ExperimentalWarning)


class BOOptimizer:
    def __init__(self, pspec: dict, names: list[str], n_init: int,
                 seed: int = 42) -> None:
        self.all_names = names
        self.fixed = {
            n: float(pspec[n]["baseline"])
            for n in names
            if pspec[n].get("fixed", False)
        }
        self.names = [n for n in names if n not in self.fixed]
        self.n_init = n_init
        self.dists = {
            n: FloatDistribution(float(pspec[n]["min"]), float(pspec[n]["max"]))
            for n in self.names
        }
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        self.study = optuna.create_study(
            direction="minimize",
            sampler=optuna.samplers.TPESampler(
                n_startup_trials=n_init, seed=seed,
                multivariate=True, group=True, constant_liar=True),
        )

    # -- seeding ----------------------------------------------------------- #
    def seed_lhs(self, n: int, seed: int = 42) -> None:
        """Queue n Latin-Hypercube points as the initial design. ask() returns
        these (FIFO) before the sampler proposes anything."""
        lows = [self.dists[n_].low for n_ in self.names]
        highs = [self.dists[n_].high for n_ in self.names]
        sampler = qmc.LatinHypercube(d=len(self.names), seed=seed)
        pts = qmc.scale(sampler.random(n=n), lows, highs)
        for row in pts:
            self.study.enqueue_trial(
                {n_: float(v) for n_, v in zip(self.names, row)}, skip_if_exists=False)

    # -- ask / tell -------------------------------------------------------- #
    def ask(self) -> tuple[dict, Trial]:
        trial = self.study.ask(self.dists)
        params = dict(self.fixed)
        params.update({n: float(trial.params[n]) for n in self.names})
        return {n: params[n] for n in self.all_names}, trial

    def tell(self, trial: Trial, value: float) -> None:
        self.study.tell(trial, value)

    def replay(self, params: dict, value: float) -> None:
        """Re-add a previously completed trial (used on resume) so the surrogate
        learns from sims run in an earlier job."""
        self.study.add_trial(create_trial(
            params={n: float(params[n]) for n in self.names},
            distributions=self.dists, value=float(value),
            state=TrialState.COMPLETE))

    # -- introspection ----------------------------------------------------- #
    def n_completed(self) -> int:
        return len([t for t in self.study.trials
                    if t.state == TrialState.COMPLETE])

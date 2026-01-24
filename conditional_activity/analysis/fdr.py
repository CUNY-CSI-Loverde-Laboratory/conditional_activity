"""
FDR correction for CONDACT Conditional Activity matrices (CONDACT-direction, >=11 transitions)

Assumes CONDACT already produced:
    - transition_times.pkl
    - Conditional_Activity_matrix.npy
    - Dihedral_Angle.csv (optional; only used to infer total_time if not passed)

Computes:
    - p-value matrix via circular-shift permutation test
    - q-value matrix via Benjamini–Hochberg FDR
    - significance mask

IMPORTANT (CONDACT exchange-time direction):
    τx(X|Y) = time until X transitions AFTER Y transitions, weighted by Y waiting times
    CA[X][Y] = -ln( τx(X|Y) / τp(X) )

TRANSITION FILTER:
    Only residues with at least 11 transitions are used.
    If tX or tY has < 11 transitions, that CA[X][Y] is NOT tested (p=NaN).

Author: Augustine C. Onyema
"""

import numpy as np
import pickle
import os
import csv


class CONDACT_FDR:
    def __init__(
        self,
        transition_times_file="transition_times.pkl",
        CA_matrix_file="Conditional_Activity_matrix.npy",
        dihedral_csv="Dihedral_Angle.csv",
        total_time=None,             # recommended: pass tau directly
        n_perm=200,
        q=0.05,
        CA_cutoff=0.0,
        random_seed=0,
        upper_triangle_only=True,
        min_transitions=11,          # <-- enforce >=11 transitions
    ):
        self.transition_times_file = transition_times_file
        self.CA_matrix_file = CA_matrix_file
        self.dihedral_csv = dihedral_csv
        self.total_time = total_time
        self.n_perm = int(n_perm)
        self.q = float(q)
        self.CA_cutoff = float(CA_cutoff)
        self.upper_triangle_only = bool(upper_triangle_only)
        self.min_transitions = int(min_transitions)
        self.rng = np.random.default_rng(random_seed)

        # ---- LOAD TRANSITION TIMES ----
        with open(self.transition_times_file, "rb") as f:
            self.transition_times = pickle.load(f)

        # ---- LOAD CA MATRIX ----
        self.CA = np.load(self.CA_matrix_file)
        self.labels = list(self.transition_times.keys())
        self.n = len(self.labels)

        if self.CA.shape != (self.n, self.n):
            raise ValueError(
                f"CA matrix shape {self.CA.shape} does not match number of labels {self.n}."
            )

        # ---- SET τ (OBSERVATION TIME) ----
        if self.total_time is None:
            self.total_time = self._read_total_time_from_dihedral_csv(self.dihedral_csv)

        if not np.isfinite(self.total_time) or self.total_time <= 0:
            raise ValueError(f"Invalid total_time (tau): {self.total_time}")

        print(f"Observation_time (tau) = {self.total_time}")
        print(f"Min transitions required per residue = {self.min_transitions}")

        # Precompute which residues pass the transition-count filter
        self.good_residue = np.zeros(self.n, dtype=bool)
        for k, lab in enumerate(self.labels):
            t = self.transition_times.get(lab, [])
            self.good_residue[k] = (t is not None) and (len(t) >= self.min_transitions)

    # ------------------------------------------------------------------
    # Read tau from last time(ps) in Dihedral_Angle.csv
    # ------------------------------------------------------------------
    def _read_total_time_from_dihedral_csv(self, csv_path):
        if not os.path.exists(csv_path):
            raise FileNotFoundError(
                f"Could not find {csv_path}. "
                f"Provide total_time explicitly or ensure Dihedral_Angle.csv exists."
            )

        chunk_size = 92
        with open(csv_path, "rb") as f:
            f.seek(0, os.SEEK_END)
            end = f.tell()
            pos = max(0, end - chunk_size)
            f.seek(pos)
            data = f.read(end - pos)

            tries = 0
            while data.count(b"\n") < 2 and pos > 0 and tries < 10:
                tries += 1
                new_pos = max(0, pos - chunk_size)
                f.seek(new_pos)
                more = f.read(pos - new_pos)
                data = more + data
                pos = new_pos

        text = data.decode("utf-8", errors="ignore")
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]

        for ln in reversed(lines):
            if ln.lower().startswith("frame"):
                continue
            parts = [p.strip() for p in ln.split(",")]
            if len(parts) < 2:
                continue
            try:
                return float(parts[1])
            except ValueError:
                continue

        raise ValueError(f"Could not parse total_time from the end of {csv_path}.")

    # ------------------------------------------------------------------
    # Timing utilities (NO wrap-around; transitions only)
    # ------------------------------------------------------------------

    def _tau_persistence(self, tX):
        """
        τp(X) = (1 / (2τ)) * Σ (Δt_X)^2
        Uses observed consecutive transitions only.
        """
        tX = np.asarray(tX, float)
        if tX.size < self.min_transitions:
            return np.nan
        tX = np.sort(tX)
        wX = np.diff(tX)
        if wX.size == 0:
            return np.nan
        return np.sum(wX * wX) / (2.0 * self.total_time)

    def _tau_XY(self, tX, tY):
        """
        τx(X|Y) (CONDACT-direction), requires both tX and tY >= min_transitions.

        τx(X|Y) = (1/τ) * Σ_i [ (tX_next - tY[i]) * (tY[i] - tY[i-1]) ], i>=1
        Stops when a Y transition has no later X transition (CONDACT-style break).
        """
        tX = np.asarray(tX, float)
        tY = np.asarray(tY, float)

        if tX.size < self.min_transitions or tY.size < self.min_transitions:
            return np.nan

        tX = np.sort(tX)
        tY = np.sort(tY)

        qY = tY[1:]          # Y transition times excluding first
        wY = np.diff(tY)     # weights: Y waiting times

        idx = np.searchsorted(tX, qY, side="right")
        valid = idx < tX.size
        if not np.any(valid):
            return np.nan

        # CONDACT-style break at first missing next-X
        first_invalid = np.where(~valid)[0]
        if first_invalid.size > 0:
            cut = first_invalid[0]
            qY = qY[:cut]
            wY = wY[:cut]
            idx = idx[:cut]

        if qY.size == 0:
            return np.nan

        WX_after_Y = tX[idx] - qY
        tx = np.sum(WX_after_Y * wY) / self.total_time
        return tx if (np.isfinite(tx) and tx > 0) else np.nan

    def _CA_from_transitions(self, tX, tY):
        """
        CA[X][Y] = -ln( τx(X|Y) / τp(X) )
        """
        tp = self._tau_persistence(tX)
        tx = self._tau_XY(tX, tY)
        if (not np.isfinite(tp)) or (not np.isfinite(tx)) or tp <= 0 or tx <= 0:
            return np.nan
        val = -np.log(tx / tp)
        return val if np.isfinite(val) else np.nan

    # ------------------------------------------------------------------
    # P-values via circular shift null (shuffle Y relative to X)
    # ------------------------------------------------------------------

    def _p_value(self, i, j):
        """
        Only tests if BOTH residues meet the >=11 transition criterion.
        """
        A_obs = self.CA[i, j]
        if (not np.isfinite(A_obs)) or (A_obs < self.CA_cutoff):
            return np.nan

        if (not self.good_residue[i]) or (not self.good_residue[j]):
            return np.nan

        tX = self.transition_times[self.labels[i]]
        tY = self.transition_times[self.labels[j]]

        # double safety
        if len(tX) < self.min_transitions or len(tY) < self.min_transitions:
            return np.nan

        A_null = []
        for _ in range(self.n_perm):
            shift = self.rng.uniform(0.0, self.total_time)
            tY_shift = (np.asarray(tY, float) + shift) % self.total_time
            tY_shift.sort()
            A_null.append(self._CA_from_transitions(tX, tY_shift))

        A_null = np.asarray(A_null, float)
        A_null = A_null[np.isfinite(A_null)]
        if A_null.size == 0:
            return np.nan

        return (1.0 + np.sum(A_null >= A_obs)) / (1.0 + A_null.size)

    # ------------------------------------------------------------------
    # BH-FDR (vector)
    # ------------------------------------------------------------------

    def _bh_fdr_vector(self, p):
        p = np.asarray(p, float)
        m = p.size
        if m == 0:
            return np.zeros(0, dtype=bool), np.zeros(0, dtype=float)

        order = np.argsort(p)
        p_sorted = p[order]

        thresh = self.q * (np.arange(1, m + 1) / m)
        passed = p_sorted <= thresh

        reject_sorted = np.zeros(m, dtype=bool)
        if np.any(passed):
            kmax = np.max(np.where(passed)[0])
            reject_sorted[:kmax + 1] = True

        q_sorted = (m / np.arange(1, m + 1)) * p_sorted
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]

        reject = np.zeros(m, dtype=bool)
        qvals = np.zeros(m, dtype=float)
        reject[order] = reject_sorted
        qvals[order] = q_sorted
        return reject, qvals

    # ------------------------------------------------------------------
    # Public run
    # ------------------------------------------------------------------

    def run(self):
        """
        Directed FDR: treat i->j and j->i as separate hypotheses.

        Compute p-values ONLY for entries where:
            - CA is finite
            - off-diagonal
            - CA >= CA_cutoff
            - BOTH residues have >= 11 transitions
        Then BH-FDR over ALL directed tested entries.

        Output:
            - CA_all_tested_edges.csv: ALL tested directed edges (significant or not)
            - CA_significant_edges.csv: ONLY significant directed edges
        """
        n = self.n
        CA = self.CA

        # Build test mask: CA finite & cutoff & off-diagonal
        mask_test = np.isfinite(CA) & (~np.eye(n, dtype=bool)) & (CA >= self.CA_cutoff)

        # Enforce transition threshold on both i and j
        good = self.good_residue
        mask_test &= (good[:, None] & good[None, :])

        # IMPORTANT: directed testing => do NOT restrict to upper triangle
        # self.upper_triangle_only is ignored by design here.

        ii, jj = np.where(mask_test)
        m_test = len(ii)

        pval_matrix = np.full((n, n), np.nan, dtype=float)
        qval_matrix = np.full((n, n), np.nan, dtype=float)
        sig_mask = np.zeros((n, n), dtype=bool)

        all_path = "CA_all_tested_edges.csv"
        sig_path = "CA_significant_edges.csv"

        # Always write headers (even if empty)
        if m_test == 0:
            with open(all_path, "w", newline="") as f:
                csv.writer(f).writerow(
                    ["i", "j", "Residue_i", "Residue_j", "Conditional_Activity", "p_value", "q_value", "Significant_FDR"]
                )
            with open(sig_path, "w", newline="") as f:
                csv.writer(f).writerow(
                    ["i", "j", "Residue_i", "Residue_j", "Conditional_Activity", "p_value", "q_value"]
                )
            print("No CA entries met the cutoff + >=11 transitions criterion for testing.")
            return pval_matrix, qval_matrix, sig_mask

        # ---- Compute p-values for ALL directed tested pairs ----
        p_test = np.full(m_test, np.nan, dtype=float)
        for k, (i, j) in enumerate(zip(ii, jj)):
            p = self._p_value(i, j)
            p_test[k] = p
            pval_matrix[i, j] = p

        # ---- BH-FDR over finite tested p-values (DIRECTED) ----
        finite = np.isfinite(p_test)
        reject_test = np.zeros(m_test, dtype=bool)
        q_test = np.full(m_test, np.nan, dtype=float)

        if np.any(finite):
            reject_f, q_f = self._bh_fdr_vector(p_test[finite])
            reject_test[finite] = reject_f
            q_test[finite] = q_f

        # Put q-values and significance back into matrices
        for k, (i, j) in enumerate(zip(ii, jj)):
            qval_matrix[i, j] = q_test[k]
            sig_mask[i, j] = reject_test[k]

        # ---- Write ALL tested directed edges ----
        with open(all_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["i", "j", "Residue_i", "Residue_j", "Conditional_Activity", "p_value", "q_value", "Significant_FDR"])
            for k, (i, j) in enumerate(zip(ii, jj)):
                w.writerow([
                    i, j,
                    self.labels[i], self.labels[j],
                    CA[i, j],
                    p_test[k],
                    q_test[k],
                    bool(reject_test[k])
                ])

        # ---- Write ONLY significant directed edges ----
        with open(sig_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["i", "j", "Residue_i", "Residue_j", "Conditional_Activity", "p_value", "q_value"])
            for k, (i, j) in enumerate(zip(ii, jj)):
                if reject_test[k]:
                    w.writerow([
                        i, j,
                        self.labels[i], self.labels[j],
                        CA[i, j],
                        p_test[k],
                        q_test[k]
                    ])

        print(f"Tested directed edges (CA >= {self.CA_cutoff} and >=11 transitions): {m_test}")
        print(f"Finite p-values: {int(np.sum(finite))}")
        print(f"Significant at FDR q <= {self.q}: {int(np.sum(reject_test))}")

        return pval_matrix, qval_matrix, sig_mask

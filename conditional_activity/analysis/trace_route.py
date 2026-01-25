import numpy as np

class CONDACT_CommPath:
    """
    Trace directed communication paths from a CONDACT Conditional Activity (CA) matrix.

    Direction convention (matches your CONDACT):
        CA[X, Y] represents directed influence/edge:  Y  ->  X

    Inputs
    ------
    CA   : (N,N) array
    dist : (N,N) array (Å)
    labels : optional list[str] of length N
        e.g. ["GLU 732", "ARG 42", ...]

    Features
    --------
    - Build directed adjacency under distance cutoff and CA cutoff.
    - Start/target can be:
        * index (int)
        * residue label (str) like "GLU 732" (requires labels)
        * start=None with start_mode="global_max" (auto-start from strongest CA edge)
    - Greedy / beam tracing: walk Y->X, then X becomes new Y (optionally multiple paths).
    - DFS enumeration of many simple paths (no cycles) with scoring.
    """

    def __init__(self, CA, dist, labels=None):
        self.CA = np.asarray(CA, float)
        self.dist = np.asarray(dist, float)

        if self.CA.ndim != 2 or self.dist.ndim != 2:
            raise ValueError("CA and dist must be 2D arrays.")
        if self.CA.shape != self.dist.shape:
            raise ValueError("CA and dist must have the same shape.")
        if self.CA.shape[0] != self.CA.shape[1]:
            raise ValueError("CA and dist must be square NxN matrices.")

        self.N = self.CA.shape[0]
        self.labels = list(labels) if labels is not None else None
        self.label_to_index = self._make_label_to_index(self.labels)

    # ------------------------------------------------------------
    # Label/index helpers
    # ------------------------------------------------------------

    @staticmethod
    def _normalize_label(s: str) -> str:
        return " ".join(str(s).strip().split())

    @classmethod
    def _make_label_to_index(cls, labels):
        if labels is None:
            return {}
        return {cls._normalize_label(lab): i for i, lab in enumerate(labels)}

    def _as_index(self, node, name="node"):
        """
        node can be:
          - int
          - str label (requires labels)
          - None
        """
        if node is None:
            return None
        if isinstance(node, (int, np.integer)):
            idx = int(node)
        elif isinstance(node, str):
            if self.labels is None:
                raise ValueError(f"{name} given as label string but self.labels is None.")
            key = self._normalize_label(node)
            if key not in self.label_to_index:
                raise ValueError(f"{name} label '{node}' not found in labels.")
            idx = int(self.label_to_index[key])
        else:
            raise TypeError(f"{name} must be int, str, or None. Got {type(node)}")

        if not (0 <= idx < self.N):
            raise ValueError(f"{name} index {idx} out of range 0..{self.N-1}.")
        return idx

    def _fmt_node(self, idx, show_index=False):
        if self.labels is None:
            return str(idx)
        return f"{self.labels[idx]}({idx})" if show_index else self.labels[idx]

    # ------------------------------------------------------------
    # Adjacency
    # ------------------------------------------------------------

    def build_adjacency(self, cutoff=5.0, min_ca=0.0, top_k=None, finite_dist_only=True):
        """
        Build directed adjacency lists:
            neighbors[y] = [(x, CA[x,y]), ...] sorted by CA desc
        """
        CA = self.CA
        dist = self.dist
        N = self.N

        M = CA.copy()
        M[~np.isfinite(M)] = -np.inf
        np.fill_diagonal(M, -np.inf)

        dist_ok = np.isfinite(dist) if finite_dist_only else np.ones_like(dist, dtype=bool)
        allowed = dist_ok & (dist <= cutoff) & np.isfinite(M) & (M > min_ca)
        np.fill_diagonal(allowed, False)

        neighbors = [[] for _ in range(N)]
        for y in range(N):
            xs = np.where(allowed[:, y])[0]  # X such that Y->X allowed
            if xs.size == 0:
                continue
            vals = M[xs, y]
            order = np.argsort(vals)[::-1]
            xs = xs[order]
            vals = vals[order]
            if top_k is not None:
                xs = xs[:top_k]
                vals = vals[:top_k]
            neighbors[y] = [(int(x), float(v)) for x, v in zip(xs, vals)]
        return neighbors

    # ------------------------------------------------------------
    # Auto-start
    # ------------------------------------------------------------

    def choose_start_from_global_max(self, cutoff=None, min_ca=0.0):
        """
        Pick start node as the source (Y) of the globally strongest edge (Y->X).
        Returns: (start_y, (y,x,val))
        """
        CA = self.CA
        M = CA.copy()
        M[~np.isfinite(M)] = -np.inf
        np.fill_diagonal(M, -np.inf)

        mask = np.isfinite(M) & (M > min_ca)
        if cutoff is not None:
            d = self.dist
            mask &= np.isfinite(d) & (d <= cutoff)

        if not np.any(mask):
            raise ValueError("No edges available under the provided filters.")

        # CA[X,Y] is edge Y->X; argmax over (X,Y)
        x, y = np.unravel_index(np.nanargmax(np.where(mask, M, -np.inf)), M.shape)
        return int(y), (int(y), int(x), float(M[x, y]))

    # ------------------------------------------------------------
    # Greedy / Beam tracing (walk: Y -> X, then Y := X)
    # ------------------------------------------------------------

    def trace_greedy(
        self,
        start=None,
        target=None,
        cutoff=5.0,
        min_ca=0.0,
        top_k=None,
        max_steps=None,
        branch_k=1,
        stop_if_no_out=True,
        score_mode="sum",
    ):
        """
        Greedy (branch_k=1) or beam search (branch_k>1) path tracing.

        Returns list of path dicts sorted by score desc:
            {"nodes":[...], "edges":[(y,x,val),...], "score":float, "reached_target":bool}
        """
        if max_steps is None:
            max_steps = self.N - 1

        if start is None:
            start, _ = self.choose_start_from_global_max(cutoff=cutoff, min_ca=min_ca)
        start = self._as_index(start, "start")
        target = self._as_index(target, "target")

        neighbors = self.build_adjacency(cutoff=cutoff, min_ca=min_ca, top_k=top_k)

        def score_update(curr_score, edge_val):
            if score_mode == "sum":
                return curr_score + edge_val
            if score_mode == "product":
                return curr_score * edge_val
            if score_mode == "min":
                return min(curr_score, edge_val)
            raise ValueError("score_mode must be 'sum', 'product', or 'min'.")

        # beam: (score, nodes, edges, visited_set)
        if score_mode == "product":
            beam = [(1.0, [start], [], {start})]
        elif score_mode == "min":
            beam = [(float("inf"), [start], [], {start})]
        else:
            beam = [(0.0, [start], [], {start})]

        for _ in range(int(max_steps)):
            new_beam = []

            for sc, nodes, edges, visited in beam:
                cur = nodes[-1]

                if target is not None and cur == target:
                    new_beam.append((sc, nodes, edges, visited))
                    continue

                outs = neighbors[cur]
                if not outs:
                    if not stop_if_no_out:
                        new_beam.append((sc, nodes, edges, visited))
                    continue

                for nxt, v in outs:
                    if nxt in visited:
                        continue
                    new_sc = score_update(sc, v)
                    new_beam.append((new_sc, nodes + [nxt], edges + [(cur, nxt, v)], visited | {nxt}))
                    if branch_k == 1:
                        # neighbors are sorted high->low; pick the best and stop
                        break

            if not new_beam:
                break

            new_beam.sort(key=lambda t: t[0], reverse=True)
            beam = new_beam[:max(1, int(branch_k))]

            if target is not None and all(p[1][-1] == target for p in beam):
                break

        out = []
        for sc, nodes, edges, _ in beam:
            out.append({
                "nodes": nodes,
                "edges": edges,
                "score": float(sc),
                "reached_target": (target is None) or (nodes[-1] == target)
            })
        out.sort(key=lambda d: d["score"], reverse=True)
        return out

    # ------------------------------------------------------------
    # Enumerate many paths (DFS)
    # ------------------------------------------------------------

    def enumerate_paths(
        self,
        start=None,
        target=None,
        cutoff=5.0,
        min_ca=0.0,
        max_depth=6,
        max_paths=50000,
        top_k=None,
        score="sum",
        start_mode="given",   # "given" | "global_max"
    ):
        """
        Enumerate many directed simple paths (no cycles), using DFS.

        Returns list of dicts: {"nodes":[...], "edges":[(y,x,val),...], "score":float}
        """
        if start is None:
            if start_mode == "global_max":
                start, _ = self.choose_start_from_global_max(cutoff=cutoff, min_ca=min_ca)
            else:
                raise ValueError("start is None. Use start_mode='global_max' or provide start explicitly.")
        start = self._as_index(start, "start")
        target = self._as_index(target, "target")

        neighbors = self.build_adjacency(cutoff=cutoff, min_ca=min_ca, top_k=top_k)

        def path_score(vals):
            if not vals:
                return 0.0
            if score == "sum":
                return float(np.sum(vals))
            if score == "product":
                p = 1.0
                for v in vals:
                    p *= v
                return float(p)
            if score == "min":
                return float(np.min(vals))
            raise ValueError("score must be 'sum', 'product', or 'min'.")

        results = []
        stack = [(start, [start], [], [], {start})]  # (current, nodes, edges, vals, visited)

        while stack and len(results) < int(max_paths):
            cur, nodes, edges, vals, visited = stack.pop()

            if target is not None:
                if cur == target:
                    results.append({"nodes": nodes, "edges": edges, "score": path_score(vals)})
                    continue
            else:
                if edges:
                    results.append({"nodes": nodes, "edges": edges, "score": path_score(vals)})

            if len(edges) >= int(max_depth):
                continue

            for nxt, v in neighbors[cur]:
                if nxt in visited:
                    continue
                stack.append((
                    nxt,
                    nodes + [nxt],
                    edges + [(cur, nxt, v)],
                    vals + [v],
                    visited | {nxt}
                ))

        results.sort(key=lambda d: d["score"], reverse=True)
        return results

    # ------------------------------------------------------------
    # Printing
    # ------------------------------------------------------------

    def print_paths(self, paths, top=20, ndigits=3, show_indices=False):
        for i, p in enumerate(paths[:top], 1):
            node_str = " → ".join(self._fmt_node(j, show_index=show_indices) for j in p["nodes"])
            edge_str = " → ".join(f"{v:.{ndigits}f}" for (_, _, v) in p["edges"])
            print(f"{i:>3}. score={p['score']:.{ndigits}f} | {node_str}")
            if p["edges"]:
                print(f"     CA: {edge_str}")
            if "reached_target" in p:
                print(f"     reached_target={p['reached_target']}")


# ===========================
# Example usage
# ===========================
"""
comm = CONDACT_CommPath(CA, dist, labels=labels)

# 1) enumerate paths between residue labels
paths = comm.enumerate_paths(start="GLU 732", target="ASP 3",
                            cutoff=8.0, min_ca=0.0,
                            max_depth=6, top_k=10, score="sum")
comm.print_paths(paths, top=10)

# 2) greedy single path
g = comm.trace_greedy(start="GLU 732", target="ASP 3",
                      cutoff=8.0, min_ca=0.0,
                      top_k=20, branch_k=1, max_steps=10)
comm.print_paths(g)

# 3) beam search (multiple plausible paths)
g = comm.trace_greedy(start="GLU 732", target="ASP 3",
                      cutoff=8.0, min_ca=0.0,
                      top_k=20, branch_k=10, max_steps=10)
comm.print_paths(g, top=10)

# 4) auto-start from global strongest edge
paths = comm.enumerate_paths(start=None, target=None,
                            cutoff=8.0, min_ca=0.0,
                            max_depth=5, top_k=10,
                            start_mode="global_max")
comm.print_paths(paths, top=10)
"""

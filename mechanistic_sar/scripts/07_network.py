#!/usr/bin/env python3
"""
Step 7: Network Analysis - Build residue contact networks from PDB,
compute shortest paths, betweenness centrality, community structure
(simple label propagation), pocket-to-gate communication scores.
Save to features/network/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
NETWORK = FEATURES / "network"
NETWORK.mkdir(parents=True, exist_ok=True)

CACHE = RESULTS / ".cache"
CACHE.mkdir(parents=True, exist_ok=True)


def _load_cache(name):
    p = CACHE / f"{name}.json"
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return None


def _save_cache(name, data):
    with open(CACHE / f"{name}.json", "w") as f:
        json.dump(data, f)


def parse_ca_atoms(pdb_path):
    """Extract C-alpha atoms from PDB."""
    ca_atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and "CA" in line[12:16]:
                try:
                    chain = line[21].strip()
                    res_name = line[17:20].strip()
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    ca_atoms.append({
                        "chain": chain,
                        "res_name": res_name,
                        "res_seq": res_seq,
                        "x": x, "y": y, "z": z,
                        "key": f"{chain}_{res_seq}",
                    })
                except (ValueError, IndexError):
                    continue
    return ca_atoms


def build_adjacency_matrix(ca_atoms, cutoff=8.0):
    """Build adjacency matrix from C-alpha distances."""
    n = len(ca_atoms)
    adj = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dx = ca_atoms[i]["x"] - ca_atoms[j]["x"]
            dy = ca_atoms[i]["y"] - ca_atoms[j]["y"]
            dz = ca_atoms[i]["z"] - ca_atoms[j]["z"]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if dist <= cutoff:
                adj[i, j] = 1.0
                adj[j, i] = 1.0
    return adj


def build_weighted_adjacency(ca_atoms, cutoff=12.0):
    """Build weighted adjacency matrix (inverse distance)."""
    n = len(ca_atoms)
    adj = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dx = ca_atoms[i]["x"] - ca_atoms[j]["x"]
            dy = ca_atoms[i]["y"] - ca_atoms[j]["y"]
            dz = ca_atoms[i]["z"] - ca_atoms[j]["z"]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if dist <= cutoff and dist > 0.1:
                adj[i, j] = 1.0 / dist
                adj[j, i] = 1.0 / dist
    return adj


def floyd_warshall(adj):
    """Compute all-pairs shortest paths using Floyd-Warshall."""
    n = adj.shape[0]
    # Initialize distance matrix
    dist = np.full((n, n), float("inf"))
    np.fill_diagonal(dist, 0)

    for i in range(n):
        for j in range(n):
            if adj[i, j] > 0:
                dist[i, j] = 1.0 / adj[i, j] if adj[i, j] > 0 else float("inf")
            elif i == j:
                dist[i, j] = 0.0

    # Floyd-Warshall
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i, k] + dist[k, j] < dist[i, j]:
                    dist[i, j] = dist[i, k] + dist[k, j]

    return dist


def compute_betweenness_centrality(adj):
    """
    Compute betweenness centrality using Brandes' algorithm (simplified).
    """
    n = adj.shape[0]
    centrality = np.zeros(n)

    for s in range(n):
        # BFS from s
        stack = []
        predecessors = [[] for _ in range(n)]
        sigma = np.zeros(n)
        sigma[s] = 1.0
        dist = np.full(n, -1.0)
        dist[s] = 0.0
        queue = [s]

        while queue:
            v = queue.pop(0)
            stack.append(v)
            for w in range(n):
                if adj[v, w] == 0:
                    continue
                if dist[w] < 0:
                    dist[w] = dist[v] + 1
                    queue.append(w)
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    predecessors[w].append(v)

        # Back-propagation
        delta = np.zeros(n)
        while stack:
            w = stack.pop()
            for v in predecessors[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
            if w != s:
                centrality[w] += delta[w]

    # Normalize
    if n > 2:
        centrality /= (n - 1) * (n - 2)

    return centrality


def label_propagation_community(adj, max_iter=100, seed=42):
    """
    Simple label propagation for community detection.
    """
    n = adj.shape[0]
    rng = np.random.RandomState(seed)
    labels = np.arange(n)  # Each node starts with unique label

    for _ in range(max_iter):
        order = rng.permutation(n)
        changed = False
        for i in order:
            # Get neighbors
            neighbors = np.where(adj[i] > 0)[0]
            if len(neighbors) == 0:
                continue
            # Count labels among neighbors
            neighbor_labels = labels[neighbors]
            unique_labels, counts = np.unique(neighbor_labels, return_counts=True)
            # Pick most common (random tiebreak)
            max_count = counts.max()
            candidates = unique_labels[counts == max_count]
            new_label = rng.choice(candidates)
            if new_label != labels[i]:
                labels[i] = new_label
                changed = True
        if not changed:
            break

    # Relabel to 0, 1, 2, ...
    unique = np.unique(labels)
    mapping = {old: new for new, old in enumerate(unique)}
    labels = np.array([mapping[l] for l in labels])

    return labels


def compute_communication_score(dist_matrix, node_i, node_j):
    """
    Compute communication score between two nodes:
    inverse of shortest path length.
    """
    d = dist_matrix[node_i, node_j]
    if d == float("inf") or d == 0:
        return 0.0
    return 1.0 / d


def find_node_index(ca_atoms, chain, res_seq):
    """Find index of a node by chain and residue number."""
    for i, a in enumerate(ca_atoms):
        if a["chain"] == chain and a["res_seq"] == res_seq:
            return i
    return None


def main():
    print("=" * 70)
    print("NETWORK ANALYSIS - Mechanistic SAR Analysis")
    print("=" * 70)

    # Find PDB files
    pdb_files = []
    for root_dir in [
        PROJECT / "phase08_boltz_affinity" / "boltz2_sar_outputs",
        PROJECT / "phase08_boltz_affinity" / "boltz2_inactive_outputs",
    ]:
        if root_dir.exists():
            for p in root_dir.glob("*_model_*.pdb"):
                pdb_files.append(p)

    ref = PROJECT / "09_docking" / "receptor_2to3.pdb"
    if ref.exists():
        pdb_files.insert(0, ref)

    print(f"\n[1/4] Found {len(pdb_files)} PDB files")

    print("\n[2/4] Computing network features...")
    network_features = {}
    processed = 0

    # Reference pocket and gate residues
    pocket_residues_alpha9 = [149, 188, 151, 150, 152, 93]
    gate_residues_alpha9 = [251, 255]

    for pdb_path in sorted(pdb_files):
        try:
            ca_atoms = parse_ca_atoms(str(pdb_path))
            if len(ca_atoms) < 10:
                continue

            n = len(ca_atoms)
            adj = build_adjacency_matrix(ca_atoms, cutoff=8.0)
            weighted_adj = build_weighted_adjacency(ca_atoms, cutoff=12.0)

            # Number of edges
            n_edges = int(np.sum(adj) / 2)

            # Average degree
            degrees = np.sum(adj, axis=1)
            mean_degree = float(np.mean(degrees)) if n > 0 else 0.0

            # Betweenness centrality
            centrality = compute_betweenness_centrality(adj)
            mean_centrality = float(np.mean(centrality))
            max_centrality = float(np.max(centrality))

            # Community detection
            communities = label_propagation_community(adj, max_iter=50)
            n_communities = len(np.unique(communities))

            # Community sizes
            comm_sizes = {}
            for c in communities:
                comm_sizes[int(c)] = comm_sizes.get(int(c), 0) + 1

            # Network efficiency
            dist_matrix = floyd_warshall(weighted_adj)
            finite_dists = dist_matrix[dist_matrix < float("inf")]
            efficiency = float(np.mean(1.0 / finite_dists)) if len(finite_dists) > 0 else 0.0

            # Pocket-to-gate communication score
            pocket_nodes = []
            gate_nodes = []
            for i, a in enumerate(ca_atoms):
                if a["chain"] == "A":
                    if a["res_seq"] in pocket_residues_alpha9:
                        pocket_nodes.append(i)
                    if a["res_seq"] in gate_residues_alpha9:
                        gate_nodes.append(i)

            # Also try chain B
            for i, a in enumerate(ca_atoms):
                if a["chain"] == "B":
                    if a["res_seq"] in pocket_residues_alpha9:
                        pocket_nodes.append(i)
                    if a["res_seq"] in gate_residues_alpha9:
                        gate_nodes.append(i)

            pocket_gate_score = 0.0
            if pocket_nodes and gate_nodes:
                scores = []
                for pi in pocket_nodes:
                    for gi in gate_nodes:
                        s = compute_communication_score(dist_matrix, pi, gi)
                        scores.append(s)
                pocket_gate_score = float(np.mean(scores)) if scores else 0.0

            # Path changes (compare to shortest paths)
            # Mean shortest path
            finite_mask = dist_matrix < float("inf")
            if np.any(finite_mask):
                mean_path = float(np.mean(dist_matrix[finite_mask]))
            else:
                mean_path = float("inf")

            name = pdb_path.stem
            network_features[name] = {
                "pdb_file": str(pdb_path.relative_to(PROJECT)),
                "n_nodes": n,
                "n_edges": n_edges,
                "mean_degree": mean_degree,
                "max_degree": float(np.max(degrees)) if n > 0 else 0.0,
                "mean_betweenness": mean_centrality,
                "max_betweenness": max_centrality,
                "n_communities": n_communities,
                "largest_community": max(comm_sizes.values()) if comm_sizes else 0,
                "network_efficiency": efficiency,
                "pocket_gate_communication": pocket_gate_score,
                "mean_shortest_path": mean_path,
                "community_sizes": comm_sizes,
            }

            processed += 1
            print(f"  Processed {pdb_path.name}: {n} nodes, {n_edges} edges, "
                  f"{n_communities} communities")

        except Exception as e:
            print(f"  WARNING: Error processing {pdb_path.name}: {e}")

    print(f"\n  Processed: {processed} structures")

    print("\n[3/4] Saving network features...")
    # Save JSON
    with open(NETWORK / "network_features.json", "w") as f:
        json.dump(network_features, f, indent=2)

    # Save CSV
    if network_features:
        names = sorted(network_features.keys())
        csv_keys = [k for k in network_features[names[0]].keys()
                    if k != "community_sizes"]
        with open(NETWORK / "network_features.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=csv_keys)
            writer.writeheader()
            for name in names:
                writer.writerow({k: network_features[name][k] for k in csv_keys})

    print("\n[4/4] Summary...")
    if network_features:
        effs = [v["network_efficiency"] for v in network_features.values()]
        print(f"  Network efficiency range: [{min(effs):.6f}, {max(effs):.6f}]")
        comms = [v["n_communities"] for v in network_features.values()]
        print(f"  Communities range: [{min(comms)}, {max(comms)}]")

    print(f"\nResults saved to: {NETWORK}")
    print("  network_features.json")
    print("  network_features.csv")


if __name__ == "__main__":
    main()

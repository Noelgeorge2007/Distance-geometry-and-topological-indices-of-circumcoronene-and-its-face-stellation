# Data Dictionary

All CSV files have a header row. `n` is the positive integer size parameter. `family` is `Cor` for circumcoronene H_n or `St` for face-stellated T_n. Graphs are simple, connected, and unweighted.

## indices.csv

One row per graph. `vertices`, `edges`, and `diameter` are the vertex count, edge count, and maximum shortest-path distance.

In the definitions below, d(u,v) is shortest-path distance, deg(u) is vertex degree, and ecc(u) is vertex eccentricity. Pair sums are over distinct unordered vertex pairs; edge sums count each undirected edge once.

| Column | Definition |
| --- | --- |
| `M1` | First Zagreb index: sum over vertices of deg(u)^2 |
| `M2` | Second Zagreb index: sum over edges uv of deg(u) deg(v) |
| `F` | Forgotten index: sum over vertices of deg(u)^3 |
| `R` | Randic index: sum over edges of 1/sqrt(deg(u) deg(v)) |
| `SC` | Sum-connectivity index: sum over edges of 1/sqrt(deg(u)+deg(v)) |
| `GA` | Geometric-arithmetic index: sum over edges of 2 sqrt(deg(u) deg(v))/(deg(u)+deg(v)) |
| `ABC` | Atom-bond connectivity index: sum over edges of sqrt((deg(u)+deg(v)-2)/(deg(u) deg(v))) |
| `W` | Wiener index: sum over pairs of d(u,v) |
| `WW` | Hyper-Wiener index: one half of the pair sum of d(u,v)+d(u,v)^2 |
| `S` | Schultz index: sum over pairs of (deg(u)+deg(v)) d(u,v) |
| `Sstar` | Modified Schultz index: sum over pairs of deg(u) deg(v) d(u,v) |
| `xi` | Eccentric connectivity index: sum over vertices of deg(u) ecc(u) |
| `zeta` | Total eccentricity: sum over vertices of ecc(u) |
| `H` | Harary index: sum over pairs of 1/d(u,v) |
| `AD` | Average distance: 2W/(vertices (vertices-1)) |
| `log10_pi` | Base-10 logarithm of the product of distances over all unordered vertex pairs |
| `Sz` | Szeged index: sum over edges uv of n_u n_v |
| `PI` | Vertex PI index: sum over edges uv of n_u+n_v |
| `Mo` | Mostar index: sum over edges uv of abs(n_u-n_v) |

For an edge uv, n_u counts vertices strictly closer to u than to v, and n_v counts those strictly closer to v. Equidistant vertices contribute to neither count. `PI` is the vertex version, not edge PI.

`H` and `AD` are exact rational strings (an integer or numerator/denominator). `R`, `SC`, `GA`, `ABC`, and `log10_pi` are floating-point values. Other descriptor columns are integers. The distance-product descriptor is stored logarithmically to avoid a very large integer; its exact value is recoverable as the product of distance^pairs from the corresponding distance-partition rows.

## distance_partitions.csv

Columns: `n`, `family`, `distance`, `pairs`. `pairs` is the number of distinct unordered vertex pairs at the given positive distance, including zero counts if present in the histogram. For each graph, the sum of `pairs` is vertices*(vertices-1)/2; the sum of `distance*pairs` is W.

## verification.csv

Columns: `n`, `family`, `index`, `formula`, `bfs`, `absolute_error`. `formula` is the evaluated expression or grouped degree calculation. The historical column name `bfs` denotes the direct graph calculation: BFS distances for distance-based descriptors, or direct vertex/edge sums for degree-based descriptors. It does not mean every entry itself uses BFS. The table contains 405 selected formula comparisons, not 19 polynomial checks for every graph.

## old_vertex_comparison.csv

One row per n. `old_pairs` counts unordered pairs of original H_n vertices. `shortened_pairs` counts those whose distance decreases in T_n. `maximum_shortening` is the largest decrease. `W_old_in_St` is the sum of T_n distances restricted to original vertices, not the Wiener index of the full T_n graph.

## verification_summary.json

Records the largest n, graph and descriptor counts, checked coordinate pairs, horizontal-edge partition checks, formula comparisons, symbolic certificates, maximum degree-index error, runtime versions, and completion status. `PASS` applies to the included checks and tested range; it is not a claim of experimental validation or proof of every expression for all n.

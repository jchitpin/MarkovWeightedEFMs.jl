# ACHMC (for AEFMs; under the hood)

This section explains the pre-processing steps under the hood of
`preprocess_all_for_atomic_chmc()` using the same multispecies reaction
network from the previous tutorial.

![Toy multispecies network](../assets/toy-network-1-achmc.png)

```@setup required
# Bodge until I figure out how to install Python and RXMapper on Git workflows
using MarkovWeightedEFMs
using SparseArrays

S = [#
  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 # Glc
  0 -1  0  0  0 -1  0  0  0  0  0  0  1  0  0  0 # ATP
  0  1 -1  0 -1  0  0  0  0  0  0  0  0  0  0  0 # G6P
  0  1  0  0  0  1  0  0  0  0  0  0  0 -1  0  0 # ADP
  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0 # 6PG
  0  0  0  0  1 -1  1  0  0  0  0  0  0  0  0  0 # F6P
  0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1  0 # Pi
  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  1 # H2O
  0  0  0  0  0  1 -1 -1  0  0  0  0  0  0  0  0 # FDP
  0  0  0  0  0  0  0  1 -1  1 -1  0  0  0  0  0 # G3P
  0  0  0  0  0  0  0  1  1 -1  0 -1  0  0  0  0 # DHAP
]

v = [10, 10, 3, 3, 7, 8, 1, 7, 1, 1, 7, 7, 18, 18, 1, 1]

mets = [#
  "Glc",
  "ATP",
  "G6P",
  "ADP",
  "6PG",
  "F6P",
  "Pi",
  "H2O",
  "FDP",
  "G3P",
  "DHAP"
]

rxns = [#
  "Source Glc",
  "Hexokinase",
  "G6P dehydrogenase",
  "Sink 6PG",
  "Phosphoglucose isomerase",
  "6-phosphofructo-1-kinase",
  "Fructose 1,6-bisphosphatase",
  "Fructose-bisphosphate aldolase",
  "Triose phosphate isomerase",
  "Triose phosphate isomerase",
  "Sink G3P",
  "Sink DHAP",
  "Source ATP",
  "Sink ADP",
  "Sink Pi",
  "Source H2O"
]

smiles = [#
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)OP(=O)(O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
  "O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O",
  "C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)OP(=O)(O)O",
  "[O-]P(=O)([O-])[O-]",
  "O",
  "C(C1C(C(C(O1)(COP(=O)(O)O)O)O)O)OP(=O)(O)O",
  "C([C@H](C=O)O)OP(=O)(O)O",
  "C(C(=O)COP(=O)(O)O)O"
]

atom = :C # carbon atom type for AEFMs

errors = CHMCAtomicErrorSummary(#
    0.0,     #absolute_flux_error
    Vector{Int64}(), #reactions_duplicated
    Vector{Int64}(), #reactions_with_zero_flux
    Vector{Int64}(), #reactions_with_negative_flux
    Vector{Int64}(), #reactions_with_non_integer_stoichiometries
    Vector{Int64}([1, 13, 16]), #reactions_unimolecular_source_stoichiometry_one
    Vector{Int64}(), #reactions_unimolecular_source_stoichiometry_not_one
    Vector{Int64}(), #reactions_multimolecular_source_stoichiometry_one
    Vector{Int64}(), #reactions_multimolecular_source_stoichiometry_not_one
    Vector{Int64}([4, 11, 12, 14, 15]), #reactions_unimolecular_sink_stoichiometry_one
    Vector{Int64}(), #reactions_unimolecular_sink_stoichiometry_not_one
    Vector{Int64}(), #reactions_multimolecular_sink_stoichiometry_one
    Vector{Int64}(), #reactions_multimolecular_sink_stoichiometry_not_one
    Vector{Int64}(), #reactions_empty
    Vector{Int64}() #unused_metabolites
)

logs = (#
    dropped_row_pseudometabolites = Int64[],
    dropped_cols_pseudometabolites = Int64[],
    dropped_cols_rxnmapper_limit = Int64[]
)
smiles = [
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"
    "O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"
    "O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O"
    "O=C(CO)[C@@H](O)[C@H](O)[C@H](O)COP(=O)(O)O"
    "O=P([O-])([O-])[O-]"
    "O"
    "O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O"
    "O=C[C@H](O)COP(=O)(O)O"
    "O=C(CO)COP(=O)(O)O"
]

rs = [#
    ">>OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O.Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O>>O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O.Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O>>O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O",
    "O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O>>",
    "O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O>>O=C(CO)[C@@H](O)[C@H](O)[C@H](O)COP(=O)(O)O",
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O.O=C(CO)[C@@H](O)[C@H](O)[C@H](O)COP(=O)(O)O>>Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O.O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O",
    "O.O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O>>O=C(CO)[C@@H](O)[C@H](O)[C@H](O)COP(=O)(O)O.O=P([O-])([O-])[O-]",
    "O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O>>O=C[C@H](O)COP(=O)(O)O.O=C(CO)COP(=O)(O)O",
    "O=C[C@H](O)COP(=O)(O)O>>O=C(CO)COP(=O)(O)O",
    "O=C(CO)COP(=O)(O)O>>O=C[C@H](O)COP(=O)(O)O",
    "O=C[C@H](O)COP(=O)(O)O>>",
    "O=C(CO)COP(=O)(O)O>>",
    ">>Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O>>",
    "O=P([O-])([O-])[O-]>>",
    ">>O"
]

ms = [#
    "",
    "[OH:5][CH2:6][C@H:7]1[O:8][CH:9]([OH:10])[C@H:11]([OH:12])[C@@H:13]([OH:14])[C@@H:15]1[OH:16].[NH2:17][c:18]1[n:19][cH:20][n:21][c:22]2[c:23]1[n:24][cH:25][n:26]2[C@@H:27]1[O:28][C@H:29]([CH2:30][O:31][P:32](=[O:33])([OH:34])[O:35][P:36](=[O:37])([OH:38])[O:39][P:2](=[O:1])([OH:3])[OH:4])[C@@H:40]([OH:41])[C@H:42]1[OH:43]>>[O:1]=[P:2]([OH:3])([OH:4])[O:5][CH2:6][C@H:7]1[O:8][CH:9]([OH:10])[C@H:11]([OH:12])[C@@H:13]([OH:14])[C@@H:15]1[OH:16].[NH2:17][c:18]1[n:19][cH:20][n:21][c:22]2[c:23]1[n:24][cH:25][n:26]2[C@@H:27]1[O:28][C@H:29]([CH2:30][O:31][P:32](=[O:33])([OH:34])[O:35][P:36](=[O:37])([OH:38])[OH:39])[C@@H:40]([OH:41])[C@H:42]1[OH:43]",
    "[O:8]=[P:7]([OH:9])([OH:10])[O:6][CH2:5][C@H:4]1[O:3][CH:2]([OH:1])[C@H:15]([OH:16])[C@@H:13]([OH:14])[C@@H:11]1[OH:12]>>[O:1]=[C:2]1[O:3][C@H:4]([CH2:5][O:6][P:7](=[O:8])([OH:9])[OH:10])[C@@H:11]([OH:12])[C@H:13]([OH:14])[C@H:15]1[OH:16]",
    "",
    "[O:14]=[P:13]([OH:15])([OH:16])[O:12][CH2:11][C@H:9]1[O:10][CH:3]([OH:4])[C@H:2]([OH:1])[C@@H:5]([OH:6])[C@@H:7]1[OH:8]>>[O:1]=[C:2]([CH2:3][OH:4])[C@@H:5]([OH:6])[C@H:7]([OH:8])[C@H:9]([OH:10])[CH2:11][O:12][P:13](=[O:14])([OH:15])[OH:16]",
    "[NH2:1][c:2]1[n:3][cH:4][n:5][c:6]2[c:7]1[n:8][cH:9][n:10]2[C@@H:11]1[O:12][C@H:13]([CH2:14][O:15][P:16](=[O:17])([OH:18])[O:19][P:20](=[O:21])([OH:22])[O:23][P:29](=[O:28])([OH:30])[OH:31])[C@@H:24]([OH:25])[C@H:26]1[OH:27].[O:35]=[C:34]([CH2:33][OH:32])[C@@H:46]([OH:47])[C@H:44]([OH:45])[C@H:36]([OH:37])[CH2:38][O:39][P:40](=[O:41])([OH:42])[OH:43]>>[NH2:1][c:2]1[n:3][cH:4][n:5][c:6]2[c:7]1[n:8][cH:9][n:10]2[C@@H:11]1[O:12][C@H:13]([CH2:14][O:15][P:16](=[O:17])([OH:18])[O:19][P:20](=[O:21])([OH:22])[OH:23])[C@@H:24]([OH:25])[C@H:26]1[OH:27].[O:28]=[P:29]([OH:30])([OH:31])[O:32][CH2:33][CH:34]1[O:35][C:36]([OH:37])([CH2:38][O:39][P:40](=[O:41])([OH:42])[OH:43])[CH:44]([OH:45])[CH:46]1[OH:47]",
    "[OH2:20].[O:17]=[P:18]([OH:19])([OH:21])[O:4][CH2:3][CH:2]1[O:1][C:9]([OH:10])([CH2:11][O:12][P:13](=[O:14])([OH:15])[OH:16])[CH:7]([OH:8])[CH:5]1[OH:6]>>[O:1]=[C:2]([CH2:3][OH:4])[C@@H:5]([OH:6])[C@H:7]([OH:8])[C@H:9]([OH:10])[CH2:11][O:12][P:13](=[O:14])([OH:15])[OH:16].[O:17]=[P:18]([O-:19])([O-:20])[O-:21]",
    "[O:18]=[P:17]([OH:19])([OH:20])[O:1][CH2:2][CH:3]1[O:4][C:12]([OH:11])([CH2:13][O:14][P:7](=[O:8])([OH:9])[OH:10])[CH:15]([OH:16])[CH:5]1[OH:6]>>[O:1]=[CH:2][C@H:3]([OH:4])[CH2:5][O:6][P:7](=[O:8])([OH:9])[OH:10].[O:11]=[C:12]([CH2:13][OH:14])[CH2:15][O:16][P:17](=[O:18])([OH:19])[OH:20]",
    "[O:4]=[CH:3][C@H:2]([OH:1])[CH2:5][O:6][P:7](=[O:8])([OH:9])[OH:10]>>[O:1]=[C:2]([CH2:3][OH:4])[CH2:5][O:6][P:7](=[O:8])([OH:9])[OH:10]",
    "[O:1]=[C:2]([CH2:3][OH:4])[CH2:5][O:6][P:7](=[O:8])([OH:9])[OH:10]>>[O:4]=[CH:3][C@H:2]([OH:1])[CH2:5][O:6][P:7](=[O:8])([OH:9])[OH:10]",
    "",
    "",
    "",
    "",
    "",
    ""
]

atom_max = [6, 10, 6, 10, 6, 6, 0, 0, 6, 3, 3]

ks = [#
    (3, 2, 1, 3),
    (2, 2, 1, 6),
    (1, 2, 1, 2),
    (2, 5, 1, 2),
    (3, 3, 1, 5),
    (3, 4, 1, 5),
    (3, 6, 1, 3),
    (2, 6, 1, 6),
    (1, 6, 1, 2),
    (2, 7, 1, 6),
    (11, 3, 1, 10),
    (2, 9, 1, 2),
    (3, 1, 1, 5),
    (3, 2, 1, 5),
    (3, 5, 1, 3),
    (2, 5, 1, 6),
    (9, 3, 1, 7),
    (1, 5, 1, 2),
    (9, 4, 1, 7),
    (11, 1, 1, 10),
    (6, 3, 1, 6),
    (9, 3, 1, 8),
    (6, 4, 1, 6),
    (9, 4, 1, 8),
    (2, 3, 1, 2),
    (11, 2, 1, 10),
    (2, 4, 1, 2),
    (3, 6, 1, 5),
    (9, 1, 1, 7),
    (6, 1, 1, 6),
    (9, 1, 1, 8),
    (2, 9, 1, 6),
    (9, 2, 1, 7),
    (2, 1, 1, 2),
    (6, 2, 1, 6),
    (9, 2, 1, 8),
    (2, 8, 1, 2),
    (2, 10, 1, 2),
    (2, 2, 1, 2),
    (3, 5, 1, 5),
    (9, 6, 1, 7),
    (10, 3, 1, 9),
    (6, 6, 1, 6),
    (3, 3, 1, 3),
    (2, 3, 1, 6),
    (9, 6, 1, 8),
    (1, 3, 1, 2),
    (3, 4, 1, 3),
    (2, 4, 1, 6),
    (1, 4, 1, 2),
    (2, 6, 1, 2),
    (2, 7, 1, 2),
    (10, 1, 1, 9),
    (3, 1, 1, 3),
    (2, 1, 1, 6),
    (9, 5, 1, 7),
    (1, 1, 1, 2),
    (6, 5, 1, 6),
    (2, 8, 1, 6),
    (9, 5, 1, 8),
    (2, 10, 1, 6),
    (10, 2, 1, 9)
]

vs = [#
    (5, 2),
    (4, 2),
    (3, 2),
    (4, 5),
    (6, 2),
    (6, 1),
    (5, 4),
    (4, 6),
    (3, 6),
    (4, 7),
    (10, 3),
    (4, 9),
    (6, 6),
    (6, 5),
    (5, 5),
    (4, 5),
    (6, 5),
    (3, 5),
    (6, 6),
    (10, 2),
    (9, 6),
    (11, 1),
    (9, 5),
    (11, 2),
    (4, 3),
    (10, 1),
    (4, 4),
    (6, 4),
    (6, 2),
    (9, 2),
    (10, 1),
    (4, 9),
    (6, 1),
    (4, 1),
    (9, 1),
    (10, 2),
    (4, 8),
    (4, 10),
    (4, 2),
    (6, 3),
    (6, 3),
    (11, 3),
    (9, 4),
    (5, 1),
    (4, 3),
    (10, 3),
    (3, 3),
    (5, 6),
    (4, 4),
    (3, 4),
    (4, 6),
    (4, 7),
    (11, 2),
    (5, 3),
    (4, 1),
    (6, 4),
    (3, 1),
    (9, 3),
    (4, 8),
    (11, 3),
    (4, 10),
    (11, 1)
]

D_C = Dict(zip(ks, vs))

src_mets = [1, 2, 8]
max_src_mets_carbon = [6, 10, 0]

I = (src_mets[1], 1, atom)

ks = Int16.([5, 4, 6, 7, 2, 8, 3, 1])
vs = [(Int16(6), Int16(6)), (Int16(0), Int16(0)), (Int16(9), Int16(4)), (Int16(11), Int16(2)), (Int16(3), Int16(1)), (Int16(10), Int16(1)), (Int16(5), Int16(3)), (Int16(1), Int16(1))]
dmc = Dict(zip(ks, vs))

ks = [#
    Int16.([1, 2, 5, 6, 7]),
    Int16.([1, 2, 5, 6]),
    Int16.([1, 2, 5, 6, 7, 4]),
    Int16.([1]),
    Int16.([1, 2, 3]),
    Int16.([1, 2, 3, 4]),
    Int16.([1, 2]),
    Int16.([1, 2, 5]),
    Int16.([1, 2, 5, 6, 7, 8]),
    Int16.([1, 2, 5, 6, 7, 8, 4])
]
vs = [#
    (id = 7, children = Int16.([8])),
    (id = 6, children = Int16.([7])),
    (id = 8, children = Int16.([])),
    (id = 1, children = Int16.([2])),
    (id = 3, children = Int16.([])),
    (id = 4, children = Int16.([])),
    (id = 2, children = Int16.([3, 5])),
    (id = 5, children = Int16.([6])),
    (id = 9, children = Int16.([])),
    (id = 10, children = Int16.([]))
]
dchmc = Dict(zip(ks, vs))

T = SparseMatrixCSC([
    0 1.0 0 0 0 0 0 0 0 0
    0 0 0.3 0 0.7 0 0 0 0 0
    0 0 0 1.0 0 0 0 0 0 0
    1.0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 1.0 0 0 0 0
    0 0 0 0 0.125 0 0.875 0 0 0
    0 0 0 0 0 0 0 0.875 0.125 0
    1.0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0.125 0 0 0.875
    1.0 0 0 0 0 0 0 0 0 0
])

R = [#
    (i = 1, j = 2, k = Int16(2)),
    (i = 2, j = 3, k = Int16(3)),
    (i = 2, j = 5, k = Int16(5)),
    (i = 3, j = 4, k = Int16(4)),
    (i = 4, j = 1, k = Int16(4)),
    (i = 5, j = 6, k = Int16(6)),
    (i = 6, j = 5, k = Int16(7)),
    (i = 6, j = 7, k = Int16(8)),
    (i = 7, j = 8, k = Int16(12)),
    (i = 8, j = 1, k = Int16(12)),
    (i = 7, j = 9, k = Int16(10)),
    (i = 9, j = 10, k = Int16(11)),
    (i = 10, j = 1, k = Int16(11)),
    (i = 9, j = 7, k = Int16(9))
]

res = CHMCAtomicSummary(#
    (1, 1, :C),
    [#
        (EFM = [6, 5, 6], Closures = [(6, 5)]),
        (EFM = [4, 1, 2, 5, 6, 7, 4], Closures = [(8, 1)]),
        (EFM = [4, 1, 2, 3, 4], Closures = [(4, 1)]),
        (EFM = [8, 7, 8], Closures = [(9, 7)]),
        (EFM = [4, 1, 2, 5, 6, 7, 8, 4], Closures = [(10, 1)])
    ],
    [0.09000000000000001, 0.56, 0.27, 0.01, 0.07],
    [1.0, 6.222222222222222, 3.0, 0.11111111111111112, 0.7777777777777778],
    dmc,
    dchmc,
    T,
    R
)

T = SparseMatrixCSC([
    0 1.0 0 0 0 0 0 0 0 0
    0 0 0.5 0 0.5 0 0 0 0 0
    0 0 0 1.0 0 0 0 0 0 0
    1.0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 1.0 0 0 0 0
    0 0 0 0 0.5 0 0.5 0 0 0
    0 0 0 0 0 0 0 0.5 0.5 0
    1.0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0.5 0 0 0.5
    1.0 0 0 0 0 0 0 0 0 0
])

res_enum = CHMCAtomicSummary(#
    (1, 1, :C),
    [#
        (EFM = [6, 5, 6], Closures = [(6, 5)]),
        (EFM = [4, 1, 2, 5, 6, 7, 4], Closures = [(8, 1)]),
        (EFM = [4, 1, 2, 3, 4], Closures = [(4, 1)]),
        (EFM = [8, 7, 8], Closures = [(9, 7)]),
        (EFM = [4, 1, 2, 5, 6, 7, 8, 4], Closures = [(10, 1)])
    ],
    nothing,
    nothing,
    dmc,
    dchmc,
    T,
    R
)

efm_seq_1 = ["FDP", "F6P", "FDP"]
efm_seq_2 = ["Glc", "G6P", "F6P", "FDP", "DHAP"]
```

## Inputs

```julia
using MarkovWeightedEFMs
S = [#
  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 # Glc
  0 -1  0  0  0 -1  0  0  0  0  0  0  1  0  0  0 # ATP
  0  1 -1  0 -1  0  0  0  0  0  0  0  0  0  0  0 # G6P
  0  1  0  0  0  1  0  0  0  0  0  0  0 -1  0  0 # ADP
  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0 # 6PG
  0  0  0  0  1 -1  1  0  0  0  0  0  0  0  0  0 # F6P
  0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1  0 # Pi
  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  1 # H2O
  0  0  0  0  0  1 -1 -1  0  0  0  0  0  0  0  0 # FDP
  0  0  0  0  0  0  0  1 -1  1 -1  0  0  0  0  0 # G3P
  0  0  0  0  0  0  0  1  1 -1  0 -1  0  0  0  0 # DHAP
]

v = [10, 10, 3, 3, 7, 8, 1, 7, 1, 1, 7, 7, 18, 18, 1, 1]

mets = [#
  "Glc",
  "ATP",
  "G6P",
  "ADP",
  "6PG",
  "F6P",
  "Pi",
  "H2O",
  "FDP",
  "G3P",
  "DHAP"
]

rxns = [#
  "Source Glc",
  "Hexokinase",
  "G6P dehydrogenase",
  "Sink 6PG",
  "Phosphoglucose isomerase",
  "6-phosphofructo-1-kinase",
  "Fructose 1,6-bisphosphatase",
  "Fructose-bisphosphate aldolase",
  "Triose phosphate isomerase",
  "Triose phosphate isomerase",
  "Sink G3P",
  "Sink DHAP",
  "Source ATP",
  "Sink ADP",
  "Sink Pi",
  "Source H2O"
]

smiles = [#
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)OP(=O)(O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
  "O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O",
  "C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)OP(=O)(O)O",
  "[O-]P(=O)([O-])[O-]",
  "O",
  "C(C1C(C(C(O1)(COP(=O)(O)O)O)O)O)OP(=O)(O)O",
  "C([C@H](C=O)O)OP(=O)(O)O",
  "C(C(=O)COP(=O)(O)O)O"
]

atom = :C # carbon atom type for AEFMs
```

We can check that the flux vector satisfies the steady state requirements.

```@example required
all(S * v .== 0) # should evaluate as true
```

## Pre-processing steps

### Checking network structure

The following functions check for issues with the inputs. The first function
`find_atomic_chmc_input_errors` identifies possible problems with the
stoichiometry matrix and flux vector.

```julia
# Confirm there are no issues with stoichiometry matrix
errors = find_atomic_chmc_input_errors(S, v)
```

```@example required
print(errors) # summary of errors associated with S/v
```

### Correcting problems in network structure

Any problems, except for the steady state flux requirement, can be
addressed via `correct_atomic_chmc_input_errors`. 

```julia
# S and v have no errors so the inputs are returned
correct_atomic_chmc_input_errors(errors, S, mets, rxns)
# S, mets, rxns = correct_atomic_chmc_input_errors(errors, S, mets, rxns) # otherwise
```

### Identifying unmappable reactions

The next function `correct_atomic_chmc_input_smiles` checks and fixes problems
relating to the SMILES strings. These problems are caused by RXNMapper being
unable to map atoms in reactions with pseudometabolites or pseudoreactions with
non-integer stoichiometries (e.g. biomass reaction). RXNMapper also has
a character limit on reaction SMILES strings. These unmappable reactions are
removed and the flux is balanced with unimolecular flux entering/exiting the
associated reaction substrates/products.

```julia
# Correct issues associated with RXNMapper character limit,
# pseudometabolites and pseudoreactions
S, v, mets, rxns, smiles, logs = correct_atomic_chmc_input_smiles(S, v, mets, rxns, smiles)
```

At this point, the SMILES strings (matching the updated `mets` if there were
errors in the initial inputs) should be canonicalized. `S` is also converted
to a `Matrix{Int16}` which is a requirement for subsequent functions.

```julia
smiles = canonicalize_smiles(smiles)
```
```@example required
smiles
```

### Atom mapping reactions

The reaction SMILES strings `rs` are next constructed from the metabolite
SMILES and the atom mapping is performed via RXNMapper and stored in `ms`.

```julia
# Construct atom traced SMILES strings
rs, ms = map_reaction_strings(S, smiles, rxns, false)
```
```@example required
rs
```
```@example required
ms
```

### Identifying all source metabolite-atom positions

The following code extracts the source metabolite indices in `mets` and
computes the total number of carbon atoms of interest.

```julia
# Total number of atom type across all metabolites
atom_max = get_max_atoms(smiles, atom)

# Identify source metabolite indices and copies of atom
src_mets = get_source_metabolites(S)

# Number of carbon atoms in each source metabolite
max_src_mets_carbon = atom_max[src_mets]
```

```@example required
# Source metabolites
mets[src_mets]
```

```@example required
# Carbons in each source metabolite
max_src_mets_carbon
```

### Enumerating metabolite-atom mappings across reactions

We then precompute an atom tracing dictionary mapping the (carbon) atom in
the stoichiometric copy of a substrate to its product atom position across
each reaction.

```julia
# Precompute atom tracing dictionary
D_C = precompute_atom_tracing_dictionary(S, ms, atom_max, atom) # S must be Matrix{Int16}
```
```@example required
D_C
```

## Conclusion

All of the code/functions described above are wrapped into
`preprocess_all_for_atomic_chmc()`. If this wrapper function fails, you
may need to step through these individual pre-processing functions to
identify the error.



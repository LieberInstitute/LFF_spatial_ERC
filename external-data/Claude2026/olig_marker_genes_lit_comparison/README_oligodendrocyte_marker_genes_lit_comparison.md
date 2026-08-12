# README — Oligodendrocyte Marker Gene Literature-Comparison Workbook

This document describes how `oligodendrocyte_marker_genes_lit_comparison.xlsx` was built,
what each column means, and how it should (and shouldn't) be used. It was written
alongside the workbook so the two should be kept together and updated together.

## What this file is

A citable comparison table of oligodendrocyte-lineage marker genes across the published
literature — built up incrementally across several review passes, each adding one or more
papers as a new column. The goal is to let someone scan a single gene and see, at a glance,
which papers use it, for which cell-type/stage, and how much agreement exists across sources.

Current size: **55 gene rows × 22 columns**, plus a second sheet listing full citations for
the 16 papers represented.

## Sheet 1: "Oligo Marker Genes"

### Column order and meaning

| Column | Purpose |
|---|---|
| `Gene` | Human gene symbol only. If a mouse-specific gene name/paralog mismatch matters, that's explained in `Notes`, not here. |
| `BroadCategory` | Exactly one of: `Broad identity` (pan-lineage, not stage-specific), `Maturation identity` (OPC→COP→NFOL→myelinating→mature), `Disease relevance` (immune/DAO/MHC/stress/amyloid). Color-banded (peach / blue / green). |
| `Category` | Finer-grained, free-text stage/context label, e.g. "OPC (regional subtype, mouse)" or "COP / premyelinating." |
| `Category_short` | A deliberately small, fixed vocabulary for plotting/faceting (see below) — collapses the ~30 distinct `Category` strings down to 10 buckets. |
| `Confidence` | `Low`/`Medium`/`High` — how confident this table's own evidence is that the gene is an effective marker specifically for **human** oligodendrocytes (see methodology below). |
| `Consensus` | One-line summary of where the field agrees/disagrees on this gene. |
| `Notes` | The longest field — cross-paper biological interpretation: dual roles, species caveats, disputed function, discrepancies between papers, secondhand-citation warnings. |
| *(remaining 15 columns)* | One column per paper, in the order it was added to the table. Each cell is either a short phrase naming the subtype/cluster where the gene is high, or `NA` if that paper doesn't use the gene. The full citation lives in a **cell comment on the header** (hover/right-click the header to see it), not in the header text itself — headers are kept as code-friendly single tokens (e.g. `Siletti2023`). |

`Notes` is placed immediately before the paper columns by design, so a reader hits the
qualitative synthesis right before the raw per-paper evidence grid that supports it.

### `Category_short` bucket definitions

Ten buckets, chosen to be plot-friendly (short, mutually exclusive, few enough to color/facet by):

- **OPC** (13 genes) — precursor-stage markers, including CSPG-module genes and regional-subtype markers
- **COP** (5) — committed OPC / early-differentiating markers
- **NFOL** (4) — newly-formed-oligodendrocyte markers
- **Myelinating** (8) — active myelination-stage / structural myelin genes
- **Mature** (8) — terminal/mature-oligodendrocyte markers
- **Pan-lineage** (2) — OLIG2, SOX10 (span most/all stages)
- **Immune/DAO** (6) — disease-associated-oligodendrocyte / immune signature genes
- **MHC/Interferon** (3)
- **Stress** (2)
- **Amyloid** (4) — amyloid-processing genes discussed in the AD literature

A few genes were force-fit into the nearest bucket rather than getting a bespoke one, to keep
the vocabulary small (e.g. `ACAN`, whose `Category` text is a mouthful — "OPC /
differentiation-attempt ECM marker" — collapses to simply `OPC` in `Category_short`).

### `Confidence` methodology

This column answers: *given only what's captured in this table, how confident should you be
that this gene is a good marker for human oligodendrocytes?* It is **not** an independent
literature review — it's a derived score computed from two signals already present in the
row:

1. **Breadth** — how many of the 15 paper columns have a non-`NA` value for that gene.
2. **Human-model support** — whether at least one of those citations comes from a
   human-derived dataset. `Jakel2019`, `Sadick2022`, and `Siletti2023` count as confirmed-human
   sources (their data is human tissue/atlases). `Zhou2020`, `Falcao2018`, and `Kenigsbuch2022`
   count as *partial*-human (each integrates some human data per its own title/abstract, but
   none of the three was reviewed directly this session — flagged accordingly). Cells whose
   own text explicitly says "mouse" are **not** counted as human support even if they appear
   in an otherwise-human paper's column (e.g. a mouse-reanalysis finding reported inside the
   Siletti et al. 2023 paper doesn't count as human evidence for that specific claim).

Rule applied:

- **High** = human-confirmed AND cited in ≥3 papers, OR partial-human AND cited in ≥5 papers
- **Medium** = solid support on only one axis (e.g. a single human dataset, or broad
  mouse-only citation across ≥3 papers)
- **Low** = single-source and/or mouse-only support, or a disqualifying caveat

Current distribution: **11 High, 21 Medium, 23 Low**.

Two manual judgment calls layered on top of the mechanical rule:

- **`CD44`** was forced to `Low` even though it technically qualifies for `Medium` (one
  human-dataset citation, from Siletti et al. 2023) — because that citation itself is
  ambiguous: the paper shows CD44 in a figure panel alongside OPALIN/RBFOX1 but never states
  which of the two mature oligodendrocyte types has higher expression. A marker whose
  direction is unconfirmed shouldn't score better than a plain single-source gene.
- **`SOX10` is a known limitation, not a real signal.** It scores `Low` in this table purely
  because only 2 of 15 paper columns happen to be populated for it (`Hughes2021`,
  `Emery2024` — both general reviews, neither flagged as a "human dataset" column), even
  though SOX10 is one of the most universally-used human oligodendrocyte-lineage markers in
  the broader field. This table's `Confidence` score is only as good as its own per-paper
  cell curation, which was not built to be an exhaustive marker panel for every paper —
  several papers' columns were populated selectively for the genes most relevant to this
  table's running narrative, not systematically for every gene they mention. **Treat `Low`/
  `Medium` here as "under-supported within this specific comparison," not as "weak marker in
  general," and cross-check against the wider literature before drawing a conclusion from
  this column alone.**

Full methodology text is also attached as a cell comment on the `Confidence` header itself,
so it travels with the file even if this README is separated from it.

## Sheet 2: "References"

One row per paper, with columns `ColumnCode | Citation | Title | Journal / Details | DOI |
Source basis in this workbook`. The last column states plainly whether the paper was
**reviewed directly** (PDF uploaded and read this session) or whether its marker genes were
**taken secondhand** from another paper's summary table — currently, nine of the papers
(`Marques2016`, `Zeisel2018`, `Falcao2018`, `Zhou2020`, `Lee2021`, `Sadick2022`,
`Kenigsbuch2022`, `Pandey2022`, `Kaya2022`) were pulled from Table 2 of a mini-review
(Valihrach et al. 2022) rather than read directly, and are flagged "verify against the
primary paper before citing." The remaining seven (`Hughes2021`, `Jakel2019`, `Kedia2025`,
`Siletti2023`, `Mironova2026`, `Emery2024`, plus Valihrach et al. 2022 itself as the
"(review source)" row) were reviewed directly from uploaded PDFs.

One DOI was corrected during this process: Mironova et al. 2026's DOI was initially entered
from memory and flagged `NOT VERIFIED`; once the PDF was uploaded and reviewed directly, the
DOI was confirmed against the paper's own "Cite this article as" line and the flag removed.

## Build history

The workbook was built incrementally, each pass adding one or more newly-reviewed papers:

1. **Starting point**: an existing workbook (37 gene rows) covering `Hughes2021`, `Jakel2019`,
   and `Kedia2025` (reviewed directly), plus nine papers pulled secondhand from Valihrach et
   al. 2022.
2. **+ DOI column** added to the References sheet (external edit), with Mironova2026's DOI
   flagged as unverified pending direct review.
3. **+ Siletti et al. 2023 and Mironova et al. 2026** (both reviewed directly from uploaded
   PDFs): updated 5 existing rows (`CSPG4`, `BCAN`, `BCAS1`, `OPALIN`, `RBFOX1`) with new
   per-paper evidence and cross-paper notes, and added 9 new rows (`VCAN`, `PTPRZ1`, `ACAN`,
   `HES5`, `IRX5`, `NELL1`, `GPC5`, `CD44`, `ASPA`). Mironova2026's DOI was confirmed directly
   from the reviewed PDF and the `NOT VERIFIED` flag removed.
4. **+ Emery & Wood 2024** (reviewed directly): updated 17 existing rows with corroborating
   mechanistic context and citations, and added 9 new rows (`SOX5`, `MYRF`, `TCF7L2`,
   `ITPR2`, `NKX2-2`, `PCDH15`, `CASR`, `APOD`, `QKI`). `SOX5` was added specifically because
   the existing table had `SOX6` (independently IHC-validated in human tissue by Jäkel et al.
   2019) but not its obligate SoxD-family partner `SOX5`, which only has mechanistic
   (mouse-based) support in this table.
5. **+ `Category_short` column** — collapsed the free-text `Category` field into 10 fixed
   buckets for plotting, and moved `Notes`/`Consensus` earlier in the column order.
6. **+ `Confidence` column** — derived Low/Medium/High per gene from citation breadth and
   human-model support (methodology above), and did a final column reorder so `Notes` sits
   immediately before the paper-column data grid.

## Formatting conventions

- Arial 10pt throughout; bold header row (white text on dark blue fill), bold `Gene` /
  `BroadCategory` / `Category_short` / `Confidence` / `Consensus` cells for scannability.
- `BroadCategory` cells are color-banded: peach = `Broad identity`, blue = `Maturation
  identity`, green = `Disease relevance`.
- Wrap text + top vertical alignment on every cell.
- Header row frozen, plus all metadata columns (`Gene` through `Notes`) frozen so the paper
  columns can be scrolled independently.
- Autofilter enabled on both sheets.
- Paper column headers are single code-friendly tokens (`AuthorYYYY`); the full citation is
  a cell comment on that header, not in the header text.

## Known limitations / things to check before citing

- **Secondhand columns aren't verified.** Nine paper columns were populated from a
  third-party summary table, not the primary papers. Treat their per-gene claims as
  provisional until checked against the original source (flagged per-row in `References`).
- **`Confidence` reflects this table's coverage, not the full literature** — see the SOX10
  caveat above. A `Low` score here can mean "genuinely weak marker" or "this table just
  hasn't captured much evidence for it yet." Check `Notes`/`Consensus` and, ideally, the
  primary papers before treating `Confidence` as a standalone verdict.
- **A few specific per-gene caveats worth remembering**: `CD44`'s Oligo1/Oligo2 direction is
  unconfirmed in the source text (figure-only); `NELL1`/`CASR`/`APOD` trace to citations
  nested inside a review rather than a primary source; one `RBFOX1` subcluster in Siletti et
  al. 2023 was flagged by the authors themselves as likely low-quality/donor-derived; the
  cluster-number labels "Oligo1"/"Oligo2" for `OPALIN`/`RBFOX1` are **not** consistent between
  Jäkel et al. 2019 and Siletti et al. 2023, even though both use the same two marker genes —
  cite the marker gene, not the cluster number, when comparing across studies.

## How to extend this table further

To add another paper: review it directly (don't rely on secondhand summaries if avoidable),
add a new single-token column at the end of the paper-column block with a header cell
comment holding the full citation, fill in `NA` for genes it doesn't discuss, add any new
genes it introduces as new rows (assign `BroadCategory`, `Category`, `Category_short`), and
recompute `Confidence` for any row whose breadth or human-support signal changed as a result.
Add a corresponding row to the `References` sheet with an accurate DOI and an honest
"reviewed directly" vs. "secondhand — verify" note in the source-basis column.

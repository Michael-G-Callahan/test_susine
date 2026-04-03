# Real Data Pipeline Reset Notes

Date: 2026-04-02

## Purpose

Document the conclusions from the April 2 discussion about the CBX8 real-data workflow, the variant-count funnel across the existing notebooks, the key architectural mistake in the historical pipeline, and the agreed plan for rebuilding the real-data case study pipeline.

---

## Executive Summary

The main issue in the historical CBX8 workflow is that the z-score export was filtered by Borzoi feasibility before LD and annotation processing. This made annotation feasibility an upstream gatekeeper for the entire downstream analysis, which is the wrong dependency structure.

The corrected logic should be:

- master analysis universe = `z ∩ LD`
- annotations are joined onto that universe
- variants without annotation should receive annotation value `0`
- annotation availability should be recorded with flags, not used as a hard inclusion criterion

This change should be implemented first on CBX8 before deciding whether to migrate from Borzoi to AlphaGenome.

---

## Historical CBX8 Funnel

These counts were reconstructed from:

- `C:\Users\mgcal\Downloads\get_z_scores.ipynb`
- `C:\Users\mgcal\Downloads\cbx8_1000g_ld.ipynb`
- `C:\Users\mgcal\Downloads\Get_annotations.ipynb`
- saved outputs in `eQTL_annotations_for_susine/output/`

### 1. Raw GTEx source

Source notebook: `get_z_scores.ipynb`

- Raw GTEx variants loaded for CBX8 in lung: `9791`

### 2. GTEx postfilter set

Source notebook: `get_z_scores.ipynb`

- After AF filter (`0.01 <= af <= 0.99`): `9466`
- After effective sample-size filter (`sample_size > 50`): `9224`

Losses:

- removed by AF filter: `325`
- removed by sample-size filter: `242`

Interpretation:

- `9224` is the first sensible "real" locus-wide GTEx set
- this is the right starting point for a rebuilt pipeline

### 3. Borzoi feasibility screen

Source notebook: `get_z_scores.ipynb`

From the `9224` postfilter GTEx variants:

- both methods feasible: `757`
- only TSS-centered feasible: `1501`
- only SNP-centered feasible: `0`
- neither feasible: `6966`

Exported with `valid_any = TRUE`:

- `2258`

Interpretation:

- this was the dominant bottleneck in the historical workflow
- Borzoi feasibility reduced the locus from `9224` to `2258`
- this step happened too early

### 4. Exported VCF and GTEx CSV

Source notebook: `get_z_scores.ipynb`

Produced:

- `CBX8_variants.vcf`: `2258` variants
- `CBX8_GTEx_z_scores.csv`: `2258` variants

Important implication:

- these are not raw GTEx-locus files
- both already inherit the Borzoi-feasible subset

### 5. Annotation input

Source notebook: `Get_annotations.ipynb`

The annotation notebook loads `CBX8_variants.vcf`, so annotation input count was:

- input to annotation pipeline: `2258`

Interpretation:

- annotation coverage was limited by the VCF exported from the z-score/Borzoi-feasibility notebook

### 6. Annotation SNP-only filter

Source notebook: `Get_annotations.ipynb`
Source code: `eQTL_annotations_for_susine/utils/data_ingestion.py`

The pipeline explicitly filters to variants with `len(ref) == 1` and `len(alt) == 1`.

- non-SNP variants removed: `148`
- SNPs retained for annotation: `2110`

Interpretation:

- this is a wrapper limitation, not a fundamental AlphaGenome/Borzoi limitation
- `grelu` itself appears capable of more than simple SNP substitution, but the current wrapper is SNP-only

### 7. Annotation output

Source notebook: `Get_annotations.ipynb`

- variants analyzed: `2110`
- variants with valid predictions: `2110`
- `CBX8_variant_effects.csv`: `2110`

Interpretation:

- once variants survived the SNP-only filter, annotation prediction itself was not the limiting step for CBX8

### 8. LD input

Source notebook: `cbx8_1000g_ld.ipynb`

The LD notebook also loads `CBX8_variants.vcf`, so LD input count was:

- input to LD matching: `2258`

Interpretation:

- LD was also limited by the Borzoi-feasible export, not by the full GTEx postfilter set

### 9. LD matching to 1000G

Source notebook: `cbx8_1000g_ld.ipynb`

- matched variants in 1000G: `2170 / 2258`
- missing from 1000G VCF: `88`

Produced:

- `CBX8_LD_variant_order.tsv`: `2170`
- `CBX8_LD_R.csv`: `2170 x 2170`

Interpretation:

- LD matching was fairly good on the exported subset
- however, this does not yet answer how much of the full `9224` GTEx set could be matched after normalization

### 10. Final harmonized set used in the R real-data pipeline

Source file: `test_susine/vignettes/real data pipeline/susine_real_case_study.Rmd`

That pipeline intersects:

- z-score IDs
- annotation IDs
- LD IDs

Using the saved outputs:

- z total: `2258`
- annotations total: `2110`
- LD total: `2170`
- all three: `2057`

Three-way overlap:

- z only: `35`
- z and annotations only: `53`
- z and LD only: `113`
- annotations only: `0`
- annotations and LD only: `0`
- LD only: `0`
- all three: `2057`

Interpretation:

- the final analyzable set for the historical CBX8 real-data pipeline was `2057`

---

## Dependency Structure of the Historical Pipeline

The pipeline actually behaved as:

`raw GTEx -> GTEx postfilter -> Borzoi feasibility filter -> shared VCF -> {annotation branch, LD branch} -> final intersection`

This means:

- the LD branch was limited by the Borzoi-feasible z-score export
- the annotation branch was also limited by that same export
- therefore the whole real-data analysis was implicitly organized around annotation feasibility

This is the core architectural mistake.

The correct dependency structure should instead be:

`raw GTEx -> GTEx postfilter -> LD matching -> master z∩LD universe -> optional annotation join`

---

## Main Conclusions from the Discussion

### 1. The current pipeline should be rebuilt

CBX8 should be rerun without using Borzoi feasibility as a hard prefilter on which variants enter the locus.

### 2. Annotation should not define the universe

The correct analysis universe is:

- variants with usable z-scores
- variants with usable LD

Annotations should be supplementary information attached to that universe.

### 3. Non-annotatable variants should be retained

For the next pipeline version:

- if a variant is outside the annotation model's usable range, assign annotation `0`
- if a variant cannot be scored for technical reasons, assign annotation `0`
- always keep a reason flag describing why the annotation is zero or missing

### 4. Indels should not be excluded by default

The next version should not automatically limit to SNPs.

At minimum:

- SNPs and indels should remain in the `z ∩ LD` universe if they can be matched and harmonized
- annotation support for indels can be implemented separately
- variants that remain unannotated can still be retained with annotation `0`

### 5. AlphaGenome is worth piloting, but not yet worth migrating wholesale

Motivation:

- likely wider usable context than the current Borzoi workflow
- likely much lower local HPC burden if the API is fast enough
- free for non-commercial use, per current docs, but API throughput and workflow friction still need real validation

The correct approach is:

- first fix the pipeline architecture on CBX8
- then test AlphaGenome on CBX8 as an annotation-engine pilot

---

## Agreed Plan Moving Forward

### Phase 1. Rebuild CBX8 with the correct architecture

Goal:

- create a new CBX8 pipeline where annotation availability does not gate variant inclusion

Target logic:

1. start from the full GTEx postfilter set (`9224` for CBX8)
2. normalize variant representation consistently
3. intersect GTEx variants with LD-reference variants
4. use `z ∩ LD` as the master locus universe
5. join annotations onto that master universe
6. set annotation to `0` for variants that are out of range or unscored
7. preserve flags explaining annotation provenance

Deliverables:

- updated CBX8 master variant table
- updated CBX8 LD object
- updated CBX8 annotation table with zero-filled annotation for unsupported variants
- updated CBX8 real-case input bundle for SuSiNE

### Phase 2. AlphaGenome pilot on CBX8

Goal:

- determine whether AlphaGenome is viable as a replacement or supplement to Borzoi for real-data annotation generation

Questions to answer:

- how easy is the API to use in practice?
- what is the scoring throughput?
- what fraction of CBX8 variants can be scored?
- how does score coverage compare to Borzoi?
- on overlapping variants, how correlated are the annotation values?

Important note:

- do not migrate the whole study to AlphaGenome until CBX8 pilot results are in hand

### Phase 3. Run one additional locus through the rebuilt pipeline

Goal:

- test whether the rebuilt architecture generalizes beyond CBX8

Suggested constraints:

- use the same `z ∩ LD` master-universe logic
- allow annotation values of `0` for unsupported variants
- do not assume AlphaGenome yet unless the CBX8 pilot looks good

### Phase 4. Benchmark difficulty-metric cost

Goal:

- determine whether exact metric computation on full real-data loci is cheap enough for broad screening

Key concern:

- once the locus size rises toward `~9K` variants, dense LD-derived metrics become substantially more expensive

Questions to answer:

- how expensive is full-LD construction on loci of this size?
- how expensive are the desired LD-based difficulty metrics?
- can exact metrics be computed comfortably for 10-20 loci?
- if not, is approximate screening needed?

### Phase 5. Decide between manual locus selection and broader screening

Decision rule:

- if metric computation is expensive, hand-select only another `2-3` loci
- if metric computation is manageable, automate screening for roughly `10-20` loci and then down-select using the dataset metrics

---

## Recommended Master-Universe Design for the Rebuilt Pipeline

For future real-data loci, the preferred structure is:

### Step A. Define the z-score universe

From GTEx:

- choose gene and tissue
- apply basic postfilters only
- produce the locus-wide z-score table

### Step B. Define the LD-compatible universe

From 1000G or another LD reference:

- normalize variant representation
- match the GTEx variants to the LD panel

### Step C. Build the master universe

Use:

- `master = GTEx_postfilter ∩ LD_reference`

This should define the variants eligible for fine-mapping.

### Step D. Join annotations

For every variant in `master`:

- if annotation is available, attach it
- if annotation is unavailable, assign `0`
- retain columns describing:
  - whether the variant was inside the annotation window
  - whether the variant type was supported
  - whether scoring was attempted
  - whether scoring succeeded
  - why a zero-filled annotation was assigned

### Step E. Compute dataset metrics on the master universe

This ensures:

- LD metrics are computed on the same variant universe that will actually be analyzed
- z-score metrics are computed on the same set
- annotation coverage becomes a property of the universe, not a filter on the universe

---

## Practical Notes on Metric Cost

At approximately `p = 9000`, dense LD is expensive but still plausible for a small number of loci.

Rough implications:

- one dense `9000 x 9000` matrix is manageable
- repeated construction and repeated quadratic summaries across many loci may become the true bottleneck

Practical guidance:

- compute LD once per locus and derive all metrics from that same object
- avoid recomputing the dense matrix repeatedly
- if screening many loci becomes too expensive, consider approximate LD-metric screening first and reserve exact LD construction for promoted loci

---

## Status of the SuSiNE Fit Stage

Current expectation from the discussion:

- the SuSiNE fitting stage itself is not the main bottleneck for this phase
- the hard part is the data preparation:
  - z-score retrieval and harmonization
  - LD construction and metric computation
  - annotation generation

Therefore the current priority is pipeline/data engineering rather than model fitting.

---

## Immediate Next Actions

1. Rebuild the CBX8 real-data preprocessing so that the master universe is `z ∩ LD`, not Borzoi-feasible variants.
2. Preserve non-annotatable variants by assigning annotation `0` instead of dropping them.
3. Audit indel handling and harmonized variant representation early.
4. Run an AlphaGenome pilot on CBX8 only after the rebuilt CBX8 master table is in place.
5. Push one more locus through the same architecture before deciding on any larger-scale screening effort.


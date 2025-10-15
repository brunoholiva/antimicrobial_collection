# Data Sources and Assay Details

This file documents the sources and assay details for all datasets used in this repository.

---

## Swanson _et al._ 2025 
**Summary:**
- **Organism:** _S. aureus_ RN4220
- **Screened:** 10,716 compounds
- **Unique molecules after deduplication:** 10,658
- **Active compounds:** 1,137 (10.7%)
- **Inactive compounds:** 9,521 (89.3%)
- **Cutoff:** Mean and standard deviation of OD600, binarized using μ − 2σ threshold
- **Assay:** OD600 after 16h at 37ºC, 50 μM in LB, biological duplicate

**Source text:**
"To create our antibacterial activity training set, we physically screened 10,716 compounds
against S. aureus RN4220 at 50 M in LB medium (Fig. 2a and Supplementary Data 1)19 131 .
Experiments were conducted in biological duplicate, with end-point growth being measured at
600 nm optical density (OD600) after 16 hours of incubation at 37ºC (Extended Data Fig. 1a). We
computed the mean  and standard deviation  OD600 value across the library and used   2 as a threshold for binarizing the OD600 values into active (low OD600) and inactive (high OD600)
136 molecules. After removing duplicate compounds based on SMILES (see Methods), we obtained
a set of 10,658 unique molecules with 1,137 active compounds (10.7%) and 9,521 (89.3%)
inactive compounds. A t-SNE visualization of our training dataset compounds and 1,007
molecules from ChEMBL20 with known antibacterial activity (Supplementary Data 2) shows that
our active compounds cover both known and novel antibiotic chemical space (Fig. 2b)"

---

## Stokes _et al._ 2020
**Summary:**
- **Organism:** _E. coli_ BW25113
- **Screened:** 2,560 molecules (1,760 FDA-approved + 800 natural products)
- **Unique molecules after deduplication:** 2,335
- **Active compounds:** 120 (growth inhibitory)
- **Inactive compounds:** 2,215
- **Cutoff:** 80% growth inhibition
- **Assay:** Growth inhibition screen

**Source text:**
"We screened for
growth inhibition against E. coli BW25113 (Zampieri _et al._,
2017) using a widely available US Food and Drug Administration
(FDA)-approved drug library consisting of 1,760 molecules of
diverse structure and function. To further increase chemical
diversity, we included an additional 800 natural products
isolated from plant, animal, and microbial sources, resulting in
a primary training set of 2,560 molecules (Figures 2A and S1A;
Table S1A)—2,335 unique compounds when deduplicated
(Figure S1B; Table S1B). Using 80% growth inhibition as a hit
cut-off, this primary screen resulted in the identification of 120
molecules with growth inhibitory activity against E. coli.
Next, all 2,335 compounds from the primary training dataset
were binarized as hit or non-hit"

---

## Swanson _et al._ 2024 
**Summary:**
- **Organism:** _A. baumannii_ ATCC 17978
- **Screened:** Three libraries (2,371 + 6,680 + 5,376 molecules)
- **Unique molecules after deduplication:** 13,524
- **Active compounds:** 470
- **Inactive compounds:** 13,054
- **Cutoff:** μ − 2σ threshold for OD600
- **Assay:** OD600 after 16h at 50 μM, 100 μl, biological duplicate

**Source text:**
"we began by physically screening
three distinct chemical libraries to use as a training dataset. Library 1
and Library 2 are collections of bioactive compounds with 2,371 and
6,680 molecules, respectively. Library 3 is a synthetic commercially
available small-molecule screening collection with 5,376 molecules.

To acquire our training dataset, we grew A. baumannii ATCC 17978 
in the presence of each chemical at 50 µM, in a volume of 100 µl, in
biological duplicate. After 16 hours of incubation, we measured the
endpoint optical density at 600 nm (OD600). Next, for each library
separately, we computed the mean μ and standard deviation σ OD600
value across the library and used μ − 2σ as a threshold for binarizing
the optical density values into active and inactive molecules (Fig. 2a,
Extended Data Fig. 1a and Supplementary Data 1)23. We then merged
the three binarized libraries and removed duplicate compounds with
conflicting activity labels (Methods). This resulted in a combined set
of 13,524 unique molecules, with 470 active compounds and 13,054"
inactive compounds.

---

## Liu _et al._ 2023
**Summary:**
- **Organism:** _A. baumannii_ ATCC 17978
- **Screened:** 7,684 small molecules (FDA approved + Broad institute)
- **Unique molecules after deduplication:** 7,684
- **Active compounds:** 480
- **Inactive compounds:** 7,204
- **Cutoff:** μ − 2σ threshold for OD600
- **Assay:** OD600 after 16h at 50 μM, 100 μl, biological duplicate

**Source text:**
"_A. baumannii_ ATCC 17978 was grown in 2 ml LB medium (Becton, Dickinson and Company) overnight at 37 °C with shaking. Cells were diluted 1/10,000 into fresh LB and 99 µl of cells was added to each well of a 96-well flat-bottom plate (Corning). Next, 1 µl of a 5 mM stock of each molecule from a collection of 7,684 small molecules (FDA-approved drugs and molecules from screening collections from the Broad Institute) was added, in duplicate, using an Agilent Bravo liquid handling system. The final concentration was 50 µM. Plates were incubated in sealed plastic bags at 37 °C for 16 h and read at 600 nm using a SpectraMax M3 plate reader. Plate data were normalized based on the interquartile mean of each plate before binarization into ‘active’ and ‘nonactive’ categories for model training. Active molecules were defined as those that resulted in growth at least 1σ below the mean growth of the entire dataset."

---

## Wong _et al_ 2024
**Summary:**
- **Organism:** _S. aureus_ RN4220
- **Screened:** 39,312 molecules
- **Active compounds:** 512 (1.3%)
- **Inactive compounds:** 38,800
- **Cutoff:** 80% normalized control growth
- **Assay:** OD600 after 16h at 50 μM, 100 μl, biological duplicate

**Source text:**
"We first screened an original set of 39,312 compounds containing most known antibiotics, natural products, and structurally diverse molecules, with molecular weights between 40 Da and 4,200 Da, for growth-inhibitory activity against a methicillin-susceptible strain, S. aureus RN4220 (Fig. 1b, Extended Data Fig. 1 and Supplementary Data 1). These compounds were screened for overnight growth-inhibitory activity in nutrient-rich medium at a final concentration of 50 μM, and their effects were binarized as active or inactive using an 80% normalized growth inhibition cut-off, resulting in a total of 512 active compounds (1.3% of all compounds)."

---

## CO-ADD
- Downloaded directly from the Community for Open Antimicrobial Drug Discovery ([single-dose dataset, downloaded September 9, 2025](https://www.co-add.org/))
- Not in the repo due to size constraints.

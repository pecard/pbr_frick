# Parametrisation for the Leslie-matrix comparison

This note documents the literature basis for the Leslie/Lefkovitch matrix
built to cross-check the Frick et al. (2026) PBR/lambda_max approach
against a simpler, standard demographic method, using the same adult
survival (s) and age at first breeding (alpha) already established for
the Parti-coloured bat analysis.

**Access limitation (2026-09):** outbound web access in this session could
only reach search-engine result snippets (`WebSearch`), not the full text
of the papers themselves (`WebFetch` was blocked for every domain tried,
including non-paywalled ones) -- so the numbers below are as reported in
search snippets/abstracts, not independently checked against the full
papers. Treat them as a reasonable starting point, not a final,
fully-verified citation; confirm the exact figures directly if this
parametrisation is to be relied upon beyond an internal, exploratory
comparison.

## What the Leslie matrix needs beyond the PBR inputs

The Niel & Lebreton demographic-invariant approximation used by Frick et
al. (and hence our PBR analysis) only requires adult survival (s) and age
at first breeding (alpha) -- it is deliberately designed to avoid needing
an explicit juvenile survival rate or a full fecundity schedule. A Leslie
matrix has no such shortcut: it needs a first-year (juvenile) survival
rate and a female fecundity value to be specified explicitly.

## Juvenile (first-year) survival

No species-specific figure for *Vespertilio murinus* was found. As a
proxy, two studies of *Myotis lucifugus* (little brown bat, also a
temperate vespertilionid, well studied demographically) were used, in the
same spirit as the Frick et al. framework's own use of *Pipistrellus
pipistrellus* as a cross-species density proxy where species-specific data
are unavailable:

| Source | Adult female survival | First-year/juvenile survival |
|---|---|---|
| Frick, W.F. et al. (2010). Influence of climate and reproductive timing on demography of little brown myotis *Myotis lucifugus*. *Journal of Animal Ecology*. | 0.63-0.90 | 0.23-0.46 |
| Schorr, R.A. et al. (2021). Population dynamics of little brown bats (*Myotis lucifugus*) at summer roosts: apparent survival, fidelity, abundance, and the influence of winter conditions. *Ecology and Evolution*. | (not extracted) | 0.45 (SE 0.06) and 0.71 (SE 0.09) at two separate roosts |

Combined, these give a first-year survival range of roughly **0.23-0.71**,
consistent with the qualitative pattern confirmed independently by
Sendor & Simon (2003, *Pipistrellus pipistrellus*, Germany: adult survival
~0.80 +/- 0.05, first-year survival explicitly lower) and by Lentini et
al. (2015, global microbat synthesis, already cited in
`alpha_first_breeding_evidence.md`: juveniles have lower apparent survival
than adults across species). No single narrower "best" value stood out
from the available snippets, so the Leslie matrix comparison treats
juvenile survival as its own explored range (0.23-0.71), the same way
`s` and `alpha` are treated as ranges rather than fixed points elsewhere
in this analysis, with 0.45 (the more precisely reported Schorr et al.
estimate) as a central working value.

## Fecundity

*Vespertilio murinus* litter size is 2, one litter per year (AnAge
database; already used in `alpha_first_breeding_evidence.md`). Assuming
an even sex ratio at birth, this gives **1 female pup per breeding female
per year** -- the fecundity value used once an individual reaches age
`alpha`.

## Matrix structure

Pre-breeding-census, female-only Leslie matrix, annual time step, ages
`0, 1, 2, ..., max_age`:

- Survival age 0 -> 1: the juvenile rate above (explored range).
- Survival age >= 1: the adult rate `s` (same range as the PBR analysis:
  0.70-0.95).
- Fecundity at age `a`: 0 for `a < alpha`; 1 female pup per female per
  year (accounting for that pup's own survival to the next census, per
  the standard pre-breeding-census convention) for `a >= alpha`. Since
  `alpha` can be non-integer (e.g. 1.5), the age class spanning `alpha`
  gets a fecundity pro-rated by the fraction of that year already at or
  past maturity, and full fecundity applies from the next whole age class
  onward.
- `max_age` is set high enough (15 years) that the dominant eigenvalue is
  insensitive to further truncation, given survival and fecundity are
  both constant from age `alpha` onward.

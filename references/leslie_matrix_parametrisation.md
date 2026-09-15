# Parametrisation for the Leslie-matrix comparison

This note documents the literature basis for the Leslie/Lefkovitch matrix
built to cross-check the Frick et al. (2026) PBR/lambda_max approach
against a simpler, standard demographic method, using the same adult
survival (s) and age at first breeding (alpha) already established for
the Parti-coloured bat analysis.

**Access limitation (2026-09):** outbound web access in this session could
only reach search-engine result snippets (`WebSearch`), not the full text
of the papers themselves (`WebFetch` was blocked for every domain tried,
including non-paywalled ones). The species-specific parameters below (Safi
2006 survival rates) were located by Paulo via a separate search engine
and could only be partially corroborated here: Safi (2006), "Die
Zweifarbfledermaus in der Schweiz: Status und Grundlagen zum Schutz"
(Haupt Verlag), is confirmed as a real monograph on *V. murinus* in
Switzerland, cited by several other *V. murinus* papers, but the exact
survival figures themselves could not be independently checked against
the full text. Zhigalin & Moskvitina's fecundity figures were
independently corroborated (see below) -- with one correction to Paulo's
citation: the paper is from **2017**, not 2007.

## What the Leslie matrix needs beyond the PBR inputs

The Niel & Lebreton demographic-invariant approximation used by Frick et
al. (and hence our PBR analysis) only requires adult survival (s) and age
at first breeding (alpha) -- it is deliberately designed to avoid needing
an explicit juvenile survival rate or a full fecundity schedule. A Leslie
matrix has no such shortcut: it needs a first-year (juvenile) survival
rate and a female fecundity value to be specified explicitly.

## Primary parametrisation: species-specific (Safi 2006; Zhigalin & Moskvitina 2017)

| Parameter | Class | Value | Source |
|---|---|---|---|
| Survival | juvenile female | 0.62 | Safi (2006) |
| Survival | juvenile male | 0.62 | assumed = juvenile female |
| Survival | adult female | 0.76 | Safi (2006) |
| Survival | adult male | 0.42 | Safi (2006) |
| Fraction reproducing | adult female | 0.87 | Safi (2006) |
| Litter size | suburban colonies | 1.8 | Zhigalin & Moskvitina (2017) |
| Litter size | urban colonies | 2.7-2.9 | Zhigalin & Moskvitina (2017) |

Independently confirmed via search: Zhigalin & Moskvitina sampled 144
individuals across 2 urban and 2 suburban colonies (Tomsk, Russia) and
found significantly larger litters in urban colonies (2.7-2.9 pups/female)
than suburban ones (1.8 pups/female) -- both the values and the
urban/suburban contrast match.

**Important limitation, flagged by Paulo and not resolved here:** Safi's
survival rates were apparently computed as simple return rates (animals
marked in year *t* recaptured in *t+1*), not from a modern open-population
mark-recapture model (e.g. Cormack-Jolly-Seber), and no confidence
intervals were reported. Treat 0.62/0.76/0.42/0.87 as central/proxy
empirical estimates, not precisely known rates -- hence the sensitivity
ranges below.

This female-only Leslie/stage matrix uses adult **female** survival
(0.76) and the juvenile rate (0.62, assumed equal by sex per the source).
Female offspring per breeding female per year, assuming a 1:1 sex ratio:

- Suburban: `F = 0.87 x (1.8 x 0.5) = 0.783`
- Urban: `F = 0.87 x (2.9 x 0.5) = 1.262` (upper bound)

Sensitivity scenarios (not published species estimates) explored
alongside the point values:

- `S_juvenile = 0.45-0.70`
- `S_adult = 0.65-0.85`

## Secondary cross-check: Myotis lucifugus proxy

Before the species-specific data above were located, juvenile survival
was proxied from two *Myotis lucifugus* studies (also a well-studied
temperate vespertilionid), in the same spirit as Frick et al.'s own use
of *Pipistrellus pipistrellus* as a cross-species density proxy:

| Source | Adult female survival | First-year/juvenile survival |
|---|---|---|
| Frick, W.F. et al. (2010). Influence of climate and reproductive timing on demography of little brown myotis *Myotis lucifugus*. *Journal of Animal Ecology*. | 0.63-0.90 | 0.23-0.46 |
| Schorr, R.A. et al. (2021). Population dynamics of little brown bats (*Myotis lucifugus*) at summer roosts. *Ecology and Evolution*. | (not extracted) | 0.45 (SE 0.06) and 0.71 (SE 0.09) at two separate roosts |

This gives a juvenile survival range of roughly 0.23-0.71 -- consistent
with, and slightly wider than, the Safi-based 0.45-0.70 sensitivity range
above -- and remains a useful independent cross-check given the
uncertainty around the *V. murinus*-specific figures.

## Matrix structure and the age-at-first-breeding question

Pre-breeding-census, female-only Leslie matrix, annual time step, ages
`0, 1, ..., max_age`. Survival age 0->1 is the juvenile rate; survival
age >= 1 is the adult rate; fecundity at age `a` is 0 below `alpha` and
`F` (female fecundity, see above) from `alpha` onward, pro-rated for the
single age class a non-integer `alpha` falls within. `max_age` = 15 is
high enough that the dominant eigenvalue is insensitive to further
truncation.

The Safi-based data structure (juveniles transitioning directly to a
breeding-capable adult class, with 87% of adult females reproducing) is
closer to an **early-maturing** structure (age at first breeding ~1 year)
than to the alpha = 1.5-3.5 range established in
`alpha_first_breeding_evidence.md` from the general bat life-history
literature and Frick et al.'s own discussion. Three alpha structures are
therefore compared explicitly in `R/leslie_vespertilio_murinus.R`:
**early** (alpha = 1), **intermediate** (alpha = 2), **delayed**
(alpha = 3), alongside the continuous 1.5-3.5 range used elsewhere in
this analysis -- to see how much this structural choice actually matters
for the resulting lambda, rather than assuming either range is correct.

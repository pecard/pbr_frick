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

## Family/order-level envelope (Paulo, 2026-10): widening the plausible
## ranges beyond Safi alone

Paulo's point (2026-10, re-read R/carrying_capacity_reference.R Section 5):
the elasticity/breakeven analysis there checks each required vital-rate
value against a plausible range, but that range (`s_range`,
`leslie_boundary_range_*` in `inputs/pbrSettings_BSH_DGY.R`) was built
mostly as ad hoc plausibility bounds, not from an actual literature
envelope beyond Safi (2006) and the single *Myotis lucifugus* cross-check
above. With full scientific reasonableness we can widen that envelope
using congeneric, confamilial, and order-level vital-rate data -- the same
logic already used for litter size (Zhigalin & Moskvitina) and for the
*M. lucifugus* juvenile-survival cross-check, just extended to all four
vital rates and to more Vespertilionidae genera.

**Same access limitation as above:** `WebFetch` was blocked for every
domain tried this session (pmc.ncbi.nlm.nih.gov, cdnsciencepub.com,
nature.com), so the figures below come from `WebSearch` result snippets
only, not the papers' full text or tables -- treat them as a documented
secondary synthesis, re-verifiable by Paulo against the primary sources
once direct access is available, not as independently confirmed numbers.

| Source | Species / scope | Adult survival | Juvenile/1st-yr survival | Breeding fraction | Litter size |
|---|---|---|---|---|---|
| Safi (2006) -- this note's baseline | *V. murinus* | 0.76 (female) | 0.62 | 0.87 | -- |
| Zhigalin & Moskvitina (2017) | *V. murinus* | -- | -- | -- | 1.8 (suburban) / 2.7-2.9 (urban) |
| Lentini et al. (2015), *Biology Letters*, "A global synthesis of survival estimates for microbats" -- 193 estimates, 44 species, 7 families, ~70% Vespertilionidae | Microbats, order-wide | 0.774, 95% CI [0.617, 0.890] (adult female, summer estimates -- the highest of the categories compared) | lower than adult (qualitative: juveniles < adults, males < females) | -- | -- (qualitative: species with more young/year have lower apparent survival) |
| Sendor & Simon (2003), *J. Animal Ecology* 72:308-320 (Marburg Castle, Germany; 15,839 captures) | *Pipistrellus pipistrellus* | 0.799 +/- 0.051 (spring: 0.892 +/- 0.028) | -- | -- | -- |
| Frick et al. (2010), *J. Animal Ecology* | *Myotis lucifugus* | 0.63-0.90 | 0.23-0.46 | -- | -- |
| Schorr et al. (2021), *Ecology and Evolution* | *Myotis lucifugus* | -- | 0.45 (SE 0.06) / 0.71 (SE 0.09), two roosts | -- | -- |
| Netherlands woodland-community study (via ScienceDirect, age/sex/climate survival comparison) | *Myotis daubentonii* | ~0.80 | ~0.50 | -- | -- |
| O'Shea, Ellison, Neubaum, Neubaum, Reynolds & Bowen (2010), *J. Mammalogy* 91:418-428 (Colorado) | *Eptesicus fuscus* | 0.79, 95% CI [0.77, 0.81] | 0.67, 95% CI [0.61, 0.73] (weaned females) | 0.64, 95% CI [0.53, 0.73] (1-yr-olds); 0.95, 95% CI [0.94, 0.96] (older) | 1.11, 95% CI [1.05, 1.17] (112 pregnant females, radiography) |

**Reading this table against the current plausible ranges:**

- **S_adult:** the order-wide meta-analytic mean (0.774) sits almost
  exactly on Safi's own point estimate (0.76) -- strong independent
  corroboration of the baseline itself. But every confamilial estimate
  found (Pipistrellus 0.799, *M. lucifugus* up to 0.90, *M. daubentonii*
  ~0.80, *E. fuscus* 0.79) and the order-level 95% CI upper bound (0.890)
  cluster well BELOW the current ad hoc upper bound of 0.95 used in
  `s_range`/`leslie_boundary_range_s_adult`. The breakeven value Section 5
  of `R/carrying_capacity_reference.R` found necessary to reach
  lambda=1.20 (S_adult=0.949) sits above every one of these confamilial/
  order-level figures, not just above the ad hoc 0.95 ceiling -- the
  broader envelope makes that route look LESS plausible than the narrower
  Safi-only range already suggested, not more.
- **S_juv:** confamilial range 0.23 (*M. lucifugus*, Frick et al.) to 0.71
  (*E. fuscus*, O'Shea et al.) comfortably contains Safi's 0.62 and is
  consistent with, if slightly wider at the bottom than, the current
  0.40-0.80 bound -- no change indicated.
- **p_breed:** previously had NO independent check at all. *E. fuscus*
  (O'Shea et al. 2010) gives 0.64 (first-time breeders) to 0.95 (older,
  experienced females), bracketing Safi's 0.87 on both sides -- the
  current 0.70-0.98 range is well supported; the lower bound could
  reasonably extend to ~0.64 to cover first-time breeders specifically.
- **Litter size:** *E. fuscus*'s 1.11 (essentially always a singleton,
  O'Shea et al. 2010) is well below V. murinus's own suburban/urban range
  (1.8-2.9). This is a genuine family-level contrast, not noise: most
  Vespertilionidae are single-pup species, and V. murinus's own
  documented fecundity is already on the high side for the family. It
  argues for anchoring litter size to the species-specific
  Zhigalin & Moskvitina figures (as this note already does) rather than
  widening the upper bound further using confamilial data, which would
  pull the range DOWN, not up.

**Net effect on the carrying-capacity elasticity analysis:** this broader
envelope does not open up new routes to lambda=1.20/1.24 that the
Safi-only range had missed -- if anything, for S_adult it narrows the
already-implausible breakeven value further out of range. The one finding
worth carrying forward precisely BECAUSE of this check: V. murinus's own
survival estimate already sits at the order-wide meta-analytic mean, so
Safi's figures are not an outlier needing correction, and the
carrying-capacity H under THAT baseline (not an inflated-survival
scenario) is the number that most deserves a literature-grounded
uncertainty interval, not a point value -- see
`R/carrying_capacity_reference.R` Section 7 for a joint (not one-at-a-time)
Monte Carlo over this combined envelope.

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

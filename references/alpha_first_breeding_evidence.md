# Evidence base for age at first breeding (alpha) — Vespertilio murinus

For PBR, `alpha` (age at first breeding, in the Niel & Lebreton lambda_max
approximation) should represent **effective age at first reproduction**
(the age at which a female actually recruits offspring into the
population), not merely physiological sexual maturity. The two can differ:
an animal can be physiologically capable of mating well before it
successfully raises young to independence.

**Access limitation (2026-09):** as with `references/leslie_matrix_parametrisation.md`,
outbound access in this session could only reach search-engine result
snippets (`WebSearch`); `WebFetch` remained blocked for `ratpenats.org` and
its sister domain `batmonitoring.org` on every attempt, including a
re-check on 2026-09-16. The Ratpenats claim below ("many females mate in
their first autumn") was cross-checked on that date across three
independent `WebSearch` queries, including one scoped to the original
Catalan-language version of the same species page
(ratpenats.org/especies/vespertilio-murinus, vs. the English
ratpenats.org/en/species/vespertilio-murinus/ cited in the table) -- all
three returned the identical sentence verbatim, the strongest
corroboration achievable without direct page access.

## Sources reviewed

| Source | Species / group | Relevant information | Implication for alpha |
|---|---|---|---|
| [AnAge database](https://genomics.senescence.info/species/entry.php?species=Vespertilio_murinus) | *Vespertilio murinus* | Female and male sexual maturity listed at **730 days** (~2 years); litter size 2; one litter per year. | Supports alpha ~= 2 as a defensible value. |
| Haensel 2010, *Nyctalus* | *Vespertilio murinus* | Reports that some females raise young already in their first year of life; the author's own data suggest most or all females may reach sexual maturity in the year of birth and raise young as one-year-olds, though this is stated cautiously. | Suggests alpha = 1-2 may be plausible; alpha = 2 is not an optimistic assumption. |
| [Ratpenats species account](https://www.ratpenats.org/en/species/vespertilio-murinus/) (original Catalan version: [ratpenats.org/especies/vespertilio-murinus](https://www.ratpenats.org/especies/vespertilio-murinus)) | *Vespertilio murinus* | States that many females mate in their first autumn. | Supports early sexual maturity, but mating is not the same as confirmed breeding. |
| Racey & Entwistle synthesis, cited in Frick et al.'s own discussion ([related: Life-history and Reproductive Strategies of Bats](https://www.researchgate.net/publication/279523560_Life-history_and_Reproductive_Strategies_of_Bats)) | Bats broadly | The general view for most bats is that age at first breeding is around 1 year, though first-year breeding probability can be lower; Frick et al. use alpha = 1.5 as a general benchmark and suggest alpha = 2 for some families (e.g. Pteropodidae). | Supports alpha = 1.5-2 as plausible for many bat species. |
| [Cryan et al. 2012](https://tethys.pnnl.gov/sites/default/files/publications/Cryan-et-al-2012.pdf) | Migratory tree bats killed at wind turbines | Found ovarian development indicating sexual maturity even in first-year females for several migratory tree-roosting species. | Useful analogue: some wind-vulnerable migratory bats mature very early. |
| [Komar et al. 2020](https://journals.biologists.com/jeb/article/223/8/jeb214825/223805/Food-restriction-delays-seasonal-sexual-maturation) | Male *Vespertilio murinus* | Food restriction delayed seasonal sexual maturation in males; well-fed males matured earlier. | Shows maturation timing is condition-dependent; male data are less directly useful for PBR (which uses female recruitment). |

## Interpretation

- **alpha = 2** is a defensible, reasonably conservative central value for
  PBR purposes, particularly under the "effective female recruitment"
  reading of alpha rather than physiological readiness. It is directly
  supported by the AnAge sexual-maturity estimate and is not contradicted
  by Haensel (2010).
- **alpha = 1.5** (the paper's own general bat benchmark, also implied by
  Haensel's more cautious lower estimate and by the Cryan et al. analogue
  in other migratory species) is a plausible *lower* bound, but rests on
  thinner, more cautiously stated evidence specific to this species
  (a single study noting "some" or "most/all" first-year breeding, not a
  confirmed rate).
- **alpha = 3-4** has no direct empirical support for *V. murinus* in the
  sources reviewed. It remains useful as a precautionary sensitivity
  scenario representing ecological constraints on *effective* recruitment
  (migration costs, poor juvenile survival, low first-year breeding
  probability) but should not be treated as the best-supported biological
  estimate. Narrowing the upper end from 4 to 3.5 years keeps this
  precautionary scenario while trimming the part of the original 2-4
  range with the least support.

## Consequence for the sensitivity/Monte Carlo range

This evidence base supports widening the alpha range used in
`R/pbr_frick_bsh_dgy_thresholds.R` and `R/pbr_frick_response_surfaces.R`
from **2-4 years** to **1.5-3.5 years**, shifting the whole plausible range
about half a year younger (better reflecting effective recruitment rather
than a purely precautionary maturity assumption) while trimming the least
-supported part of the old upper tail. See
`R/pbr_frick_bsh_dgy_thresholds.R` for the resulting comparison between the
two ranges.

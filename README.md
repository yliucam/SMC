D&C-melding
================

[`dc_melding`](./dc_melding) contains the running examples for *Chained
Markov melding using divide and conquer sequential Monte Carlo*
(<https://doi.org/10.48550/arXiv.2605.22301>).

## **Owls example**

[`data`](./owls/data) contains the data used in the owls example.

Use `recap_run.R` to draw particles from the capture-recapture model
$p_1$.

Use `fecundity_run.R` to draw particles from the fecundity model $p_3$.

Then, use `owls_run.R` to update those particles above in $p_2$.

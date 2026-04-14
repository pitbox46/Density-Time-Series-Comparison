# Todo

## Known Issues

Changing the density of `dens_grid` affects the choice of k and subsequent results.
Higher density values at the tails leads to some odd overfitting type problem at the tails
where we get lower Wasserstein distances, but it all comes from the tails.
I suspect this issue is with the Wasserstein distance calculation.
Similarly, the WAR model visually outperforms the others, but doesn't always
come on top in terms of Wasserstein Distance.

LQD consistently exhibits larger tails than all other methods.
This could be due to precision issues similar to above.
UPDATE: This is due to the tail ends not being properly evaluated.
They don't descend to zero.

I am using a KNN bandwidth selection. This was chosen for log-normal data,
but may not be optimal for other data types.

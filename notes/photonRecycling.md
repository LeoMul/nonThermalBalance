# Introducing the problem

For a fraction of element $X$ in $X^{i+}$ denoted by $f_i$, the dynamics of the ionization looks something like:

$$
    \frac{d f_i}{dt} = f_{i-1}\mathcal{I}_{i \to i+1} +f_{i+1}\mathcal{R}_{i+1 \to i} -f_i \mathcal{I}_{i \to i+1} -f_i\mathcal{R}_{i \to i-1},
$$

for $0\leq i \leq Z-1$. Under the condition $\frac{d f_i}{dt} = 0$ and performing cumulative sums (or inverting a tri-diagonal matrix) gives balance of each of the individual channels,

$$
f_i \mathcal{I}_{i \to i+1} = f_{i+1}\mathcal{R}_{i+1 \to i}, 
\implies f_{i+1}   = \frac{\mathcal{I}_{i \to i+1}}{\mathcal{R}_{i+1 \to i}} f_i,
$$

which combined with the constraint $\sum_i f_i = 1$ gives the balance. Suppose now we wanted to account for the impact of recombination photons internally photoionizing the plasma. We model the rate of photoionization from a recombination photon as something like,

$$
P_{ij} f_{j+1} \mathcal{R}_{j+1 \to j},
$$

where the photon emitted from $X^{(j+1)+} \to X^{j+} + h\nu_j$ has probability $P_{ij}$ of resulting in a photoionization of $X^{i+}$. The ionization dynamics look something now like

$$
    \frac{d f_i}{dt} =   f_{i-1}\mathcal{I}_{i \to i+1} + f_{i+1} (1-P_{i,i})\mathcal{R}_{i+1 \to i} + \sum_{j>{i-1}}P_{i-1,j}\mathcal{R}_{j+1 \to j}f_{j+1} \\ -  f_i \mathcal{I}_{i \to i+1} -  f_i (1-P_{i-1,i-1})\mathcal{R}_{i \to i-1} - \sum_{j>i}P_{ij} \mathcal{R}_{j+1 \to j}f_{j+1}   ,
$$

where the gain and loss terms from adjacent recombinations must now be modulated by one minus the probability of self-photoionization. There is an additionally an extra gain and loss term from the photoionization by photons non-adjacent recombinations. 

Again, we can take the partial sums in the steady state,

$$
\sum_{k=0}^{i} \frac{d f_k}{dt} = 0,
$$

to arrive at the new set of equations

$$
f_i \mathcal{I}_{i \to i+1} = f_{i+1} (1-P_{ii})\mathcal{R}_{i+1 \to i} - \sum_{j>i}P_{ij} \mathcal{R}_{j+1 \to j}f_{j+1},
$$

which can be implicitly inverted as 

$$
f_{i+1} = \frac{f_i \mathcal{I}_{i \to i+1} + \sum_{j>i}P_{ij} \mathcal{R}_{j+1 \to j}f_{j+1}}{(1-P_{ii})\mathcal{R}_{i+1 \to i}}.
$$

Which can be solved iteratively. For example, first solve the rate equations with $P_{ij}=0$ to obtain some initial ionization balance $f^{(0)}_i$ and calculate the probabilities $P^{(0)}_{ij}$. Then at iteration $k$, 

$$
f^{(k+1)}_{i+1} = \frac{f^{(k+1)}_i \mathcal{I}_{i \to i+1} + \sum_{j>i}P^{(k)}_{ij} \mathcal{R}_{j+1 \to j}f^{(k)}_{j+1}}{(1-P^{(k)}_{ii})\mathcal{R}_{i+1 \to i}}.
$$

What remains is to calculate the probabilities $P_{ij}$.

# How Axelrod handles it

From Axelrod's Thesis, the probability of the recombination $X^{(j+1)+} \to X^{j+} + h\nu_j$ resulting in a photoionization of $X^{i+}$ is

$$
P_{ij} = \frac{f_i \sigma_{i}(h\nu_j)}{\sum_{k=j}^Z  f_k\sigma_k(h\nu_j)}  (1 - e^{-\tau_j}) \times \phi_R
$$

with an 'optical depth' of the recombination transition to $j$

$$
\tau_j = nR\sum_{k=j}^Z  \sigma_k(h\nu_j)
$$

where (it seems like) Axelrod is assuming that all recombination occurs through the ground state (probably mostly or approximately true) - so there is no need to integrate over the spectrum of recombination photons? I think he introduces this coefficient $\phi_R$ as a way to get around that. Note that the above sum is probably a typo - and is probably supposed to go from 0 to $j$, as one notices that $\sigma_k(h\nu_j) =0$ for $k>j$ and the above summations would be empty.  Probably Axelrod had everything summed from 0 to $z$ and eliminated the wrong part of the summation, and the correct formulae should be:

$$
P_{ij} = \frac{f_i \sigma_{i}(h\nu_j)}{\sum_{k=0}^j  f_k\sigma_k(h\nu_j)}  (1 - e^{-\tau_j}) \times \phi_R
$$

with an 'optical depth' of the recombination transition to $j$

$$
\tau_j = nR\sum_{k=0}^j  \sigma_k(h\nu_j)
$$

which gives rough agreement with what is in his thesis around the discussion of Eq. 4.30. 


Note that in the optically thick limit of $\tau_j >> 1$, we have 

$$
\sum_{i=0}^z P_{ij} = \phi_R \sum_i \frac{f_i \sigma_{i}(h\nu_j)}{\sum_{k=0}^j  f_k\sigma_k(h\nu_j)}(1 - e^{-\tau_j}) \phi_R  = (1 - e^{-\tau_j}) \phi_R\sum_{i=0}^z \frac{f_i \sigma_{i}(h\nu_j)}{\sum_{k=0}^z  f_k\sigma_k(h\nu_j)}=(1 - e^{-\tau_j}) \phi_R \to_{\tau_j\to\infty} \phi_R,
$$

which is to say that the recombination event $X^{(j+1)+} \to X^{j+} + h\nu_j$ has probability $(1 - e^{-\tau_j})\phi_R$ of photoionizing _something_.

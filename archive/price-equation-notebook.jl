### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ a46ae414-7e4b-11ee-1937-35ac76f5d627
md"""

## price-equation-notebook.jl

by Rohan Maddamsetti.

This is a derivation of the Price equation and consequences, following derivation in Chapter 6 of "Evolutionary Theory" by Sean Rice

"""

# ╔═╡ 71709273-312d-462e-8883-a7397741e8f1
md"""

The following more or less comes straight from Chapter 6 of Rice's book.  

### Price's theorem.  
\
Consider a trait ϕ, which can be any property we can assign a numerical value to.  

We define the following terms: \
$N$ = Population size \
$\phi_i$ = Phenotype of individual $i$ $(0 < i \leq N)$ \
$\overline{\phi}$ = Mean phenotype in the population \
$\delta_{i,j}$ = Difference between phenotype of the jth descendant of individual i and i's phenotype \
$\overline{\delta_i}$ = Difference between the mean value of $\phi$ among i's descendants and $\phi_i$ \
$W_i$ = Number of descendants of individual $i$ \
$\overline{W}$ = Mean number of descendants per individual.
"""

# ╔═╡ 9a96014d-346d-430c-bfc2-ed0af7727ea9
md"""
Using this notation, the phenotypic value of descendant j of individual i in the current generation is $\phi_i + \delta_{i,j}$. \
The mean phenotype of the descendants, $\overline{\phi'}$, is then: \
\
$\overline{\phi'} = \frac{\sum_{i=1}^N \sum_{j=1}^{W_{i}} (\phi_i + \delta_{i,j})}{\sum_{i=1}^N W_i}$. \
\
The numerator adds up the phenotypic values of each of the descendants by looking in term at each of the N individuals in the current generation, and then summing over the phenotypic values of each of their $W_i$ descendants. The denominator is simply the total number of descendants.
"""

# ╔═╡ 3e6201d5-51fa-4929-8b87-826974c31e81
md"""
We then use the following identities, which follow from the definition of the arithmetic mean:  \
\
$\sum_{j=1}^{W_{i}} \phi_i = W_{i}\phi_{i}$ \
$\sum_{j=1}^{W_{i}} \delta_{i,j} = W_{i}\overline{\delta_{i}}$ \
$\sum_{j=1}^N W_{i} = N\overline{W_{i}}$ \
"""

# ╔═╡ cd46948d-045a-43c5-b667-a11d655ef21b
md"""
These allow us to rewrite the equation as: \
$\overline{\phi'} = \frac{1}{N\overline{W}}[\sum_{i=1}^N W_{i}\phi_{i} + \sum_{i=1}^N W_{i}\overline{\delta_{i}}]$, so, \
$\overline{\phi'} = \frac{1}{\overline{W}}[E(W\phi) + E(W\overline{\delta})]$
"""

# ╔═╡ 46b60f17-8d9b-4207-80fe-0f91de85b0c0
md"""
Now, using the fact that $Cov(x,y) = E(xy) - \overline{x}\cdot\overline{y}$, \
we can substitute for $E(W\phi)$ to get: \
$\overline{\phi'} = \frac{1}{\overline{W}}[Cov(W_, \phi) + \overline{W}\cdot\overline{\phi} + E(W\overline{\delta})]$, so, \

$\overline{\phi'} = \frac{1}{\overline{W}}[Cov(W_, \phi) + E(W\overline{\delta})] + \overline{\phi}$. \

"""

# ╔═╡ 8be9d012-0bad-4415-b455-e0b1c358a94c
md""" 
Finally, subtracting $\overline{\phi}$ from both sides yields: \
\
$\Delta\overline{\phi} = \frac{1}{\overline{W}} [Cov(W, \phi) + E(W\overline{\delta})]$ \
\
This is *Price's theorem* (Price 1970).
"""

# ╔═╡ ba95d92e-3fd5-4d47-a286-80fe94d93e97
md""" 

**Remarks:** \
The $\frac{1}{\overline{W}}Cov(W, \phi)$ term represents the change due to differential survival and reproduction, encompassing selection and genetic drift. \
\
The $\frac{1}{\overline{W}}E(W\overline{\delta})$ term represents the change due to processes involved in reproduction, such as recombination, regression toward the mean phenotype, or selection at a lower level of organization (i.e. plasmids or other genetic elements that bias their own transmission into daughter cells, potentially at the expense of other, competing genetic elements, as in meiotic drive or CRISPR drives) \
\
Note that although we used absolute fitness values in our derivation of Price's theorem, we could just as well use relative fitness $w$, since any constant that multiplies all fitness values will appear in the numerator and denominator, and will thus cancel out.

"""

# ╔═╡ 808e94f1-f946-4c02-8bb7-700d698cfc52
md"""
### Consequences of Price's theorem for incompatible plasmid systems.

Suppose $\overline\delta = 0$. This means that offspring do not systematically deviate from their parents in fitness, which implies no transmission bias in the plasmids, and implies that the environment is not changing (e.g. the antibiotic concentration has not changed).

Then, $E(W\overline{\delta}) = 0$, so \
\
$\Delta\overline{\phi} = \frac{1}{\overline{W}}Cov(W, \phi)$, so \
\
$\Delta\overline{\phi} = \frac{1}{\overline{W}} \beta_{w,\phi} \cdot Var(\phi)$,
\
\
as the regression coefficient $\beta_{w,\phi}$ is defined as $\frac{Cov(W, \phi)}{Var(\phi)}$

"""

# ╔═╡ 0bcf8cc7-7d78-40de-a71f-f8157a62f10f
md"""
Let $\phi = W$. Then, \
\
$\Delta\overline{W} = \frac{1}{\overline{W}}Cov(W, W) = \frac{1}{\overline{W}}Var(W)$.
\
This result is known as *Fisher's fundamental theorem of natural selection*.
\
\
**More importantly:**
\
Let $\phi = TCN$, where $TCN$ is the transposon copy number. Then, \
\
$\Delta\overline{TCN} = \frac{1}{\overline{W}} \beta_{W,TCN} \cdot Var(TCN)$.
\

This result gives the change in TCN over one generation, and shows that it depends on the slope of the linear regression of fitness on TCN, and the variance of TCN in the population. This has a nice geometric interpretation (Figure 1).
"""

# ╔═╡ Cell order:
# ╠═a46ae414-7e4b-11ee-1937-35ac76f5d627
# ╠═71709273-312d-462e-8883-a7397741e8f1
# ╠═9a96014d-346d-430c-bfc2-ed0af7727ea9
# ╠═3e6201d5-51fa-4929-8b87-826974c31e81
# ╠═cd46948d-045a-43c5-b667-a11d655ef21b
# ╠═46b60f17-8d9b-4207-80fe-0f91de85b0c0
# ╠═8be9d012-0bad-4415-b455-e0b1c358a94c
# ╠═ba95d92e-3fd5-4d47-a286-80fe94d93e97
# ╠═808e94f1-f946-4c02-8bb7-700d698cfc52
# ╠═0bcf8cc7-7d78-40de-a71f-f8157a62f10f

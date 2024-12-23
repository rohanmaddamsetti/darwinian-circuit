### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 498a2906-9a95-47ef-b3e3-a5ef24ee0a77
begin
	using Plots
	using PlutoUI
	using Statistics
	using LinearAlgebra
	using DifferentialEquations
	using Distributions
end

# ╔═╡ 60f6f4ca-a684-11ee-26bd-6b50dd5567a6
md"""
# darwinian-circuit-quasispecies-model-v0.5.jl 

by Rohan Maddamsetti  

Julia version 1.10.  

## Programmed population control by diversity-maintaining genetic elements

##### Abstract

Extrachromosomal DNAs (ecDNAs) amplify gene expression and generate phenotypic diversity in cancers and bacteria. A diverse population of ecDNA variants can be maintained within a single cell; however, the implications of this observation are understudied. Here, we use bacterial ecDNAs (plasmids) with a special property— balancing selection maintains coexisting variants within a single cell— to program population-level gene expression dynamics in response to pulses of antibiotic. We show that population-level gene expression dynamics are determined by a simple quantitative rule grounded in evolutionary theory.\
\
When plasmid copy number is increased, stochastic segregation during cell division generates more diverse daughter cells, which increases the speed at which selection can tune population-level gene expression in response to environmental change.\
\
When ecDNA copy number is sufficiently high to maintain diverse ecDNAs in single cells, tunable evolutionary dynamics emerge in otherwise clonal populations. Our work demonstrates a universal quantitative principle governing the speed at which ecDNAs allow populations to adapt to fluctuating environments. 


"""

# ╔═╡ fea3c300-7aaf-4da8-bed7-6e9426b4ac3d
md"""

## TODO: 

### A) double check the key logical steps in my tentative argument.

###### 1) does increasing PCN increase variance in fitness? Yes, but extremely small effect.
###### 2) does increasing PCN increase tetA-copy-number-fitness covariance? Yes, but extremely small effect.

###### 3) does increasing variance in fitness speed adaptation?
###### 4) does increasing tetA-fitness covariance speed adaptation?

### B) fix scaling factor bugs. In fact, scaling makes the fit worse??

Note that the fit is really good when covariance is positive, but not good when covariance is negative... this may be a clue as to the bug.

## DEBUGGING: 
1) in the constant [Tet] model, fitness variance divided by mean fitness does not exactly equal rate of change of mean fitness.

2) in the pulse [Tet] model, tetA-fitness covariance divided by mean fitness does not exactly equal rate of change of mean fitness. 

"""

# ╔═╡ 772351cf-cddd-45d8-9f34-f3ef6ccb2a66
md"""
## **Model description**

I built a simple evolutionary model to examine how plasmid copy number affects tetA-GFP copy number dynamics in the tetA-transposon system.


##### Quasispecies model formulation.

We use the quasispecies model (see "Unifying Evolutionary Dynamics" by Karen Page and Martin Nowak in Journal of Theoretical Biology, 2002). This model has a physical interpretation as a model of continuous culture in a idealized turbidostat, such that the dilution rate is equal to the bulk growth rate (i.e. mean population fitness). The idea here is that the dilution rate increases as the population adapts, such that the total population density is conserved in the model. Physically, this corresponds to a turbidostat in which the dilution rate instantaneously increases as the population increases in mean growth rate through evolutionary adaptation, such that the (optical) density of the culture stays constant.

**Model Assumptions**  

The population is modeled as a distribution over tetA copy numbers, ranging from 1 (found on the chromosome) to $n+1$, where $n$ is the maximum plasmid copy number. We therefore represent the population as a vector 

$\mathbf{x}(t) = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ ... \\ x_{n+1} \end{pmatrix}$  

The sum of the $n+1$ entries of $\mathbf{x}(t)$ gives the total population size $N(t) = \sum_{i=1}^{n+1} \mathbf{x}_i(t)$.  

"""

# ╔═╡ d3b956b4-885a-447c-bd00-de5fa18fdcad
md"""
*Growth dynamics*  

For a given tetracycline concentration, each tetA copy number class has a growth rate, or *fitness* $r_i$.  Let $g_i$ be be the log-fitness of tetA copy number class $i$, such that $g_i = log(r_i)$. We assume that there is some optimal tetA copy number given some tetracycline concentration. We assume that $r_i$ is non-negative, with a single peak at the optimum tetA copy number for a given $[Tet]$ concentration, so a natural choice is to define the following quadratic log-fitness function:  

$g_i = g_{max} - \frac{(i - \alpha*[Tet])^2}{2\sigma^2}$  

Here, $g_{max}$ is the maximum log-growth-rate, $i$ is the number of tetA copies in this strain, [Tet] is the antibiotic concentration, $\sigma$ is a free parameter that determines the width of the quadratic function, and $\alpha$ is a unit-conversion factor with units of "tetA gene copies/[Tet]". Without loss of generality we assume $\alpha = 1$, so:

$g_i = g_{max} - \frac{(i - [Tet])^2}{2\sigma^2}$ 

We assume fitness is non-negative, so we pass this $g_i$ through an exponential to get a Gaussian fitness function:

$r_i = e^{g_i}$ 

Then, the growth rates (fitnesses) for each subpopulation is a vector:

$\mathbf{r} = \begin{pmatrix} r_1 \\ r_2 \\ r_3 \\ ... \\ r_{n+1} \end{pmatrix}$  


We then define the diagonal matrix $\mathbf{G} = \begin{bmatrix}
r_1 & 0 & 0 & ... & 0 \\
0 & r_2 & 0 & ... & 0 \\
0 & 0 & r_3 & ... & 0 \\
... & ... & ... & ... & ... \\
0 & 0 & 0 & ... & r_{n+1} \\
\end{bmatrix}$.  

This represents how each subpopulation of $\mathbf{x}$ grows based on $\mathbf{r}$.

And the average population growth rate (mean population fitness) of the whole population is $\overline{r}(t) = \mathbf{r} \cdot \mathbf{x}(t)$.

"""

# ╔═╡ 624831f0-b669-44e9-ae78-ab0e09716ea4
md"""
*Phenotypic switching dynamics*

Subpopulations switch phenotypes based on plasmid segregration during cell division.

We assume that the total plasmid copy number is fixed per cell. However, the cells of each subpopulation can gain or lose plasmids containing the tetA transposon, based on sampling the plasmids of the parental cell with replacement. This assumption leads to a stochastic switching matrix based on binomial probabilities.  

Since the first index represents a state in which _none_ of the plasmids have the transposon,  we use a change of variable from zero-based indexing to one-based indexing, so that we can use the binomial formula (zero-based indexing) and then use a change-of-variable to use one-based indexing for the matrix.  

We define the entries of the stochastic transition matrix $\mathbf{S}$ with zero-based indexing:  

$\mathbf{S}_{ab} = \binom{n}{a} \left( \frac{b}{n} \right) ^a \left( \frac{n-b}{n} \right) ^{n-a}$  

where $a$ is the zero-based row-index (representing an offspring with $a$ copies of the transposon) and $b$ is the zero-based column-index (representing a parent with $b$ copies of the transposon). Recall that $n$ is the plasmid copy number, so $0 \leq a \leq n$ and $0 \leq b \leq n$.

We then use the following change-of-variables to define $\mathbf{S}$ with one-based indexing. The matrix is exactly the same-- only the indexing to refer to entries is different:

$i = a + 1$  
$j = b + 1$


(6-dimensional case shown):  

$\mathbf{S} = \begin{bmatrix}
1 & \binom{5}{0} \left( \frac{1}{5} \right) ^0 \left( \frac{4}{5} \right) ^5 & \binom{5}{0} \left( \frac{2}{5} \right) ^0 \left( \frac{3}{5} \right) ^5 & \binom{5}{0} \left( \frac{3}{5} \right) ^0 \left( \frac{2}{5} \right) ^5 & \binom{5}{0} \left( \frac{4}{5} \right) ^0 \left( \frac{1}{5} \right) ^5 & 0 \\
0 & \binom{5}{1} \left( \frac{1}{5} \right) ^1 \left( \frac{4}{5} \right) ^4 & \binom{5}{1} \left( \frac{2}{5} \right) ^1 \left( \frac{3}{5} \right) ^4 & \binom{5}{1} \left( \frac{3}{5} \right) ^1 \left( \frac{2}{5} \right) ^4 & \binom{5}{1} \left( \frac{4}{5} \right) ^1 \left( \frac{1}{5} \right) ^4 & 0 \\
0 & \binom{5}{2} \left( \frac{1}{5} \right) ^2 \left( \frac{4}{5} \right) ^3 & \binom{5}{2} \left( \frac{2}{5} \right) ^2 \left( \frac{3}{5} \right) ^3 & \binom{5}{2} \left( \frac{3}{5} \right) ^2 \left( \frac{2}{5} \right) ^3 & \binom{5}{2} \left( \frac{4}{5} \right) ^2 \left( \frac{1}{5} \right) ^3 & 0 \\
0 & \binom{5}{3} \left( \frac{1}{5} \right) ^3 \left( \frac{4}{5} \right) ^2 & \binom{5}{3} \left( \frac{2}{5} \right) ^3 \left( \frac{3}{5} \right) ^2 & \binom{5}{3} \left( \frac{3}{5} \right) ^3 \left( \frac{2}{5} \right) ^2 & \binom{5}{3} \left( \frac{4}{5} \right) ^3 \left( \frac{1}{5} \right) ^2 & 0 \\
0 & \binom{5}{4} \left( \frac{1}{5} \right) ^4 \left( \frac{4}{5} \right) ^1 & \binom{5}{4} \left( \frac{2}{5} \right) ^4 \left( \frac{3}{5} \right) ^1 & \binom{5}{4} \left( \frac{3}{5} \right) ^4 \left( \frac{2}{5} \right) ^1 & \binom{5}{4} \left( \frac{4}{5} \right) ^4 \left( \frac{1}{5} \right) ^1 & 0 \\
0 & \binom{5}{5} \left( \frac{1}{5} \right) ^5 \left( \frac{4}{5} \right) ^0 & \binom{5}{5} \left( \frac{2}{5} \right) ^5 \left( \frac{3}{5} \right) ^0 & \binom{5}{5} \left( \frac{3}{5} \right) ^5 \left( \frac{2}{5} \right) ^0 & \binom{5}{5} \left( \frac{4}{5} \right) ^5 \left( \frac{1}{5} \right) ^0 & 1 \\
\end{bmatrix}$  

"""

# ╔═╡ 37d2e248-1f6d-4c39-8cc8-8080a853cb43
md"""
Mathematically, the subpopulation in which all plasmids have the transposon is an absorbing state of the Markov chain.  

However, the subpopulation in which _none_ of the plasmids have the transposon is *not* an absorbing state, because transposons can jump from the chromosome to a plasmid. To model this, we assume a transposition rate $\eta$, and assume that $\eta$ is small compared to the binomial transition probabilities so that it can be ignored in the other entries of the matrix.

$\mathbf{S} = \begin{bmatrix}
1- \eta & \binom{5}{0} \left( \frac{1}{5} \right) ^0 \left( \frac{4}{5} \right) ^5 & \binom{5}{0} \left( \frac{2}{5} \right) ^0 \left( \frac{3}{5} \right) ^5 & \binom{5}{0} \left( \frac{3}{5} \right) ^0 \left( \frac{2}{5} \right) ^5 & \binom{5}{0} \left( \frac{4}{5} \right) ^0 \left( \frac{1}{5} \right) ^5 & 0 \\
\eta & \binom{5}{1} \left( \frac{1}{5} \right) ^1 \left( \frac{4}{5} \right) ^4 & \binom{5}{1} \left( \frac{2}{5} \right) ^1 \left( \frac{3}{5} \right) ^4 & \binom{5}{1} \left( \frac{3}{5} \right) ^1 \left( \frac{2}{5} \right) ^4 & \binom{5}{1} \left( \frac{4}{5} \right) ^1 \left( \frac{1}{5} \right) ^4 & 0 \\
0 & \binom{5}{2} \left( \frac{1}{5} \right) ^2 \left( \frac{4}{5} \right) ^3 & \binom{5}{2} \left( \frac{2}{5} \right) ^2 \left( \frac{3}{5} \right) ^3 & \binom{5}{2} \left( \frac{3}{5} \right) ^2 \left( \frac{2}{5} \right) ^3 & \binom{5}{2} \left( \frac{4}{5} \right) ^2 \left( \frac{1}{5} \right) ^3 & 0 \\
0 & \binom{5}{3} \left( \frac{1}{5} \right) ^3 \left( \frac{4}{5} \right) ^2 & \binom{5}{3} \left( \frac{2}{5} \right) ^3 \left( \frac{3}{5} \right) ^2 & \binom{5}{3} \left( \frac{3}{5} \right) ^3 \left( \frac{2}{5} \right) ^2 & \binom{5}{3} \left( \frac{4}{5} \right) ^3 \left( \frac{1}{5} \right) ^2 & 0 \\
0 & \binom{5}{4} \left( \frac{1}{5} \right) ^4 \left( \frac{4}{5} \right) ^1 & \binom{5}{4} \left( \frac{2}{5} \right) ^4 \left( \frac{3}{5} \right) ^1 & \binom{5}{4} \left( \frac{3}{5} \right) ^4 \left( \frac{2}{5} \right) ^1 & \binom{5}{4} \left( \frac{4}{5} \right) ^4 \left( \frac{1}{5} \right) ^1 & 0 \\
0 & \binom{5}{5} \left( \frac{1}{5} \right) ^5 \left( \frac{4}{5} \right) ^0 & \binom{5}{5} \left( \frac{2}{5} \right) ^5 \left( \frac{3}{5} \right) ^0 & \binom{5}{5} \left( \frac{3}{5} \right) ^5 \left( \frac{2}{5} \right) ^0 & \binom{5}{5} \left( \frac{4}{5} \right) ^5 \left( \frac{1}{5} \right) ^0 & 1 \\
\end{bmatrix}$  

"""

# ╔═╡ c270d434-d426-47f7-aec2-1423f3354fa6
md"""
*Dilution dynamics*

We assume that cells are diluted out at a rate equal to the bulk population growth rate, which is the average population growth rate (mean population fitness) of the whole population $\overline{r}(t) = \mathbf{r} \cdot \mathbf{x}(t)$.

*Full dynamics*  

The growth and stochastic switching dynamics are combined into a matrix $\mathbf{A} = \mathbf{S}\mathbf{G}$, where $\mathbf{S}$ is a stochastic matrix.

Including dilution, the full dynamics are modeled by the following matrix system of ODEs:  

$\frac{d\mathbf{x}}{dt} = \mathbf{A}\mathbf{x}(t) - \overline{r}(t) \mathbf{x}(t)$ Note that $\overline{r}(t)$ is a scalar and not a vector.


This results in this form of the matrix $\mathbf{A}$ (6-dimensional case shown):

$\mathbf{A} = \begin{bmatrix}
r_1(1- \eta) & \binom{5}{0} \left( \frac{1}{5} \right) ^0 \left( \frac{4}{5} \right) ^5 & \binom{5}{0} \left( \frac{2}{5} \right) ^0 \left( \frac{3}{5} \right) ^5 & \binom{5}{0} \left( \frac{3}{5} \right) ^0 \left( \frac{2}{5} \right) ^5 & \binom{5}{0} \left( \frac{4}{5} \right) ^0 \left( \frac{1}{5} \right) ^5 & 0 \\
\eta & r_2 \binom{5}{1} \left( \frac{1}{5} \right) ^1 \left( \frac{4}{5} \right) ^4 & \binom{5}{1} \left( \frac{2}{5} \right) ^1 \left( \frac{3}{5} \right) ^4 & \binom{5}{1} \left( \frac{3}{5} \right) ^1 \left( \frac{2}{5} \right) ^4 & \binom{5}{1} \left( \frac{4}{5} \right) ^1 \left( \frac{1}{5} \right) ^4 & 0 \\
0 & \binom{5}{2} \left( \frac{1}{5} \right) ^2 \left( \frac{4}{5} \right) ^3 & r_3\binom{5}{2} \left( \frac{2}{5} \right) ^2 \left( \frac{3}{5} \right) ^3 & \binom{5}{2} \left( \frac{3}{5} \right) ^2 \left( \frac{2}{5} \right) ^3 & \binom{5}{2} \left( \frac{4}{5} \right) ^2 \left( \frac{1}{5} \right) ^3 & 0 \\
0 & \binom{5}{3} \left( \frac{1}{5} \right) ^3 \left( \frac{4}{5} \right) ^2 & \binom{5}{3} \left( \frac{2}{5} \right) ^3 \left( \frac{3}{5} \right) ^2 & r_4\binom{5}{3} \left( \frac{3}{5} \right) ^3 \left( \frac{2}{5} \right) ^2 & \binom{5}{3} \left( \frac{4}{5} \right) ^3 \left( \frac{1}{5} \right) ^2 & 0 \\
0 & \binom{5}{4} \left( \frac{1}{5} \right) ^4 \left( \frac{4}{5} \right) ^1 & \binom{5}{4} \left( \frac{2}{5} \right) ^4 \left( \frac{3}{5} \right) ^1 & \binom{5}{4} \left( \frac{3}{5} \right) ^4 \left( \frac{2}{5} \right) ^1 & r_5 \binom{5}{4} \left( \frac{4}{5} \right) ^4 \left( \frac{1}{5} \right) ^1 & 0 \\
0 & \binom{5}{5} \left( \frac{1}{5} \right) ^5 \left( \frac{4}{5} \right) ^0 & \binom{5}{5} \left( \frac{2}{5} \right) ^5 \left( \frac{3}{5} \right) ^0 & \binom{5}{5} \left( \frac{3}{5} \right) ^5 \left( \frac{2}{5} \right) ^0 & \binom{5}{5} \left( \frac{4}{5} \right) ^5 \left( \frac{1}{5} \right) ^0 & r_6 \\
\end{bmatrix}$  
"""

# ╔═╡ a7b38e19-1f54-4ce1-9977-79e9968ce84f
md"""
We can generalize this model to the case where the plasmid copy number limits tetA copy number. In this case, higher tetA copy numbers beyond plasmid copy number cannot be reached.

Suppose the maximum plasmid copy number is 4. Then the truncated matrix $\mathbf{A'}$ is (for the 6-dimensional case):   

$\mathbf{A'} = \begin{bmatrix}
r_1(1- \eta) & \binom{5}{0} \left( \frac{1}{5} \right) ^0 \left( \frac{4}{5} \right) ^5 & \binom{5}{0} \left( \frac{2}{5} \right) ^0 \left( \frac{3}{5} \right) ^5 & 0 & 0 & 0 \\
\eta & r_2 \binom{5}{1} \left( \frac{1}{5} \right) ^1 \left( \frac{4}{5} \right) ^4 & \binom{5}{1} \left( \frac{2}{5} \right) ^1 \left( \frac{3}{5} \right) ^4 & 0 & 0 & 0 \\
0 & \binom{5}{2} \left( \frac{1}{5} \right) ^2 \left( \frac{4}{5} \right) ^3 & r_3\binom{5}{2} \left( \frac{2}{5} \right) ^2 \left( \frac{3}{5} \right) ^3 & 0 & 0 & 0 \\
0 & \binom{5}{3} \left( \frac{1}{5} \right) ^3 \left( \frac{4}{5} \right) ^2 & \binom{5}{3} \left( \frac{2}{5} \right) ^3 \left( \frac{3}{5} \right) ^2 & r_4 & 0 & 0 \\
0 & \binom{5}{4} \left( \frac{1}{5} \right) ^4 \left( \frac{4}{5} \right) ^1 & \binom{5}{4} \left( \frac{2}{5} \right) ^4 \left( \frac{3}{5} \right) ^1 & 0 & r_5 & 0 \\
0 & \binom{5}{5} \left( \frac{1}{5} \right) ^5 \left( \frac{4}{5} \right) ^0 & \binom{5}{5} \left( \frac{2}{5} \right) ^5 \left( \frac{3}{5} \right) ^0 & 0 & 0 & r_6 \\
\end{bmatrix}$  

"""

# ╔═╡ 623f1230-77ff-4eeb-865f-01db8f673c01
md"""

### Notes on Price's theorem.  
\
The following more or less comes straight from Chapter 6 of Sean Rice's book "Evolutionary Theory".  
\
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

# ╔═╡ 5d83343b-e40a-4193-8283-dd33e453855a
md"""
Using this notation, the phenotypic value of descendant j of individual i in the current generation is $\phi_i + \delta_{i,j}$. \
The mean phenotype of the descendants, $\overline{\phi'}$, is then: \
\
$\overline{\phi'} = \frac{\sum_{i=1}^N \sum_{j=1}^{W_{i}} (\phi_i + \delta_{i,j})}{\sum_{i=1}^N W_i}$. \
\
The numerator adds up the phenotypic values of each of the descendants by looking in term at each of the N individuals in the current generation, and then summing over the phenotypic values of each of their $W_i$ descendants. The denominator is simply the total number of descendants.
"""

# ╔═╡ f5a21fa6-f7b4-42f5-8ad4-58bb813b98b6
md"""
We then use the following identities, which follow from the definition of the arithmetic mean:  \
\
$\sum_{j=1}^{W_{i}} \phi_i = W_{i}\phi_{i}$ \
$\sum_{j=1}^{W_{i}} \delta_{i,j} = W_{i}\overline{\delta_{i}}$ \
$\sum_{j=1}^N W_{i} = N\overline{W_{i}}$ \
"""

# ╔═╡ 853b7a17-bc10-4e1c-8e73-2a0e80592863
md"""
These allow us to rewrite the equation as: \
$\overline{\phi'} = \frac{1}{N\overline{W}}[\sum_{i=1}^N W_{i}\phi_{i} + \sum_{i=1}^N W_{i}\overline{\delta_{i}}]$, so, \
$\overline{\phi'} = \frac{1}{\overline{W}}[E(W\phi) + E(W\overline{\delta})]$
"""

# ╔═╡ 0dcf8ff4-bff7-448c-be66-92d3555d9efb
md"""
Now, using the fact that $Cov(x,y) = E(xy) - \overline{x}\cdot\overline{y}$, \
we can substitute for $E(W\phi)$ to get: \
$\overline{\phi'} = \frac{1}{\overline{W}}[Cov(W_, \phi) + \overline{W}\cdot\overline{\phi} + E(W\overline{\delta})]$, so, \

$\overline{\phi'} = \frac{1}{\overline{W}}[Cov(W_, \phi) + E(W\overline{\delta})] + \overline{\phi}$. \

"""

# ╔═╡ b1e5c35b-35ae-48c2-babf-2058a91f1ea8
md""" 
Finally, subtracting $\overline{\phi}$ from both sides yields: \
\
$\Delta\overline{\phi} = \frac{1}{\overline{W}} [Cov(W, \phi) + E(W\overline{\delta})]$ \
\
This is *Price's theorem* (Price 1970).
"""

# ╔═╡ 41de1814-a88d-4c4e-901b-a69a0a6ca655
md""" 
**Remarks:** \
The $\frac{1}{\overline{W}}Cov(W, \phi)$ term represents the change due to differential survival and reproduction, encompassing selection and genetic drift. \
\
The $\frac{1}{\overline{W}}E(W\overline{\delta})$ term represents the change due to processes involved in reproduction, such as recombination, regression toward the mean phenotype, or selection at a lower level of organization (i.e. plasmids or other genetic elements that bias their own transmission into daughter cells, potentially at the expense of other, competing genetic elements, as in meiotic drive or CRISPR drives). \
\
Note that although we used absolute fitness values in our derivation of Price's theorem, we could just as well use relative fitness $w$, since any constant that multiplies all fitness values will appear in the numerator and denominator, and will thus cancel out.
"""

# ╔═╡ 22ac1b86-a2b6-4562-adb3-61718a1c9a59
md"""
### Consequences of Price's theorem for diversity-maintaining genetic elements (DGEs).

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

# ╔═╡ 60e1b0f7-b91f-4121-80f8-673ce117dccf
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

# ╔═╡ 6d43ce5f-5e34-4d4e-bcd1-b876c5744be8
begin
	## Define global constants.
	SIGMA = 10.0 ## set the width of the fitness function
	R_MAX = 1.0 ## set the max growth rate.
	TIMESPAN = 300.0
	η₀ = 0.001 ## set the transposition rate η₀.
end

# ╔═╡ ddc8420f-1212-459f-91b1-be62db8032dd
function raw_quadratic_fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	fitness is a function of tetA copies and antibiotic concentration.
	for simplicity, define a quadratic function whose peak is shifted
	left or right based on Tet_conc.
	"""
	fitness = r_max - ((Tet_conc - tetA_copy_number)/2σ)^2
	return fitness
end

# ╔═╡ 9e38c88c-fb6a-4623-b3aa-c7c487ff919f
function gaussian_fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	fitness is a function of tetA copies and antibiotic concentration.
	for simplicity, define a gaussian function whose peak is shifted
	left or right based on Tet_conc.
	"""
	fitness = exp(raw_quadratic_fitness_function(tetA_copy_number, Tet_conc, σ, r_max))
	return fitness
end

# ╔═╡ 6dbe79d8-f3f5-45d7-946a-7454542df74e
function fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	syntactic sugar: many downstream functions depend on this function,
	so simpler to change the *definition* of this function than to update the *name* of this function across many different blocks of code.
	"""
	return gaussian_fitness_function(tetA_copy_number, Tet_conc, σ, r_max)
end

# ╔═╡ ce9769d3-a234-4718-b3b1-2ee6fda2d733
function calc_mean_fitness(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	## get the fitness for each fitness class (defined by tetA copy number).
	raw_fitness_vec = fitness_function.(tetA_classes, Tet_conc)
	## calculate mean fitness in the population,
	## by weighting fitness for each class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_fitness = sum(raw_fitness_vec .* frequency_vec)
	return mean_fitness
end

# ╔═╡ 29daeb98-190b-4972-9a7b-1407ed3c4b1a
function calc_relative_fitness(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	fitness_classes = collect(1:length(pop_vec))
	## get the fitness for each fitness class (defined by tetA copy number).
	raw_fitness_vec = fitness_function.(fitness_classes, Tet_conc)
	## calculate mean fitness in the population,
	## by weighting fitness for each class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_fitness = sum(raw_fitness_vec .* frequency_vec)
	## now calculate relative fitness, based on mean fitness.
	relative_fitness_vec = raw_fitness_vec/mean_fitness
	return relative_fitness_vec
end

# ╔═╡ 4c7767df-bcc4-4fb7-9d66-c807ac8ce6e1
function calc_fitness_variance(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	## get the fitness for each fitness class (defined by tetA copy number).
	raw_fitness_vec = fitness_function.(tetA_classes, Tet_conc)
	## calculate mean fitness in the population,
	## by weighting fitness for each class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_fitness = sum(raw_fitness_vec .* frequency_vec)
	raw_squared_deviation_vec = [(x - mean_fitness)^2 for x in raw_fitness_vec]
	fitness_variance = sum(raw_squared_deviation_vec .* frequency_vec)
	return fitness_variance
end

# ╔═╡ 6a470e0e-de4a-4398-893a-d2edbdf44b9f
function calc_mean_tetA_copy_number(pop_vec)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	## calculate mean tetA copy number in the population,
	## by weighting each tetA copy number class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_tetA_copy_number = sum(tetA_classes .* frequency_vec)
	return mean_tetA_copy_number
end

# ╔═╡ c04868af-914e-406d-a82c-6e75b189600f
function calc_tetA_copy_number_variance(pop_vec)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	
	## calculate mean copy number in the population,
	## by weighting each tetA copy number class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)

	mean_copy_number = sum(frequency_vec .* tetA_classes)
	copy_num_squared_deviation_vec = [(x - mean_copy_number)^2 for x in tetA_classes]

	copy_num_variance = sum(copy_num_squared_deviation_vec .* frequency_vec)	
	return copy_num_variance
end

# ╔═╡ 3b3e97a8-ccff-4ab6-94d6-e6d2d7c92323
function calc_tetA_copy_number_fitness_covariance(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	## get the fitness for each fitness class (defined by tetA copy number).
	raw_fitness_vec = fitness_function.(tetA_classes, Tet_conc)
	## calculate mean fitness in the population,
	## by weighting fitness for each class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_fitness = sum(raw_fitness_vec .* frequency_vec)
	
	## calculate covariance by multiplying fitness deviations
	## with tetA copy number deviations.
	mean_tetA_copy_number = sum(tetA_classes .* frequency_vec)
	tetA_deviation_vec = [(x - mean_tetA_copy_number) for x in tetA_classes]
	raw_fitness_deviation_vec = [(x - mean_fitness) for x in raw_fitness_vec]
	tetA_fitness_covariance = sum(tetA_deviation_vec .* raw_fitness_deviation_vec .* frequency_vec)
	return tetA_fitness_covariance
end

# ╔═╡ 492dbbd5-28d3-45da-99ca-b5f7d5e3c31e
function ZeroBasedHypergeometricSwitchingEntry(a, b, plasmid_copy_num)
	""" Input parameters:
	    a: number of tetA transposons in offspring.
	    b: number of tetA transposons in parent.
	    plasmid_copy_num: total plasmid copy number. """
	## Hypergeometric(s, f, n)  
	## Hypergeometric distribution for a population with s successes and f failures, and a sequence of n trials.	
	## We add a factor of 2 to encode the assumption that the plasmid copy number
	## doubles before cell division-- otherwise copy number gets trapped at 1 copy on the plasmid (since sampling without replacement).
	return pdf(Hypergeometric(2b, 2(plasmid_copy_num - b), plasmid_copy_num), a)
end

# ╔═╡ 8ecaeeb3-2a93-48f5-ba21-49899d20c538
function OneBasedHypergeometricSwitchingEntry(i, j, plasmid_copy_num)
	""" Input parameters:
	    i: row-index of hypergeometric switching matrix: 1 + number of tetA transposons in offspring.
	    j: column-index of hypergeometric switching matrix: 1 + number of tetA transposons in parent.
	    plasmid_copy_num: plasmid copy number. 

	Use ZeroBasedHypergeometricSwitchingEntry() with this change-of-variables:
	i = a + 1 ## a is number of tetA transposons in offspring.
	j = b + 1 ## b is number of tetA transposons in parent.
	
	"""
	return ZeroBasedHypergeometricSwitchingEntry(i-1,j-1,plasmid_copy_num)
end

# ╔═╡ 7715798e-aaaa-46a9-bd16-14b866cf39f6
function ZeroBasedBinomialSwitchingEntry(a, b, plasmid_copy_num)
	""" Input parameters:
	    a: number of tetA transposons in offspring.
	    b: number of tetA transposons in parent.
	    plasmid_copy_num: total plasmid copy number. """
	p = b/plasmid_copy_num ## probability of sampling a plasmid with the tetA transposon.
	## 'a' is the number of times that the tetA-transposon plasmid picked from parent into offspring.
	return pdf(Binomial(plasmid_copy_num, p), a)
end

# ╔═╡ 3e8936df-d80e-4b66-b035-a1007c8f5dca
function OneBasedBinomialSwitchingEntry(i, j, plasmid_copy_num)
	""" Input parameters:
	    i: row-index of binomial switching matrix: 1 + number of tetA transposons in offspring.
	    j: column-index of binomial switching matrix: 1 + number of tetA transposons in parent.
	    plasmid_copy_num: plasmid copy number. 

	Use ZeroBasedBinomialSwitchingEntry() with this change-of-variables:
	i = a + 1 ## a is number of tetA transposons in offspring.
	j = b + 1 ## b is number of tetA transposons in parent.
	
	"""
	return ZeroBasedBinomialSwitchingEntry(i-1,j-1,plasmid_copy_num)
end

# ╔═╡ 5e33d199-29f5-493e-9863-f02b61711e26
md""" ##### sometimes the PCN and TET_CONC variables do not update right when the slider is changed-- we have to double check the actual values of the PCN and TET_CONC variables before each slider.."""

# ╔═╡ a2d049bb-4d27-469d-8626-d0c0fedbced3
PCNSlider = @bind PCN Slider(1:100, default=30, show_value=true)

# ╔═╡ b5fcc845-0d6b-41ab-b27f-05f67ff05325
function SwitchingBinomialMatrix(plasmid_copy_num=PCN, η=η₀)
	""" We set a final absorbing state in the Markov chain,
	based on maximum plasmid copy number (like 3 for SC101, 2000 for pUC)."""
	
	## max transposon copy number is plasmid copy number + one chromosomal copy.
	max_TCN = plasmid_copy_num + 1

	## initialize as an m x m identity matrix with ones on the diagonal.
	switching_matrix = Matrix{BigFloat}(I, max_TCN, max_TCN)

	## update binomial entries (overwrite ones and zeros as needed).
	for i in 1:max_TCN
		for j in 1:max_TCN
			if (i <= max_TCN) && (j <= max_TCN)
				switching_matrix[i,j] = OneBasedBinomialSwitchingEntry(i,j,plasmid_copy_num)
			end
		end
	end

	## now add in transpositions from chromosome to plasmid in the zero-plasmid state.
	switching_matrix[1,1] = 1 - η
	switching_matrix[2,1] = η
	
	return(switching_matrix)
end

# ╔═╡ ad412b59-fe6c-499b-9de3-d0bd59c34288
function SwitchingHypergeometricMatrix(plasmid_copy_num=PCN, η=η₀)
	""" We set a final absorbing state in the Markov chain,
	based on maximum plasmid copy number (like 3 for SC101, 2000 for pUC)."""
	
	## max transposon copy number is plasmid copy number + one chromosomal copy.
	max_TCN = plasmid_copy_num + 1

	## initialize as an m x m identity matrix with ones on the diagonal.
	switching_matrix = Matrix{BigFloat}(I, max_TCN, max_TCN)

	## update binomial entries (overwrite ones and zeros as needed).
	for i in 1:max_TCN
		for j in 1:max_TCN
			if (i <= max_TCN) && (j <= max_TCN)
				switching_matrix[i,j] = OneBasedHypergeometricSwitchingEntry(i,j,plasmid_copy_num)
			end
		end
	end

	## now add in transpositions from chromosome to plasmid in the zero-plasmid state.
	switching_matrix[1,1] = 1 - η
	switching_matrix[2,1] = η
	
	return(switching_matrix)
end

# ╔═╡ a40dc8e6-606c-4eed-b48d-a78837497861
function SwitchingTridiagonalMatrix(plasmid_copy_num=PCN)
	""" We set a final absorbing state in the Markov chain,
	based on maximum plasmid copy number (like 3 for SC101, 2000 for pUC)."""
	## max transposon copy number is plasmid copy number + one chromosomal copy.
	max_TCN = plasmid_copy_num + 1
	n = max_TCN ## syntactic sugar
	K = 0.001
	superdiag = vcat([K for i in 2:n-1], 0)
	diag = vcat(-K, [-2K for i in 1:n-2], 1.0)
	subdiag = [K for i in 1:n-1]
	## have to add an identity matrix to get a stochastic mutation matrix
	switching_matrix = Tridiagonal(subdiag, diag, superdiag) + I(n)
	return(switching_matrix)
end

# ╔═╡ ee5cd57e-cacb-44a3-bc17-95c382c285ee
PCN

# ╔═╡ 9424e197-9e5e-482e-973f-d57d66a2dfcb
MAX_TCN = PCN + 1 ## max transposon copy number is tied to plasmid_copy_number

# ╔═╡ c84145c4-277e-4d5f-8582-b1fec3715780
TetConcSlider = @bind TET_CONC Slider(0:50, default=35, show_value=true)

# ╔═╡ f6840ba3-7cdd-46a6-8ba5-43bddd99a20f
function SelectionDiagonalMatrix(plasmid_copy_num=PCN, Tet_conc=TET_CONC)
	
	## max transposon copy number is plasmid copy number + one chromosomal copy.
	max_TCN = plasmid_copy_num + 1
	
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:max_TCN)
	
	## get the fitness for each fitness class (defined by tetA copy number).
	## use BigFloats to get better numerical precision.
	growth_rate_vec = big.(fitness_function.(tetA_classes, Tet_conc))
	diagonal_growth_matrix = Diagonal(growth_rate_vec)
	return(diagonal_growth_matrix)
end

# ╔═╡ cd57af70-9193-41e5-b963-2af7b9896df5
function MutSelMatrix(plasmid_copy_num=PCN, Tet_conc=TET_CONC)
	mutsel_matrix = SwitchingBinomialMatrix(plasmid_copy_num) * SelectionDiagonalMatrix(plasmid_copy_num, Tet_conc)
	return(mutsel_matrix)
end

# ╔═╡ fe1990d5-7f2f-4782-ac12-fc8f5309d557
function MutSelMatrix2(plasmid_copy_num=PCN, Tet_conc=TET_CONC)
	mutsel_matrix = SwitchingHypergeometricMatrix(plasmid_copy_num) * SelectionDiagonalMatrix(plasmid_copy_num, Tet_conc)
	return(mutsel_matrix)
end

# ╔═╡ cc6cf6cd-2995-4069-b5cf-e38b8c2781eb
function MutSelMatrix3(plasmid_copy_num=PCN, Tet_conc=TET_CONC)
	mutsel_matrix = SwitchingTridiagonalMatrix(plasmid_copy_num) * SelectionDiagonalMatrix(plasmid_copy_num, Tet_conc)
	return(mutsel_matrix)
end

# ╔═╡ 1409a6ad-ab49-4413-b3d5-7ab42d699728
TET_CONC

# ╔═╡ 459c43a7-88c6-4481-a713-3ce62631d41b
begin
	xvec = collect(1:MAX_TCN)
	fitnesses = fitness_function.(xvec, TET_CONC)
	plot(xvec, fitnesses, label="Fitness", xlabel="tetA copy number", ylabel="Fitness")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ 67ce4118-efe4-4117-ba17-48000b046a52
md""" ### let's examine the eigenvalues and eigenvectors of the MutSel matrix, as PCN is varied. """

# ╔═╡ a88246ff-7f53-4ea2-b122-0a315b740448
md""" See page 35 of Nowak's book, Evolutionary Dynamics to see discussion of eigenanalysis and analytical solution of the quasispecies equation, and **most importantly** the original 1988 quasispecies paper by Eigen, McCaskill, and Schuster. \
\
In addition, note that the Kussell and Leibler (2005) Science paper also does an eigenanalysis of the growth-switching matrix in their paper, and has a time-lag interpretation of their core result on Lyapunov exponents (see endnotes 23, 24, 25 in their paper).

"""

# ╔═╡ dc7a7934-2360-418f-99bc-b1d69ff7a86b
function get_eigenvalues(A)
	## Compute eigenvalues and eigenvectors
	eigen_result = eigen(A)

	# Extract eigenvalues
	λ = eigen_result.values
	return λ
end

# ╔═╡ 5e97161e-a5f0-4b8b-ad59-18527b6425c0
function get_eigenvectors(A)
	## Compute eigenvalues and eigenvectors
	eigen_result = eigen(A)

	# Extract eigenvectors
	V = eigen_result.vectors
	return V
end

# ╔═╡ ce669101-0716-466b-8846-bf5d338cd16c
function StudyEigenBehavior(A)

	# Compute eigenvalues and eigenvectors
	eigen_result = eigen(A)

	# Extract eigenvalues and eigenvectors
	λ = eigen_result.values
	V = eigen_result.vectors

	# Diagonalize the matrix
	D = diagm(λ)

	# Verify the diagonalization A = V * D * V^(-1)
	A_reconstructed = V * D * inv(V)

	# Display the results
	println("Original matrix A:")
	println(A)
	println("\nEigenvalues λ:")
	println(λ)
	println("\nEigenvectors V:")
	println(V)
	println("\nDiagonal matrix D:")
	println(D)
	println("\nReconstructed matrix A = V * D * V^(-1):")
	println(A_reconstructed)
	println("check A against A_reconstructed")
	
	return(λ, V, D, A_reconstructed)
end

# ╔═╡ e7536229-5c27-44d9-a58a-5fa090b35554
begin
	testmatrix1 = MutSelMatrix(30, 15)
	testmatrix1_eigenvals = get_eigenvalues(testmatrix1)
	testmatrix1_eigenvals
end

# ╔═╡ 1db8400a-f38f-4aae-be7c-b689f3d9407f
begin
		testmatrix1_eigenvectors = get_eigenvectors(testmatrix1)
		testmatrix1_eigenvectors[:,31]
end

# ╔═╡ cc8af80a-c3d2-40fd-bf24-215f0bd2e1b3
begin
	testmatrix2 = MutSelMatrix(20, 15)
	testmatrix2_eigenvals = get_eigenvalues(testmatrix2)
	testmatrix2_eigenvals
end

# ╔═╡ b3a236a6-68c6-4d02-b9a1-c2cc234748c5
begin
	testmatrix3 = MutSelMatrix(15, 15)
	testmatrix3_eigenvals = get_eigenvalues(testmatrix3)
	testmatrix3_eigenvals
end

# ╔═╡ 9162780a-5c20-4d89-9984-9bf2b9550203
md""" #### Let's define the quasispecies equation here. """

# ╔═╡ 4a679082-b16c-4e91-bb8f-04ffdab18026
function normalize_vector!(u)
    total = sum(big.(u))
    if total != 0
        u .= u / total
    end
end

# ╔═╡ 6654f511-9709-40e8-abdd-ea44c1874fa2
function quasispecies_odefunc(du, u, p, t)

	## pcn and tet_conc needs to be passed in as parameters.
	pcn, tet_conc = p
	
	## Define the ODE system
	A = MutSelMatrix(pcn, tet_conc)

	## Enforce positivity constraint
    u .= max.(u, 0.0)

	## Normalize the vector to ensure it sums to one
	normalize_vector!(u)
	
	## subtract the mean population growth during this dt interval.
	Au = A*u ## sugar to avoid recomputation
	du .= Au - sum(Au) * u
end

# ╔═╡ 46d75d31-4d10-4247-9daf-4561d1b9a9fa
begin
	## parameters used across time course simulations.
	initial_pop_vec = zeros(BigFloat, MAX_TCN)
	## initialize the population with one cell with 1 tetA copy.
	initial_pop_vec[1] = big"1.0"

	## Define the time span
	tspan = (0.0, TIMESPAN)  # Replace with your desired time span
end

# ╔═╡ 548e51f4-8c76-46e9-91a9-3ab4fdbb3fe7
function check_stochastic_matrix(my_matrix)
	## this function is for debugging the stochastic switching matrix.
	num_columns = size(my_matrix, 2)
	for col in 1:num_columns
        column_sum = sum(my_matrix[:, col])
        println(column_sum)
        if isapprox(column_sum, 1.0)
            println("Sum of column $col is 1.0")
        else
            println("Sum of column $col is not equal to 1.0")
        end
    end
end

# ╔═╡ 959627bd-489a-4e66-aac9-63de0df9fc52
md""" ## model population dynamics under constant [Tet]."""

# ╔═╡ 61d1fa00-b844-47ee-bd54-5dce24f243d0
begin
	## run the simulation for constant [Tet] concentration.
	
	## Create an ODEProblem
	prob = ODEProblem(quasispecies_odefunc, initial_pop_vec, tspan, (PCN, TET_CONC))

	## Solve the ODE system
	sol1 = solve(prob, Tsit5())
	
	## Extract the solution
	final_t1 = sol1.t
	final_x1 = sol1.u
end

# ╔═╡ 9b0089e8-7865-4a70-bc1d-b8731f2222a8
begin
	## Extract the subpopulation timecourses.
	## Example for subpopulation 1: x1 = [u[1] for u in final_x]
	## Example for subpopulation 2: x2 = [u[2] for u in final_x]

	# Initialize an empty vector to store vectors
	subpopulation_time_courses = Vector{Vector{BigFloat}}()
	num_subpopulations = length(initial_pop_vec)
	for i in 1:num_subpopulations
		cur_time_course = [u[i] for u in final_x1]
		push!(subpopulation_time_courses, cur_time_course)
	end
end

# ╔═╡ 3d1d138a-a5c4-4961-9550-85b00b50c31a
let
	## plot the first component over time
	x_1 = subpopulation_time_courses[1]
	timecourse_plot = plot(final_t1, x_1, label="1", xlabel="Time", ylabel="Subpopulation size")
	
	## plot the remaining components over time.
	for i in 2:num_subpopulations
		x_i = subpopulation_time_courses[i]
		plot!(timecourse_plot, final_t1, x_i, label=nothing)
	end

	# Add vertical dashed lines every 25 timesteps
	for i in 25:25:TIMESPAN
    	vline!([i], line=:dash, color=:black, label=nothing)
	end
	
	timecourse_plot
end

# ╔═╡ c8768ded-7dc9-4b1a-b2d9-b29ff3aa36d2
let
	## plot the second component over time
	x_2 = subpopulation_time_courses[2]
	timecourse_plot = plot(final_t1, x_2, label="2", xlabel="Time", ylabel="Subpopulation size")
	
	## plot the remaining components over time.
	for i in 3:(num_subpopulations)
		x_i = subpopulation_time_courses[i]
		plot!(timecourse_plot, final_t1, x_i, label=nothing)
	end
	timecourse_plot
end

# ╔═╡ 0d7e2067-8878-4753-a954-f534423504a5
let
	## plot the total population size over time.
	total_population_time_course = reduce(+, subpopulation_time_courses)
	
	plot(final_t1, total_population_time_course, xlabel="Time", ylabel="Total population size")
end

# ╔═╡ c4f7163e-b1a0-4a2c-bbb3-fb228ba229bd
begin ## let's make a matrix saving the population dynamics.

	result_matrix = [] ## Initialize an empty matrix

	for i in 1:length(final_x1)
		total_N = sum(final_x1[i])
		cur_frequency_vec = final_x1[i]/total_N
		##concatenate the state vector horizontally into the result_matrix
		result_matrix = push!(result_matrix, cur_frequency_vec)
	end
end

# ╔═╡ 0876298b-0efc-42ce-a123-c3924f0f76d2
final_t1

# ╔═╡ 387f22bd-8ea6-4ac0-847c-d1dd3f33389e
time_step_slider = @bind cur_timestep Slider(1:length(final_x1), default=1, show_value=true)

# ╔═╡ 912d93e7-c073-46b3-9287-a8a646381492
let
	xvec = collect(1:length(final_x1[1]))
	plot(xvec, result_matrix[cur_timestep], label="", xlabel="tetA copy number", ylabel="Current population distribution")#, ylims=[0,1])
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ 3d85da6c-a48f-46c3-b7e5-a3125c5a985b
md"""
##### calculate the response time in the constant [Tet] population.

We define the response time as the time for the system to hit equilibrium.

"""

# ╔═╡ 1f7f2b7d-5972-49eb-9b51-52f0a98f0d73
function calc_response_time(my_sol)
	response_time = -1
	prev_vec = my_sol(0)
	## IMPORTANT: t is the 'real time' in the timespan-- 
	## and NOT an index for the my_sol.t object.
	for t in 1:TIMESPAN
		cur_vec = my_sol(t)
		diff_vec = cur_vec - prev_vec
		diff_vec_length = sqrt(diff_vec' * diff_vec)
		if diff_vec_length < 1e-4
			response_time = t
			break
		end
		prev_vec = cur_vec
	end
	return response_time
end

# ╔═╡ 3218d18f-f054-488d-84ea-c15cd2c0116f
md""" ### vary plasmid copy number (PCN) and calculate response times."""

# ╔═╡ 6cf74de2-7695-46c7-bc9b-358d7adcd4fb
begin
	## PCN == 20
	initial_pop_vec2 = zeros(BigFloat, 20+1)
	initial_pop_vec2[1] = big"1.0"
	prob2 = ODEProblem(quasispecies_odefunc, initial_pop_vec2, tspan, (20, TET_CONC))
	sol2 = solve(prob2, Tsit5())

	## PCN == 15
	initial_pop_vec3 = zeros(BigFloat, 15+1)
	initial_pop_vec3[1] = big"1.0"
	prob3 = ODEProblem(quasispecies_odefunc, initial_pop_vec3, tspan, (15, TET_CONC))
	sol3 = solve(prob3, Tsit5())

	## PCN == 10
	initial_pop_vec4 = zeros(BigFloat, 10+1)
	initial_pop_vec4[1] = big"1.0"
	prob4 = ODEProblem(quasispecies_odefunc, initial_pop_vec4, tspan, (10, TET_CONC))
	sol4 = solve(prob4, Tsit5())

	## PCN == 5
	initial_pop_vec5 = zeros(BigFloat, 5+1)
	initial_pop_vec5[1] = big"1.0"
	prob5 = ODEProblem(quasispecies_odefunc, initial_pop_vec5, tspan, (5, TET_CONC))
	sol5 = solve(prob5, Tsit5())

	## PCN == 1
	initial_pop_vec6 = zeros(BigFloat, 1+1)
	initial_pop_vec6[1] = big"1.0"
	prob6 = ODEProblem(quasispecies_odefunc, initial_pop_vec6, tspan, (1, TET_CONC))
	sol6 = solve(prob6, Tsit5())
end

# ╔═╡ e66f0c7b-a312-47a7-b4f9-39089d58d66f
md""" #### Response time and eigenvalue result.

Set PCN to 30 and TET_CONC to 35. When [Tet] > maximum PCN in this model, then increasing PCN decreases the response time (time to the stationary distribution) and this behavior corresponds to the maximum eigenvalue for the matrix for each PCN. That is, increasing PCN monotonically increases the maximum eigenvalue for the value matrix $\mathbf{A}$.

However, when maximum PCN > [Tet], this monotonic behavior does not happen. Why? I am not sure. My guess is the following: \
\
From Price's theorem, the rate of change of tetA copy number depends on _both_ the regression slope as well as variance in tetA copy number. When maximum PCN > [Tet], the regression slope is *lower* for the highest PCN matrix, compared to ones with somewhat lower PCNs. Yuanchi made this important observation!
\
The simplest interpretation is that when [Tet] > maximum PCN, having higher PCN allows for "mutation" to the most fit offspring, which then rapidly takes over the population. In this situation, fitness and response time have a clear inverse relationship (higher fitness --> lower response time).

""" 

# ╔═╡ d9795311-2021-4de2-9f10-8d4e9b4cb3bd
calc_response_time(sol1)

# ╔═╡ 351198e1-8ae5-45c0-8ae2-1ea1a0108e55
calc_response_time(sol2)

# ╔═╡ 8b7544e0-0ba5-4fb6-8e7c-d3868f6ab31e
calc_response_time(sol3)

# ╔═╡ 28a629d3-7216-4d16-850a-8bb884fef420
calc_response_time(sol4)

# ╔═╡ 6f97795b-ffbd-415f-b16c-e1628779e518
calc_response_time(sol5)

# ╔═╡ b5b6ce5b-2bca-4e82-8197-043b380ed606
calc_response_time(sol6)

# ╔═╡ 79f6428b-0158-4ad1-a4c4-64868eedb29e
md""" ##### cross-check these response times with the maximum eigenvalue for each matrix. """

# ╔═╡ 8c0fe267-209d-46e9-971c-a612d62d86ab
begin
	sol1matrix = MutSelMatrix(PCN, TET_CONC)
	sol1matrix_eigenvals = [real(x) for x in get_eigenvalues(sol1matrix)]
	maximum(sol1matrix_eigenvals)
end

# ╔═╡ a16a2bac-7987-4479-ab3e-5ab2a665739e
begin
	sol2matrix = MutSelMatrix(20, TET_CONC)
	sol2matrix_eigenvals = [real(x) for x in get_eigenvalues(sol2matrix)]
	maximum(sol2matrix_eigenvals)
end

# ╔═╡ 3ac756a1-4103-45ec-aa4c-8532ed98975d
begin
	sol3matrix = MutSelMatrix(15, TET_CONC)
	sol3matrix_eigenvals = [real(x) for x in get_eigenvalues(sol3matrix)]
	maximum(sol3matrix_eigenvals)
end

# ╔═╡ a9469715-32cf-4d36-b83e-97a98a099055
begin
	sol4matrix = MutSelMatrix(10, TET_CONC)
	sol4matrix_eigenvals = [real(x) for x in get_eigenvalues(sol4matrix)]
	maximum(sol4matrix_eigenvals)
end

# ╔═╡ fe9abb1e-0591-4ec9-930d-8313c1d51fb6
begin
	sol5matrix = MutSelMatrix(5, TET_CONC)
	sol5matrix_eigenvals = [real(x) for x in get_eigenvalues(sol5matrix)]
	maximum(sol5matrix_eigenvals)
end

# ╔═╡ 11e85c84-60ec-4f66-9312-0032492a9773
begin
	sol6matrix = MutSelMatrix(1, TET_CONC)
	sol6matrix_eigenvals = [real(x) for x in get_eigenvalues(sol6matrix)]
	maximum(sol6matrix_eigenvals)
end

# ╔═╡ f7fe6bf0-86c5-4eb5-875b-be935719ac18
md""" ##### calculate mean fitness in the constant [Tet] population."""

# ╔═╡ 7a8ef3b1-2d49-4bbe-b722-16f9c2f21432
begin
mean_fitness_vec1 = [] 
	for i in 1:length(result_matrix)
		cur_pop_vec = result_matrix[i]
		cur_mean_fitness = calc_mean_fitness(cur_pop_vec, TET_CONC)
		append!(mean_fitness_vec1, cur_mean_fitness)
	end
end

# ╔═╡ 30df8e1a-28c3-49ed-8193-d6a594ce90d9
let
	plot(final_t1, mean_fitness_vec1, label="Mean fitness", xlabel="Time", ylabel="Mean fitness")
end

# ╔═╡ a2929307-3f1b-4184-8de4-d00b331daf55
md""" ##### calculate mean copy number in the constant [Tet] population."""

# ╔═╡ e0267c4e-9ea7-4943-9ef0-eeb74018073d
begin
	mean_copy_num_vec1 = []
	for i in 1:length(final_x1)
		cur_pop_vec = result_matrix[i]
		cur_mean_copy_number = calc_mean_tetA_copy_number(cur_pop_vec)
		append!(mean_copy_num_vec1, cur_mean_copy_number)
	end
end

# ╔═╡ 2c19c0fe-8cf7-4dfa-a795-22dd56c5b886
let
	plot(final_t1, mean_copy_num_vec1, label="Mean tetA copy number", xlabel="Time", ylabel="Mean tetA copy number")
end

# ╔═╡ a239233a-7677-4312-9209-d64570195948
md""" ##### calculate fitness variance in the constant [Tet] population."""

# ╔═╡ dd8312c1-509e-468b-a7d0-234e52ffc8a1
begin
fitness_variance_vec1 = [] 
	for i in 1:length(final_x1)
		cur_pop_vec = result_matrix[i]
		cur_fitness_variance = calc_fitness_variance(cur_pop_vec, TET_CONC)
		append!(fitness_variance_vec1, cur_fitness_variance)
	end
end

# ╔═╡ dc03903c-4c3d-40a9-93a5-ecced3018aa1
md""" ###### scale fitness variance by mean fitness to examine Fisher's fundamental theorem of natural selection. """

# ╔═╡ d9be6f05-3951-44db-837c-564cd29334d3
scaled_fitness_variance_vec1 = fitness_variance_vec1 .* [BigFloat(1.0)/BigFloat(x) for x in mean_fitness_vec1]

# ╔═╡ e7f0bb97-1d37-465a-bf5d-666fdb42ac28
let
	plot(final_t1, fitness_variance_vec1, label="Fitness variance", xlabel="Time", ylabel="Fitness variance")
end

# ╔═╡ 2ac10b20-110c-44af-9a17-3df97e62e577
let
	plot(final_t1, fitness_variance_vec1, label="Fitness variance")
	plot!(final_t1, mean_fitness_vec1, label="Mean fitness", xlabel="Time", ylabel="Fitness mean and variance")
end

# ╔═╡ cc32bc71-427a-49c4-a1bf-2b4e6dcaa801
md""" ##### calculate rate of mean fitness increase in the constant [Tet] population."""

# ╔═╡ 33819244-5197-4571-865f-7e1bc9d88b98
## append a zero to the front of the difference vector.
d_mean_fitness_vec1 = [0; diff(mean_fitness_vec1)]

# ╔═╡ 4ea73ac5-3a0f-4a07-80af-447b2c0198d7
let
	plot(final_t1, scaled_fitness_variance_vec1, label="Scaled Fitness variance")
	plot!(final_t1, d_mean_fitness_vec1, label="Rate of change of mean fitness", xlabel="Time", ylabel="Fitness variance and\nrate of mean fitness change")
end

# ╔═╡ 4231b6a6-a7f4-4960-9d92-707939c49258
md""" ##### calculate rate of copy number change in the constant [Tet] population."""

# ╔═╡ 1f0a45e5-c475-4f6d-8f01-ab512ce333f9
## append a zero to the front of the difference vector.
d_mean_copy_num_vec1 = [0; diff(mean_copy_num_vec1)]

# ╔═╡ 44e93e91-1886-4f40-9578-1228dabbfa79
md""" ##### calculate tetA copy number--fitness covariance in the constant [Tet] population. """ 

# ╔═╡ f4df8237-23a9-4631-873c-2f537c415cd1
begin
	copy_num_covariance_vec1 = []
	for i in 1:length(final_x1)
		cur_pop_vec = result_matrix[i]
		cur_copy_num_covariance = calc_tetA_copy_number_fitness_covariance(cur_pop_vec, TET_CONC)
		
		append!(copy_num_covariance_vec1, cur_copy_num_covariance)
	end
end

# ╔═╡ aa098647-6390-4858-83f0-3f52e6776420
scaled_copy_num_covariance_vec1 = copy_num_covariance_vec1 .* [BigFloat(1.0)/BigFloat(x) for x in mean_fitness_vec1]

# ╔═╡ 4724d83c-201f-453f-aeac-52202712d97a
let
	plot(final_t1, copy_num_covariance_vec1, label="tetA copy number fitness covariance")
	plot!(final_t1, d_mean_copy_num_vec1, label="Rate of change of mean tetA copy number", xlabel="Time", ylabel="tetA copy number derivative and variance")
end

# ╔═╡ 472637ca-befd-49b6-9c5e-c829e7e84ede
let
	plot(final_t1, scaled_copy_num_covariance_vec1, label="scaled tetA copy number fitness covariance")
	plot!(final_t1, d_mean_copy_num_vec1, label="Rate of change of mean tetA copy number", xlabel="Time", ylabel="tetA copy number derivative and variance")
end

# ╔═╡ 1f84ef10-2776-415b-ad60-6c2a5eb5b568
d_mean_copy_num_vec1

# ╔═╡ fef4d315-5440-40db-9e35-122187b55567
mean_copy_num_vec1

# ╔═╡ a4fcf0b1-9625-4be3-a71b-0a4966a017e4
md""" ## model the Tet pulse conditions of the Darwin experiment."""

# ╔═╡ 208f6604-701a-47f9-b934-1e4e546ee722
## We model [Tet] pulses over time by dividing the current time by 20, and
## setting [Tet] on if in the first half.
TetPulseFunction = t -> mod(t,100) < 50 ? TET_CONC : 0

# ╔═╡ 15129817-9644-4d23-9c3f-c46fbc6bc267
function pulse_quasispecies_odefunc(du, u, p, t)
	## Define the ODE system.
	## p is a TetPulseFunction of time that is passed in as a parameter.
	cur_tet_conc = p(t)
	A = MutSelMatrix(PCN, cur_tet_conc)

	## set subpopulations smaller than ϵ to zero, to model extinction.
	##ϵ = 1e-7
	## use a one-line lambda function.
	##set_values_less_than_ϵ_to_zero = x -> ifelse(x < ϵ, 0.0, x)
	##u .= set_values_less_than_ϵ_to_zero.(u)
	
	## Enforce positivity constraint
    u .= max.(u, 0.0)

	## Normalize the vector to ensure it sums to one
	normalize_vector!(u)
	
	## This is more stable-- the numerical error in the fitness calculation
	## is large enough to cause errors.
	Au = A*u ## sugar to avoid recomputation
    du .= Au - sum(Au) * u
end

# ╔═╡ 0b5b688b-c508-49e0-bb2a-edb4fe371d8d
begin
	## run the simulation for pulses of [Tet] concentration.
	
	## Create an ODEProblem
	pulse_prob = ODEProblem(pulse_quasispecies_odefunc, initial_pop_vec, tspan, TetPulseFunction)

	## Solve the ODE system
	pulse_sol = solve(pulse_prob, Tsit5())

	## Extract the solution
	pulse_final_t = pulse_sol.t
	pulse_final_x = pulse_sol.u
end

# ╔═╡ 99fc3172-c63c-4289-92ea-3f7f75e8f2bd
begin ## let's make a matrix saving the population dynamics.

	pulse_result_matrix = [] ## Initialize an empty matrix

	for i in 1:length(pulse_final_x)
		total_N = sum(pulse_final_x[i])
		cur_frequency_vec = pulse_final_x[i]/total_N
		##concatenate the state vector horizontally into the pulse_result_matrix
		pulse_result_matrix = push!(pulse_result_matrix, cur_frequency_vec)
	end
end

# ╔═╡ ffc9c8b9-7e67-4b62-b87d-b89a47d8ef08
pulse_time_step_slider = @bind pulse_cur_timestep Slider(1:length(pulse_final_x), default=1, show_value=true)

# ╔═╡ c6fe6a20-0e83-4392-8d9e-eb1ae5c47409
let
	pulse_xvec = collect(1:length(pulse_final_x[1]))
	plot(pulse_xvec, pulse_result_matrix[pulse_cur_timestep], label="", xlabel="tetA copy number", ylabel="Current population distribution") ##, ylims=[0,1])
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ 4e920578-22b6-46fa-94e1-809f80e86425
pulse_result_matrix[30]

# ╔═╡ 71943ad7-d763-48b8-96f2-dcc91e86a666
md""" ### plot the antibiotic pulse regime over time."""

# ╔═╡ 53e1105a-62de-4fc0-b9eb-60172596ec3f
begin
	tet_pulse_vec = [TetPulseFunction(t) for t in pulse_final_t]
end

# ╔═╡ 59b79d42-900e-4946-a593-08200503cf34
plot(pulse_final_t, tet_pulse_vec,label="Tet pulse regime", xlabel="Time", ylabel="Tet concentration")

# ╔═╡ bfbba79b-4a26-4ef5-910f-b001e82ac1e6
md""" ##### calculate mean fitness in the population."""

# ╔═╡ 966f754c-2467-4cda-b1fe-8cb28ecfc0e7
begin
pulse_mean_fitness_vec = [] 
	for i in 1:length(pulse_result_matrix)
		cur_pop_vec = pulse_result_matrix[i]
		cur_tet_conc = TetPulseFunction(pulse_final_t[i])
		cur_mean_fitness = calc_mean_fitness(cur_pop_vec, cur_tet_conc)
		append!(pulse_mean_fitness_vec, cur_mean_fitness)
	end
end

# ╔═╡ dc7b5cbe-8c1f-48f8-960e-82d5be823c28
let
	plot(pulse_final_t, pulse_mean_fitness_vec, label="Mean fitness", xlabel="Time", ylabel="Mean fitness")
end

# ╔═╡ 445c9f91-c742-4491-8728-04b7831fdfa5
md""" ##### calculate mean copy number in the population."""

# ╔═╡ f33ff779-5455-4afd-841d-2eb81238e399
begin
	pulse_mean_copy_num_vec = []
	for i in 1:length(pulse_final_x)
		cur_pop_vec = pulse_result_matrix[i]
		cur_mean_copy_number = calc_mean_tetA_copy_number(cur_pop_vec)
		append!(pulse_mean_copy_num_vec, cur_mean_copy_number)
	end
end

# ╔═╡ 7fee8b2c-4cfe-4ae2-94fc-2660de749e74
let
	plot(pulse_final_t, pulse_mean_copy_num_vec, label="Mean tetA copy number", xlabel="Time", ylabel="Mean tetA copy number")
end

# ╔═╡ eec10cb1-3880-4a5e-8996-abb2bb48a12f
md""" ##### calculate fitness variance in the population."""

# ╔═╡ fa8acd76-9abe-4c7f-bc7f-fcb35ea2ebd7
begin
pulse_fitness_variance_vec = [] 
	for i in 1:length(pulse_final_x)
		cur_pop_vec = pulse_result_matrix[i]
		cur_tet_conc = tet_pulse_vec[i]
		cur_fitness_variance = calc_fitness_variance(cur_pop_vec, cur_tet_conc)
		append!(pulse_fitness_variance_vec, cur_fitness_variance)
	end
end

# ╔═╡ e2ec8e48-fe32-4e56-869b-194a5468c9f3
let
	plot(pulse_final_t, pulse_fitness_variance_vec, label="Fitness variance", xlabel="Time", ylabel="Fitness variance")
end

# ╔═╡ 590b3b14-0ca0-4a3c-8212-ab2f441b914f
let
	plot(pulse_final_t, pulse_fitness_variance_vec, label="Fitness variance")
	plot!(pulse_final_t, pulse_mean_fitness_vec, label="Mean fitness", xlabel="Time", ylabel="Fitness mean and variance")
end

# ╔═╡ 3f5c1fa4-8e08-47a7-9cae-ac91a8e0df9e
md""" ##### calculate rate of mean fitness increase in the population."""

# ╔═╡ 6b2ec37f-cbe5-4a22-91eb-05baf160d452
## append a zero to the front of the difference vector.
d_pulse_mean_fitness_vec = [0; diff(pulse_mean_fitness_vec)]

# ╔═╡ 6234576b-51e6-4967-adf1-5859185d3a48
let
	plot(pulse_final_t, pulse_fitness_variance_vec, label="Fitness variance")
	plot!(pulse_final_t, d_pulse_mean_fitness_vec, label="Rate of change of mean fitness", xlabel="Time", ylabel="Fitness variance\nand rate of mean fitness change")
end

# ╔═╡ 0e0b6318-a434-4b11-b719-78bb4a39557f
md""" ##### calculate rate of copy number change in the population."""

# ╔═╡ 19fd0381-bf9a-4e39-aeee-2187b6275140
## append a zero to the front of the difference vector.
d_pulse_mean_copy_num_vec = [0; diff(pulse_mean_copy_num_vec)]

# ╔═╡ 1e347047-233d-4162-8cf3-3378b2460f7a
md""" ##### calculate tetA copy number--fitness covariance. """ 

# ╔═╡ d19c78c2-7099-405f-94b9-cb091aa63df1
begin
	copy_num_covariance_vec = []
	for i in 1:length(pulse_final_x)
		cur_pop_vec = pulse_result_matrix[i]
		cur_Tet_conc = tet_pulse_vec[i]
		## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
		## with a minimum of 1 copy (on the chromosome).
		cur_copy_num_covariance = calc_tetA_copy_number_fitness_covariance(cur_pop_vec, cur_Tet_conc)
		
		append!(copy_num_covariance_vec, cur_copy_num_covariance)
	end
end

# ╔═╡ c9670300-5cee-42c0-bc65-e5aba726297b
md""" ##### and calculate tetA copy number--fitness covariance scaled by mean fitness. """ 

# ╔═╡ 9b9b54a3-e4e4-467c-9b01-c280315d0ee4
scaled_copy_num_covariance_vec = copy_num_covariance_vec .* [big(1.0)/big(x) for x in pulse_mean_fitness_vec]

# ╔═╡ c9e51c52-0766-4d72-93cb-39ce14c7edeb
test_scaled_copy_num_covariance_vec = copy_num_covariance_vec .* pulse_mean_fitness_vec

# ╔═╡ 1e80925c-b770-4a6e-8332-3131b1ff004b
pulse_mean_fitness_vec

# ╔═╡ 44638f03-9a8f-45aa-8120-2964e7f728ef
[big(1.0)/big(x) for x in pulse_mean_fitness_vec]

# ╔═╡ fa74a3d2-6f13-4286-9165-b3d34f534d4a
md"""
##### debugging plot: comparing fitness-scaled covariance to just covariance.
"""

# ╔═╡ d41de4e9-66a7-4e20-9dfd-7ffc1a430300
let
	plot(pulse_final_t, copy_num_covariance_vec, label="tetA copy number-fitness covariance", xlabel="Time", ylabel="tetA copy number--fitness covariance")
	plot!(pulse_final_t, scaled_copy_num_covariance_vec, label="scaled-tetA copy number-fitness covariance", legend=:left)	
end

# ╔═╡ d3fa42c0-deb3-400d-9b62-ca7b32749955
let
	plot(pulse_final_t, copy_num_covariance_vec, label="tetA copy number fitness covariance")
	plot!(pulse_final_t, pulse_mean_copy_num_vec, label="mean tetA copy number", xlabel="Time", ylabel="tetA copy number mean and variance")
end

# ╔═╡ 7d1fda00-8528-4186-9a8d-80ea6ab9071d
let
	plot(pulse_final_t, copy_num_covariance_vec, label="tetA copy number fitness covariance")
	plot!(pulse_final_t, d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="Time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:right)
end

# ╔═╡ b618efad-30b6-42c3-ba78-4add20b925fc
let
	plot(pulse_final_t, scaled_copy_num_covariance_vec, label="scaled tetA copy number fitness covariance")
	plot!(pulse_final_t, d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="Time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:left)
end

# ╔═╡ 7a3dc053-90f5-4229-be06-4b503bb17f16
md""" ##### Result. I tested my hypothesis that when plasmid copy number is increased, stochastic segregation during cell division generates more diverse daughter cells. In some ways yes, but there is some important subtlety here.
\

Consider all possible clonal parental populations $\mathbf{x}$, with 1 in the clonal state, and zeros everywhere else. By multiplying $\mathbf{x}$ with $\mathbf{S}$, we "select" the column of $\mathbf{S}$ that corresponds to that particular parental clone. The variance of this column of $\mathbf{S}$ represents the variance in offspring per dt. \
\
When PCN increases, the *absolute* number of possible offspring states increases. But this means that the probability of generating any particular offspring state decreases as $1/(PCN+1)$, since the probability distribution over offspring states has to sum to one. \
\
Since each column vector of $\mathbf{S}$ sums to one ($\mathbf{S}$ is a stochastic matrix), the mean of each offspring distribution vector has to be $1/(PCN+1)$, where $PCN+1$ is the number of states (columns and rows) in the S matrix (i.e. PCN+1 = maximum tetA copy number). \
\
This means that the variance in offspring per dt *decreases* as PCN increases-- the fraction of new mutants in each state decreases, as the number of states increases. \
\
The relative variance $Var(\mathbf{S_{-b}})/Mean(\mathbf{S_{-b}})$ does increase a *little bit* but this is a pretty small effect. \
\
Altogether, this numerical result suggests the following:\
1) higher PCN only increases the speed of adaptation when its most fit mutant cannot be generated by lower PCN.
2) in some cases, lower PCN may even adapt faster than higher PCN, in the specific case that it generates a larger initial fraction of the most fit mutant in the given environment.

"""

# ╔═╡ 8868caec-6e6f-4c9b-b20c-0ac8db38847c
testmatrix1

# ╔═╡ 520679b5-ced1-4a9f-bc41-fdf31a7c40ce
testmatrix2

# ╔═╡ 4538ec89-e65b-4c55-a719-cdd80b44e478
testmatrix3

# ╔═╡ c57f83dd-b672-425d-9b32-cc418ea7c4b1
Smatrix1 = SwitchingBinomialMatrix(PCN)

# ╔═╡ 434db880-3825-4f3d-b736-af8109bedf4d
Smatrix2 = SwitchingBinomialMatrix(20)

# ╔═╡ 45452600-600b-4ccc-a6e6-7a92febdd80b
Smatrix3 = SwitchingBinomialMatrix(15)

# ╔═╡ fa4a59d9-079a-4feb-a74a-4630f21b6442
Smatrix4 = SwitchingBinomialMatrix(10)

# ╔═╡ 641eff1b-401c-4dad-88ac-b2aa0a3052dc
Smatrix5 = SwitchingBinomialMatrix(5)

# ╔═╡ 175db452-1d36-4a3c-9c24-1da14d7d769a
Smatrix6 = SwitchingBinomialMatrix(1)

# ╔═╡ a2d1ffca-4981-4db2-a6eb-7b6aad2f2b75
function calculateMatrixColumnVariances(matrix)
	rows, cols = size(matrix)
	column_vec_variances = BigFloat[]

	for j in 1:cols
		cur_col = matrix[:, j]
		cur_col_variance = var(cur_col)
		push!(column_vec_variances, cur_col_variance)
	end
	
	return column_vec_variances
end

# ╔═╡ a9118a54-a13e-4bc6-9b7d-1f7b9d79fd19
calculateMatrixColumnVariances(testmatrix1)

# ╔═╡ 070a672b-7955-47c1-b381-b507f44bb149
calculateMatrixColumnVariances(testmatrix2)

# ╔═╡ 4e4bda01-cdb4-4f44-8827-bb398fb12fbe
calculateMatrixColumnVariances(testmatrix3)

# ╔═╡ 75a17150-1e2a-4914-8f0f-e8f1c41c0546
calculateMatrixColumnVariances(Smatrix1)

# ╔═╡ f7aac374-5012-44c4-9999-ae86b7d4b6c9
calculateMatrixColumnVariances(Smatrix2)

# ╔═╡ f4c57c28-dab0-4d63-9d93-972e14b96a80
calculateMatrixColumnVariances(Smatrix3)

# ╔═╡ 7b2e7e05-9a98-4442-86c1-03621f3bbc77
calculateMatrixColumnVariances(Smatrix4)

# ╔═╡ b870c7b6-449b-40b4-badd-338fa4ce8a12
calculateMatrixColumnVariances(Smatrix5)

# ╔═╡ 77fea351-fe99-4d57-9602-b69ac915afec
calculateMatrixColumnVariances(Smatrix6)

# ╔═╡ 14e01b5e-aad2-48d2-ac9c-10a8bb595fd6
md""" since each column sums to one, each column has the same mean: 1/k, where k is the number of rows, or 1/(n+1), where n is the plasmid copy number. """

# ╔═╡ a9c5c21a-69de-4fa6-ae89-33e25b74ce1d
md""" Calculate var(x)/mean(x) for each column vector for each Smatrix."""

# ╔═╡ b8dffaff-5788-4cb4-bf40-1c07d505fe3f
101*calculateMatrixColumnVariances(Smatrix1)

# ╔═╡ 34e72ef5-08b4-4cb5-8164-90ef49bf29cc
21*calculateMatrixColumnVariances(Smatrix2)

# ╔═╡ c0209008-ab91-4c60-9449-d708a97f230c
16*calculateMatrixColumnVariances(Smatrix3)

# ╔═╡ 6c43017a-966a-4d43-a5ff-8ef55ce5ccc7
11*calculateMatrixColumnVariances(Smatrix4)

# ╔═╡ c0140abb-4ee8-4f34-b26b-d202030790ac
6*calculateMatrixColumnVariances(Smatrix5)

# ╔═╡ d43486f7-5d63-4a9d-9779-89c1b2469cf7
2*calculateMatrixColumnVariances(Smatrix6)

# ╔═╡ 08728379-9a61-4177-8d80-3eff252b8232
md""" ##### increasing PCN increases fitness variance among offspring, but the magnitude of the effect is extremely small."""

# ╔═╡ 144dbd76-005e-4a76-ba77-db881e965e04
function calculateOffspringFitnessVariances(matrix, Tet_conc=TET_CONC)
	rows, cols = size(matrix)
	column_vec_fitness_variances = BigFloat[]

	for j in 1:cols
		cur_col = matrix[:, j]
		cur_col_fitness_variance = calc_fitness_variance(cur_col, Tet_conc)
		push!(column_vec_fitness_variances, cur_col_fitness_variance)
	end
	
	return column_vec_fitness_variances
end

# ╔═╡ e5e20978-e24a-49e3-bc02-00463ccec77b
function calculateOffspringTetAFitnessCovariances(matrix, Tet_conc=TET_CONC)
	rows, cols = size(matrix)
	column_vec_tetA_fitness_covariances = BigFloat[]

	for j in 1:cols
		cur_col = matrix[:, j]
		cur_col_tetA_fitness_covariance = calc_tetA_copy_number_fitness_covariance(cur_col, Tet_conc)
		push!(column_vec_tetA_fitness_covariances, cur_col_tetA_fitness_covariance)
	end
	
	return column_vec_tetA_fitness_covariances
end

# ╔═╡ db5cb755-af4c-4a29-b76d-3cd71284e8a1
calculateOffspringFitnessVariances(testmatrix1)

# ╔═╡ 51b3e6bd-82dc-411f-a6f7-c9b27602f185
calculateOffspringFitnessVariances(testmatrix2)

# ╔═╡ dad60661-2c76-40d7-997f-fd1139a2aea5
calculateOffspringFitnessVariances(testmatrix3)

# ╔═╡ 4f5f3fe2-08d0-4b0f-b9e2-05171e7a8456
calculateOffspringFitnessVariances(Smatrix1)

# ╔═╡ e455e16f-e2d4-4984-9dcb-329d02ba5d95
calculateOffspringFitnessVariances(Smatrix2)

# ╔═╡ 84c35c75-3f3a-4979-a9c1-705f5ae0cc8f
calculateOffspringFitnessVariances(Smatrix3)

# ╔═╡ bbb8b8aa-5fd4-4408-991c-e9f281d57647
calculateOffspringFitnessVariances(Smatrix4)

# ╔═╡ 4aa9e4f8-f49e-4ca9-89d2-d167bfd7d5af
calculateOffspringFitnessVariances(Smatrix5)

# ╔═╡ 140c83dd-7112-4651-aa68-2dd441b239c0
calculateOffspringFitnessVariances(Smatrix6)

# ╔═╡ fa246234-6f8e-4808-af59-5e5cf163718c
md""" ##### increasing PCN increases tetA-copy-number-fitness covariance among offspring, but the magnitude of the effect is extremely small."""

# ╔═╡ 26d9bac8-0e64-40c9-ac6f-42df7e2595ba
calculateOffspringTetAFitnessCovariances(testmatrix1)

# ╔═╡ 1193b90a-0528-4877-9abb-3670b608be86
calculateOffspringTetAFitnessCovariances(testmatrix2)

# ╔═╡ 2dc593ca-1f50-4714-a5c7-7d7de04e01d3
calculateOffspringTetAFitnessCovariances(testmatrix3)

# ╔═╡ a689a276-5ef2-47f1-9577-31a4513cb398
calculateOffspringTetAFitnessCovariances(Smatrix1)

# ╔═╡ 2e9f116c-97f0-4d0a-9ac0-054a58eb6a9d
calculateOffspringTetAFitnessCovariances(Smatrix2)

# ╔═╡ e65fb06f-b80b-4c3b-aadc-47cf38a608cd
calculateOffspringTetAFitnessCovariances(Smatrix3)

# ╔═╡ 4f3640f8-8224-4acb-a3df-94a71bb62347
calculateOffspringTetAFitnessCovariances(Smatrix4)

# ╔═╡ 65389d85-99f1-4682-b35b-8eff20cfee1c
calculateOffspringTetAFitnessCovariances(Smatrix5)

# ╔═╡ c54a479b-97e0-4f39-89a1-74cfb38f659a
calculateOffspringTetAFitnessCovariances(Smatrix6)

# ╔═╡ 2109d5d2-a710-4f47-bd01-1c613f16dd83
md""" ##### TODO: check the assumption that generating more diverse daughter cells increases the speed at which selection can tune population-level gene expression in response to environmental change."""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DifferentialEquations = "~7.12.0"
Distributions = "~0.25.104"
Plots = "~1.39.0"
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "036588975bb5fc9a27cddd53733874a3766ca07a"

[[deps.ADTypes]]
git-tree-sha1 = "41c37aa88889c171f1300ceac1313c06e891d245"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.6"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "bbec08a37f8722786d87bedf84eae19c020c4efa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "b08a4043e1c14096ef8efe4dd97e07de5cacf240"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.4.5"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "27baf04c642465b4289179f29bb7127f0673d4f1"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.4.0"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "ForwardDiff", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools", "Tricks", "TruncatedStacktraces", "UnPack"]
git-tree-sha1 = "dd234c9a030350d5ff4c45761d6cad0cfb358cb9"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.6.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"
    BoundaryValueDiffEqOrdinaryDiffEqExt = "OrdinaryDiffEq"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"
    OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "886826d76ea9e72b35fcd000e535588f7b60f21d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "8cfa272e8bdedfa88b6aefbbca7c19f1befac519"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "9534eac1eb01f9fbc9319799cabfdb6ec9c56e42"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.45.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces"]
git-tree-sha1 = "f45adda33e3eb24ae1428dfa9b57ded8248ced99"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.145.4"

    [deps.DiffEqBase.extensions]
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "cf334da651a6e42c50e1477d6ab978f1b8be3057"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.36.1"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "319377c927a4aa1f491228b2ac23f3554a3497c6"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.20.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "8864b6a953eeba7890d23258aca468d90ca73fd6"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.12.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9242eec9b7e2e14f9952e8ea1c7e31a50501d587"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.104"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
deps = ["Adapt"]
git-tree-sha1 = "2efe862de93cd87f620ad6ac9c9e3f83f1b2841b"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.6.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "602e4585bcbd5a25bc06f514724593d13ff9e862"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.25.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "178316d87f883f0702e79d9c83a8049484c9f619"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "a6e756a880fc419c8b41592010aebe6a5ce09136"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.8"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "b12f05108e405dadcc2aff0008db7f831374e051"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "73d1214fec245096717847c62d389a5d2ac86504"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.22.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9a68d75d466ccc1218d0552a8e1631151c569545"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "abbbb9ec3afd783a7cbd82ef01dcd088ea051398"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "UnPack"]
git-tree-sha1 = "c451feb97251965a9fe40bacd62551a72cc5902c"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.10.1"
weakdeps = ["FastBroadcast"]

    [deps.JumpProcesses.extensions]
    JumpProcessFastBroadcastExt = "FastBroadcast"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "884c2968c2e8e7e6bf5956af88cb46aa745c854b"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.1"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "8a6837ec02fe5fb3def1abc907bb802ef11a0729"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays"]
git-tree-sha1 = "9cfca23ab83b0dfac93cb1a1ef3331ab9fe596a5"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "1.8.3"
weakdeps = ["StaticArrays"]

    [deps.LazyArrays.extensions]
    LazyArraysStaticArraysExt = "StaticArrays"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "97dc499678d50d989f1a74170840808641ce9880"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.22.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveEnzymeExt = ["Enzyme", "EnzymeCore"]
    LinearSolveFastAlmostBandedMatricesExt = ["FastAlmostBandedMatrices"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "0f5648fbae0d015e3abe5867bca2b362f67a5894"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.166"

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.LoopVectorization.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "b211c553c199c111d998ecdaf7623d1b89b69f93"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.12"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "78f6e33434939b0ac9ba1df81e6d005ee85a7396"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "2.1.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "a85c6a98c9e5a2a7046bc1bb89f28a3241e1de4d"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "EnumX", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "72b036b728461272ae1b1c3f7096cb4c319d8793"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.4.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSymbolicsExt = "Symbolics"
    NonlinearSolveZygoteExt = "Zygote"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc6e1927ac521b659af340e0ca45828a3ffc748f"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.12+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "01f85d9269b13fedc61e63cc72ee2213565f7a72"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.8"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "7857772d839bf0773f5db9e456149c832ba5facb"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.65.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "fca25670784a1ae44546bcb17288218310af2778"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.9"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "01ac95fca7daabe77a9cb705862bd87016af9ddb"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.13"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9ebcd48c498668c7fa0e97a9cae873fbee7bfee1"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "c860e84651f58ce240dd79e5d9e055d55234c35a"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.2"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "3a81c52c1fa8522e359c63b2c52d9625ea1e0992"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.3.2"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "8bc86c78c7d8e2a5fe559e3721c0f9c9e303b2ed"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.21"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "2c3706caa9adab5031f30937699c46e159ae477f"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.13.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "8d672bd91dc432fb286b6d4bcf1a5dc417e932a3"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.2.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolvePolyesterForwardDiffExt = "PolyesterForwardDiff"

    [deps.SimpleNonlinearSolve.weakdeps]
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "c281e11db4eacb36a292a054bac83c5a0aca2a26"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.15.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "fba11dbe2562eecdfcac49a05246af09ee64d055"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.8.1"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "a735fd5053724cf4de31c81b4e2cc429db844be5"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.0.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "753219de57ac7aab0feb88871d3c51e0eb5e3b03"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.64.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "d6415f66f3d89c615929af907fdc6a3e17af0d8c"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.2"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "ded52f017fe7faa3d004427f10ecce4c0491c16a"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.23.1"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "ba4d38faeb62de7ef47155ed321dce40a549c305"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.2+0"

[[deps.SymbolicIndexingInterface]]
git-tree-sha1 = "be414bfd80c2c91197823890c66ef4b74f5bf5fe"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TranscodingStreams]]
git-tree-sha1 = "1fbeaaca45801b4ba17c251dd8603ef24801dd84"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.2"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "fadebab77bf3ae041f77346dd1c290173da5a443"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.20"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "7209df901e6ed7489fe9b7aa3e46fb788e15db85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.65"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522b8414d40c4cbbab8dee346ac3a09f9768f25d"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.5+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═60f6f4ca-a684-11ee-26bd-6b50dd5567a6
# ╠═fea3c300-7aaf-4da8-bed7-6e9426b4ac3d
# ╠═498a2906-9a95-47ef-b3e3-a5ef24ee0a77
# ╟─772351cf-cddd-45d8-9f34-f3ef6ccb2a66
# ╟─d3b956b4-885a-447c-bd00-de5fa18fdcad
# ╟─624831f0-b669-44e9-ae78-ab0e09716ea4
# ╟─37d2e248-1f6d-4c39-8cc8-8080a853cb43
# ╟─c270d434-d426-47f7-aec2-1423f3354fa6
# ╟─a7b38e19-1f54-4ce1-9977-79e9968ce84f
# ╟─623f1230-77ff-4eeb-865f-01db8f673c01
# ╟─5d83343b-e40a-4193-8283-dd33e453855a
# ╟─f5a21fa6-f7b4-42f5-8ad4-58bb813b98b6
# ╟─853b7a17-bc10-4e1c-8e73-2a0e80592863
# ╟─0dcf8ff4-bff7-448c-be66-92d3555d9efb
# ╟─b1e5c35b-35ae-48c2-babf-2058a91f1ea8
# ╟─41de1814-a88d-4c4e-901b-a69a0a6ca655
# ╟─22ac1b86-a2b6-4562-adb3-61718a1c9a59
# ╟─60e1b0f7-b91f-4121-80f8-673ce117dccf
# ╠═6d43ce5f-5e34-4d4e-bcd1-b876c5744be8
# ╠═ddc8420f-1212-459f-91b1-be62db8032dd
# ╠═9e38c88c-fb6a-4623-b3aa-c7c487ff919f
# ╠═6dbe79d8-f3f5-45d7-946a-7454542df74e
# ╠═ce9769d3-a234-4718-b3b1-2ee6fda2d733
# ╠═29daeb98-190b-4972-9a7b-1407ed3c4b1a
# ╠═4c7767df-bcc4-4fb7-9d66-c807ac8ce6e1
# ╠═6a470e0e-de4a-4398-893a-d2edbdf44b9f
# ╠═c04868af-914e-406d-a82c-6e75b189600f
# ╠═3b3e97a8-ccff-4ab6-94d6-e6d2d7c92323
# ╠═492dbbd5-28d3-45da-99ca-b5f7d5e3c31e
# ╠═8ecaeeb3-2a93-48f5-ba21-49899d20c538
# ╠═7715798e-aaaa-46a9-bd16-14b866cf39f6
# ╠═3e8936df-d80e-4b66-b035-a1007c8f5dca
# ╠═b5fcc845-0d6b-41ab-b27f-05f67ff05325
# ╠═ad412b59-fe6c-499b-9de3-d0bd59c34288
# ╠═f6840ba3-7cdd-46a6-8ba5-43bddd99a20f
# ╠═a40dc8e6-606c-4eed-b48d-a78837497861
# ╠═cd57af70-9193-41e5-b963-2af7b9896df5
# ╠═fe1990d5-7f2f-4782-ac12-fc8f5309d557
# ╠═cc6cf6cd-2995-4069-b5cf-e38b8c2781eb
# ╠═5e33d199-29f5-493e-9863-f02b61711e26
# ╠═a2d049bb-4d27-469d-8626-d0c0fedbced3
# ╠═ee5cd57e-cacb-44a3-bc17-95c382c285ee
# ╠═9424e197-9e5e-482e-973f-d57d66a2dfcb
# ╠═c84145c4-277e-4d5f-8582-b1fec3715780
# ╠═1409a6ad-ab49-4413-b3d5-7ab42d699728
# ╠═459c43a7-88c6-4481-a713-3ce62631d41b
# ╟─67ce4118-efe4-4117-ba17-48000b046a52
# ╟─a88246ff-7f53-4ea2-b122-0a315b740448
# ╠═dc7a7934-2360-418f-99bc-b1d69ff7a86b
# ╠═5e97161e-a5f0-4b8b-ad59-18527b6425c0
# ╠═ce669101-0716-466b-8846-bf5d338cd16c
# ╠═e7536229-5c27-44d9-a58a-5fa090b35554
# ╠═1db8400a-f38f-4aae-be7c-b689f3d9407f
# ╠═cc8af80a-c3d2-40fd-bf24-215f0bd2e1b3
# ╠═b3a236a6-68c6-4d02-b9a1-c2cc234748c5
# ╠═9162780a-5c20-4d89-9984-9bf2b9550203
# ╠═4a679082-b16c-4e91-bb8f-04ffdab18026
# ╠═6654f511-9709-40e8-abdd-ea44c1874fa2
# ╠═46d75d31-4d10-4247-9daf-4561d1b9a9fa
# ╠═548e51f4-8c76-46e9-91a9-3ab4fdbb3fe7
# ╠═959627bd-489a-4e66-aac9-63de0df9fc52
# ╠═61d1fa00-b844-47ee-bd54-5dce24f243d0
# ╠═9b0089e8-7865-4a70-bc1d-b8731f2222a8
# ╠═3d1d138a-a5c4-4961-9550-85b00b50c31a
# ╠═c8768ded-7dc9-4b1a-b2d9-b29ff3aa36d2
# ╠═0d7e2067-8878-4753-a954-f534423504a5
# ╠═c4f7163e-b1a0-4a2c-bbb3-fb228ba229bd
# ╠═0876298b-0efc-42ce-a123-c3924f0f76d2
# ╠═387f22bd-8ea6-4ac0-847c-d1dd3f33389e
# ╠═912d93e7-c073-46b3-9287-a8a646381492
# ╠═3d85da6c-a48f-46c3-b7e5-a3125c5a985b
# ╠═1f7f2b7d-5972-49eb-9b51-52f0a98f0d73
# ╠═3218d18f-f054-488d-84ea-c15cd2c0116f
# ╠═6cf74de2-7695-46c7-bc9b-358d7adcd4fb
# ╠═e66f0c7b-a312-47a7-b4f9-39089d58d66f
# ╠═d9795311-2021-4de2-9f10-8d4e9b4cb3bd
# ╠═351198e1-8ae5-45c0-8ae2-1ea1a0108e55
# ╠═8b7544e0-0ba5-4fb6-8e7c-d3868f6ab31e
# ╠═28a629d3-7216-4d16-850a-8bb884fef420
# ╠═6f97795b-ffbd-415f-b16c-e1628779e518
# ╠═b5b6ce5b-2bca-4e82-8197-043b380ed606
# ╠═79f6428b-0158-4ad1-a4c4-64868eedb29e
# ╠═8c0fe267-209d-46e9-971c-a612d62d86ab
# ╠═a16a2bac-7987-4479-ab3e-5ab2a665739e
# ╠═3ac756a1-4103-45ec-aa4c-8532ed98975d
# ╠═a9469715-32cf-4d36-b83e-97a98a099055
# ╠═fe9abb1e-0591-4ec9-930d-8313c1d51fb6
# ╠═11e85c84-60ec-4f66-9312-0032492a9773
# ╠═f7fe6bf0-86c5-4eb5-875b-be935719ac18
# ╠═7a8ef3b1-2d49-4bbe-b722-16f9c2f21432
# ╠═30df8e1a-28c3-49ed-8193-d6a594ce90d9
# ╠═a2929307-3f1b-4184-8de4-d00b331daf55
# ╠═e0267c4e-9ea7-4943-9ef0-eeb74018073d
# ╠═2c19c0fe-8cf7-4dfa-a795-22dd56c5b886
# ╠═a239233a-7677-4312-9209-d64570195948
# ╠═dd8312c1-509e-468b-a7d0-234e52ffc8a1
# ╠═dc03903c-4c3d-40a9-93a5-ecced3018aa1
# ╠═d9be6f05-3951-44db-837c-564cd29334d3
# ╠═e7f0bb97-1d37-465a-bf5d-666fdb42ac28
# ╠═2ac10b20-110c-44af-9a17-3df97e62e577
# ╠═cc32bc71-427a-49c4-a1bf-2b4e6dcaa801
# ╠═33819244-5197-4571-865f-7e1bc9d88b98
# ╠═4ea73ac5-3a0f-4a07-80af-447b2c0198d7
# ╠═4231b6a6-a7f4-4960-9d92-707939c49258
# ╠═1f0a45e5-c475-4f6d-8f01-ab512ce333f9
# ╠═44e93e91-1886-4f40-9578-1228dabbfa79
# ╠═f4df8237-23a9-4631-873c-2f537c415cd1
# ╠═aa098647-6390-4858-83f0-3f52e6776420
# ╠═4724d83c-201f-453f-aeac-52202712d97a
# ╠═472637ca-befd-49b6-9c5e-c829e7e84ede
# ╠═1f84ef10-2776-415b-ad60-6c2a5eb5b568
# ╠═fef4d315-5440-40db-9e35-122187b55567
# ╠═a4fcf0b1-9625-4be3-a71b-0a4966a017e4
# ╠═208f6604-701a-47f9-b934-1e4e546ee722
# ╠═15129817-9644-4d23-9c3f-c46fbc6bc267
# ╠═0b5b688b-c508-49e0-bb2a-edb4fe371d8d
# ╠═99fc3172-c63c-4289-92ea-3f7f75e8f2bd
# ╠═ffc9c8b9-7e67-4b62-b87d-b89a47d8ef08
# ╠═c6fe6a20-0e83-4392-8d9e-eb1ae5c47409
# ╠═4e920578-22b6-46fa-94e1-809f80e86425
# ╠═71943ad7-d763-48b8-96f2-dcc91e86a666
# ╠═53e1105a-62de-4fc0-b9eb-60172596ec3f
# ╠═59b79d42-900e-4946-a593-08200503cf34
# ╠═bfbba79b-4a26-4ef5-910f-b001e82ac1e6
# ╠═966f754c-2467-4cda-b1fe-8cb28ecfc0e7
# ╠═dc7b5cbe-8c1f-48f8-960e-82d5be823c28
# ╠═445c9f91-c742-4491-8728-04b7831fdfa5
# ╠═f33ff779-5455-4afd-841d-2eb81238e399
# ╠═7fee8b2c-4cfe-4ae2-94fc-2660de749e74
# ╠═eec10cb1-3880-4a5e-8996-abb2bb48a12f
# ╠═fa8acd76-9abe-4c7f-bc7f-fcb35ea2ebd7
# ╠═e2ec8e48-fe32-4e56-869b-194a5468c9f3
# ╠═590b3b14-0ca0-4a3c-8212-ab2f441b914f
# ╠═3f5c1fa4-8e08-47a7-9cae-ac91a8e0df9e
# ╠═6b2ec37f-cbe5-4a22-91eb-05baf160d452
# ╠═6234576b-51e6-4967-adf1-5859185d3a48
# ╠═0e0b6318-a434-4b11-b719-78bb4a39557f
# ╠═19fd0381-bf9a-4e39-aeee-2187b6275140
# ╠═1e347047-233d-4162-8cf3-3378b2460f7a
# ╠═d19c78c2-7099-405f-94b9-cb091aa63df1
# ╠═c9670300-5cee-42c0-bc65-e5aba726297b
# ╠═9b9b54a3-e4e4-467c-9b01-c280315d0ee4
# ╠═c9e51c52-0766-4d72-93cb-39ce14c7edeb
# ╠═1e80925c-b770-4a6e-8332-3131b1ff004b
# ╠═44638f03-9a8f-45aa-8120-2964e7f728ef
# ╠═fa74a3d2-6f13-4286-9165-b3d34f534d4a
# ╠═d41de4e9-66a7-4e20-9dfd-7ffc1a430300
# ╠═d3fa42c0-deb3-400d-9b62-ca7b32749955
# ╠═7d1fda00-8528-4186-9a8d-80ea6ab9071d
# ╠═b618efad-30b6-42c3-ba78-4add20b925fc
# ╠═7a3dc053-90f5-4229-be06-4b503bb17f16
# ╠═8868caec-6e6f-4c9b-b20c-0ac8db38847c
# ╠═520679b5-ced1-4a9f-bc41-fdf31a7c40ce
# ╠═4538ec89-e65b-4c55-a719-cdd80b44e478
# ╠═c57f83dd-b672-425d-9b32-cc418ea7c4b1
# ╠═434db880-3825-4f3d-b736-af8109bedf4d
# ╠═45452600-600b-4ccc-a6e6-7a92febdd80b
# ╠═fa4a59d9-079a-4feb-a74a-4630f21b6442
# ╠═641eff1b-401c-4dad-88ac-b2aa0a3052dc
# ╠═175db452-1d36-4a3c-9c24-1da14d7d769a
# ╠═a2d1ffca-4981-4db2-a6eb-7b6aad2f2b75
# ╠═a9118a54-a13e-4bc6-9b7d-1f7b9d79fd19
# ╠═070a672b-7955-47c1-b381-b507f44bb149
# ╠═4e4bda01-cdb4-4f44-8827-bb398fb12fbe
# ╠═75a17150-1e2a-4914-8f0f-e8f1c41c0546
# ╠═f7aac374-5012-44c4-9999-ae86b7d4b6c9
# ╠═f4c57c28-dab0-4d63-9d93-972e14b96a80
# ╠═7b2e7e05-9a98-4442-86c1-03621f3bbc77
# ╠═b870c7b6-449b-40b4-badd-338fa4ce8a12
# ╠═77fea351-fe99-4d57-9602-b69ac915afec
# ╠═14e01b5e-aad2-48d2-ac9c-10a8bb595fd6
# ╠═a9c5c21a-69de-4fa6-ae89-33e25b74ce1d
# ╠═b8dffaff-5788-4cb4-bf40-1c07d505fe3f
# ╠═34e72ef5-08b4-4cb5-8164-90ef49bf29cc
# ╠═c0209008-ab91-4c60-9449-d708a97f230c
# ╠═6c43017a-966a-4d43-a5ff-8ef55ce5ccc7
# ╠═c0140abb-4ee8-4f34-b26b-d202030790ac
# ╠═d43486f7-5d63-4a9d-9779-89c1b2469cf7
# ╠═08728379-9a61-4177-8d80-3eff252b8232
# ╠═144dbd76-005e-4a76-ba77-db881e965e04
# ╠═e5e20978-e24a-49e3-bc02-00463ccec77b
# ╠═db5cb755-af4c-4a29-b76d-3cd71284e8a1
# ╠═51b3e6bd-82dc-411f-a6f7-c9b27602f185
# ╠═dad60661-2c76-40d7-997f-fd1139a2aea5
# ╠═4f5f3fe2-08d0-4b0f-b9e2-05171e7a8456
# ╠═e455e16f-e2d4-4984-9dcb-329d02ba5d95
# ╠═84c35c75-3f3a-4979-a9c1-705f5ae0cc8f
# ╠═bbb8b8aa-5fd4-4408-991c-e9f281d57647
# ╠═4aa9e4f8-f49e-4ca9-89d2-d167bfd7d5af
# ╠═140c83dd-7112-4651-aa68-2dd441b239c0
# ╠═fa246234-6f8e-4808-af59-5e5cf163718c
# ╠═26d9bac8-0e64-40c9-ac6f-42df7e2595ba
# ╠═1193b90a-0528-4877-9abb-3670b608be86
# ╠═2dc593ca-1f50-4714-a5c7-7d7de04e01d3
# ╠═a689a276-5ef2-47f1-9577-31a4513cb398
# ╠═2e9f116c-97f0-4d0a-9ac0-054a58eb6a9d
# ╠═e65fb06f-b80b-4c3b-aadc-47cf38a608cd
# ╠═4f3640f8-8224-4acb-a3df-94a71bb62347
# ╠═65389d85-99f1-4682-b35b-8eff20cfee1c
# ╠═c54a479b-97e0-4f39-89a1-74cfb38f659a
# ╠═2109d5d2-a710-4f47-bd01-1c613f16dd83
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

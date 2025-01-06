### A Pluto.jl notebook ###
# v0.20.0

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

# ╔═╡ 8c80a36a-372a-4fad-92bc-bb4426277f55
begin
	using Plots
	using PlutoUI
	using Images
	using Statistics
	using LinearAlgebra
	using DifferentialEquations
	using Distributions
end

# ╔═╡ 44e8209c-8d79-11ef-24db-ff1ea28d01b3
md"""
# darwinian-circuit-quasispecies-model-v0.11.jl 

by Rohan Maddamsetti  

Julia version 1.11.1.  

## Emergence of population-level feedback control by transposon-plasmid coevolution

##### Abstract

The origins of complex adaptive functions from simpler components remains poorly understood. Here, we report the spontaneous evolution of population-level feedback control in clonal bacterial strains, enabling these strains to express pulses of GFP in response to pulses of tetracycline. Balancing selection maintains diverse intracellular plasmids within single cells, due to negative feedback from the toxic effects of a TetA-GFP tetracycline resistance transposon. The intracellular balance of plasmids containing the transposon responds to antibiotic selection. The stability of the balancing selection regime and the magnitude of the GFP response increases with plasmid copy number. In the balancing selection regime, the distribution of plasmids within cells tracks antibiotic concentrations, effectively carrying out population-level Bayesian inference of the environment. The rate at which GFP expression changes in response to antibiotic is determined by the covariance between GFP expression and bacterial fitness. This work shows how mobile genetic elements allow host populations to rapidly evolve the ability to track changes in environmental stressors and describes a fundamental principle for engineering population-level gene expression with intracellular populations of mobile genetic elements. Our findings imply that novel computational abilities can emerge in biological systems, on laboratory timescales, without direct selection for such properties.

"""

# ╔═╡ d4d1f72f-7aa7-46e7-9be0-f38fe991bc0c
md""" ## In this notebook, we test the following claims with our mathematical model:"""

# ╔═╡ b706fac1-3a38-423b-bed3-db6518475c73
md"""
##### Claim (1): The speed at which selection tunes population-level gene expression in response to environmental change is determined by the covariance between gene expression and fitness in the population (Price's theorem).

We examine Claim (1) by randomly sampling 1,000 random initial conditions with random plasmid copy numbers (PCN) and [Tet] concentrations.

"""

# ╔═╡ c3d81737-fc16-44a0-b3de-bab5e0f7aab2
md""" ##### Claim (2) When ecDNA copy number is sufficiently high to maintain diverse ecDNAs in single cells, tunable evolutionary dynamics can emerge in otherwise clonal populations.
"""

# ╔═╡ 3931b5c1-51df-42f5-86e9-09ddba2d2f11
md"""
##### IMPORTANT NOTE: It is not the case that increasing PCN always results in more tunable dynamics. Here is a counterexample (PCN = 40, [Tet] = 5). In this model, [Tet] == the optimal tetA copy number with the highest fitness. Notice that this point can be seen directly in the 'allele frequency' phase diagram at the bottom of this notebook.

We set the initial configuration of the population to 100% "fitness optimum", that is, at TCN == TET_CONC. Then we examine how the distribution evolves.

Set TET_CONC == 5.

Set PCN to 40. The distribution is no longer tunable, due to stochastic loss (too close to bottom absorbing state.)

BUT-- when I set PCN to 5, the population is tunable again! The population goes to the tetA == 6 state, but can still go down to the tetA == 1 state by mutation since the population is large.
"""

# ╔═╡ f8be1737-038a-41c0-9f61-e1980b005ed2
md"""
## DEBUGGING TODO: 

The Price equation prediction is _almost_ correct in the model--there is some apparent numerical error between the Price equation prediction and actual TetA copy number change when the derivatives are very large. Are there any numerical techniques to reduce the remaining numerical error?

"""

# ╔═╡ b1f80124-822d-46e2-9386-54a0117f833d
md"""
## **Model description**

I built a mathematical model to examine how plasmid copy number affects tetA-GFP copy number dynamics in the tetA-transposon system.


##### Quasispecies model formulation.

We use the quasispecies model (see "Unifying Evolutionary Dynamics" by Karen Page and Martin Nowak in Journal of Theoretical Biology, 2002). This model has a physical interpretation as a model of continuous culture in a idealized turbidostat, such that the dilution rate is equal to the bulk growth rate (i.e. mean population fitness). The idea here is that the dilution rate increases as the population adapts, such that the total population density is conserved in the model. Physically, this corresponds to a turbidostat in which the dilution rate instantaneously increases as the population increases in mean growth rate through evolutionary adaptation, such that the (optical) density of the culture stays constant.

**Model Assumptions**  

The population is modeled as a distribution over tetA copy numbers, ranging from 1 (found on the chromosome) to $n+1$, where $n$ is the maximum plasmid copy number. We therefore represent the population as a vector 

$\mathbf{x}(t) = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ ... \\ x_{n+1} \end{pmatrix}$  

The sum of the $n+1$ entries of $\mathbf{x}(t)$ gives the total population size $N(t) = \sum_{i=1}^{n+1} \mathbf{x}_i(t)$.  

"""

# ╔═╡ 749b2bd1-4ccc-48e0-9ab4-bc701500c728
md"""
*Growth dynamics*  

For a given tetracycline concentration, each subpopulation (tetA copy number class) has a growth rate, or * malthusian fitness* $r_i$.  Let $g_i$ be be the log-fitness of tetA copy number class $i$, such that $g_i = log(r_i)$. We assume that there is some optimal tetA copy number given some tetracycline concentration. We assume that $r_i$ is non-negative, with a single peak at the optimum tetA copy number for a given $[Tet]$ concentration, so a natural choice is to define the following quadratic log-fitness function:  

$g_i = g_{max} - \frac{(i - \alpha*[Tet])^2}{2\sigma^2}$  

Here, $g_{max}$ is the maximum log-growth-rate, $i$ is the number of tetA copies in this strain, [Tet] is the antibiotic concentration, $\sigma$ is a free parameter that determines the width of the quadratic function, and $\alpha$ is a unit-conversion factor with units of "tetA gene copies/[Tet]". Without loss of generality we assume $\alpha = 1$, so:

$g_i = g_{max} - \frac{(i - [Tet])^2}{2\sigma^2}$ 

We assume fitness is non-negative, so we pass this $g_i$ through an exponential to get a Gaussian fitness function:

$r_i = e^{g_i}$ 

Then, the growth rates (malthusian fitnesses) for each subpopulation is a vector:

$\mathbf{r} = \begin{pmatrix} r_1 \\ r_2 \\ r_3 \\ ... \\ r_{n+1} \end{pmatrix}$  


We then define the diagonal matrix $\mathbf{R} = \begin{bmatrix}
r_1 & 0 & 0 & ... & 0 \\
0 & r_2 & 0 & ... & 0 \\
0 & 0 & r_3 & ... & 0 \\
... & ... & ... & ... & ... \\
0 & 0 & 0 & ... & r_{n+1} \\
\end{bmatrix}$.  

This represents how each subpopulation of $\mathbf{x}$ grows based on $\mathbf{r}$.

And the average population growth rate (mean population fitness) of the whole population is $\overline{r}(t) = \mathbf{r} \cdot \mathbf{x}(t)$.

"""

# ╔═╡ 137b0ae4-4718-49bc-8e91-60b0e9f90934
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

# ╔═╡ 35b3eeca-74a6-4480-8ee6-0c52d3f31a57
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

# ╔═╡ 53bcbf1c-a2e2-4931-a4f8-88b7c0ea01f4
md"""
*Dilution dynamics*

We assume that cells are diluted out at a rate equal to the bulk population growth rate, which is the average population growth rate (mean population fitness) of the whole population $\overline{r}(t) = \mathbf{r} \cdot \mathbf{x}(t)$.

*Full dynamics*  

The growth and stochastic switching dynamics are combined into a matrix $\mathbf{A} = \mathbf{S}\mathbf{R}$, where $\mathbf{S}$ is a stochastic matrix.

Including dilution, the full dynamics are modeled by the following matrix system of ODEs:  

$\frac{d\mathbf{x}}{dt} = \mathbf{A}\mathbf{x}(t) - \overline{r}(t) \mathbf{x}(t)$
Note that $\overline{r}(t)$ is a scalar and not a vector.


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

# ╔═╡ 4f840ddf-058f-48cc-90ea-609d6cf09ccb
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

# ╔═╡ b8e8c64d-841d-4591-9fce-f2648d5d4f53
md"""

### Notes on continuous-time Price equation.
\
the following comes from the 2023 arxiv paper "The virial theorem and the Price equation" by Steinunn Liorsdottir and Lior Pachter. (we use '$x$' instead of '$p$' to denote subpopulations)

Consider a numerical trait in $n$ subpopulations at time $t$ denoted $\mathbf{z}(t) = (z_1(t), ..., z_n(t))$.
\
\
The subpopulations have
sizes $x_1(t),...,x_n(t)$,
\
\
and have Wrightian fitness $\mathbf{w}(t) = (w_1(t),...,w_n(t))$ defined by $w_i(t) = \frac{x_i(t+∆t)}{x_i(t)}$
where $∆t$ denotes the time interval of one generation.
\
\
Let $q_i(t) = \frac{x_i(t)}{\sum_{j=1}^{n} x_j(t)}$ be the relative size
of the $i$th population.
\
\
Define the population average of fitness to be $\mathbf{\overline{w}}(t) = \sum_{i=1}^{n} q_i(t)w_i(t)$.
\
\
Note that $q(t)$ forms a probability distribution for $w(t)$ viewed as a random variable, and $\mathbb{E}(\mathbf{w}(t)) = \mathbf{\overline{w}}(t)$.
\
\
Let $∆z_i(t) = z_i(t+∆t)−z_i(t)$,
\
\
$∆\mathbf{z}(t) = \mathbf{z}(t+∆t)−\mathbf{z}(t)$,
\
\
and $\mathbf{\overline{z}}(t) = \sum_{i=1}^{n} q_i(t)z_i(t)$,
\
\
with $∆\mathbf{\overline{z}}(t) = \mathbf{\overline{z}}(t+∆t)−\mathbf{\overline{z}}(t)$.

"""

# ╔═╡ 74892533-c56c-4fbe-82d0-3de7127702f0
md""" ### Discrete-time Price equation.

Then, the discrete-time Price equation holds:  

$∆\mathbf{\overline{z}}(t) = \frac{1}{\mathbf{\overline{w}}(t)} cov(\mathbf{w}(t), \mathbf{z}(t)) + \frac{1}{\mathbf{\overline{w}}(t)}\mathbb{E}(\mathbf{w}(t) ⊙ ∆\mathbf{z}(t))$

where $\mathbb{E}(\mathbf{w}(t) ⊙ ∆\mathbf{z}(t))$ is the expected value of the Hadamard product of $\mathbf{w}(t)$ and $∆\mathbf{z}(t)$ with respect to the relative subpopulation sizes, and $cov(\mathbf{w}(t), \mathbf{z}(t)) = E(\mathbf{w} ⊙ \mathbf{z}) − E(\mathbf{w})E(\mathbf{z})$ is the covariance between the subpopulation fitnesses and trait values with respect to the relative subpopulation sizes.

See Liorsdottir and Pachter (2023) for a proof of this fact.

"""

# ╔═╡ 7cc3e6c5-f6f8-42da-8aec-38147950da8a
md""" 
**Remarks:**
\
(from Chapter 6 of Sean Rice's book, Evolutionary Theory)
\
\
The $\frac{1}{\mathbf{\overline{w}}}cov(\mathbf{w}(t),\mathbf{z}(t))$ term represents the change due to differential survival and reproduction, encompassing selection and genetic drift. \
\
The $\frac{1}{\mathbf{\overline{w}}(t)}\mathbb{E}(\mathbf{w}(t) ⊙ ∆\mathbf{z}(t))$ term represents the change due to processes involved in reproduction, such as recombination, regression toward the mean phenotype, or selection at a lower level of organization (i.e. plasmids or other genetic elements that bias their own transmission into daughter cells, potentially at the expense of other, competing genetic elements, as in meiotic drive or CRISPR drives).
"""

# ╔═╡ b9277783-9e31-44e5-82a3-a13af62c62e4
md""" 
### Continuous-time Price equation.
The discrete time Price equation has a continuous time analog. It is formulated using the Malthusian fitness $r(t) = r_1,...,r_n$ given by $r_i(t) = \frac{1}{x_i(t)} \frac{dx_i(t)}{dt} = \frac{d}{dt}ln(x_i(t))$  instead of the Wrightian fitness $w(t)$. Note that this implies that $\frac{dx_i(t)}{dt} = r_i(t)x_i(t)$.
\
\
$\frac{d}{dt}\mathbb{E}(\mathbf{z}(t)) = cov(\mathbf{r}(t), \mathbf{z}(t)) + \mathbb{E}(\frac{d\mathbf{z}(t)}{dt})$.

This is the continuum limit of the discrete-time Price equation; See Liorsdottir and Pachter (2023) for a proof of this fact.
"""

# ╔═╡ b753c51f-6682-47e5-a015-c183e221aa32
md"""
### combining continuous-time Price equation with the quasispecies equation.

By changing the notation for variables in the continuous-time Price equation derived in the Liorsdottir paper, we get equation (4) in the Page and Nowak paper (continuous-time Price equation with different notation). By combining the continuous-time Price equation with the quasispecies equation (equation 1 in the Page and Nowak paper), we get equation 5, which is the invariant we want for this quasispecies model of transposon-plasmid evolution:
\
\
$\frac{d}{dt}\mathbb{E}(\mathbf{z}(t)) = cov(\mathbf{r}(t), \mathbf{z}(t)) + \mathbb{E}(\frac{d\mathbf{z}(t)}{dt}) + \mathbb{E}(\mathbf{r}(t)⊙Δ_m \mathbf{z}(t))$.
\
\
Here, this final term $\mathbb{E}(\mathbf{r}(t) ⊙ Δ_m \mathbf{z}(t))$ describes mutation among types (in our model, due to plasmid segregation during cell division), with $Δ_m \mathbf{z}(t)$ = $\sum_{j} \mathbf{S}_{ji}(\mathbf{z}_j(t) - \mathbf{z}_i(t))$ denoting the expected change in trait value when mutating from type i.
"""

# ╔═╡ 0abeaf9e-5322-4b6d-ba8f-81efd831080b
md"""
##### Key assumption in this model:
$\mathbb{E}({\frac{dz(t)}{dt}}) = \sum_{i} x_i \frac{dz_i(t)}{dt}= 0$

It is easy to see that this condition is always true, because for all times $t$ and trait values $z_i$, $\frac{dz_i(t)}{dt}= 0$, even as $x_i(t)$ changes over time.
\
\
Here is a simple proof: consider arbitrary time $t$ and subpopulation $x_i$ with trait value $z_i$. At time $(t+ \Delta t)$, $z_i(t+ \Delta t) = z_i = i$ because each subpopulation is defined by the number of tetA copies it has, which is $z_i = i$.
\
\
Therefore, $\frac{d}{dt}\mathbb{E}(\mathbf{z}(t)) = cov(\mathbf{r}(t), \mathbf{z}(t)) + \mathbb{E}(\mathbf{r}(t)⊙Δ_m \mathbf{z}(t))$.
"""

# ╔═╡ 08535a54-414b-4a3e-a718-b356cb22cb23
md"""
### Consequences of Price's theorem for the tetA-transposon-plasmid system, and diversity-maintaining genetic elements (DGEs) more generally.


Suppose $Δ_m \mathbf{z}(t)) \approx 0$. This means that offspring do not systematically deviate from their parents in fitness, which implies no transmission bias in the plasmids. This condition holds in this model because the binomial transitions are symmetric, leading to an expected change of zero. The only asymmetric mutational pressure is the rate at which transposons jump from chromosomes to plasmid, which we assume to be very small.

Since $\mathbb{E}(\mathbf{r}(t)⊙Δ_m \mathbf{z}(t)) \approx 0$,
\
\
$\frac{d}{dt}\mathbb{E}(\mathbf{z}(t)) = Cov(\mathbf{r}(t), \mathbf{z}(t))$,
\
\
so $\frac{d}{dt}\mathbb{E}(\mathbf{z}(t)) = \beta_{\mathbf{r},\mathbf{z}} \cdot Var(\mathbf{z})$,
\
\
as the regression coefficient $\beta_{\mathbf{r},\mathbf{z}}$ is defined as $\frac{Cov(\mathbf{r}(t), \mathbf{z}(t))}{Var(\mathbf{z}(t))}$

"""

# ╔═╡ 83685c04-3a04-4d30-90e4-ab9ec7eee180
md"""
Let $\mathbf{z} = TCN$, where $TCN$ is the transposon copy number. Then, \
\
$\frac{d}{dt}\mathbb{E}(TCN(t)) =  \beta_{\mathbf{r},TCN} \cdot Var(TCN)$.
\

This result gives the instantaneous change in TCN, and shows that it depends on the slope of the linear regression of fitness on TCN, and the variance of TCN in the population. This has a nice geometric interpretation.
"""

# ╔═╡ 39b0ddde-3cb4-4630-a4ff-a3ab163c57d8
md""" ### notes on eigenvalues and eigenvectors of the MutSel matrix. """

# ╔═╡ 35adaf46-3b45-4208-aa74-e89a54c4a6d5
md""" See page 35 of Nowak's book, Evolutionary Dynamics to see discussion of eigenanalysis and analytical solution of the quasispecies equation, and **most importantly** the original 1988 quasispecies paper by Eigen, McCaskill, and Schuster. \
\
In addition, note that the Kussell and Leibler (2005) Science paper also does an eigenanalysis of the growth-switching matrix in their paper, and has a time-lag interpretation of their core result on Lyapunov exponents (see endnotes 23, 24, 25 in their paper).

"""

# ╔═╡ 248e7e49-749a-4110-a3df-12d9b873b949
load("../data/Nowak-textbook-page-35.png")

# ╔═╡ 039ffb9c-2c1b-423c-a809-93f9e62bc688
md"""
## Module imports
"""

# ╔═╡ cfd4050d-31a8-49ed-9994-076e55c02ceb
md"""
## Function definitions
"""

# ╔═╡ 934ef964-bdc4-4d5d-bd97-e4e32b4d0380
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

# ╔═╡ 892bb0dd-da94-493f-bf72-cdc53a8b0181
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

# ╔═╡ 412c944d-168d-4146-8bc8-7219db9c291a
function ZeroBasedBinomialSwitchingEntry(a, b, plasmid_copy_num)
	""" Input parameters:
	    a: number of tetA transposons in offspring.
	    b: number of tetA transposons in parent.
	    plasmid_copy_num: total plasmid copy number. """
	p = b/plasmid_copy_num ## probability of sampling a plasmid with the tetA transposon.
	## 'a' is the number of times that the tetA-transposon plasmid picked from parent into offspring.
	return pdf(Binomial(plasmid_copy_num, p), a)
end

# ╔═╡ de13afba-dec7-4a7a-a26b-85e26d36a84f
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

# ╔═╡ 1ab2e8f1-753a-4374-8b5b-1ad3de48a710
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

# ╔═╡ 348e1723-08b6-4968-896c-8459afc8bb08
function normalize_vector!(u)
    total = sum(big.(u))
    if total != 0
        u .= u / total
    end
end

# ╔═╡ f1b278a5-1ce4-47b3-b93e-969e1f9ec6bf
function sample_unit_vector(n)
    ## Generate n random numbers
    x = rand(n)
    ## Normalize the vector to make the sum equal to one
    x /= sum(x)
    return x
end

# ╔═╡ c10d54d1-eb66-4008-8f83-a54342ea8128
function CalcMeanTetACopyNumberVec(sol)
	mean_tetA_copy_num_vec = []
	for i in 1:length(sol.u)
		cur_pop_vec = sol.u[i]
		cur_mean_tetA_copy_num = calc_mean_tetA_copy_number(cur_pop_vec)
		append!(mean_tetA_copy_num_vec, cur_mean_tetA_copy_num)
	end
	return mean_tetA_copy_num_vec
end

# ╔═╡ 66ee7dae-f87d-4b12-91e0-fffdbc420478
function CalcTetACopyNumberVelocity(sol)	
	mean_tetA_copy_num_vec = CalcMeanTetACopyNumberVec(sol)

	## append a zero to the front of the difference vector.
	d_mean_tetA_copy_num_vec = [0; diff(mean_tetA_copy_num_vec)]
	## calculate the dt vec; append a one to the front to avoid division by zero.
	dt_vec = [1; diff(sol.t)]
	## now calculate the derivative vector.
	d_mean_tetA_copy_num_dt_vec = d_mean_tetA_copy_num_vec ./ dt_vec
	return d_mean_tetA_copy_num_dt_vec
end

# ╔═╡ 192ab415-b73a-49d3-ba74-99b5d91ad482
function Entropy(probability_vec)
	vec_sum = sum(probability_vec)
	@assert vec_sum ≈ 1.0 "Error: $probability_vec should sum to one but sums to $vec_sum"
	my_entropy = sum(map(p -> ifelse(p == 0, 0, -p*log(p)), probability_vec))
	return my_entropy
end

# ╔═╡ 0d55f4ed-f01c-43c6-8a58-6edf9b77e9d9
function chromosomal_tetA_initial_pop_vec(my_pcn)
	my_max_tcn = my_pcn + 1
	initial_pop_vec = zeros(BigFloat, my_max_tcn)
	## initialize the population as 100% 1 chromosomal tetA copy.
	initial_pop_vec[1] = big"1.0"
	
	return initial_pop_vec
end

# ╔═╡ f238e7a6-86f1-4898-bc03-3fd80f4e1f3e
function optimal_tetA_initial_pop_vec(my_pcn, my_optimal_tcn)
	my_max_tcn = my_pcn + 1
	initial_pop_vec = zeros(BigFloat, my_max_tcn)
	@assert my_max_tcn >= my_optimal_tcn ## make sure we don't go out of bounds.
	## initialize the population as 100%  my_optimal_tetA_copies.
	initial_pop_vec[my_optimal_tcn] = big"1.0"
	
	return initial_pop_vec
end

# ╔═╡ 8dfdd593-117b-430e-b988-2eaa70484b4c
md"""
### Let's start modeling!
"""

# ╔═╡ 9bd490b7-3f11-487e-b454-785c7d87391b
md""" default parameter settings: \
SIGMA = 10.0 \
R_MAX = 1.0 \
TIMESPAN = 400.0 \
η₀ = 0.001
"""

# ╔═╡ 28bd0eee-f54e-4ab5-aecc-275fe3e8319a
begin
	## Define global constants.
	SIGMA = 10.0 ## set the width of the fitness function
	R_MAX = 1.0 ## set the max growth rate.

	## Define the time span for solving the ODEs
	TIMESPAN = 400.0
	tspan = (0.0, TIMESPAN)
	
	η₀ = 0.001 ## set the transposition rate η₀.

	## set up plot dimensions
	width_to_height_ratio=1.5
	my_plot_height = 3*72
	my_plot_width = width_to_height_ratio * my_plot_height
	my_plot_size = (my_plot_width, my_plot_height)
	
	## Set plot default parameters:
	## set the font to Arial and increase font sizes.
	## Helvetica does not render properly in PDF.
	default(fontfamily="Arial", 
		tickfontsize=12,
		xlabelfontsize=9,
		ylabelfontsize=9,
		legendfontsize=6,
		titlefontsize=10,
		size=my_plot_size,
		grid=false
	)
	
end

# ╔═╡ 9e9cc2c2-55be-456e-95f0-224b08fbeae3
function raw_quadratic_fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	fitness is a function of tetA copies and antibiotic concentration.
	for simplicity, define a quadratic function whose peak is shifted
	left or right based on Tet_conc.
	"""
	fitness = r_max - ((tetA_copy_number - Tet_conc)/2σ)^2
	return fitness
end

# ╔═╡ ef893778-f28f-4ade-ab70-978ded8cfcd5
function gaussian_fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	fitness is a function of tetA copies and antibiotic concentration.
	for simplicity, define a gaussian function whose peak is shifted
	left or right based on Tet_conc.
	"""
	fitness = exp(raw_quadratic_fitness_function(tetA_copy_number, Tet_conc, σ, r_max))
	return fitness
end

# ╔═╡ 38fa13dc-2eca-4ac9-87fd-5d82e2e0ce41
function fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	syntactic sugar: many downstream functions depend on this function,
	so simpler to change the *definition* of this function than to update the *name* of this function across many different blocks of code.
	"""
	return gaussian_fitness_function(tetA_copy_number, Tet_conc, σ, r_max)
end

# ╔═╡ 0a576694-2971-464e-9c91-35e7e3f29fd0
function calc_mean_fitness(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	## get the fitness for each fitness class (defined by tetA copy number).
	fitness_vec = fitness_function.(tetA_classes, Tet_conc)
	## calculate mean fitness in the population,
	## by weighting fitness for each class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_fitness = sum(fitness_vec .* frequency_vec)
	return mean_fitness
end

# ╔═╡ bed5464d-47c8-46c4-9037-e73130a5b0e6
function CalcFitnessVelocity(sol, tet_conc)	
	mean_fitness_vec = [] 
	for i in 1:length(sol.u)
		cur_pop_vec = sol.u[i]
		cur_mean_fitness = calc_mean_fitness(cur_pop_vec, tet_conc)
		append!(mean_fitness_vec, cur_mean_fitness)
	end

	## append a zero to the front of the difference vector.
	d_mean_fitness_vec = [0; diff(mean_fitness_vec)]
	## calculate the dt vec; append a one to the front to avoid division by zero.
	dt_vec = [1; diff(sol.t)]
	## now calculate the derivative vector.
	d_mean_fitness_dt_vec = d_mean_fitness_vec ./ dt_vec
	
	return d_mean_fitness_dt_vec
end

# ╔═╡ bbb2c832-201f-409e-b262-ec7c74b6c490
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

# ╔═╡ b78463f4-1deb-45bd-b091-7ae15a23471b
function calc_tetA_copy_number_fitness_covariance(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))
	## get the fitness for each fitness class (defined by tetA copy number).
	fitness_vec = fitness_function.(tetA_classes, Tet_conc)
	## calculate mean fitness in the population,
	## by weighting fitness for each class by the frequency of each class.
	frequency_vec = pop_vec/sum(pop_vec)
	mean_fitness = sum(fitness_vec .* frequency_vec)
	
	## calculate covariance by multiplying fitness deviations
	## with tetA copy number deviations.
	mean_tetA_copy_number = sum(tetA_classes .* frequency_vec)
	tetA_deviation_vec = [(x - mean_tetA_copy_number) for x in tetA_classes]
	fitness_deviation_vec = [(x - mean_fitness) for x in fitness_vec]
	tetA_fitness_covariance = sum(tetA_deviation_vec .* fitness_deviation_vec .* frequency_vec)
	return tetA_fitness_covariance
end

# ╔═╡ 0b129025-f83c-4560-8b3e-a56356966db4
function CalcTetACopyNumberFitnessCovariance(sol, tet_conc_vec)
	copy_num_covariance_vec = []
	for i in 1:length(sol.u)
		cur_pop_vec = sol.u[i]
		cur_Tet_conc = tet_conc_vec[i]
		cur_copy_num_covariance = calc_tetA_copy_number_fitness_covariance(cur_pop_vec, cur_Tet_conc)	
		append!(copy_num_covariance_vec, cur_copy_num_covariance)
	end
	return copy_num_covariance_vec
end

# ╔═╡ 2f5d53c1-5b8b-40eb-a897-c1efbf1977a9
## Retrieve the current default size
println("Default plot size: ", Plots.default(:size))

# ╔═╡ 6bd86018-8a50-4b8d-a2bd-40bfbe45829b
PCNSlider = @bind PCN Slider(1:100, default=50, show_value=true)

# ╔═╡ 0cf84c01-0378-44be-b878-3c2119b468b9
function BinomialTransitionMatrix(plasmid_copy_num=PCN)
	""" Basic transition matrix modeling plasmid segregation without transposition."""
	
	## max transposon copy number is plasmid copy number + one chromosomal copy.
	max_TCN = plasmid_copy_num + 1

	## initialize as an m x m identity matrix with ones on the diagonal.
	transition_matrix = Matrix{BigFloat}(I, max_TCN, max_TCN)

	## update binomial entries (overwrite ones and zeros as needed).
	for i in 1:max_TCN
		for j in 1:max_TCN
			if (i <= max_TCN) && (j <= max_TCN)
				transition_matrix[i,j] = OneBasedBinomialSwitchingEntry(i,j,plasmid_copy_num)
			end
		end
	end
	
	return(transition_matrix)
end

# ╔═╡ 260df35d-6d1c-42be-8270-01fb62627938
function TransposonSwitchingMatrix(plasmid_copy_num=PCN, η=η₀)
	""" Take the binomial transition matrix, and add 
	transposon jumping from chromosome to plasmid."""
	switching_matrix = BinomialTransitionMatrix(plasmid_copy_num)
	## now add in transpositions from chromosome to plasmid in the zero-plasmid state.
	switching_matrix[1,1] = 1 - η
	switching_matrix[2,1] = η
	return(switching_matrix)
end

# ╔═╡ 1955042b-e29d-4c54-84c3-8d32bed550a3
function calc_Δₘp_vec(pop_vec, Tet_conc)
	""" calculate the expected change in trait value due to mutation (in this case, plasmid segregation)"""

	## length of pop_vec is max_TCN, which is pcn + 1.
	my_pcn = length(pop_vec) - 1
	## switching_matrix has dimensions (my_pcn+1) x (my_pcn+1)
	## NOTE: this code is valid because this notebook assumes the use
	## of the TransposonSwitchingMatrix throughout.
	switching_matrix = TransposonSwitchingMatrix(my_pcn)
	nrow, ncol = size(switching_matrix)
	@assert nrow == ncol == length(pop_vec) ## self-consistency check

	Δₘp_vec = zeros(ncol)
	
	for j in 1:ncol
		for i in 1:nrow
			"""
			IMPORTANT: i,j indices in this code are swapped compared to
			notation in the Page and Nowak (2002) paper.
			In Page and Nowak, i->j represents ancestor -> mutant.
			This code follows standard matrix multiplication notation,
			such that j->i represents ancestor to mutant.
			
			So, the appropriate formula to implement is:
			Δₘpⱼ = ∑ᵢ qᵢⱼ(pᵢ - pⱼ)
			where qᵢⱼ is switching_matrix[i,j]
			and pᵢ is the trait value--
			the number of TetA transposons, or i in this case,
			since there is always a minimum of one copy on the chromosome.
			"""
			Δₘp_vec[j] += switching_matrix[i,j]*(i - j)
		end
	end

	return Δₘp_vec
end

# ╔═╡ 2cb2dae4-33af-4118-87c8-f41c8aba8225
function calc_expected_trait_change_by_mutation(pop_vec, Tet_conc)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:length(pop_vec))

	## get the fitness for each fitness class (defined by tetA copy number).
	## use BigFloats to get better numerical precision.
	growth_rate_vec = big.(fitness_function.(tetA_classes, Tet_conc))

	frequency_vec = pop_vec/sum(pop_vec)

	Δₘp_vec = calc_Δₘp_vec(pop_vec, Tet_conc)

	## see the formula given in Page and Nowak (2002).
	expected_trait_change_by_mutation = sum(frequency_vec .* growth_rate_vec .* Δₘp_vec)

	return expected_trait_change_by_mutation
end

# ╔═╡ c99d3673-9f9a-4508-b257-360cc6d0733d
function calc_delta_expected_trait_change(pop_vec, Tet_conc)
	""" calculate the time derivative of the expected change in trait value (Equation 5 in Page and Nowak 2002)."""
	covariance_term = calc_tetA_copy_number_fitness_covariance(pop_vec, Tet_conc)

	""" IMPORTANT: we assume that the expected derivative of the trait value vector == 0,
	 such that we can ignore this term in the RHS of the Price equation.
	 It is easy to see that this is true in this model, because the number of
	 tetA copies is fixed per subpopulation (any change is equivalent to a
	 mutation transition to a different subpopulation).
	 This condition may not necessarily hold for an arbitrary trait,
	 say if the environment changes, causing the trait to change (plasticity)
	 or if the trait is frequency-dependent 
	(for instance, frequency-dependent fitness). See Page and Nowak (2002).
	"""
	
	mut_term = calc_expected_trait_change_by_mutation(pop_vec, Tet_conc)

	return (covariance_term + mut_term)
end

# ╔═╡ 36bab24f-30a2-489f-94e0-fa5962a49b37
function CalcPriceEquationLHS_vec(sol, tet_conc_vec)
	Price_equation_LHS_vec = []
	for i in 1:length(sol.u)
		cur_pop_vec = sol.u[i]
		cur_tet_conc = tet_conc_vec[i]
		cur_Price_LHS = calc_delta_expected_trait_change(cur_pop_vec, cur_tet_conc)
		append!(Price_equation_LHS_vec, cur_Price_LHS)
	end
	return Price_equation_LHS_vec
end

# ╔═╡ a506ff86-41c6-44ac-adf5-c3fdb368cb02
TetConcSlider = @bind TET_CONC Slider(0:50, default=25, show_value=true)

# ╔═╡ 45362f0d-76ac-44e8-b6cd-0d1824b3a3b4
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

# ╔═╡ 49e1dd14-e839-45cf-8e0a-f0badd573d67
function MutSelMatrix(plasmid_copy_num=PCN, Tet_conc=TET_CONC)
	mutsel_matrix = TransposonSwitchingMatrix(plasmid_copy_num) * SelectionDiagonalMatrix(plasmid_copy_num, Tet_conc)
	return(mutsel_matrix)
end

# ╔═╡ ea1bdd07-3d5a-49ef-a626-95d9b264d297
function quasispecies_odefunc(du, u, p, t; ϵ=1e-7)
	## ϵ is a (1/N) threshold to model bottlenecks in transfers/finite population sizes in experiments.

	## Unpack the parameters in p
	pcn, tet_conc_function = p

	# Check if tet_conc_function is indeed callable
    if !isa(tet_conc_function, Function)
        error("The second element of p must be a TetFunction")
    end

	## get the current [Tet] concentration
	cur_tet_conc = tet_conc_function(t)
	
	## Define the ODE system
	A = MutSelMatrix(pcn, cur_tet_conc)

	## Enforce positivity constraint
    u .= max.(u, 0.0)

	## set subpopulations with population size < ϵ to zero, to model extinction.
	u[u .< ϵ] .= 0

	## Normalize the vector to ensure it sums to one
	normalize_vector!(u)
	
	## subtract the mean population growth during this dt interval.
	Au = A*u ## sugar to avoid recomputation
	du .= Au - sum(Au) * u
end

# ╔═╡ 638d7f76-512c-4c11-951d-76a77341b59d
function SolveTetQuasispeciesSystem(my_initial_pop_vec, my_pcn, my_tet_function;
	my_tspan=tspan)
	""" 
	my_tet_function is a function from time to tet_conc.
	"""
	## max transposon copy number is tied to plasmid_copy_number
	max_tcn = my_pcn + 1
	
	## Create an ODEProblem
	prob = ODEProblem(quasispecies_odefunc, my_initial_pop_vec, my_tspan, (my_pcn, my_tet_function))
	## Solve the ODE system
	sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-9)
	return sol
end

# ╔═╡ 621be509-5e7e-415a-bcc6-daa7b1903834
function CalcNormalizedTopEigenvector(pcn, tet_conc)
	## max transposon copy number is tied to plasmid_copy_number
	max_tcn = pcn + 1
	
	A = MutSelMatrix(pcn, tet_conc)

	my_top_eigenvector = [real(x) for x in eigen(A).vectors[:,max_tcn]]
	## normalize eigenvector to sum to unity
	normalize_vector!(my_top_eigenvector)
	return my_top_eigenvector	
end

# ╔═╡ ea419d1b-1246-4999-a8a4-706c73e14ed8
function TetPulseFunction(t;tet_conc=TET_CONC)
	## We model [Tet] pulses over time by dividing the current time by 100, and
	## setting "[Tet] ON" if in the first half.
	pulse_cycle_time = 100
	return mod(t,pulse_cycle_time) < (pulse_cycle_time/2) ? tet_conc : 0
end

# ╔═╡ c899121e-2b9e-4519-b168-98994684a86d
function CreateTetPulseFunction(tet_conc=TET_CONC)
	## We model [Tet] pulses over time by dividing the current time by 100, and
	## setting "[Tet] ON" if in the first half.
	pulse_cycle_time = 100
	return t -> mod(t,pulse_cycle_time) < (pulse_cycle_time/2) ? tet_conc : 0
end

# ╔═╡ bd19d151-c3fc-4ff0-b2c5-4ffe07adcae0
function SolvePulseTetQuasispeciesSystem(my_initial_pop_vec, my_pcn, my_tet_conc;
	my_tspan=tspan)
	## provide pulses of [Tet].
	my_tet_function = CreateTetPulseFunction(my_tet_conc)
	return SolveTetQuasispeciesSystem(my_initial_pop_vec, my_pcn, my_tet_function; my_tspan)
end


# ╔═╡ ac741b84-00f6-4f70-89a4-07dbb24ec9cb
function ConstantTetFunction(t;tet_conc=TET_CONC)
	## The purpose of this function is to provide a common interface for downstream code that may use either TetPulseFunction(t) or ConstantTetFunction(t).
	return tet_conc
end

# ╔═╡ 5bd6b242-9a1e-44ea-b319-6b7604e57624
function CreateConstantTetFunction(tet_conc=TET_CONC)
	return t -> tet_conc
end

# ╔═╡ 162efa78-e3a7-4ebb-ad7a-395a567071b1
function SolveConstantTetQuasispeciesSystem(my_initial_pop_vec, my_pcn, my_tet_conc;
	my_tspan=tspan)
	## provide a constant [Tet] concentration.
	my_tet_function = CreateConstantTetFunction(my_tet_conc)
	return SolveTetQuasispeciesSystem(my_initial_pop_vec, my_pcn, my_tet_function; my_tspan)
end

# ╔═╡ 9241c94d-857c-4a31-93f1-cfed036cf2f2
md""" ##### sometimes the PCN and TET\_CONC variables do not update properly when the slider is changed-- double check the actual values of the PCN and TET_CONC variables..."""

# ╔═╡ 09c74702-0d19-4049-8f8e-6f91a0163dcc
begin
	MAX_TCN = PCN + 1 ## max transposon copy number is tied to plasmid_copy_number
	xvec = collect(1:MAX_TCN) ## tetA copy number classes for plot x-axes
	println("PCN is $PCN, [Tet] is $TET_CONC")
end

# ╔═╡ eb9ea298-95f2-4df7-966a-069ca8965c67
InitialCloneSlider = @bind INITIAL_CLONE_TCN Slider(1:MAX_TCN, default=1, show_value=true)

# ╔═╡ e1cd955e-a069-426b-b222-c60790c809b3
INITIAL_CLONE_TCN

# ╔═╡ c40e9cfa-f58f-459c-a994-34bab25c20ad
let
	fitnesses = fitness_function.(xvec, TET_CONC)
	plot(xvec, fitnesses, label="Fitness", xlabel="tetA copy number", ylabel="Fitness")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="[Tet]")
end

# ╔═╡ 6fbce7a0-af24-4501-92e4-fbd6b5d27b65
md""" #### Set the initial population state. """

# ╔═╡ 4fd0feac-346b-42bc-b0de-1bf0a4718928
begin
	## parameters used across time course simulations.
	initial_pop_vec = zeros(BigFloat, MAX_TCN)

	## initialize the population as a uniform distribution. 
	##initial_pop_vec .= big(1/MAX_TCN)
	
	## initialize the population as 100% INITIAL_CLONE_TCN tetA copies.
	initial_pop_vec[INITIAL_CLONE_TCN] = big"1.0"
	
	## to show how increasing PCN increases the stability of the optimal state,
	## initialize the population at the highest fitness state (TCN == TET_CONC). 
	#initial_pop_vec[TET_CONC] = big"1.0"	
end

# ╔═╡ 96241706-df79-4c86-838d-398302160e3d
md""" ## model population dynamics under constant [Tet]."""

# ╔═╡ 1d8f69c1-e2b2-492b-86a5-a9947fe9ac8a
begin
	## run the simulation for constant [Tet] concentration.
	sol = SolveConstantTetQuasispeciesSystem(initial_pop_vec, PCN, TET_CONC)

	## get the final stationary distribution in the constant [Tet] population.
	final_const_Tet_population = sol.u[end]

	## make a vector of tet concentration over time in the model solution
	## for calculating summary statistics.
	const_tet_vec = [ConstantTetFunction(t) for t in sol.t]

	## calculate mean copy number.
	mean_copy_num_vec = CalcMeanTetACopyNumberVec(sol)

	## calculate rate of copy number change.
	d_mean_copy_num_vec_dt_vec = CalcTetACopyNumberVelocity(sol)

	## calculate tetA copy number--fitness covariance.
	copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(sol, const_tet_vec)

	## calculate Price equation LHS (covariance + mutation terms).
	constant_tet_pop_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(sol, const_tet_vec)
end

# ╔═╡ 2ce5fa35-27ae-4382-a41d-fa835052412d
time_step_slider = @bind cur_timestep Slider(1:length(sol.u), default=1, show_value=true)

# ╔═╡ 6e3e1eec-03c2-4bb3-bb5a-9c27182753c7
let
	bar(xvec, sol.u[cur_timestep], label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="[Tet]")
end

# ╔═╡ a25df7df-eb57-48a2-a538-44f7b7432387
md""" ##### check whether the final stationary distribution matches the eigenvector corresponding to the dominant eigenvalue."""

# ╔═╡ af8a47bd-e823-41ae-9813-0b484f040e4f
begin
	Amatrix = MutSelMatrix(PCN, TET_CONC)
	Amatrix_top_eigenvector = [real(x) for x in eigen(Amatrix).vectors[:,MAX_TCN]]
	## normalize the stationary population distribution.
	normalized_Amatrix_top_eigenvector = Amatrix_top_eigenvector/sum(Amatrix_top_eigenvector)

	## Invert the Amatrix
	Amatrix_inv = inv(Amatrix)
	# and transform eigenvector back to the original basis
	original_basis_vector = Amatrix_inv * Amatrix_top_eigenvector
	## and normalize the stationary population distribution.
	normalized_original_basis_vector = original_basis_vector/sum(original_basis_vector)
end

# ╔═╡ d03e6d49-2c29-4101-9918-b8917fb037d8
let
	bar(xvec, final_const_Tet_population, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="[Tet]")
end

# ╔═╡ bc4a7835-dd9f-42e4-b61e-118a33c7ff96
let
	bar(xvec, normalized_Amatrix_top_eigenvector, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="[Tet]")
end

# ╔═╡ 9028f452-a9cd-42ad-94ff-f86e54462a7c
let
	bar(xvec, normalized_original_basis_vector, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="[Tet]")
end

# ╔═╡ e4e8c61d-5503-490f-993d-bb1464def4c9
md""" ##### plot mean copy number in the constant [Tet] population."""

# ╔═╡ d0c028ae-b991-4558-87e6-41b4f5bfe274
let
	plot(sol.t, mean_copy_num_vec, label="mean tetA copy number", xlabel="time", ylabel="mean tetA copy number")
end

# ╔═╡ 7bfbb45d-57fc-4c9e-b929-ac44f7f88fef
md""" ##### plot tetA copy number--fitness covariance in the constant [Tet] population. """ 

# ╔═╡ 73873c97-daf6-4133-a657-90de30803517
let
	plot(sol.t, copy_num_covariance_vec, label="tetA-fitness covariance")
	plot!(sol.t, d_mean_copy_num_vec_dt_vec, label="rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative and covariance")
end

# ╔═╡ 91d1760c-bbaf-466b-8eb7-623cb0dcd686
md""" ##### plot LHS of quasispecies continuous-time Price equation in the constant [Tet] population. """ 

# ╔═╡ 9b92d86b-9f19-4fed-be2b-691d58ffaab6
let
	plot(sol.t, constant_tet_pop_Price_equation_LHS_vec, label="Price equation LHS from Page and Nowak 2002")
	plot!(sol.t, d_mean_copy_num_vec_dt_vec, label="rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative and covariance")
end

# ╔═╡ 49830d58-03a3-4379-aeea-767c9a3eeb26
let
	plot(sol.t, (constant_tet_pop_Price_equation_LHS_vec - d_mean_copy_num_vec_dt_vec), label="error check",xlabel="time", yaxis="error")
end

# ╔═╡ 1cff27ae-4465-4c38-8b26-5cace8a833d7
constant_tet_pop_Price_equation_LHS_vec - copy_num_covariance_vec

# ╔═╡ e1046a07-ad24-45e0-b871-069640aa6df3
md""" ### Here, we show that the response time depends on the linear regression of population fitness on the landscape. In this model, the population that is best able to "hit the peak", has the fastest response time. This is a key theoretical prediction of the Price equation. 

"""

# ╔═╡ 4bddf691-83a7-4101-bb00-0fe463e5de77
md""" #### Vary PCN (5, 15, 25, 50), keep [Tet] == 10, and compare rate of changes of mean tetA copy number"""

# ╔═╡ f96766c2-bd3b-4782-a884-a2d72667d0e4
begin
	pcn5_tet10_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(5), 5, 10)
	pcn5_tet10_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn5_tet10_sol, [10 for t in pcn5_tet10_sol.t])
	pcn5_tet10_copy_num_velocity = CalcTetACopyNumberVelocity(pcn5_tet10_sol)

	
	pcn15_tet10_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(15), 15, 10)
	pcn15_tet10_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn15_tet10_sol, [10 for t in pcn15_tet10_sol.t])
	pcn15_tet10_copy_num_velocity = CalcTetACopyNumberVelocity(pcn15_tet10_sol)

	
	pcn25_tet10_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(25), 25, 10)
	pcn25_tet10_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn25_tet10_sol, [10 for t in pcn25_tet10_sol.t])
	pcn25_tet10_copy_num_velocity = CalcTetACopyNumberVelocity(pcn25_tet10_sol)

	
	pcn50_tet10_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(50), 50, 10)
	pcn50_tet10_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn50_tet10_sol, [10 for t in pcn50_tet10_sol.t])
	pcn50_tet10_copy_num_velocity = CalcTetACopyNumberVelocity(pcn50_tet10_sol)
end

# ╔═╡ d153fde0-3435-49f2-b493-2c2fc8fec41f
let
	plot(pcn5_tet10_sol.t, pcn5_tet10_copy_num_covariance, label="PCN=5, TET=10", xlabel="time", ylabel="tetA copy number-fitness covariance", ylabelfontsize=7)
	
	plot!(pcn15_tet10_sol.t, pcn15_tet10_copy_num_covariance, label="PCN=15, TET=10")
	plot!(pcn25_tet10_sol.t, pcn25_tet10_copy_num_covariance, label="PCN=25, TET=10")
	plot!(pcn50_tet10_sol.t, pcn50_tet10_copy_num_covariance, label="PCN=50, TET=10")
end

# ╔═╡ 11e9cbd5-19f2-4401-ba19-8590d3f4dca7
let
	tet10_response_plot = plot(pcn5_tet10_sol.t, pcn5_tet10_copy_num_velocity, label="PCN=5, TET=10", xlabel="time", ylabel="change in mean tetA copy number", ylabelfontsize=7)

	plot!(pcn15_tet10_sol.t, pcn15_tet10_copy_num_velocity, label="PCN=15, TET=10")
	plot!(pcn25_tet10_sol.t, pcn25_tet10_copy_num_velocity, label="PCN=25, TET=10")
	plot!(pcn50_tet10_sol.t, pcn50_tet10_copy_num_velocity, label="PCN=50, TET=10")
	
	savefig(tet10_response_plot, 
		"../results/modeling-results/tet10_response_plot.pdf")
	tet10_response_plot
end

# ╔═╡ 9a48f3c4-a5d3-4d53-8248-557c4dd86efb
md""" #### Vary PCN (5, 15, 25, 40, 50), keep [Tet] == 40, and compare rate of changes of mean tetA copy number"""

# ╔═╡ cb0d4f61-2a7d-4fe2-b2a7-8b2aafae22bd
begin
	pcn5_tet40_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(5), 5, 40)
	pcn5_tet40_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn5_tet40_sol, [40 for t in pcn5_tet40_sol.t])
	pcn5_tet40_copy_num_velocity = CalcTetACopyNumberVelocity(pcn5_tet40_sol)

	
	pcn15_tet40_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(15), 15, 40)
	pcn15_tet40_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn15_tet40_sol, [40 for t in pcn15_tet40_sol.t])
	pcn15_tet40_copy_num_velocity = CalcTetACopyNumberVelocity(pcn15_tet40_sol)

	
	pcn25_tet40_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(25), 25, 40)
	pcn25_tet40_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn25_tet40_sol, [40 for t in pcn25_tet40_sol.t])
	pcn25_tet40_copy_num_velocity = CalcTetACopyNumberVelocity(pcn25_tet40_sol)

	pcn40_tet40_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(40), 40, 40)
	pcn40_tet40_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn40_tet40_sol, [40 for t in pcn40_tet40_sol.t])
	pcn40_tet40_copy_num_velocity = CalcTetACopyNumberVelocity(pcn40_tet40_sol)
	
	pcn100_tet40_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(100), 100, 40)
	pcn100_tet40_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(pcn100_tet40_sol, [40 for t in pcn100_tet40_sol.t])
	pcn100_tet40_copy_num_velocity = CalcTetACopyNumberVelocity(pcn100_tet40_sol)
end

# ╔═╡ 86e16ddc-a9e0-4e2c-822c-abc9f19d0ccf
let
	plot(pcn5_tet40_sol.t, pcn5_tet40_copy_num_covariance, label="PCN=5, TET=40", xlabel="time", ylabel="tetA copy number-fitness covariance", ylabelfontsize=7)

	plot!(pcn15_tet40_sol.t, pcn15_tet40_copy_num_covariance, label="PCN=15, TET=40")
	plot!(pcn25_tet40_sol.t, pcn25_tet40_copy_num_covariance, label="PCN=25, TET=40")
	plot!(pcn100_tet40_sol.t, pcn100_tet40_copy_num_covariance, label="PCN=100,TET=40")
end

# ╔═╡ 933c9aef-3f09-42aa-b7c5-8f7e5c63d1c9
let
	tet40_response_plot = plot(pcn5_tet40_sol.t, pcn5_tet40_copy_num_velocity, label="PCN=5, TET=40", xlabel="time", ylabel="change of mean tetA copy number", ylabelfontsize=7)

	plot!(pcn15_tet40_sol.t, pcn15_tet40_copy_num_velocity, label="PCN=15, TET=40")
	plot!(pcn25_tet40_sol.t, pcn25_tet40_copy_num_velocity, label="PCN=25, TET=40")
	plot!(pcn100_tet40_sol.t, pcn100_tet40_copy_num_velocity, label="PCN=100, TET=40")

	savefig(tet40_response_plot, "../results/modeling-results/tet40_response_plot.pdf")
	tet40_response_plot
end

# ╔═╡ a2603db4-3fd8-4cc7-a943-8fa142961a61
md""" ## model the [Tet] pulse conditions of the Darwin experiment."""

# ╔═╡ 55f2289b-9268-4bca-9f0f-7f27e29697f8
begin 
	## run the simulation for pulses of [Tet] concentration.	
	pulse_sol = SolvePulseTetQuasispeciesSystem(initial_pop_vec, PCN, TET_CONC)
	
	## make a vector of tet concentration over time in the model solution
	## for calculating summary statistics.
	tet_pulse_vec = [TetPulseFunction(t) for t in pulse_sol.t]

	## calculate mean copy number.
	pulse_mean_copy_num_vec = CalcMeanTetACopyNumberVec(pulse_sol)
	
	## calculate rate of copy number change.
	d_pulse_mean_copy_num_vec = CalcTetACopyNumberVelocity(pulse_sol)

	## calculate tetA copy number--fitness covariance. 
	pulse_copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(pulse_sol, tet_pulse_vec)

	## calculate Price equation LHS (covariance + mutation terms).
	pulse_tet_pop_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(pulse_sol, tet_pulse_vec)
end

# ╔═╡ c83fbfaa-1f5f-40b9-a6c1-68d480c4dfe7
pulse_time_step_slider = @bind pulse_cur_timestep Slider(1:length(pulse_sol.u), default=1, show_value=true)

# ╔═╡ 348fabc0-3caa-4585-975e-688a26b7fa8a
let
	bar(xvec, pulse_sol.u[pulse_cur_timestep], label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ 4acf6a58-6ca5-43e2-8169-75f169cae851
md""" ### plot the antibiotic pulse regime over time."""

# ╔═╡ e22f600c-4468-47ff-93f5-ad004ff5ea23
plot(pulse_sol.t, tet_pulse_vec,label="Tet pulse regime", xlabel="time", ylabel="[Tet] concentration")

# ╔═╡ d20f047d-68d3-44eb-a258-3642d1197a88
md""" ##### plot mean copy number in the population."""

# ╔═╡ 362c0eff-b354-481d-bc5f-673613e8b0a6
let
	plot(pulse_sol.t, pulse_mean_copy_num_vec, label="mean tetA copy number", xlabel="time", ylabel="Mean tetA copy number")
end

# ╔═╡ 5f571b12-244c-4e7e-8774-883f0427ff06
let
	plot(pulse_sol.t, pulse_copy_num_covariance_vec, label="tetA copy number fitness covariance",alpha=0.5)
	plot!(pulse_sol.t, d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:right,alpha=0.5)
end

# ╔═╡ ab3d5352-abfa-4de3-8421-eb86472b91ea
let
	plot(pulse_sol.t, (d_pulse_mean_copy_num_vec - pulse_copy_num_covariance_vec), label="error check", xlabel="time", ylabel="error")
end

# ╔═╡ c9c7c6e9-75b3-4900-b146-010dd37f4123
let
	scatter(pulse_copy_num_covariance_vec, d_pulse_mean_copy_num_vec, xlabel="tetA copy number covariance with fitness", ylabel="tetA copy number velocity",label="")
end

# ╔═╡ cef4643c-fd32-4a75-b9ab-0276341fc0bb
let
	scatter(pulse_tet_pop_Price_equation_LHS_vec, d_pulse_mean_copy_num_vec, xlabel="Price equation LHS", ylabel="tetA copy number velocity",label="")
end

# ╔═╡ 76311c6d-8e96-4265-8bcf-fbc80a6870f3
let
	scatter(pulse_copy_num_covariance_vec, pulse_tet_pop_Price_equation_LHS_vec, xlabel="tetA copy number covariance with fitness", ylabel="Price equation LHS",label="")
end

# ╔═╡ bc0d22dc-b1d2-48d8-aa9d-b0210ab3f9a7
pulse_tet_pop_Price_equation_LHS_vec - pulse_copy_num_covariance_vec

# ╔═╡ e9d58996-30f0-43b2-b13b-9cf7edfcd9eb
numerical_comparison_matrix = hcat(d_pulse_mean_copy_num_vec, pulse_tet_pop_Price_equation_LHS_vec)


# ╔═╡ 026d3047-cca8-462f-8df3-b73bdeb4af31
for row in eachrow(numerical_comparison_matrix)
    println(join(row, " "))
	println()
end

# ╔═╡ e15bbef0-42f1-47d5-af03-ce95468cba93
md""" #### Let's test Claim 1: the speed at which selection tunes population-level gene expression in response to environmental change is determined by the covariance between gene expression and fitness in the population (Price's theorem).

We examine Claim (1) by randomly sampling 1,000 random initial conditions with random plasmid copy numbers (PCN) and [Tet] concentrations.
"""

# ╔═╡ 153e0346-5632-4b96-929c-07cac8f8cfb6
md"""

## TODO: run for N = 1000, then save the results to file. Then if these results exist on disk, load them. If they don't exist on disk, then run this (slow) code to create them. Move the code for making the N = 1000 file into a separate piece of Julia code for speed so that this notebook is not slowed down.
"""

# ╔═╡ bc46fd90-a4c5-4ec0-80f6-c96c6901ac6d
md"""
## DEBUGGING TODO: 

The Price equation prediction is _almost_ correct here, it is probably good enough to publish-- but there is some apparent numerical error between the Price equation prediction and actual TetA copy number change when the magnitude of the derivatives are very large. Are there any numerical techniques to reduce the remaining numerical error?
"""

# ╔═╡ f3109642-2e16-4e11-8b1c-19667a380321
md""" I asked ChatGPT the following question:
\
what are some tricks for reducing numerical error in julia in estimating derivatives of a time series of a quantity of interest, where the dynamics are generated by solving an ODEProblem in Julia?
"""

# ╔═╡ a24529bf-7679-45b9-ac84-a84169ec7c7e
md""" And I got the following answer:
\
\
Reducing numerical error in estimating derivatives of a time series in Julia, especially when solving an ODE problem, involves careful consideration of both numerical integration and derivative estimation. Here are some tips and tricks:
\
1. Choose an Appropriate Solver and Tolerance
ODE Solver Choice: Use a high-accuracy ODE solver suitable for your problem, such as Tsit5 for non-stiff problems or RODAS5 for stiff problems.
Adjust Tolerances: Set abstol and reltol in the ODE solver to sufficiently low values (e.g., 1e-8 or smaller). Lower tolerances reduce integration error but may increase computation time.
julia
Copy code
sol = solve(prob, Tsit5(), abstol=1e-10, reltol=1e-10)
\
2. Ensure Smoothness of the Solution
Dense Output: Use the solver’s interpolation feature (e.g., sol(t)) to get a smooth representation of the solution at any point, rather than relying on discrete time points.
Avoid Overly Sparse Time Points: If you output too few points, interpolation errors can dominate.
\
3. Use Higher-Order Differencing for Derivatives
Numerical differentiation amplifies noise. Use central difference formulas with higher-order accuracy to reduce truncation error. For example:
julia
Copy code
function central_difference(y, t, dt)
    return (y(t + dt) - y(t - dt)) / (2 * dt)
end
Choose dt small enough to minimize truncation error without amplifying round-off error.
\
4. Symbolic or Automatic Differentiation
If the equations of motion are available symbolically or programmatically, use automatic differentiation (via tools like ForwardDiff.jl or Zygote.jl) to compute derivatives directly, avoiding finite difference errors.
\
5. Post-Processing with Smoothing
Smoothing Splines: Fit a smoothing spline to your time series and differentiate the spline to estimate derivatives.
julia
Copy code
using Dierckx
spline = Spline1D(t, y, k=4, s=0.01)  # Adjust smoothing parameter `s`
dy_dt = spline'(t)
Filter Noise: Apply a low-pass filter to the time series before differentiation to reduce high-frequency noise.
\
6. Normalize and Scale Variables
Scale the state variables and time to reduce the magnitude of values. This can mitigate round-off errors during computation.
\
7. Adaptive Time-Stepping
Let the ODE solver determine the step size dynamically to balance error across the integration. Avoid manually fixing the step size unless absolutely necessary.
\
8. Estimate Derivatives During Integration
Some solvers (like DifferentialEquations.jl) allow the computation of derivative terms (e.g., state velocities) directly during integration, reducing the need for post hoc numerical differentiation.
\
9. Validate with Analytic Solutions
If possible, validate your approach using a problem with a known analytical solution. This helps you tune parameters and understand the error characteristics.
\
10. Refine at Key Points
If specific time intervals are of interest, refine the solver’s output near these points using saveat or dense output for higher resolution.
By carefully combining these strategies, you can minimize numerical errors in both solving the ODE and estimating derivatives from the resulting time series.
"""

# ╔═╡ 38cd0bb1-137f-401b-8279-cc3a28ead6f3
begin
	
	sampled_solution_covariance_velocity_tuples = []

	## sample N times.
	N = 10
	## N = 100
	for i in 1:N

		## sample a random [Tet] concentration between 0 and 50
		my_Tet = rand(0:50)
		
		## sample a random plasmid copy number between 1 and 100
		## We do this to avoid the following error:
		## "Binomial: the condition zero(p) <= p <= one(p) is not satisfied."
		my_PCN = rand(1:100)
		## then the initial condition has my_PCN + 1 entries.
		my_max_TCN = my_PCN + 1
		println("PCN is $my_PCN and [Tet] is $my_Tet.")
		
		## using a random initial condition that sums to one.
		my_random_initial_vec = big.(sample_unit_vector(my_max_TCN))
		my_sol = SolveConstantTetQuasispeciesSystem(my_random_initial_vec, my_PCN, my_Tet)

		## using an initial condition of one tetA copy on the chromosome.
		##my_sol = SolveConstantTetQuasispeciesSystem(chromosomal_tetA_initial_pop_vec(my_PCN), my_PCN, my_Tet)
		
		## calculate a vector of constant [Tet] over time for covariance
		## calculations
		my_tet_vec = [my_Tet for t in my_sol.t]
		
		my_copy_num_covariance = CalcTetACopyNumberFitnessCovariance(my_sol, my_tet_vec)
		my_copy_num_velocity = CalcTetACopyNumberVelocity(my_sol)

		## append the (sol, covariance, copy_num_velocity) result
		my_result_tuple = (my_sol, my_copy_num_covariance, my_copy_num_velocity)
		push!(sampled_solution_covariance_velocity_tuples, my_result_tuple)
	end

end

# ╔═╡ ab7f4870-6967-4910-b272-ed75fb8fbb54
md""" ### let's look at the first random simulation as a sanity check."""

# ╔═╡ 506fb600-7043-4562-8470-ffd2fd48d073
begin
	test_result_tuple = sampled_solution_covariance_velocity_tuples[1]
	test_sol, test_copy_num_covariance, test_copy_num_velocity = test_result_tuple
end

# ╔═╡ 47445b3c-bf48-485a-b291-79437595b45a
let
	plot(label="")
	plot!(test_sol.t, test_copy_num_covariance, xlabel="time", ylabel="tetA copy number covariance")
	plot!(test_sol.t, test_copy_num_velocity, xlabel="time", ylabel="tetA copy number velocity")
end

# ╔═╡ adfa712e-80ca-405e-9854-1d0024c3d6f1
md""" 
##### NOTE: right now the first (velocity, covariance) point is off, when starting from a random initial vector, because the first value of the difference vector is always zero. Perhaps a better numerical estimate can fix this? Otherwise, OK to skip this value. 
""" 

# ╔═╡ 5f9855dd-004f-4dc7-9673-ca5f81fb441b
let
	scatter(legend=false) # Initialize an empty plot
	scatter!(test_copy_num_velocity[1:end], test_copy_num_covariance[1:end], xlabel="tetA copy number velocity ", ylabel="tetA copy number-fitness covariance", aspect_ratio=1)
end

# ╔═╡ b1697098-07e0-4577-b144-9ecc3879422b
let
	scatter(legend=false) # Initialize an empty plot
	scatter!(test_copy_num_velocity[2:end], test_copy_num_covariance[2:end], xlabel="tetA copy number velocity ", ylabel="tetA copy number-fitness covariance", aspect_ratio=1)
end

# ╔═╡ d530e10c-63ff-49bb-a7b1-ea4b363738d0
begin
	my_big_scatterplot = scatter(legend=false)  # Initialize an empty plot
	
	for my_result_tuple in sampled_solution_covariance_velocity_tuples
		
		my_sol, my_copy_num_covariance, my_copy_num_velocity = my_result_tuple
			
		scatter!(my_copy_num_covariance[2:end], my_copy_num_velocity[2:end], xlabel="tetA copy number-fitness covariance", ylabel="tetA copy number velocity",  aspect_ratio=1,markersize=3)
	end
	
end

# ╔═╡ 1f030b53-5503-490e-abf0-20101555d3b1
my_big_scatterplot

# ╔═╡ 798394c1-79a3-4886-b42d-c6f20c73c572
savefig(my_big_scatterplot, "../results/modeling-results/velocity-covariance-scatterplot.pdf")

# ╔═╡ 138d4602-86bf-4923-b648-b6d273f6adac
md"""  #### let's test Claim 2: when ecDNA copy number is sufficiently high to maintain diverse ecDNAs in single cells, tunable evolutionary dynamics can emerge in otherwise clonal populations.

"""

# ╔═╡ ba17fa2d-7243-45aa-86d2-c740029e544b
md"""

We directly show the stability of the internal equilibrium by varying PCN and the fitness optimum (by varying [Tet] concentration) in three ways:
\
\
1) examining the phase diagram in terms of allele frequency by dividing mean TetA copy number by TCN (TCN = PCN+1). The nice thing about the eigenvector result is that we can immediately get the stationary distribution without having to simulate the differential equations. So, we can show how the stationary distribution changes as a function of the matrix, which is a function of PCN and [Tet] concentration. We can make a phase diagram in which we vary PCN and [Tet] systematically, and show how does the stability of the internal equilibrium (in allele frequency space) varies.
\
When there is an internal allele frequency equilibrium, then the stationary state has a distribution of cells with varying TetA transposon copy number.
\
\
2) examining the stationary distribution for representative PCN and [Tet] cases, based on the phase diagram.
\
\
3) show pulse dynamics for these representative PCN and [Tet] cases, starting from the stationary distribution in the constant [Tet] scenario.

#### I put these panels together into a figure demonstrating Claim 2.
##### Panel Set A: phase diagrams.
##### Panel Set B: stationary distributions for representative PCN and [Tet] cases, based on the phase diagram.
##### Panel Set C: pulse dynamics for the corresponding PCN and [Tet] parameters.

#### Summary of findings

When PCN == 5, and TET_CONC == 4, and initial population to 100% TetA == 4, **Fitness declines!!** This result shows that the "optimal state" at the top of the landscape is not stable. This is a very nice result, that demonstrates the phase transition to tunable dynamics. Also note that the stationary distribution depends on which absorbing state (tetA == 1 or tetA == PCN+1) has higher fitness.

When PCN == 10 and TET_CONC == 5, the distribution is completely found at the boundaries in the constant [Tet] simulation. The distribution is not tunable in the pulsed [Tet] simulation. 

When PCN == 16 and TET_CONC == 8, it is almost flat in the middle, and most of the  distribution is at each boundary. The distribution weakly tunable in the pulsed [Tet] simulation. 

When PCN == 20 and TET_CONC == 10, there is a small bump in the middle, and some distribution at each boundary. The distribution is somewhat tunable in the pulsed [Tet] simulation. 

When PCN == 40 and TET_CONC == 20, there is a nice bell curve, and no distribution at the boundaries. The distribution is highly tunable in the pulsed [Tet] simulation.


"""

# ╔═╡ 83f4511d-ffe2-4677-9940-3cbed562845f
md"""  #### let's first make a phase diagram of how the stationary distribution depends on PCN and [Tet] concentration.  

The stationary distribution is a function of two parameters, PCN, [Tet]. To draw a heatmap, we need some summary statistic of the distribution. Let's try entropy,  variance, and 
\
\
$equilibrium\ allele\ frequency = \frac{mean\ tetA\ copy\ number}{PCN+1}$.

"""

# ╔═╡ 1fe4e134-5347-4b18-8926-db04bcd104b3
begin ## make eigenspecies_stationary_distribution_matrix.
	PCN_range_max = 40
	Tet_conc_range_max = 30 

	## Initialize the matrix. 
	## Each row is a plasmid copy number, each column is a [Tet] concentration.
	## Each entry holds the stationary distribution given (PCN, [Tet]).
	eigenspecies_stationary_distribution_matrix = Matrix{Vector{Float64}}(undef, PCN_range_max, Tet_conc_range_max)

	for pcn in 1:PCN_range_max
    	for tet in 1:Tet_conc_range_max
			## The first argument is PCN, the second is [Tet] concentration
        	eigenspecies_stationary_distribution_matrix[pcn, tet] = 	CalcNormalizedTopEigenvector(pcn, tet)
    	end
	end
end

# ╔═╡ e3f29f0d-2197-4f72-b84b-098e4372f424
## calculate the fitness variance of each stationary distribution in the matrix
final_fitness_variance_eigenspecies_matrix = [calc_fitness_variance(eigenspecies_stationary_distribution_matrix[pcn, tet], tet) for pcn in 1:PCN_range_max, tet in 1:Tet_conc_range_max]

# ╔═╡ e3ab744f-2068-4b86-b87d-85a347137409
stationary_distribution_fitness_variance_map = heatmap(final_fitness_variance_eigenspecies_matrix', ylabel="[Tet] concentration", xlabel="plasmid copy number", title="fitness variance")

# ╔═╡ c0bbcd2b-09c7-4422-91e6-d78b7550cc43
savefig(stationary_distribution_fitness_variance_map, "../results/modeling-results/stationary_distribution_fitness_variance.pdf")

# ╔═╡ 6982529f-7aab-4f1d-bffe-50b0642e0112
## calculate the tetA copy number variance of each stationary distribution in the matrix
final_tetA_copy_variance_eigenspecies_matrix = [calc_tetA_copy_number_variance(eigenspecies_stationary_distribution_matrix[pcn, tet]) for pcn in 1:PCN_range_max, tet in 1:Tet_conc_range_max]

# ╔═╡ d1cdb20e-df8f-45be-a9a9-3c445189877d
stationary_distribution_tetA_variance_map = heatmap(final_tetA_copy_variance_eigenspecies_matrix', ylabel="[Tet] concentration", xlabel="plasmid copy number", title="tetA copy number variance")

# ╔═╡ 648fa994-88ed-470c-9f85-c5a9713a49b0
savefig(stationary_distribution_tetA_variance_map, "../results/modeling-results/stationary_distribution_tetA_variance.pdf")

# ╔═╡ 27bc3701-8120-4be6-9dac-22f68687dbe5
## calculate the entropy of each stationary distribution in the matrix.
final_entropy_eigenspecies_matrix = map(Entropy, eigenspecies_stationary_distribution_matrix)

# ╔═╡ 0d52cfa4-142c-4780-a2b0-03a7f2b4e43d
stationary_distribution_entropy_map = heatmap(final_entropy_eigenspecies_matrix', ylabel="[Tet] concentration", xlabel="plasmid copy number", title="Shannon entropy")

# ╔═╡ d5c3f4e7-e57b-47df-acb5-af76ed7f6748
savefig(stationary_distribution_entropy_map, "../results/modeling-results/stationary_distribution_entropy.pdf")

# ╔═╡ e90eeda2-9b71-43b5-8ad8-68d8ad66298d
md"""  #### let's now make a phase diagram showing the stationary distribution in terms of _allele frequency_ (tetA copy number/TCN), where TCN = PCN+1. Note that the optimal number of TetA copies == [Tet] concentration in this model, for easy interpretation of the model results.

"""

# ╔═╡ ac632907-7a80-48fe-b6b7-d796a48d8032
## calculate the mean allele frequency (tetA copy number/(PCN+1)) of each stationary distribution in the matrix
final_allele_frequency_eigenspecies_matrix = [calc_mean_tetA_copy_number(eigenspecies_stationary_distribution_matrix[pcn, tet])/(pcn+1) for pcn in 1:PCN_range_max, tet in 1:Tet_conc_range_max]

# ╔═╡ 2a991b7a-b2bc-42bb-8f0d-6ffc3408ada9
stationary_distribution_allele_frequency_map = heatmap(final_allele_frequency_eigenspecies_matrix', ylabel="[Tet] concentration", xlabel="plasmid copy number", title="equilibrium TetA allele frequency")

# ╔═╡ 135f8788-0aa6-4e9f-a706-10eafb638842
savefig(stationary_distribution_allele_frequency_map, "../results/modeling-results/stationary_distribution_equilibrium_tetA_allele_frequency.pdf")

# ╔═╡ f1eb2fd8-5604-4ff8-9f88-421e72d970ab
md""" ### make a phase diagram showing how mean population fitness deviates from the optimal population fitness.
"""

# ╔═╡ bc5318e4-728c-41da-af1d-f2b08bf5fa23
## calculate how mean population fitness deviates from the optimal population fitness stationary distribution in the matrix: 
final_mean_fitness_deviation_eigenspecies_matrix = [calc_mean_fitness(eigenspecies_stationary_distribution_matrix[pcn, tet],tet)/fitness_function(tet, tet) for pcn in 1:PCN_range_max, tet in 1:Tet_conc_range_max]

# ╔═╡ a31fbb70-9b83-4f33-a612-f2b2c115a5f4
let
	final_mean_fitness_deviation_map = heatmap(final_mean_fitness_deviation_eigenspecies_matrix', ylabel="[Tet] concentration", xlabel="plasmid copy number", title="fitness deviation of stationary distribution")

	savefig(final_mean_fitness_deviation_map, "../results/modeling-results/final_mean_fitness_deviation_map.pdf")
	final_mean_fitness_deviation_map
end

# ╔═╡ 74279b9d-4de8-4b03-9d6b-0985d7e48ede
md""" Notice how the fitness deviation heat map and equilibrium allele frequency phase diagram shows that, when PCN == 5 and TET_CONC == 4 and given an initial population to 100% TetA == 4, **fitness declines!!** This result shows that the "optimal state" at the top of the landscape is not stable. This is a very nice result, that demonstrates the phase transition to tunable dynamics. Also note that the stationary distribution depends on which absorbing state (tetA == 1 or tetA == PCN+1) has higher fitness."""

# ╔═╡ 6ccd78c9-6e5c-453c-9a98-57c775ec61d9
md""" #### Draw stationary distributions for representative PCN & [Tet] values, based on these phase diagrams.
"""

# ╔═╡ 3067c892-5fdb-4cbf-aa04-110dccfde19c
begin
	final_pop_pcn5_tet4 = CalcNormalizedTopEigenvector(5, 4)
	final_pop_pcn10_tet5 = CalcNormalizedTopEigenvector(10, 5)
	final_pop_pcn16_tet8 = CalcNormalizedTopEigenvector(16, 8)
	final_pop_pcn20_tet10 = CalcNormalizedTopEigenvector(20, 10)
	final_pop_pcn40_tet20 = CalcNormalizedTopEigenvector(40, 20)
end

# ╔═╡ ccdce86d-e8a1-4f36-8f23-2296aef613db
let
	my_xvec = collect(1:(5+1))
	final_pop_pcn5_tet4_fig = bar(my_xvec, final_pop_pcn5_tet4, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([4], linestyle=:dash, label="[Tet]")
	savefig(final_pop_pcn5_tet4_fig, "../results/modeling-results/final_pop_pcn5_tet4.pdf")
	final_pop_pcn5_tet4_fig
end

# ╔═╡ 5cea363f-de01-41be-8e28-1335d657a7cd
let
	my_xvec = collect(1:(10+1))
	final_pop_pcn10_tet5_fig = bar(my_xvec, final_pop_pcn10_tet5, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([5], linestyle=:dash, label="[Tet]")
	savefig(final_pop_pcn10_tet5_fig, "../results/modeling-results/final_pop_pcn10_tet5.pdf")
	final_pop_pcn10_tet5_fig
end

# ╔═╡ c9455bc3-b3d7-4934-9a4c-cc138aa1e57a
savefig("../results/modeling-results/final_pop_pcn10_tet5.pdf")

# ╔═╡ f1323fcd-ea2f-4165-a5fa-e0888f7cba6b
let
	my_xvec = collect(1:(16+1))
	final_pop_pcn16_tet8_fig = bar(my_xvec, final_pop_pcn16_tet8, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([8], linestyle=:dash, label="[Tet]")
	savefig(final_pop_pcn16_tet8_fig, "../results/modeling-results/final_pop_pcn16_tet8.pdf")
	final_pop_pcn16_tet8_fig
end

# ╔═╡ 562ecba7-53c6-4f7b-8eee-e3a4999ea22f
let
	my_xvec = collect(1:(20+1))
	final_pop_pcn20_tet10_fig = bar(my_xvec, final_pop_pcn20_tet10, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([10], linestyle=:dash, label="[Tet]")
	savefig(final_pop_pcn20_tet10_fig, "../results/modeling-results/final_pop_pcn20_tet10.pdf")
	final_pop_pcn20_tet10_fig
end

# ╔═╡ 25faf55d-12b6-445f-98c2-d14e2fa6ac0d
let
	my_xvec = collect(1:(40+1))
	final_pop_pcn40_tet20_fig = bar(my_xvec, final_pop_pcn40_tet20, label="", xlabel="tetA copy number", ylabel="subpopulation frequency")
	# Add a vertical dashed line at x = TET_CONC
	vline!([20], linestyle=:dash, label="[Tet]")
	savefig(final_pop_pcn40_tet20_fig, "../results/modeling-results/final_pop_pcn40_tet20.pdf")
	final_pop_pcn40_tet20_fig
end

# ╔═╡ b1a2b45e-2a63-4dff-b3f8-546a7e203791
md"""
When PCN == 5, and TET_CONC == 4, and initial population to 100% TetA == 4, **Fitness declines!!** This result shows that the "optimal state" at the top of the landscape is not stable. This is a very nice result, that demonstrates the phase transition to tunable dynamics. Also note that the stationary distribution depends on which absorbing state (tetA == 1 or tetA == PCN+1) has higher fitness.

When PCN == 10 and TET_CONC == 5, the distribution is completely found at the boundaries in the constant [Tet] simulation. The distribution is not tunable in the pulsed [Tet] simulation. 

When PCN == 16 and TET_CONC == 8, it is almost flat in the middle, and most of the  distribution is at each boundary. The distribution is weakly tunable in the pulsed [Tet] simulation. 

When PCN == 20 and TET_CONC == 10, there is a small bump in the middle, and some distribution at each boundary. The distribution is somewhat tunable in the pulsed [Tet] simulation. 

When PCN == 40 and TET_CONC == 20, there is a nice bell curve, and no distribution at the boundaries. The distribution is highly tunable in the pulsed [Tet] simulation.
"""

# ╔═╡ 11e0bfe0-c9c5-41bc-8c17-0eed797691b1
md""" 
#### Draw pulse dynamics for these same representative PCN & [Tet] values.
"""

# ╔═╡ 950aff23-d57e-45ce-a57b-2f3c3571e619
md"""
For these time courses, set the initial configuration of the population to the stationary distribution under the corresponding constant [Tet] scenario. We can see whether the populations are tunable based on the dynamics under pulse [Tet] conditions.
"""

# ╔═╡ 8648e178-59f5-49db-b47a-b00b703537ca
begin
	## solve the ODEs for the dynamics
	pcn5_tet4_sol = SolvePulseTetQuasispeciesSystem(final_pop_pcn5_tet4, 5, 4)
	pcn10_tet5_sol = SolvePulseTetQuasispeciesSystem(final_pop_pcn10_tet5, 10, 5)
	pcn16_tet8_sol = SolvePulseTetQuasispeciesSystem(final_pop_pcn16_tet8, 16, 8)
	pcn20_tet10_sol = SolvePulseTetQuasispeciesSystem(final_pop_pcn20_tet10, 20, 10)
	pcn40_tet20_sol = SolvePulseTetQuasispeciesSystem(final_pop_pcn40_tet20, 40, 20)
	
	## get the mean tetA copy number dynamics from the ODE solutions
	pcn5_tet4_mean_copy_num_vec = CalcMeanTetACopyNumberVec(pcn5_tet4_sol)
	pcn10_tet5_mean_copy_num_vec = CalcMeanTetACopyNumberVec(pcn10_tet5_sol)
	pcn16_tet8_mean_copy_num_vec = CalcMeanTetACopyNumberVec(pcn16_tet8_sol)
	pcn20_tet10_mean_copy_num_vec = CalcMeanTetACopyNumberVec(pcn20_tet10_sol)
	pcn40_tet20_mean_copy_num_vec = CalcMeanTetACopyNumberVec(pcn40_tet20_sol)

	## make vectors of tet concentration over time in the model solutions
	## for calculating summary statistics.
	pcn5_tet4_pulse_vec = [TetPulseFunction(t, tet_conc=4) for t in pcn5_tet4_sol.t]
	pcn10_tet5_pulse_vec = [TetPulseFunction(t, tet_conc=5) for t in pcn10_tet5_sol.t]
	pcn16_tet8_pulse_vec = [TetPulseFunction(t, tet_conc=8) for t in pcn16_tet8_sol.t]
	pcn20_tet10_pulse_vec = [TetPulseFunction(t, tet_conc=10) for t in pcn20_tet10_sol.t]
	pcn40_tet20_pulse_vec = [TetPulseFunction(t, tet_conc=20) for t in pcn40_tet20_sol.t]

	## calculate rates of copy number change from the ODE solutions.
	pcn5_tet4_d_pulse_mean_copy_num_vec = CalcTetACopyNumberVelocity(pcn5_tet4_sol)
	pcn10_tet5_d_pulse_mean_copy_num_vec = CalcTetACopyNumberVelocity(pcn10_tet5_sol)
	pcn16_tet8_d_pulse_mean_copy_num_vec = CalcTetACopyNumberVelocity(pcn16_tet8_sol)
	pcn20_tet10_d_pulse_mean_copy_num_vec = CalcTetACopyNumberVelocity(pcn20_tet10_sol)
	pcn40_tet20_d_pulse_mean_copy_num_vec = CalcTetACopyNumberVelocity(pcn40_tet20_sol)

	### calculate Price equation LHS from the ODE solutions.
	pcn5_tet4_pulse_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(pcn5_tet4_sol, pcn5_tet4_pulse_vec)
	
	pcn10_tet5_pulse_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(pcn10_tet5_sol, pcn10_tet5_pulse_vec)
	
	pcn16_tet8_pulse_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(pcn16_tet8_sol, pcn16_tet8_pulse_vec)
	
	pcn20_tet10_pulse_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(pcn20_tet10_sol, pcn20_tet10_pulse_vec)
	
	pcn40_tet20_pulse_Price_equation_LHS_vec = CalcPriceEquationLHS_vec(pcn40_tet20_sol, pcn40_tet20_pulse_vec)
	
	### calculate tetA copy number--fitness covariance from the ODE solutions.
	pcn5_tet4_pulse_copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(pcn5_tet4_sol, pcn5_tet4_pulse_vec)

	pcn10_tet5_pulse_copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(pcn10_tet5_sol, pcn10_tet5_pulse_vec)

	pcn16_tet8_pulse_copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(pcn16_tet8_sol, pcn16_tet8_pulse_vec)

	pcn20_tet10_pulse_copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(pcn20_tet10_sol, pcn20_tet10_pulse_vec)

	pcn40_tet20_pulse_copy_num_covariance_vec = CalcTetACopyNumberFitnessCovariance(pcn40_tet20_sol, pcn40_tet20_pulse_vec)


	
end

# ╔═╡ 7ba197fe-33fc-4a49-a52f-1055fd53030e
let
	pcn5_tet4_pulsefig = plot(pcn5_tet4_sol.t, pcn5_tet4_mean_copy_num_vec, label="Mean tetA copy number", xlabel="time", ylabel="mean tetA copy number")
	# Add a horizontal dashed line at y = TET_CONC
	hline!([4], linestyle=:dash, label="[Tet] = 4")
	savefig(pcn5_tet4_pulsefig, "../results/modeling-results/pcn5_tet4_pulsefig.pdf")
	pcn5_tet4_pulsefig
end

# ╔═╡ 2905fb49-99a6-41e2-90a6-7f8fa36d4cbb
let
	pcn10_tet5_pulsefig = plot(pcn10_tet5_sol.t, pcn10_tet5_mean_copy_num_vec, label="Mean tetA copy number", xlabel="time", ylabel="mean tetA copy number")
		# Add a horizontal dashed line at y = TET_CONC
	hline!([5], linestyle=:dash, label="[Tet] = 5")
	savefig(pcn10_tet5_pulsefig, "../results/modeling-results/pcn10_tet5_pulsefig.pdf")
	pcn10_tet5_pulsefig
end

# ╔═╡ f2702a90-936a-45b3-8162-ff07936f4c80
let
	pcn16_tet8_pulsefig = plot(pcn16_tet8_sol.t, pcn16_tet8_mean_copy_num_vec, label="Mean tetA copy number", xlabel="time", ylabel="mean tetA copy number")
	# Add a horizontal dashed line at y = TET_CONC
	hline!([8], linestyle=:dash, label="[Tet] = 8")
	savefig(pcn16_tet8_pulsefig, "../results/modeling-results/pcn16_tet8_pulsefig.pdf")
	pcn16_tet8_pulsefig
end

# ╔═╡ 489d6a31-7b62-42ec-8b54-82cd47420b68
let
	pcn20_tet10_pulsefig = plot(pcn20_tet10_sol.t, pcn20_tet10_mean_copy_num_vec, label="Mean tetA copy number", xlabel="time", ylabel="mean tetA copy number")
		# Add a horizontal dashed line at y = TET_CONC
	hline!([10], linestyle=:dash, label="[Tet] = 10")
	savefig(pcn20_tet10_pulsefig, "../results/modeling-results/pcn20_tet10_pulsefig.pdf")
	pcn20_tet10_pulsefig
end

# ╔═╡ b3b9c85d-9090-4c8f-b8a4-fc2696288cca
let
	pcn40_tet20_pulsefig = plot(pcn40_tet20_sol.t, pcn40_tet20_mean_copy_num_vec, label="Mean tetA copy number", xlabel="time", ylabel="mean tetA copy number")
	# Add a horizontal dashed line at y = TET_CONC
	hline!([20], linestyle=:dash, label="[Tet] = 20")
	savefig(pcn40_tet20_pulsefig, "../results/modeling-results/pcn40_tet20_pulsefig.pdf")
	pcn40_tet20_pulsefig
end

# ╔═╡ 60a88dd5-243e-433a-aa91-1a29526e6e8c
md""" 
#### Draw Price theorem predictions for these same representative PCN & [Tet] values.
"""

# ╔═╡ bef56080-368b-4d64-af08-740394dfe5d2
let
	pcn5_tet4_price_eq_fig = plot(pcn5_tet4_sol.t, 
		pcn5_tet4_pulse_Price_equation_LHS_vec,  label="Price equation prediction",alpha=0.5)
	
	plot!(pcn5_tet4_sol.t, pcn5_tet4_pulse_copy_num_covariance_vec, label="tetA copy number fitness covariance",alpha=0.5)
	
	plot!(pcn5_tet4_sol.t, pcn5_tet4_d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:right,alpha=0.5)

	savefig(pcn5_tet4_price_eq_fig, "../results/modeling-results/pcn5_tet4_price_eq_fig.pdf")
	pcn5_tet4_price_eq_fig
end

# ╔═╡ 58535a41-226d-4193-9c8a-70d4c3bbd1a0
let
	pcn10_tet5_price_eq_fig = plot(pcn10_tet5_sol.t, 
		pcn10_tet5_pulse_Price_equation_LHS_vec,  label="Price equation prediction",alpha=0.5)
	
	plot!(pcn10_tet5_sol.t, pcn10_tet5_pulse_copy_num_covariance_vec, label="tetA copy number fitness covariance",alpha=0.5)
	
	plot!(pcn10_tet5_sol.t, pcn10_tet5_d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:right,alpha=0.5)
	
	savefig(pcn10_tet5_price_eq_fig, "../results/modeling-results/pcn10_tet5_price_eq_fig.pdf")
	pcn10_tet5_price_eq_fig
end

# ╔═╡ 126fe9a6-3fe1-4805-8156-9b112cdc6132
let
	pcn16_tet8_price_eq_fig = plot(pcn16_tet8_sol.t, 
		pcn16_tet8_pulse_Price_equation_LHS_vec,  label="Price equation prediction",alpha=0.5)
		
	plot!(pcn16_tet8_sol.t, pcn16_tet8_pulse_copy_num_covariance_vec, label="tetA copy number fitness covariance",alpha=0.5)
	
	plot!(pcn16_tet8_sol.t, pcn16_tet8_d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:right,alpha=0.5)

	savefig(pcn16_tet8_price_eq_fig, "../results/modeling-results/pcn16_tet8_price_eq_fig.pdf")
	pcn16_tet8_price_eq_fig
end

# ╔═╡ 5768920c-15b3-4be2-bab5-03c0d4d3c9e3
let
	pcn20_tet10_price_eq_fig = plot(pcn20_tet10_sol.t, 
		pcn20_tet10_pulse_Price_equation_LHS_vec,  label="Price equation prediction",alpha=0.5)
	
	plot!(pcn20_tet10_sol.t, pcn20_tet10_pulse_copy_num_covariance_vec, label="tetA copy number fitness covariance",alpha=0.5)
	
	plot!(pcn20_tet10_sol.t, pcn20_tet10_d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:none,alpha=0.5)

	savefig(pcn20_tet10_price_eq_fig, "../results/modeling-results/pcn20_tet10_price_eq_fig.pdf")
	pcn20_tet10_price_eq_fig
end

# ╔═╡ ee8b810c-4b6e-4b6a-8c0a-7304392217f8
let
	pcn40_tet20_price_eq_fig = plot(pcn40_tet20_sol.t, 
		pcn40_tet20_pulse_Price_equation_LHS_vec,  label="Price equation prediction",alpha=0.5)
		
	plot!(pcn40_tet20_sol.t, pcn40_tet20_pulse_copy_num_covariance_vec, label="tetA copy number fitness covariance",alpha=0.5)
	
	plot!(pcn40_tet20_sol.t, pcn40_tet20_d_pulse_mean_copy_num_vec, label="Rate of change of mean tetA copy number", xlabel="time", ylabel="tetA copy number derivative\nand covariance with fitness", legend=:none,alpha=0.5)

	savefig(pcn40_tet20_price_eq_fig, "../results/modeling-results/pcn40_tet20_price_eq_fig.pdf")
	pcn40_tet20_price_eq_fig
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DifferentialEquations = "~7.14.0"
Distributions = "~0.25.112"
Images = "~0.26.1"
Plots = "~1.40.8"
PlutoUI = "~0.7.60"
Statistics = "~1.11.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "16d06f37a028b528329e766485637f4b055ca375"

[[deps.ADTypes]]
git-tree-sha1 = "eea5d80188827b35333801ef97a40c2ed653b081"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.9.0"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0dd7edaff278e346eb0ca07a7e75c9438408a3ce"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.10.3"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "a2c85f53ddcb15b4099da59867868bd40f005579"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.7.5"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "8dba4ea86fce8ed94977ed9e56cf5c6d67c494e5"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.11.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "9ebb045901e9bbf58767a9f34ff89831ed711aae"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "f9d7112bfff8a19a3a4ea4e03a8e6a91fe8456bf"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "066f60231c1b0ae2905ffd2651e207accd91f627"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.48.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "TruncatedStacktraces"]
git-tree-sha1 = "9c36ac68cf6c59a1c2569aadd7a62c47ed5c8eb5"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.158.1"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NonlinearSolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "19dbd44d18bbfdfcf5e56c99cea9b0ed23df350a"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "3.9.1"
weakdeps = ["OrdinaryDiffEq", "OrdinaryDiffEqCore", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ab1e6515ce15f01316a9825b02729fefa51726bd"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.23.0"

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
git-tree-sha1 = "d851f2ca05f3cec9988f081b047a778a58b48aaf"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.14.0"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "16611777adf4a818f18c33210895947814180964"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.16"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = "Enzyme"
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = "ForwardDiff"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = "PolyesterForwardDiff"
    DifferentiationInterfaceReverseDiffExt = "ReverseDiff"
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "d7477ecdafb813ddee2ae727afa94e9dcb5f3fb0"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.112"

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

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "9c3a42611e525352e9ad5e4134ddca5c692ff209"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.4"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "3f03d94c71126b6cfe20d3cbcc41c5cd27e1c419"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.4"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cbf5edddb61a43669710cbc2241bc08b36d9e660"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.4"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "62ca0547a14c57e98154423419d8a342dca75ca9"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.4"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "b10bdafd1647f57ace3885143936749d61638c3b"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

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
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

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
git-tree-sha1 = "64d8e93700c7a3f28f717d265382d52fac9fa1c1"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.12"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee28ddcd5517d54e417182fec3886e7412d3926f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f31929b9e67066bee48eec8b03c0df47d31a74b3"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.8+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1dc470db8b1131cfc7fb4c115de89fe391b9e780"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "33485b4e40d1df46c806498c73ea32dc17475c59"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.1"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "3447781d4c80dbe6d71d239f7cfb1f8049d4c84f"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.6"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d65554bad8b16d9562050c67e7223abf91eaba2f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.13+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "6f0a801136cb9c229aebea0df296cdcd471dbcd1"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.5"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "44664eea5408828c03e5addb84fa4f916132fc26"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.8.1"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "e0884bdf01bbbb111aea77c348368a86fb4b5ab6"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.1"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "12fdd617c7fe25dc4a6cc804d657cc4b2230302b"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "be8e690c3973443bec584db3346ddc904d4884eb"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "Requires", "TranscodingStreams"]
git-tree-sha1 = "a0746c21bdc986d0dc293efa6b1faee112c37c28"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.53"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "6e0fb267f2c869df59b6bf9f8cb0f2794f2d85e9"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.13.7"
weakdeps = ["FastBroadcast"]

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "71dd823bf0a5aaa8264e1143b96a0baa631f4ce0"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.7"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "360f6039babd6e4d6364eff0d4fc9120834a2d9a"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.2.1"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "dc82b3c3640a4362f09e4d7594b4387a820857e4"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.3"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "8941ad4bdd83768359801982e143008349b1a827"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.35.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
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
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "fa7fd067dca76cadd880f1ca937b4f387975a9f5"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.16.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8084c25a250e00ae427a379a5b607e7aed96a2dd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.171"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "16a726dba99685d9e94c8d0a8f655383121fc608"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.0.1"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "1130dbe1d5276cb656f6e1094ce97466ed700e5a"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

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

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "3cebfc94a0754cc329ebc3bab1e6c89621e791ad"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.20"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "4d8944f32db2b07a2bdf8477e878bcb9c9ea2308"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.15.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
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
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "f4cb457ffac5f5cf695699f82c537073958a6a6c"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.2+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "cd892f12371c287dc50d6ad3af075b088b6f2d48"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.89.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "8e3c5978d0531a961f70d2f2730d1d16ed3bbd12"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.1.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b4498d40bf35da0b6d22652ff2e9d8820590b3c6"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.1.2"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "TruncatedStacktraces"]
git-tree-sha1 = "f4b7b11e2e4c0f4a7fe15d2edcec1e1ce2917d67"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.7.1"
weakdeps = ["EnzymeCore"]

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"

[[deps.OrdinaryDiffEqDefault]]
deps = ["DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport"]
git-tree-sha1 = "c8223e487d58bef28a3535b33ddf8ffdb44f46fb"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.1.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "e63ec633b1efa99e3caa2e26a01faaa88ba6cef9"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.1.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "f63938b8e9e5d3a05815defb3ebdbdcf61ec0a74"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.1.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["DiffEqBase", "FastBroadcast", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "fea595528a160ed5cade9eee217a9691b1d97714"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.1.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLOperators"]
git-tree-sha1 = "795221c662698851328cb7787965ab4a180d9468"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.1.1"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a7cc74d3433db98e59dc3d58bc28174c6c290adf"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.1.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "925a91583d1ab84f1f0fea121be1abf1179c5926"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.1.1"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "103e017ff186ac39d731904045781c9bacfca2b0"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.1.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "9f8f52aad2399d7714b400ff9d203254b0a89c4a"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.1.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "0f81a77ede3da0dc714ea61e81c76b25db4ab87a"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.1.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "d4bb32e09d6b68ce2eb45fb81001eab46f60717a"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.2.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "590561f3af623d5485d070b4d7044f8854535f5a"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.2.1"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "a2a4119f3e35f7982f78e17beea7b12485d179e9"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.2.1"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "ef44754f10e0dfb9bb55ded382afed44cd94ab57"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.1.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "StaticArrays"]
git-tree-sha1 = "a8b7f8107c477e07c6a6c00d1d66cac68b801bbc"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.1.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport"]
git-tree-sha1 = "da525d277962a1b76102c79f30cb0c31e13fe5b9"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.1.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "332f9d17d0229218f66a73492162267359ba85e9"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.1.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "41c09d9c20877546490f907d8dffdd52690dd65f"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.1.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "96b47cdd12cb4ce8f70d701b49f855271a462bd4"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.2.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "f6683803a58de600ab7a26d2f49411c9923e9721"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.1.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "7dbe4ac56f930df5e9abd003cedb54e25cbbea86"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.2.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "348fd6def9a88518715425025eadd58517017325"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.1.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "1b0d894c880e25f7d0b022d7257638cf8ce5b311"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.1.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "4e8b8c8b81df3df17e2eb4603115db3b30a88235"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.1.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "96552f7d4619fabab4038a29ed37dd55e9eb513a"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.1.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "81d7841e73e385b9925d5c8e4427f2adcdda55db"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.1.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

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
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "650a022b2ce86c7dcfbdecf00f78afeeb20e5655"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

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
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6d38fea02d983051776a856b7df75b30cf9a3c1f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.16"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "1a9cfb2dc2c2f1bd63f1906d72af39a79b49b736"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.11"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "b034171b93aebc81b3e1890a036d13a9c4a9e3e0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

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
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "98ca7c29edd6fc79cd74c61accb7010a4e7aee33"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.6.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "50ed64cd5ad79b0bef71fdb6a11d10c3448bfef0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.56.1"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "991d2a8900e687e2c693d587daa739c8fda01366"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.0"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "e39c5f217f9aca640c8e27ab21acf557a3967db5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.10"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "44021f3efc023be3871195d8ad98b865001a2fa1"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.12.3"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "4b33e0e081a825dbfaf314decf58fa47e53d6acb"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "35d346f6aa78814859f1f69cd4f41c24108afa00"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.7"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerDataInterpolationsExt = "DataInterpolations"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "b906758c107b049b6b71599b9f928d9b14e5554a"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.23.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DataStructures", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "ccc32032d8f6790ad198c99fb8ef20d8100a0de4"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.7"
weakdeps = ["Colors"]

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsColorsExt = "Colors"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "87d51a3ee9a4b0d2fe054bdd3fc2436258db2603"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.1.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "920acf6ae36c86f23969fea1d317e040dbfccf53"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.4.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "b49030a24847e8c3e85eed3359825297ffb70dfa"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.69.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "e87efb31e5360cb223a151c2398903dc2faeb32b"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.26.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "4bc96df5d71515b1cb86dd626915f06f4c0d46f5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.33"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

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
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "38f139cc4abf345dd4f22286ec000728d5e8e097"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.2"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3a6f063d690135f5c1ba351412c82bae4d1402bf"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.25"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "be986ad9dac14888ba338c2554dcfec6939e1393"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.1"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

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
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "e7f5b81c65eb858bed630fe006837b935518aca5"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.70"

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
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

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
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "7dfa0fd9c783d3d0cc43ea1af53d69ba45c447df"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

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
# ╟─44e8209c-8d79-11ef-24db-ff1ea28d01b3
# ╟─d4d1f72f-7aa7-46e7-9be0-f38fe991bc0c
# ╟─b706fac1-3a38-423b-bed3-db6518475c73
# ╟─c3d81737-fc16-44a0-b3de-bab5e0f7aab2
# ╟─3931b5c1-51df-42f5-86e9-09ddba2d2f11
# ╟─f8be1737-038a-41c0-9f61-e1980b005ed2
# ╟─b1f80124-822d-46e2-9386-54a0117f833d
# ╟─749b2bd1-4ccc-48e0-9ab4-bc701500c728
# ╟─137b0ae4-4718-49bc-8e91-60b0e9f90934
# ╟─35b3eeca-74a6-4480-8ee6-0c52d3f31a57
# ╟─53bcbf1c-a2e2-4931-a4f8-88b7c0ea01f4
# ╟─4f840ddf-058f-48cc-90ea-609d6cf09ccb
# ╟─b8e8c64d-841d-4591-9fce-f2648d5d4f53
# ╟─74892533-c56c-4fbe-82d0-3de7127702f0
# ╟─7cc3e6c5-f6f8-42da-8aec-38147950da8a
# ╟─b9277783-9e31-44e5-82a3-a13af62c62e4
# ╟─b753c51f-6682-47e5-a015-c183e221aa32
# ╟─0abeaf9e-5322-4b6d-ba8f-81efd831080b
# ╟─08535a54-414b-4a3e-a718-b356cb22cb23
# ╟─83685c04-3a04-4d30-90e4-ab9ec7eee180
# ╟─39b0ddde-3cb4-4630-a4ff-a3ab163c57d8
# ╟─35adaf46-3b45-4208-aa74-e89a54c4a6d5
# ╟─248e7e49-749a-4110-a3df-12d9b873b949
# ╟─039ffb9c-2c1b-423c-a809-93f9e62bc688
# ╠═8c80a36a-372a-4fad-92bc-bb4426277f55
# ╟─cfd4050d-31a8-49ed-9994-076e55c02ceb
# ╠═9e9cc2c2-55be-456e-95f0-224b08fbeae3
# ╠═ef893778-f28f-4ade-ab70-978ded8cfcd5
# ╠═38fa13dc-2eca-4ac9-87fd-5d82e2e0ce41
# ╠═0a576694-2971-464e-9c91-35e7e3f29fd0
# ╠═bbb2c832-201f-409e-b262-ec7c74b6c490
# ╠═934ef964-bdc4-4d5d-bd97-e4e32b4d0380
# ╠═892bb0dd-da94-493f-bf72-cdc53a8b0181
# ╠═b78463f4-1deb-45bd-b091-7ae15a23471b
# ╠═412c944d-168d-4146-8bc8-7219db9c291a
# ╠═de13afba-dec7-4a7a-a26b-85e26d36a84f
# ╠═0cf84c01-0378-44be-b878-3c2119b468b9
# ╠═260df35d-6d1c-42be-8270-01fb62627938
# ╠═45362f0d-76ac-44e8-b6cd-0d1824b3a3b4
# ╠═49e1dd14-e839-45cf-8e0a-f0badd573d67
# ╠═1ab2e8f1-753a-4374-8b5b-1ad3de48a710
# ╠═1955042b-e29d-4c54-84c3-8d32bed550a3
# ╠═2cb2dae4-33af-4118-87c8-f41c8aba8225
# ╠═c99d3673-9f9a-4508-b257-360cc6d0733d
# ╠═348e1723-08b6-4968-896c-8459afc8bb08
# ╠═ea1bdd07-3d5a-49ef-a626-95d9b264d297
# ╠═f1b278a5-1ce4-47b3-b93e-969e1f9ec6bf
# ╠═638d7f76-512c-4c11-951d-76a77341b59d
# ╠═162efa78-e3a7-4ebb-ad7a-395a567071b1
# ╠═bd19d151-c3fc-4ff0-b2c5-4ffe07adcae0
# ╠═c10d54d1-eb66-4008-8f83-a54342ea8128
# ╠═0b129025-f83c-4560-8b3e-a56356966db4
# ╠═36bab24f-30a2-489f-94e0-fa5962a49b37
# ╠═66ee7dae-f87d-4b12-91e0-fffdbc420478
# ╠═bed5464d-47c8-46c4-9037-e73130a5b0e6
# ╠═621be509-5e7e-415a-bcc6-daa7b1903834
# ╠═192ab415-b73a-49d3-ba74-99b5d91ad482
# ╠═ea419d1b-1246-4999-a8a4-706c73e14ed8
# ╠═c899121e-2b9e-4519-b168-98994684a86d
# ╠═ac741b84-00f6-4f70-89a4-07dbb24ec9cb
# ╠═5bd6b242-9a1e-44ea-b319-6b7604e57624
# ╠═0d55f4ed-f01c-43c6-8a58-6edf9b77e9d9
# ╠═f238e7a6-86f1-4898-bc03-3fd80f4e1f3e
# ╟─8dfdd593-117b-430e-b988-2eaa70484b4c
# ╟─9bd490b7-3f11-487e-b454-785c7d87391b
# ╠═28bd0eee-f54e-4ab5-aecc-275fe3e8319a
# ╠═2f5d53c1-5b8b-40eb-a897-c1efbf1977a9
# ╠═6bd86018-8a50-4b8d-a2bd-40bfbe45829b
# ╠═a506ff86-41c6-44ac-adf5-c3fdb368cb02
# ╟─9241c94d-857c-4a31-93f1-cfed036cf2f2
# ╠═09c74702-0d19-4049-8f8e-6f91a0163dcc
# ╠═eb9ea298-95f2-4df7-966a-069ca8965c67
# ╠═e1cd955e-a069-426b-b222-c60790c809b3
# ╠═c40e9cfa-f58f-459c-a994-34bab25c20ad
# ╟─6fbce7a0-af24-4501-92e4-fbd6b5d27b65
# ╠═4fd0feac-346b-42bc-b0de-1bf0a4718928
# ╟─96241706-df79-4c86-838d-398302160e3d
# ╠═1d8f69c1-e2b2-492b-86a5-a9947fe9ac8a
# ╠═2ce5fa35-27ae-4382-a41d-fa835052412d
# ╠═6e3e1eec-03c2-4bb3-bb5a-9c27182753c7
# ╠═a25df7df-eb57-48a2-a538-44f7b7432387
# ╠═af8a47bd-e823-41ae-9813-0b484f040e4f
# ╠═d03e6d49-2c29-4101-9918-b8917fb037d8
# ╠═bc4a7835-dd9f-42e4-b61e-118a33c7ff96
# ╠═9028f452-a9cd-42ad-94ff-f86e54462a7c
# ╟─e4e8c61d-5503-490f-993d-bb1464def4c9
# ╠═d0c028ae-b991-4558-87e6-41b4f5bfe274
# ╟─7bfbb45d-57fc-4c9e-b929-ac44f7f88fef
# ╠═73873c97-daf6-4133-a657-90de30803517
# ╟─91d1760c-bbaf-466b-8eb7-623cb0dcd686
# ╠═9b92d86b-9f19-4fed-be2b-691d58ffaab6
# ╠═49830d58-03a3-4379-aeea-767c9a3eeb26
# ╠═1cff27ae-4465-4c38-8b26-5cace8a833d7
# ╟─e1046a07-ad24-45e0-b871-069640aa6df3
# ╟─4bddf691-83a7-4101-bb00-0fe463e5de77
# ╠═f96766c2-bd3b-4782-a884-a2d72667d0e4
# ╠═d153fde0-3435-49f2-b493-2c2fc8fec41f
# ╠═11e9cbd5-19f2-4401-ba19-8590d3f4dca7
# ╟─9a48f3c4-a5d3-4d53-8248-557c4dd86efb
# ╠═cb0d4f61-2a7d-4fe2-b2a7-8b2aafae22bd
# ╠═86e16ddc-a9e0-4e2c-822c-abc9f19d0ccf
# ╠═933c9aef-3f09-42aa-b7c5-8f7e5c63d1c9
# ╟─a2603db4-3fd8-4cc7-a943-8fa142961a61
# ╠═55f2289b-9268-4bca-9f0f-7f27e29697f8
# ╠═c83fbfaa-1f5f-40b9-a6c1-68d480c4dfe7
# ╠═348fabc0-3caa-4585-975e-688a26b7fa8a
# ╟─4acf6a58-6ca5-43e2-8169-75f169cae851
# ╠═e22f600c-4468-47ff-93f5-ad004ff5ea23
# ╟─d20f047d-68d3-44eb-a258-3642d1197a88
# ╠═362c0eff-b354-481d-bc5f-673613e8b0a6
# ╠═5f571b12-244c-4e7e-8774-883f0427ff06
# ╠═ab3d5352-abfa-4de3-8421-eb86472b91ea
# ╠═c9c7c6e9-75b3-4900-b146-010dd37f4123
# ╠═cef4643c-fd32-4a75-b9ab-0276341fc0bb
# ╠═76311c6d-8e96-4265-8bcf-fbc80a6870f3
# ╠═bc0d22dc-b1d2-48d8-aa9d-b0210ab3f9a7
# ╠═e9d58996-30f0-43b2-b13b-9cf7edfcd9eb
# ╠═026d3047-cca8-462f-8df3-b73bdeb4af31
# ╟─e15bbef0-42f1-47d5-af03-ce95468cba93
# ╠═153e0346-5632-4b96-929c-07cac8f8cfb6
# ╟─bc46fd90-a4c5-4ec0-80f6-c96c6901ac6d
# ╟─f3109642-2e16-4e11-8b1c-19667a380321
# ╟─a24529bf-7679-45b9-ac84-a84169ec7c7e
# ╠═38cd0bb1-137f-401b-8279-cc3a28ead6f3
# ╠═ab7f4870-6967-4910-b272-ed75fb8fbb54
# ╠═506fb600-7043-4562-8470-ffd2fd48d073
# ╠═47445b3c-bf48-485a-b291-79437595b45a
# ╟─adfa712e-80ca-405e-9854-1d0024c3d6f1
# ╠═5f9855dd-004f-4dc7-9673-ca5f81fb441b
# ╠═b1697098-07e0-4577-b144-9ecc3879422b
# ╠═d530e10c-63ff-49bb-a7b1-ea4b363738d0
# ╠═1f030b53-5503-490e-abf0-20101555d3b1
# ╠═798394c1-79a3-4886-b42d-c6f20c73c572
# ╟─138d4602-86bf-4923-b648-b6d273f6adac
# ╟─ba17fa2d-7243-45aa-86d2-c740029e544b
# ╟─83f4511d-ffe2-4677-9940-3cbed562845f
# ╠═1fe4e134-5347-4b18-8926-db04bcd104b3
# ╠═e3f29f0d-2197-4f72-b84b-098e4372f424
# ╠═e3ab744f-2068-4b86-b87d-85a347137409
# ╠═c0bbcd2b-09c7-4422-91e6-d78b7550cc43
# ╠═6982529f-7aab-4f1d-bffe-50b0642e0112
# ╠═d1cdb20e-df8f-45be-a9a9-3c445189877d
# ╠═648fa994-88ed-470c-9f85-c5a9713a49b0
# ╠═27bc3701-8120-4be6-9dac-22f68687dbe5
# ╠═0d52cfa4-142c-4780-a2b0-03a7f2b4e43d
# ╠═d5c3f4e7-e57b-47df-acb5-af76ed7f6748
# ╠═e90eeda2-9b71-43b5-8ad8-68d8ad66298d
# ╠═ac632907-7a80-48fe-b6b7-d796a48d8032
# ╠═2a991b7a-b2bc-42bb-8f0d-6ffc3408ada9
# ╠═135f8788-0aa6-4e9f-a706-10eafb638842
# ╠═f1eb2fd8-5604-4ff8-9f88-421e72d970ab
# ╠═bc5318e4-728c-41da-af1d-f2b08bf5fa23
# ╠═a31fbb70-9b83-4f33-a612-f2b2c115a5f4
# ╟─74279b9d-4de8-4b03-9d6b-0985d7e48ede
# ╠═6ccd78c9-6e5c-453c-9a98-57c775ec61d9
# ╠═3067c892-5fdb-4cbf-aa04-110dccfde19c
# ╠═ccdce86d-e8a1-4f36-8f23-2296aef613db
# ╠═5cea363f-de01-41be-8e28-1335d657a7cd
# ╠═c9455bc3-b3d7-4934-9a4c-cc138aa1e57a
# ╠═f1323fcd-ea2f-4165-a5fa-e0888f7cba6b
# ╠═562ecba7-53c6-4f7b-8eee-e3a4999ea22f
# ╠═25faf55d-12b6-445f-98c2-d14e2fa6ac0d
# ╟─b1a2b45e-2a63-4dff-b3f8-546a7e203791
# ╠═11e0bfe0-c9c5-41bc-8c17-0eed797691b1
# ╠═950aff23-d57e-45ce-a57b-2f3c3571e619
# ╠═8648e178-59f5-49db-b47a-b00b703537ca
# ╠═7ba197fe-33fc-4a49-a52f-1055fd53030e
# ╠═2905fb49-99a6-41e2-90a6-7f8fa36d4cbb
# ╠═f2702a90-936a-45b3-8162-ff07936f4c80
# ╠═489d6a31-7b62-42ec-8b54-82cd47420b68
# ╠═b3b9c85d-9090-4c8f-b8a4-fc2696288cca
# ╠═60a88dd5-243e-433a-aa91-1a29526e6e8c
# ╠═bef56080-368b-4d64-af08-740394dfe5d2
# ╠═58535a41-226d-4193-9c8a-70d4c3bbd1a0
# ╠═126fe9a6-3fe1-4805-8156-9b112cdc6132
# ╠═5768920c-15b3-4be2-bab5-03c0d4d3c9e3
# ╠═ee8b810c-4b6e-4b6a-8c0a-7304392217f8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

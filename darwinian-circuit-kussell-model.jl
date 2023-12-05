### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ 96e88c96-ffa3-440b-a5c4-df380f1ee881
begin
	using FileIO
	using Plots
	using PlutoUI
	using Statistics
	using LinearAlgebra
	using DifferentialEquations
end

# ╔═╡ bb669570-7cb6-11ee-00f9-95e5d1f0bf9d
md"""
# darwinian-circuit-kussell-model.jl 

by Rohan Maddamsetti  

Julia version 1.9.  

## Tentative manuscript title: Programmed gene expression control by diversity-maintaining synthetic genetic elements

##### We compare fitness variance to the rate of mean fitness change, and compare tetA copy number-fitness covariance to the rate of tetA copy number change.

"""

# ╔═╡ 1a7794a3-ee9b-4b0c-8ae0-3930a63f5e54
##load("../results/diagrams/Darwin-Markov-model-figure.jpg")

# ╔═╡ a66ee148-7cf2-43a4-9b14-407eb91da8ea
md""" 

### Notes from Lingchong
High copy number plasmids are expensive/burdensome resistance/cost. Despite this, many studies have shown widespread persistence of multicopy plasmids in diverse environments (CITATIONS).

1) we need previous literature on prevalence of high copy number plasmids in nature, and 2) work on its why high copy plasmids exist in nature.

Others have demonstrated that high copy plasmids can cause heterogeneity.

We can resolve this apparent paradox in a new way. Here we propose that high-copy plasmids can enable intracellular heterogeneity in plasmid identity, and this configuration allows very fast adaptation to different environments, to switch between highly resistant and highly fast growing states extremely quickly.

This story depends on speed of adaptation.

Critical data is on how copy number affects the timescale of response (both gain as well as the loss of the saturated state). Settling down into the fast growing state: hypothesis is that the high-copy plasmid will be the best, due to the ghost effect with the low copy number plasmid.

Also: with the high copy number plasmid and costly selection, the plasmid will never saturated. Also purifying selection will cause the rapid going back to the baseline.

Also can emphasize that adaptation in either direction does not depend on de novo genetic mutation (this is supporting evidence for this story, not main point).


"""

# ╔═╡ 6d0b7c31-69c0-4f9c-9f7e-1cb43502f58f
##load("../results/diagrams/Darwin-Markov-model-figure.pdf")

# ╔═╡ 676cbaa5-449b-4111-9b7f-2e58420d8f35
md"""
## **Model description**

I built a simple evolutionary model to examine how plasmid copy number affects tetA-GFP copy number dynamics in the tetA-transposon system.

This model uses the standard population dynamics framework described by Edo Kussell and Stan Liebler in "Phenotypic Diversity, Population Growth, and Information in Fluctuating Environments" (Science 2005) and "Individual histories and selection in heterogeneous populations" (PNAS 2010).  

**Model Assumptions**  

The population is modeled as a distribution over tetA copy numbers, ranging from 1 (found on the chromosome) to $n$ (maximum plasmid copy number). We therefore represent the population as a vector 

$\mathbf{x}(t) = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ ... \\ x_n \end{pmatrix}$  

The sum of the $n$ entries of $\mathbf{x}(t)$ gives the total population size $N(t) = \sum_{i=1}^{n} \mathbf{x}_i(t)$.  

*Growth dynamics*  

For a given tetracycline concentration, each tetA copy number class has a growth rate, or *fitness*.  Let $r_i$ be be the fitness of tetA copy number class $i$. We assume that there is some optimal tetA copy number given some tetracycline concentration. We assume that $r_i$ can be positive or negative, with a single peak at the optimum tetA copy number for a given $[Tet]$ concentration, so a natural choice is to define the following quadratic fitness function:  

$r_i = r_{max} - \frac{(i - \alpha*[Tet])^2}{2\sigma^2}$  

Here, $r_{max}$ is the maximum growth rate, $i$ is the number of tetA copies in this strain, [Tet] is the antibiotic concentration, $\sigma$ is a free parameter that determines the width of the quadratic function, and $\alpha$ is a conversion factor with units of "tetA gene copies/[Tet]". Without loss of generality we assume $\alpha = 1$, so:

$r_i = r_{max} - \frac{(i - [Tet])^2}{2\sigma^2}$  

Then, the growth rates (fitnesses) for each subpopulation is a vector:

$\mathbf{r} = \begin{pmatrix} r_1 \\ r_2 \\ r_3 \\ ... \\ r_n \end{pmatrix}$  


We then define the diagonal matrix $\mathbf{G} = \begin{bmatrix}
r_1 & 0 & 0 & ... & 0 \\
0 & r_2 & 0 & ... & 0 \\
0 & 0 & r_3 & ... & 0 \\
... & ... & ... & ... & ... \\
0 & 0 & 0 & ... & r_n \\
\end{bmatrix}$.  

This represents how each subpopulation grows based on $\mathbf{r}$.

And the average population growth rate (mean population fitness) of the whole population is $\overline{r}(t) = \mathbf{r} \cdot \mathbf{x}(t)$.

*Phenotypic switching dynamics*

Following the Kussell and Leibler papers cited above, we assume that individuals can switch phenotypes according to a set of rates $s_{ij}$, which denotes the rate of switching from phenotype $j$ to $i$.  

We define $s_{jj} = -\sum_{i \ne j} s_{ij}$, so that the diagonal elements are the total rate of switching out of each phenotype.

Subpopulations switch phenotypes based on plasmid segregration during cell division.

We assume that the total plasmid copy number is fixed per cell. However, the cells of each subpopulation can gain or lose one plasmid containing the tetA transposon, unless all plasmids in the cell have the tetA transposon. In that case, cell division always results in daughter cells in which all plasmids have the tetA transposon. Mathematically, the subpopulation in which all plasmids have the transposon is an absorbing state of the Markov chain.  

We also assume that gain and loss rates are equal and assume detailed balance. Without loss of generality we let gain rate = loss rate = $k$. This results in the following tridiagonal matrix (6-dimensional case shown):  

$\mathbf{S} = \begin{bmatrix}
-k & k & 0 & 0 & 0 & 0 \\
k & -2k & k & 0 & 0 & 0 \\
0 & k & -2k & k & 0 & 0 \\
0 & 0 & k & -2k & k & 0 \\
0 & 0 & 0 & k & -2k & 0 \\
0 & 0 & 0 & 0 & k & 0 \\
\end{bmatrix}$.  

*Full dynamics*  

These dynamics are combined into a matrix $\mathbf{A} = \mathbf{G} + \mathbf{S}$,

so the full dynamics are modeled by the following matrix system of ODEs:  

$\frac{d\mathbf{x}}{dt} = \mathbf{A}\mathbf{x}(t)$.


This results in this form of the matrix $\mathbf{A}$ (6-dimensional case shown):

$\mathbf{A} = \begin{bmatrix}
r_1-k & k & 0 & 0 & 0 & 0 \\
k & r_2-2k & k & 0 & 0 & 0 \\
0 & k & r_3-2k & k & 0 & 0 \\
0 & 0 & k & r_4-2k & k & 0 \\
0 & 0 & 0 & k & r_5-2k & 0 \\
0 & 0 & 0 & 0 & k & r_6 \\
\end{bmatrix}$.  
 

We can generalize this model to the case where the plasmid copy number limits tetA copy number. In this case, higher tetA copy numbers beyond plasmid copy number cannot be reached.

Suppose the maximum plasmid copy number is 4. Then the truncated matrix $\mathbf{A'}$ is (for the 6-dimensional case):   

$\mathbf{A'} = \begin{bmatrix}
r_1-k & k & 0 & 0 & 0 & 0 \\
k & r_2-2k & k & 0 & 0 & 0 \\
0 & k & r_3-2k & 0 & 0 & 0 \\
0 & 0 & k & r_4 & 0 & 0 \\
0 & 0 & 0 & 0 & r_5 & 0 \\
0 & 0 & 0 & 0 & 0 & r_6 \\
\end{bmatrix}$.  

"""

# ╔═╡ 43ee91c7-8061-463d-ab97-e585cd7df0a9
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

# ╔═╡ dfb082a8-d421-4ea3-99ec-564ca067ed38
md"""
Using this notation, the phenotypic value of descendant j of individual i in the current generation is $\phi_i + \delta_{i,j}$. \
The mean phenotype of the descendants, $\overline{\phi'}$, is then: \
\
$\overline{\phi'} = \frac{\sum_{i=1}^N \sum_{j=1}^{W_{i}} (\phi_i + \delta_{i,j})}{\sum_{i=1}^N W_i}$. \
\
The numerator adds up the phenotypic values of each of the descendants by looking in term at each of the N individuals in the current generation, and then summing over the phenotypic values of each of their $W_i$ descendants. The denominator is simply the total number of descendants.
"""

# ╔═╡ cfe80654-c904-4a58-9cdf-31449da87828
md"""
We then use the following identities, which follow from the definition of the arithmetic mean:  \
\
$\sum_{j=1}^{W_{i}} \phi_i = W_{i}\phi_{i}$ \
$\sum_{j=1}^{W_{i}} \delta_{i,j} = W_{i}\overline{\delta_{i}}$ \
$\sum_{j=1}^N W_{i} = N\overline{W_{i}}$ \
"""

# ╔═╡ 1439fde6-4bc0-4359-86f7-6eb9d072d1c8
md"""
These allow us to rewrite the equation as: \
$\overline{\phi'} = \frac{1}{N\overline{W}}[\sum_{i=1}^N W_{i}\phi_{i} + \sum_{i=1}^N W_{i}\overline{\delta_{i}}]$, so, \
$\overline{\phi'} = \frac{1}{\overline{W}}[E(W\phi) + E(W\overline{\delta})]$
"""

# ╔═╡ 04d78bdc-bbd4-4e15-bc89-c73b6ab06a9a
md"""
Now, using the fact that $Cov(x,y) = E(xy) - \overline{x}\cdot\overline{y}$, \
we can substitute for $E(W\phi)$ to get: \
$\overline{\phi'} = \frac{1}{\overline{W}}[Cov(W_, \phi) + \overline{W}\cdot\overline{\phi} + E(W\overline{\delta})]$, so, \

$\overline{\phi'} = \frac{1}{\overline{W}}[Cov(W_, \phi) + E(W\overline{\delta})] + \overline{\phi}$. \

"""

# ╔═╡ 612c1242-d078-48a2-a7c5-f870c1b3672e
md""" 
Finally, subtracting $\overline{\phi}$ from both sides yields: \
\
$\Delta\overline{\phi} = \frac{1}{\overline{W}} [Cov(W, \phi) + E(W\overline{\delta})]$ \
\
This is *Price's theorem* (Price 1970).
"""

# ╔═╡ ba20d844-8719-43da-8831-02a61eeb247a
md""" 
**Remarks:** \
The $\frac{1}{\overline{W}}Cov(W, \phi)$ term represents the change due to differential survival and reproduction, encompassing selection and genetic drift. \
\
The $\frac{1}{\overline{W}}E(W\overline{\delta})$ term represents the change due to processes involved in reproduction, such as recombination, regression toward the mean phenotype, or selection at a lower level of organization (i.e. plasmids or other genetic elements that bias their own transmission into daughter cells, potentially at the expense of other, competing genetic elements, as in meiotic drive or CRISPR drives) \
\
Note that although we used absolute fitness values in our derivation of Price's theorem, we could just as well use relative fitness $w$, since any constant that multiplies all fitness values will appear in the numerator and denominator, and will thus cancel out.
"""

# ╔═╡ 30362c2d-d4d1-46c9-829b-d4f190ba3e51
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

# ╔═╡ bb7e3db9-6adb-4ab1-906d-98ca68065875
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

# ╔═╡ de6a80c0-138f-4b5f-8c25-e7c29c18d3e3
begin
	## Define global constants.
	MAX_TCN = 100 ## maximum attainable transposon copy number (TCN)
	MAX_PCN = 80 ## maximum attainable plasmid copy number (TCN)
	SIGMA = 15 ## set the width of the fitness function
	R_MAX = 1 ## set the max growth rate.
	K = 1/3 ## set the switching rate.
	TIMESPAN = 100.0
end

# ╔═╡ 80dd3108-c7aa-4f8b-9170-f61766bb18bc
function fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	fitness is a function of tetA copies and antibiotic concentration.
	for simplicity, define a quadratic function whose peak is shifted
	left or right based on Tet_conc.
	"""
	fitness = r_max - ((Tet_conc - tetA_copy_number)^2 / 2σ^2)
	return fitness
end

# ╔═╡ 773c3a96-841d-4802-ba55-995fb38e2d29
fitness_function(20.0, 0.0)

# ╔═╡ 97efc2ea-1814-4326-a8e8-a566d923feca
function gaussian_fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX)
	"""
	fitness is a function of tetA copies and antibiotic concentration.
	for simplicity, define a gaussian function whose peak is shifted
	left or right based on Tet_conc.
	"""
	fitness = exp(fitness_function(tetA_copy_number, Tet_conc, σ=SIGMA, r_max=R_MAX))
	return fitness
end

# ╔═╡ d3a7d259-da19-4cb0-b087-b5eba6c8c0d7
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

# ╔═╡ 67f9c177-02dd-4ea1-a477-adf77618db4f
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

# ╔═╡ e007814b-d92a-4189-ae4c-c05638ea6531
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

# ╔═╡ eb0e9882-80b4-4612-a466-8cfee96cd7b9
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

# ╔═╡ 5475908e-154e-49ec-b694-4efe759fa55b
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

# ╔═╡ 02af596b-6273-42b5-ae15-8a8b5588c12a
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

# ╔═╡ d4fda4d5-bddd-44f9-bd4a-4f3c4790cb44
function SwitchingTridiagonalMatrix(max_tetA_copy_number=MAX_TCN, plasmid_max_copies=MAX_PCN)
	""" We set a final absorbing state in the Markov chain,
	based on maximum plasmid copy number (like 3 for SC101, 2000 for pUC)."""
	n = max_tetA_copy_number
	superdiag = vcat([K for i in 2:n-1], 0)
	diag = vcat(-K, [-2K for i in 1:n-2], 1.0)
	subdiag = [K for i in 1:n-1]
	
	switching_matrix = Tridiagonal(subdiag, diag, superdiag)
	return(switching_matrix)
end

# ╔═╡ 023a8f11-e921-47a5-bef2-5b2c8a50efe6
function SwitchingTridiagonalMatrixWithAbsorbingState(max_tetA_copy_number=MAX_TCN, plasmid_max_copies=MAX_PCN)
	""" We set a final absorbing state in the Markov chain,
	based on maximum plasmid copy number (like 3 for SC101, 2000 for pUC)."""
	n = max_tetA_copy_number
	@assert plasmid_max_copies <= n "ERROR: plasmid_max_copies out of bounds"
	k = plasmid_max_copies ## shorthand notation.
	superdiag = [i < k ? K : 0.0 for i in 2:n]
	diag = vcat(-K, [-2K for i in 2:k-1], [0.0 for i in k:n])
	subdiag = [i < k ? K : 0.0 for i in 1:n-1]
	
	switching_matrix = Tridiagonal(subdiag, diag, superdiag)
	return(switching_matrix)
end

# ╔═╡ 8051d7a7-67fe-4c0c-b55d-0e0d1f424723
TetConcSlider = @bind TET_CONC Slider(0:100, default=20, show_value=true)

# ╔═╡ 71f13821-0c57-4320-b347-88963ed9be39
function SelectionDiagonalMatrix(max_tetA_copy_number=MAX_TCN, Tet_conc=TET_CONC)
	## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
	## with a minimum of 1 copy (on the chromosome).
	tetA_classes = collect(1:max_tetA_copy_number)
	## get the fitness for each fitness class (defined by tetA copy number).
	growth_rate_vec = fitness_function.(tetA_classes, Tet_conc)
	diagonal_growth_matrix = Diagonal(growth_rate_vec)
	return(diagonal_growth_matrix)
end

# ╔═╡ d0bd89fe-234c-426d-9d72-feb708c19e58
function MutSelMatrix1(max_tetA_copy_number=MAX_TCN, plasmid_max_copies=MAX_PCN, Tet_conc=TET_CONC)
	mutsel_matrix = SelectionDiagonalMatrix(max_tetA_copy_number, Tet_conc) + SwitchingTridiagonalMatrix(max_tetA_copy_number, plasmid_max_copies) 
	return(mutsel_matrix)
end

# ╔═╡ 4b630f32-06e1-40a0-af3f-59febaef5b93
function MutSelMatrix2(max_tetA_copy_number=MAX_TCN, plasmid_max_copies=MAX_PCN, Tet_conc=TET_CONC)
	mutsel_matrix = SelectionDiagonalMatrix(max_tetA_copy_number, Tet_conc) + SwitchingTridiagonalMatrixWithAbsorbingState(max_tetA_copy_number, plasmid_max_copies)  
	return(mutsel_matrix)
end

# ╔═╡ e603fa00-4f48-4c8c-9382-6435374bf20a
begin
	xvec = collect(1:MAX_TCN)
	fitnesses = fitness_function.(xvec, TET_CONC)
	plot(xvec, fitnesses, label="Fitness", xlabel="tetA copy number", ylabel="Fitness")
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ 21ab9391-bd06-4629-ab3b-88487713d6fe
function odefunc1(du, u, p, t)
	## Define the ODE system
	A = MutSelMatrix1()
    du .= A * u
end

# ╔═╡ 9e2ff9ae-76ee-4dfe-8ac0-296fa9dd160a
function odefunc2(du, u, p, t)
	## Define the ODE system
	A = MutSelMatrix2()
    du .= A * u
end

# ╔═╡ a8267ba1-a248-46bb-a69d-6494210bc7fa
begin
	## parameters used across time course simulations.
	initial_pop_vec = zeros(MAX_TCN)
	## initialize the population with one cell with 1 tetA copy.
	initial_pop_vec[1] = 1
	## Define the time span
	tspan = (0.0, TIMESPAN)  # Replace with your desired time span
end

# ╔═╡ 48af8690-f8ac-4a99-bb25-672fc4fdd40f
md""" ## model population dynamics under constant [Tet]."""

# ╔═╡ 73758488-d007-4af5-a245-508cc232af7e
begin
	## run the simulation for constant [Tet] concentration.
	
	## Create an ODEProblem
	prob = ODEProblem(odefunc1, initial_pop_vec, tspan)

	## Solve the ODE system
	sol = solve(prob, Tsit5())

	## Extract the solution
	final_t1 = sol.t
	final_x1 = sol.u
end

# ╔═╡ 9897d389-bd11-40bd-ab8a-3d2595458eec
begin
	## Extract the subpopulation timecourses.
	## Example for subpopulation 1: x1 = [u[1] for u in final_x]
	## Example for subpopulation 2: x2 = [u[2] for u in final_x]

	# Initialize an empty vector to store vectors
	subpopulation_time_courses = Vector{Vector{Float64}}()
	num_subpopulations = length(initial_pop_vec)
	for i in 1:num_subpopulations
		cur_time_course = [u[i] for u in final_x1]
		push!(subpopulation_time_courses, cur_time_course)
	end
end

# ╔═╡ 73f12bf5-dce2-4fbe-8aa0-9686874b69a0
let
	## plot the first component over time
	x_1 = subpopulation_time_courses[1]
	timecourse_plot = plot(final_t1, x_1, label="1", xlabel="Time", ylabel="Subpopulation size")
	
	## plot the remaining components over time.
	for i in 2:num_subpopulations
		x_i = subpopulation_time_courses[i]
		plot!(timecourse_plot, final_t1, x_i, label=string(i))
	end
	timecourse_plot
end

# ╔═╡ 1462bf95-6e84-453e-bd60-5068cdc48b99
begin ## let's make a matrix saving the population dynamics.

	result_matrix = [] ## Initialize an empty matrix

	for i in 1:length(final_x1)
		total_N = sum(final_x1[i])
		cur_frequency_vec = final_x1[i]/total_N
		##concatenate the state vector horizontally into the result_matrix
		result_matrix = push!(result_matrix, cur_frequency_vec)
	end
end

# ╔═╡ a783d983-03d1-42cc-a3b5-fdd755734dfc
time_step_slider = @bind cur_timestep Slider(1:length(final_x1), default=1, show_value=true )

# ╔═╡ 965b6b02-2866-4350-adb2-85c38784542f
let
	xvec = collect(1:length(final_x1[1]))
	plot(xvec, result_matrix[cur_timestep], label="", xlabel="tetA copy number", ylabel="Current population distribution", ylims=[0,1])
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ 922f75f1-3697-44a8-a0b1-039754718023
md""" ##### calculate mean fitness in the constant [Tet] population."""

# ╔═╡ f8f2c684-b939-49b4-9aba-ef844c63d18e
begin
mean_fitness_vec1 = [] 
	for i in 1:length(result_matrix)
		cur_pop_vec = result_matrix[i]
		cur_mean_fitness = calc_mean_fitness(cur_pop_vec, TET_CONC)
		append!(mean_fitness_vec1, cur_mean_fitness)
	end
end

# ╔═╡ 5ff4baaa-37f6-4559-a461-2d71804b39ff
md""" ##### calculate mean copy number in the constant [Tet] population."""

# ╔═╡ e08cfdd5-e5bf-4e35-8586-60215c45c29e
begin
	mean_copy_num_vec1 = []
	for i in 1:length(final_x1)
		cur_pop_vec = result_matrix[i]
		cur_mean_copy_number = calc_mean_tetA_copy_number(cur_pop_vec)
		append!(mean_copy_num_vec1, cur_mean_copy_number)
	end
end

# ╔═╡ d16a8dd4-3244-4487-a2ae-39f1e1c25e5e
let
	plot(final_t1, mean_copy_num_vec1, label="Mean tetA copy number", xlabel="Time", ylabel="Mean tetA copy number")
end

# ╔═╡ 74320ce9-74e1-4e34-8e96-e29220ef97d7
md""" ##### calculate fitness variance in the constant [Tet] population."""

# ╔═╡ e4c1987c-79fe-4787-b555-810ce216676e
begin
fitness_variance_vec1 = [] 
	for i in 1:length(final_x1)
		cur_pop_vec = result_matrix[i]
		cur_fitness_variance = calc_fitness_variance(cur_pop_vec, TET_CONC)
		append!(fitness_variance_vec1, cur_fitness_variance)
	end
end

# ╔═╡ eab8e464-dd8a-48d5-bfde-52a11d3b6d58
let
	plot(final_t1, fitness_variance_vec1, label="Fitness variance", xlabel="Time", ylabel="Fitness variance")
end

# ╔═╡ ca4a41bb-8a43-4f2e-b8c0-c68d34d0d868
let
	plot(final_t1, fitness_variance_vec1, label="Fitness variance")
	plot!(final_t1, mean_fitness_vec1, label="Mean fitness", xlabel="Time", ylabel="Fitness mean and variance")
end

# ╔═╡ 86a1e012-58ff-4c1e-9fd9-d849e572c7ee
md""" ##### calculate rate of mean fitness increase in the constant [Tet] population."""

# ╔═╡ fcb28d74-cbd1-4037-b2d5-71176c2c936d
## append a zero to the front of the difference vector.
d_mean_fitness_vec1 = [0; diff(mean_fitness_vec1)]

# ╔═╡ 4627a6b9-7607-4941-bd4b-ce48a72152a2
let
	plot(final_t1, fitness_variance_vec1, label="Fitness variance")
	plot!(final_t1, d_mean_fitness_vec1, label="Rate of change of mean fitness", xlabel="Time", ylabel="Fitness variance and rate of mean fitness change")
end

# ╔═╡ 3e1e527f-8f62-40b1-8b7e-02e4815c3d61
md""" ##### calculate rate of copy number change in the constant [Tet] population."""

# ╔═╡ c596d4a7-03e4-4169-af09-3b8175fd2674
## append a zero to the front of the difference vector.
d_mean_copy_num_vec1 = [0; diff(mean_copy_num_vec1)]

# ╔═╡ bc920a63-df48-490d-ab8b-ef5e31c83168
md""" ##### calculate tetA copy number--fitness covariance in the constant [Tet] population. """ 

# ╔═╡ 3f7f2cdf-6826-4bc9-ad85-393de803e2d0
begin
	copy_num_covariance_vec1 = []
	for i in 1:length(final_x1)
		cur_pop_vec = result_matrix[i]
		cur_copy_num_covariance = calc_tetA_copy_number_fitness_covariance(cur_pop_vec, TET_CONC)
		
		append!(copy_num_covariance_vec1, cur_copy_num_covariance)
	end
end

# ╔═╡ 78a55b31-9a12-48d8-ae19-368a9aa2090b
let
	plot(final_t1, copy_num_covariance_vec1, label="tetA copy number fitness covariance")
	plot!(final_t1, d_mean_copy_num_vec1, label="Rate of change of mean tetA copy number", xlabel="Time", ylabel="tetA copy number derivative and variance")
end

# ╔═╡ 205a5b06-c3af-4241-a248-0fb59658c402
d_mean_copy_num_vec1

# ╔═╡ 801861a2-416b-4f08-be22-3d3cac5e066a
mean_copy_num_vec1

# ╔═╡ d1dd08d6-6216-4433-9ce3-f83ea63b62dc
md""" ## model the Tet pulse conditions of the Darwin experiment."""

# ╔═╡ 0c1b9281-e87e-4017-8ffa-5b0bde454c1b
## We model [Tet] pulses over time by dividing the current time by 100, and
## setting [Tet] on if in the first half.
TetPulseFunction = t -> t < 80 ? TET_CONC : 0

# ╔═╡ ade78925-2e39-4df8-a361-06c246dd508b
function odefunc3(du, u, p, t)
	## Define the ODE system.

	## p is a TetPulseFunction of time that is passed in as a parameter.
	cur_tet_conc = p(t)
	A = MutSelMatrix1(MAX_TCN, MAX_PCN, cur_tet_conc)
    du .= A * u
end

# ╔═╡ c938a144-9f81-45e3-a012-e52134d49c7d
begin
	## run the simulation for pulses of [Tet] concentration.
	
	## Create an ODEProblem
	prob2 = ODEProblem(odefunc3, initial_pop_vec, tspan, TetPulseFunction)

	## Solve the ODE system
	sol2 = solve(prob2, Tsit5())

	## Extract the solution
	final_t2 = sol2.t
	final_x2 = sol2.u
end

# ╔═╡ 846f1cb1-cb33-4052-8fcd-9a4b3b63f5cd
begin ## let's make a matrix saving the population dynamics.

	pulse_result_matrix = [] ## Initialize an empty matrix

	for i in 1:length(final_x2)
		total_N = sum(final_x2[i])
		cur_frequency_vec = final_x2[i]/total_N
		##concatenate the state vector horizontally into the pulse_result_matrix
		pulse_result_matrix = push!(pulse_result_matrix, cur_frequency_vec)
	end
end

# ╔═╡ f987af28-9868-4ff4-8734-8ea10a4458c7
pulse_time_step_slider = @bind pulse_cur_timestep Slider(1:length(final_x2), default=1, show_value=true)

# ╔═╡ 4f2f25ab-e5c7-4d7a-a2c9-9d366b523c9b
let
	xvec2 = collect(1:length(final_x2[1]))
	plot(xvec2, pulse_result_matrix[pulse_cur_timestep], label="", xlabel="tetA copy number", ylabel="Current population distribution", ylims=[0,1])
	# Add a vertical dashed line at x = TET_CONC
	vline!([TET_CONC], linestyle=:dash, label="TET_CONC")
end

# ╔═╡ e392d94e-5e8c-46d3-8b9a-9fcbb0d9d1b5
md""" ### plot the antibiotic pulse regime over time."""

# ╔═╡ 4d47b937-4a39-41b0-8a70-21a4b97c529e
let
	tvec = collect(1:TIMESPAN)
	plot(tvec, [TetPulseFunction(t) for t in tvec],label="Tet pulse regime", xlabel="Time", ylabel="Tet concentration")
end

# ╔═╡ 98d9548b-fa2e-4cf6-8fa0-53ae2e4d5c99
md""" ##### calculate mean fitness in the population."""

# ╔═╡ f0a8f82b-e82b-4d8e-b460-3d9ab494cb22
begin
mean_fitness_vec2 = [] 
	for i in 1:length(pulse_result_matrix)
		cur_pop_vec = pulse_result_matrix[i]
		cur_tet_conc = TetPulseFunction(final_t2[i])
		cur_mean_fitness = calc_mean_fitness(cur_pop_vec, cur_tet_conc)
		append!(mean_fitness_vec2, cur_mean_fitness)
	end
end

# ╔═╡ cbc3a1b7-4114-4f4c-a9cc-ab9f52ec4453
let
	plot(final_t2, mean_fitness_vec2, label="Mean fitness", xlabel="Time", ylabel="Mean fitness")
end

# ╔═╡ 39fd8d0b-60e1-485e-920d-9fc2d4402bcd
md""" ##### calculate mean copy number in the population."""

# ╔═╡ 2eea090c-d8e7-4efa-a575-96a9a3daf9e2
begin
	mean_copy_num_vec2 = []
	for i in 1:length(final_x2)
		cur_pop_vec = pulse_result_matrix[i]
		cur_mean_copy_number = calc_mean_tetA_copy_number(cur_pop_vec)
		append!(mean_copy_num_vec2, cur_mean_copy_number)
	end
end

# ╔═╡ 5d868baf-6e1e-4cbe-ad62-67344f6d0fe7
let
	plot(final_t2, mean_copy_num_vec2, label="Mean tetA copy number", xlabel="Time", ylabel="Mean tetA copy number")
end

# ╔═╡ cd2b0606-26ac-4a20-9f7d-dd5cd7bedfcb
md""" ##### calculate fitness variance in the population."""

# ╔═╡ 234d1cc2-6be6-4712-a62d-ef8d6b3f4f14
begin
fitness_variance_vec2 = [] 
	for i in 1:length(final_x2)
		cur_pop_vec = pulse_result_matrix[i]
		cur_tet_conc = TET_CONC##tet_pulse_vec[i]
		cur_fitness_variance = calc_fitness_variance(cur_pop_vec, cur_tet_conc)
		append!(fitness_variance_vec2, cur_fitness_variance)
	end
end

# ╔═╡ 92272b01-1866-4125-b334-e8a3f1de2f4e
let
	plot(final_t2, fitness_variance_vec2, label="Fitness variance", xlabel="Time", ylabel="Fitness variance")
end

# ╔═╡ 2c878d03-cf1f-437b-a50c-a31776fac4c9
let
	plot(final_t2, fitness_variance_vec2, label="Fitness variance")
	plot!(final_t2, mean_fitness_vec2, label="Mean fitness", xlabel="Time", ylabel="Fitness mean and variance")
end

# ╔═╡ d1379ccf-979c-4fc5-a568-272aab2ab025
md""" ##### calculate rate of mean fitness increase in the population."""

# ╔═╡ 9fe8cd10-bd4b-4b4c-8534-59fad9a89e89
## append a zero to the front of the difference vector.
d_mean_fitness_vec2 = [0; diff(mean_fitness_vec2)]

# ╔═╡ 2503b3f3-e53a-4e43-a5b7-f2c2bfff1654
let
	plot(final_t2, fitness_variance_vec2, label="Fitness variance")
	plot!(final_t2, d_mean_fitness_vec2, label="Rate of change of mean fitness", xlabel="Time", ylabel="Fitness variance and rate of mean fitness change")
end

# ╔═╡ 5fa6f2d7-fa09-41c8-b08c-561907be01ae
md""" ##### calculate rate of copy number change in the population."""

# ╔═╡ 47c37b1b-2692-4e37-84cd-ae43d3195dda
## append a zero to the front of the difference vector.
d_mean_copy_num_vec2 = [0; diff(mean_copy_num_vec2)]

# ╔═╡ 1b2fa624-7afc-48bb-b924-e3638e8862c3
md""" ##### calculate tetA copy number--fitness covariance. """ 

# ╔═╡ 7442c9af-41d9-48fe-bb15-ec1264a89592
begin
	copy_num_covariance_vec = []
	for i in 1:length(final_x2)
		cur_pop_vec = pulse_result_matrix[i]
		cur_Tet_conc = TET_CONC
		## KEY MODELING ASSUMPTION: the index is the number of tetA copies,
		## with a minimum of 1 copy (on the chromosome).
		cur_copy_num_covariance = calc_tetA_copy_number_fitness_covariance(cur_pop_vec, cur_Tet_conc)
		
		append!(copy_num_covariance_vec, cur_copy_num_covariance)
	end
end

# ╔═╡ a5981dc1-08c1-4006-a481-de3c55a2bc51
let
	plot(final_t2, copy_num_covariance_vec, label="tetA copy number-fitness covariance", xlabel="Time", ylabel="tetA copy number--fitness covariance")
end

# ╔═╡ 4a832714-9009-45a3-a77c-75819fe8960f
let
	plot(final_t2, copy_num_covariance_vec, label="tetA copy number fitness covariance")
	plot!(final_t2, mean_copy_num_vec2, label="mean tetA copy number", xlabel="Time", ylabel="tetA copy number mean and variance")
end

# ╔═╡ 75d40406-5867-42d6-8964-3776cc3ac8f0
let
	plot(final_t2, copy_num_covariance_vec, label="tetA copy number fitness covariance")
	plot!(final_t2, d_mean_copy_num_vec2, label="Rate of change of mean tetA copy number", xlabel="Time", ylabel="tetA copy number derivative and variance")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DifferentialEquations = "~7.11.0"
FileIO = "~1.16.1"
Plots = "~1.39.0"
PlutoUI = "~0.7.53"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "eba41f5071cdf2306a22f0a782dfd483cc97cabd"

[[deps.ADTypes]]
git-tree-sha1 = "332e5d7baeff8497b923b730b994fa480601efc7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.5"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Test"]
git-tree-sha1 = "a7055b939deae2455aa8a67491e034f735dd08d3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.33"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "02f731463748db57cc2ebfbd9fbc9ce8280d3433"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.1"
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
git-tree-sha1 = "247efbccf92448be332d154d6ca56b9fcdd93c31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.6.1"

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
git-tree-sha1 = "af43df5704827c8618afd36eb56fcab20d3041ee"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.4.3"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "67bcff3f50026b6fa952721525d3a04f0570d432"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.2.1"
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
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "ConcreteStructs", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools", "Tricks", "TruncatedStacktraces", "UnPack"]
git-tree-sha1 = "8a19e2457da8a7e5ae54ee9479885738d8fd926b"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.4.0"

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
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

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
git-tree-sha1 = "df712c77bb43b37ea966feb72cb2e92d51a3face"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.43.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces"]
git-tree-sha1 = "309efb205c30d43b595466283bbecf2769283e22"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.141.0"

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
git-tree-sha1 = "4e4de57a0ac47b2f20aae62f132355b058e9f0cd"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.34.0"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "57ed4597a309c5b2a10cab5f9813adcb78f92117"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.19.0"

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
git-tree-sha1 = "19a5b6314715139ddefea4108a105bb9b90dc4fb"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "5225c965635d8c21168e32a12954675e7bea1151"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.10"

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
git-tree-sha1 = "a6c00f894f24460379cb7136633cef54ac9f6f4a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.103"

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
git-tree-sha1 = "ab81396e4e7b61f5590db02fa1c17fae4f16d7ab"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.6.3"

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

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "35f0c0f345bff2c6d636f95fdb136323b5a796ef"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.7.0"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

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
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

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

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

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
git-tree-sha1 = "9fb0b890adab1c0a4a475d4210d51f228bfc250d"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.6"

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
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "3de1d557e382cad270d921fbc22351f5628e7b1f"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.8.0"
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
git-tree-sha1 = "17e462054b42dcdda73e9a9ba0c67754170c88ae"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.4"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

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

[[deps.LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

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
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

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
deps = ["ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "EnzymeCore", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "UnPack"]
git-tree-sha1 = "051943b8b8e81c548e9d099d6eb3d3ed23093c35"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.20.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveEnzymeExt = "Enzyme"
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
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

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
version = "2022.10.11"

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
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "EnumX", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LineSearches", "LinearAlgebra", "LinearSolve", "PrecompileTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "6166ccd8f79c93c636ca61ab4cd18f555932563d"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "2.8.2"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

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
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SciMLOperators", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "42d0b4515472a25a3b47be228c03ec842bdc9d49"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.59.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4e5be6bb265d33669f98eb55d2a57addd1eeb72c"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.30"

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
git-tree-sha1 = "a935806434c9d4c506ba941871b327b96d41f2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.0"

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
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

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
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "1d05623b5952aed1307bf8b43bec8b8d1ef94b6e"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.5"

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

[[deps.QuasiMonteCarlo]]
deps = ["Accessors", "ConcreteStructs", "LatticeRules", "LinearAlgebra", "Primes", "Random", "Requires", "Sobol", "StatsBase"]
git-tree-sha1 = "cc086f8485bce77b6187141e1413c3b55f9a4341"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.3.3"
weakdeps = ["Distributions"]

    [deps.QuasiMonteCarlo.extensions]
    QuasiMonteCarloDistributionsExt = "Distributions"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "552f30e847641591ba3f39fd1bed559b9deb0ef3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.1"

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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
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
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "QuasiMonteCarlo", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "164773badb9ee8c62af2ff1a7778fd4867142a07"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.9.0"

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

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "765b788339abd7d983618c09cfc0192e2b6b15fd"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.9"

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
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "69b1a53374dd14d7c165d98cb646aeb5f36f8d07"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.25"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveNNlibExt = "NNlib"

    [deps.SimpleNonlinearSolve.weakdeps]
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "07272c80c278947baca092df0a01da4a10622ad5"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.13.0"

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
git-tree-sha1 = "03fec6800a986d191f64f5c0996b59ed526eda25"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.1"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "5ef59aea6f18c25168842bded46b16662141ab87"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.7.0"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

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
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "2ca69f4be3294e4cd987d83d6019037d420d9fc1"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.16.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "7a71f1e67cbcfcd5387707e6621431d1afff62a9"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.63.2"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "e7dd250422df290cee14960c1ee144b44ac3dd77"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.1"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "71dc65a2d7decdde5500299c9b04309e0138d1b4"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.20.1"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

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

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

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
git-tree-sha1 = "242982d62ff0d1671e9029b52743062739255c7e"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.18.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

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
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

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
git-tree-sha1 = "da69178aacc095066bad1f69d2f59a60a1dd8ad1"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.0+0"

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
version = "1.2.13+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "47cf33e62e138b920039e8ff9f9841aafe1b733e"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.35.1+0"

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
version = "5.8.0+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

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
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
# ╠═bb669570-7cb6-11ee-00f9-95e5d1f0bf9d
# ╠═96e88c96-ffa3-440b-a5c4-df380f1ee881
# ╠═1a7794a3-ee9b-4b0c-8ae0-3930a63f5e54
# ╠═a66ee148-7cf2-43a4-9b14-407eb91da8ea
# ╠═6d0b7c31-69c0-4f9c-9f7e-1cb43502f58f
# ╟─676cbaa5-449b-4111-9b7f-2e58420d8f35
# ╟─43ee91c7-8061-463d-ab97-e585cd7df0a9
# ╟─dfb082a8-d421-4ea3-99ec-564ca067ed38
# ╟─cfe80654-c904-4a58-9cdf-31449da87828
# ╟─1439fde6-4bc0-4359-86f7-6eb9d072d1c8
# ╟─04d78bdc-bbd4-4e15-bc89-c73b6ab06a9a
# ╟─612c1242-d078-48a2-a7c5-f870c1b3672e
# ╟─ba20d844-8719-43da-8831-02a61eeb247a
# ╟─30362c2d-d4d1-46c9-829b-d4f190ba3e51
# ╟─bb7e3db9-6adb-4ab1-906d-98ca68065875
# ╠═de6a80c0-138f-4b5f-8c25-e7c29c18d3e3
# ╠═80dd3108-c7aa-4f8b-9170-f61766bb18bc
# ╠═773c3a96-841d-4802-ba55-995fb38e2d29
# ╠═97efc2ea-1814-4326-a8e8-a566d923feca
# ╠═d3a7d259-da19-4cb0-b087-b5eba6c8c0d7
# ╠═67f9c177-02dd-4ea1-a477-adf77618db4f
# ╠═e007814b-d92a-4189-ae4c-c05638ea6531
# ╠═eb0e9882-80b4-4612-a466-8cfee96cd7b9
# ╠═5475908e-154e-49ec-b694-4efe759fa55b
# ╠═02af596b-6273-42b5-ae15-8a8b5588c12a
# ╠═d4fda4d5-bddd-44f9-bd4a-4f3c4790cb44
# ╠═023a8f11-e921-47a5-bef2-5b2c8a50efe6
# ╠═71f13821-0c57-4320-b347-88963ed9be39
# ╠═d0bd89fe-234c-426d-9d72-feb708c19e58
# ╠═4b630f32-06e1-40a0-af3f-59febaef5b93
# ╠═8051d7a7-67fe-4c0c-b55d-0e0d1f424723
# ╠═e603fa00-4f48-4c8c-9382-6435374bf20a
# ╠═21ab9391-bd06-4629-ab3b-88487713d6fe
# ╠═9e2ff9ae-76ee-4dfe-8ac0-296fa9dd160a
# ╠═a8267ba1-a248-46bb-a69d-6494210bc7fa
# ╠═48af8690-f8ac-4a99-bb25-672fc4fdd40f
# ╠═73758488-d007-4af5-a245-508cc232af7e
# ╠═9897d389-bd11-40bd-ab8a-3d2595458eec
# ╠═73f12bf5-dce2-4fbe-8aa0-9686874b69a0
# ╠═1462bf95-6e84-453e-bd60-5068cdc48b99
# ╠═a783d983-03d1-42cc-a3b5-fdd755734dfc
# ╠═965b6b02-2866-4350-adb2-85c38784542f
# ╠═922f75f1-3697-44a8-a0b1-039754718023
# ╠═f8f2c684-b939-49b4-9aba-ef844c63d18e
# ╠═5ff4baaa-37f6-4559-a461-2d71804b39ff
# ╠═e08cfdd5-e5bf-4e35-8586-60215c45c29e
# ╠═d16a8dd4-3244-4487-a2ae-39f1e1c25e5e
# ╠═74320ce9-74e1-4e34-8e96-e29220ef97d7
# ╠═e4c1987c-79fe-4787-b555-810ce216676e
# ╠═eab8e464-dd8a-48d5-bfde-52a11d3b6d58
# ╠═ca4a41bb-8a43-4f2e-b8c0-c68d34d0d868
# ╠═86a1e012-58ff-4c1e-9fd9-d849e572c7ee
# ╠═fcb28d74-cbd1-4037-b2d5-71176c2c936d
# ╠═4627a6b9-7607-4941-bd4b-ce48a72152a2
# ╠═3e1e527f-8f62-40b1-8b7e-02e4815c3d61
# ╠═c596d4a7-03e4-4169-af09-3b8175fd2674
# ╠═bc920a63-df48-490d-ab8b-ef5e31c83168
# ╠═3f7f2cdf-6826-4bc9-ad85-393de803e2d0
# ╠═78a55b31-9a12-48d8-ae19-368a9aa2090b
# ╠═205a5b06-c3af-4241-a248-0fb59658c402
# ╠═801861a2-416b-4f08-be22-3d3cac5e066a
# ╠═d1dd08d6-6216-4433-9ce3-f83ea63b62dc
# ╠═0c1b9281-e87e-4017-8ffa-5b0bde454c1b
# ╠═ade78925-2e39-4df8-a361-06c246dd508b
# ╠═c938a144-9f81-45e3-a012-e52134d49c7d
# ╠═846f1cb1-cb33-4052-8fcd-9a4b3b63f5cd
# ╠═f987af28-9868-4ff4-8734-8ea10a4458c7
# ╠═4f2f25ab-e5c7-4d7a-a2c9-9d366b523c9b
# ╠═e392d94e-5e8c-46d3-8b9a-9fcbb0d9d1b5
# ╠═4d47b937-4a39-41b0-8a70-21a4b97c529e
# ╠═98d9548b-fa2e-4cf6-8fa0-53ae2e4d5c99
# ╠═f0a8f82b-e82b-4d8e-b460-3d9ab494cb22
# ╠═cbc3a1b7-4114-4f4c-a9cc-ab9f52ec4453
# ╠═39fd8d0b-60e1-485e-920d-9fc2d4402bcd
# ╠═2eea090c-d8e7-4efa-a575-96a9a3daf9e2
# ╠═5d868baf-6e1e-4cbe-ad62-67344f6d0fe7
# ╠═cd2b0606-26ac-4a20-9f7d-dd5cd7bedfcb
# ╠═234d1cc2-6be6-4712-a62d-ef8d6b3f4f14
# ╠═92272b01-1866-4125-b334-e8a3f1de2f4e
# ╠═2c878d03-cf1f-437b-a50c-a31776fac4c9
# ╠═d1379ccf-979c-4fc5-a568-272aab2ab025
# ╠═9fe8cd10-bd4b-4b4c-8534-59fad9a89e89
# ╠═2503b3f3-e53a-4e43-a5b7-f2c2bfff1654
# ╠═5fa6f2d7-fa09-41c8-b08c-561907be01ae
# ╠═47c37b1b-2692-4e37-84cd-ae43d3195dda
# ╠═1b2fa624-7afc-48bb-b924-e3638e8862c3
# ╠═7442c9af-41d9-48fe-bb15-ec1264a89592
# ╠═a5981dc1-08c1-4006-a481-de3c55a2bc51
# ╠═4a832714-9009-45a3-a77c-75819fe8960f
# ╠═75d40406-5867-42d6-8964-3776cc3ac8f0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

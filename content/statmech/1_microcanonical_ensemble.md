---
numbering:
  title:
    offset: 0

kernelspec:
  name: python3
  display_name: 'Python 3'
---

(ch-microcanonical)=
# The Microcanonical Ensemble

Statistical mechanics is the bridge between the microscopic world — atoms bouncing, colliding, and interacting — and the macroscopic world of temperature, pressure, and entropy. When you touch a hot cup of tea, you are not detecting the kinetic energy of any single molecule; you are sensing the collective statistical average of some $10^{23}$ molecules in frantic thermal motion.

The central question of statistical mechanics is: **given only the macroscopic constraints on a system (total energy, volume, number of particles), what can we predict about its thermodynamic behaviour?**

The answer rests on a profound and elegant idea: *we do not need to track every single particle*. Instead, we ask what fraction of all possible microscopic configurations are consistent with our macroscopic constraints, and we treat all such configurations as equally likely. This is the **fundamental postulate of statistical mechanics**, and the mathematical framework built around it is the **ensemble theory**.

In this chapter we develop the first and most fundamental ensemble — the **microcanonical ensemble** — which describes a system that is completely isolated from its surroundings. We will:

1. Define the **microcanonical ensemble** and its probability distribution over microstates.
2. Introduce the **phase space** and the volume of accessible microstates.
3. Postulate and derive properties of the **Boltzmann entropy** $S = k_B \ln \Gamma$.
4. Define **temperature**, **pressure**, and **chemical potential** through entropy.
5. Give a physical interpretation of temperature.
6. Discuss the **validity** of the statistical description.
7. Derive the **first law of thermodynamics** from statistical mechanics.
8. Discuss the **statistical interpretation of entropy** and the second law.
9. Apply the full machinery to the **classical ideal gas**, deriving the Sackur-Tetrode equation and the equations of state.
10. Resolve the **Gibbs paradox** and understand the role of indistinguishability.

(ch_mce_s_definition)=
## Definition of the Microcanonical Ensemble

(ch_mce_ss_isolated_system)=
### The Isolated System

We consider a system that is **completely isolated** from its surroundings. This means:

- It contains a **fixed number of particles** $N$.
- It has a **fixed total energy** $E$ (to within a small uncertainty $\Delta$).
- It occupies a **fixed volume** $V$.

Such a system exchanges **neither energy nor particles** with anything outside. It is described entirely by the triple $(E, V, N)$.

The microscopic state of such a system — its **microstate** — is specified by the positions and momenta of all $N$ particles:

$$\text{microstate} = \{(\mathbf{q}_1, \mathbf{p}_1), (\mathbf{q}_2, \mathbf{p}_2), \ldots, (\mathbf{q}_N, \mathbf{p}_N)\}$$

For a system in three dimensions, this is a point in a **$6N$-dimensional phase space** $\Gamma$, with coordinates $(q_1, \ldots, q_{3N}, p_1, \ldots, p_{3N})$.

(ch_mce_ss_liouville)=
### Liouville's Theorem and Stationarity

The time evolution of the phase-space density $\rho(\mathbf{q}, \mathbf{p}, t)$ is governed by **Liouville's theorem**:

$$\frac{d\rho}{dt} = \frac{\partial \rho}{\partial t} + \{\rho, H\} = 0$$

where $\{A, B\} = \sum_i \left(\frac{\partial A}{\partial q_i}\frac{\partial B}{\partial p_i} - \frac{\partial A}{\partial p_i}\frac{\partial B}{\partial q_i}\right)$ is the Poisson bracket and $H(\mathbf{q}, \mathbf{p})$ is the Hamiltonian. This equation says that phase-space density is **conserved** along trajectories — it flows like an incompressible fluid in phase space.

For a stationary (equilibrium) distribution we require $\partial_t \rho = 0$, so $\{H, \rho\} = 0$. The general solution is:

$$\rho(\mathbf{q}, \mathbf{p}) = \rho\bigl(H(\mathbf{q}, \mathbf{p})\bigr)$$

meaning the equilibrium distribution can only depend on the coordinates through the Hamiltonian (the energy). For an isolated system, this leads immediately to a **uniform distribution** over the energy surface.

(ch_mce_ss_equal_probability)=
### The Fundamental Postulate

```{important}
**Postulate of Equal A Priori Probabilities (Fundamental Postulate):** In an isolated system at equilibrium, all accessible microstates — those consistent with the macroscopic constraints $E$, $V$, $N$ — are *equally probable*.
```

This is not something we can derive from first principles; it is a postulate whose validity is confirmed by the extraordinary agreement of its consequences with experiment. It is the cornerstone of all of statistical mechanics.

(ch_mce_ss_energy_shell)=
### The Energy Shell and the Probability Distribution

In practice, we allow the energy to lie in a narrow shell $[E, E + \Delta]$ rather than at a precise value $E$ (both for mathematical convenience and because energy cannot be specified with infinite precision in any real situation).

The **microcanonical ensemble** is then defined by the probability density:

$$\boxed{\rho(\mathbf{q}, \mathbf{p}) = \begin{cases} \dfrac{1}{\Gamma(E, V, N)} & \text{if } E < H(\mathbf{q}, \mathbf{p}) < E + \Delta \\ 0 & \text{otherwise} \end{cases}}$$

where

$$\Gamma(E, V, N) = \int\!\!\!\int_{E < H(\mathbf{q},\mathbf{p}) < E+\Delta} d^{3N}q\, d^{3N}p \tag{1}$$

is the **phase-space volume** of the energy shell. The dependence on the spatial volume $V$ enters through the limits of the position integrals $dq_i$, which are bounded to the physical volume of the container.

```{note}
The distribution $\rho(\mathbf{q},\mathbf{p})$ is uniform over the energy shell and zero outside it. This is the most natural and democratic choice: no microstate within the accessible region is favoured over any other.
```

The expectation value of any classical observable $O(\mathbf{q}, \mathbf{p})$ is:

$$\langle O \rangle = \int d^{3N}q\, d^{3N}p\; \rho(\mathbf{q}, \mathbf{p})\, O(\mathbf{q}, \mathbf{p}) = \frac{1}{\Gamma(E,V,N)}\int\!\!\!\int_{E<H<E+\Delta} d^{3N}q\, d^{3N}p\; O(\mathbf{q},\mathbf{p})$$

(ch_mce_s_phase_space)=
## Volume of Phase Space

(ch_mce_ss_Phi_Omega)=
### Cumulative Phase-Space Volume and the Density of States

Let $\Phi(E, V, N)$ denote the **total phase-space volume enclosed by the energy surface $H = E$**:

$$\Phi(E, V, N) = \int\!\!\!\int_{H(\mathbf{q},\mathbf{p}) \leq E} d^{3N}q\, d^{3N}p$$

The volume of the energy shell is the difference:

$$\Gamma(E) = \Phi(E + \Delta) - \Phi(E) \approx \frac{\partial \Phi(E)}{\partial E}\,\Delta \equiv \Omega(E)\,\Delta$$

where we defined the **density of states**:

$$\boxed{\Omega(E, V, N) = \frac{\partial \Phi(E)}{\partial E} = \int d^{3N}q \int d^{3N}p\;\delta\bigl(E - H(\mathbf{q},\mathbf{p})\bigr)} \tag{2}$$

The density of states $\Omega(E)$ counts the number of microstates per unit energy at energy $E$. It is a fundamental quantity in statistical mechanics.

```{note}
**Shell vs. sphere:** In the thermodynamic limit $N \to \infty$, both $\ln \Gamma(E)$ and $\ln \Phi(E)$ give the same entropy, because the volume of a $3N$-dimensional sphere is dominated by its outermost shell. More precisely:

$$\ln \Phi(E) \sim 3N \ln R, \qquad \ln \Gamma(E) \sim (3N-1)\ln R + \ln \Delta \approx 3N \ln R$$

for a sphere of radius $R \propto \sqrt{E}$ in $3N$ dimensions. The correction is of order $\ln N$, which is negligible compared to the dominant $N$ scaling.
```

(ch_mce_ss_hypersphere)=
### Hyperspheres in High Dimensions

For later use, we compute the volume of an $n$-dimensional hypersphere of radius $R$:

$$V_n(R) = \int_{\sum_{i=1}^n x_i^2 \leq R^2} d^n x = R^n \, V_n(1)$$

To find $V_n(1)$, we use the Gaussian integral trick. Write:

$$\pi^{n/2} = \left(\int_{-\infty}^{+\infty} e^{-x^2} dx\right)^n = \int e^{-R^2} V_n(1)\, n R^{n-1}\, dR = \frac{n}{2} V_n(1)\,\Gamma\!\left(\frac{n}{2}\right)$$

where $\Gamma(z) = \int_0^\infty x^{z-1} e^{-x}\,dx$ is the Euler gamma function. This gives:

$$\boxed{V_n(1) = \frac{\pi^{n/2}}{\Gamma(n/2 + 1)}, \qquad V_n(R) = \frac{\pi^{n/2}}{\Gamma(n/2+1)} R^n} \tag{3}$$

For reference: $\Gamma(n) = (n-1)!$ for integer $n$, $\Gamma(1/2) = \sqrt{\pi}$, and $\Gamma(n+1) = n\,\Gamma(n)$.

```{example} Hypersphere volumes
- **2D** (disk): $V_2(R) = \pi R^2$, i.e. $\Gamma(2) = 1! = 1$, giving $V_2(1) = \pi$. ✓
- **3D** (ball): $V_3(R) = \frac{4}{3}\pi R^3$, i.e. $\Gamma(5/2) = \frac{3}{4}\sqrt{\pi}$, giving $V_3(1) = \frac{4\pi}{3}$. ✓
- **High $n$**: The volume is heavily concentrated in a thin shell near the surface — a striking geometric fact that has profound physical consequences.
```

(ch_mce_ss_shell_thickness)=
### How Thick Should the Energy Shell Be?

The entropy depends on $\ln \Gamma(E)$, which involves $\Delta$. Is the entropy then ambiguous?

$$\ln \Gamma = \ln\bigl(\Omega(E)\,\Delta\bigr) = \ln \Omega(E) + \ln \Delta$$

In the **thermodynamic limit** $N \to \infty$, we have $\ln \Omega(E) \propto N$ (as we shall see), while $\ln \Delta$ is at most $\mathcal{O}(\ln N)$ for any reasonable choice of $\Delta$ (e.g., $\Delta \sim k_B T$). Therefore:

$$\frac{\ln \Delta}{\ln \Omega(E)} \sim \frac{\ln N}{N} \to 0 \quad \text{as } N \to \infty$$

The entropy is **independent of $\Delta$** in the thermodynamic limit. This is a crucial consistency check of the whole framework.

(ch_mce_s_entropy)=
## Definition of Entropy

(ch_mce_ss_boltzmann_entropy)=
### The Boltzmann Postulate

Entropy cannot be obtained as an expectation value of a classical observable — it is instead a property of the **overall distribution** of microstates. Boltzmann's great insight (1872) was:

```{important}
**Boltzmann's Entropy Postulate:** The entropy is proportional to the logarithm of the number of accessible microstates:

$$\boxed{S(E, V, N) = k_B \ln \left[\frac{\Gamma(E, V, N)}{\Gamma_0(N)}\right]} \tag{4}$$

where $k_B = 1.381 \times 10^{-23}$ J/K is **Boltzmann's constant** and $\Gamma_0(N)$ is a normalisation factor that makes the argument dimensionless.
```

In the thermodynamic limit, this is equivalent to:

$$S = k_B \ln \Omega(E)\,\Delta \approx k_B \ln \Phi(E) \approx k_B \ln \Omega(E)$$

(all three expressions agree up to $\mathcal{O}(\ln N)$ corrections that are negligible for $N \sim 10^{23}$).

(ch_mce_ss_normalisation)=
### The Normalisation Factor $\Gamma_0(N)$

The normalisation constant $\Gamma_0(N)$ serves two purposes:

1. **Dimensional consistency:** $\Gamma$ has units of $(\text{position} \times \text{momentum})^{3N}$, so $\Gamma_0$ must cancel these units to make $\ln(\Gamma/\Gamma_0)$ dimensionless.

2. **Extensivity:** The $N$-dependence of $\Gamma_0$ is physically significant.

Classical mechanics cannot determine $\Gamma_0(N)$ uniquely. Quantum mechanics gives the answer:

$$\boxed{\Gamma_0(N) = h^{3N}\, N!} \tag{5}$$

where $h = 6.626 \times 10^{-34}$ m²·kg/s is Planck's constant.

- The factor $h^{3N}$ provides the natural unit of phase-space volume: each quantum state occupies a phase-space cell of volume $h^{3N}$ (from the Heisenberg uncertainty principle, $\Delta q \Delta p \sim h$).
- The factor $N!$ corrects for the **indistinguishability** of identical particles — exchanging two identical particles does not create a new microstate.

```{note}
Even in situations where classical mechanics gives an excellent description of particle dynamics (high temperature, low density), we must still account for the quantum-mechanical indistinguishability of particles. As Gibbs noted: even though one can treat the motion of molecules classically, one cannot label individual atomic particles as though they were macroscopic billiard balls. The factor $N!$ is a quantum-mechanical input into an otherwise classical theory.
```

(ch_mce_ss_entropy_expectation)=
### Entropy as an Expectation Value

The Boltzmann entropy can be rewritten as an expectation value. Since $\rho(\mathbf{q}, \mathbf{p}) = 1/\Gamma$ within the energy shell:

$$S = -k_B \ln\bigl[\Gamma_0(N)\, \rho(\mathbf{q}, \mathbf{p})\bigr] = -k_B \left\langle \ln\bigl[\Gamma_0(N)\, \rho(\mathbf{q}, \mathbf{p})\bigr] \right\rangle$$

This is precisely the **Gibbs–Shannon entropy**:

$$\boxed{S = -k_B \langle \ln[\Gamma_0 \rho] \rangle} \tag{6}$$

This form makes clear that entropy measures our **ignorance** about the microstate: the more microstates are accessible (large $\Gamma$), the less we know about which one the system is actually in, and the larger the entropy.

(ch_mce_s_temperature)=
## Definition of Temperature

(ch_mce_ss_thermal_contact)=
### Two Systems in Thermal Contact

Consider two systems, 1 and 2, characterised by $(E_1, V_1, N_1)$ and $(E_2, V_2, N_2)$ respectively, brought into **thermal contact**. They may exchange energy, but the total energy $E = E_1 + E_2$ is conserved. The volumes and particle numbers of each subsystem are kept fixed.

The combined phase-space volume factorises:

$$\frac{\Gamma(E, V, N)}{\Gamma_0(N)} = \sum_{E_1} \frac{\Gamma_1(E_1, V_1, N_1)}{\Gamma_0(N_1)} \cdot \frac{\Gamma_2(E - E_1, V_2, N_2)}{\Gamma_0(N_2)}$$

The **most probable** energy split is the one that maximises this product (equivalently, maximises the sum of logarithms, i.e., maximises the total entropy). Setting the derivative to zero:

$$\frac{d}{dE_1}\bigl[S_1(E_1) + S_2(E - E_1)\bigr] = 0 \implies \frac{\partial S_1}{\partial E_1}\bigg|_{V_1, N_1} = \frac{\partial S_2}{\partial E_2}\bigg|_{V_2, N_2}$$

(ch_mce_ss_temperature_def)=
### Definition of Temperature

The equilibrium condition above tells us that at equilibrium, a certain quantity is equalised between the two systems. We **define** the **thermodynamic temperature** as:

$$\boxed{\frac{1}{T} \equiv \left(\frac{\partial S}{\partial E}\right)_{V, N}} \tag{7}$$

so that the equilibrium condition becomes simply $T_1 = T_2$.

This definition has the correct properties:

- $T$ is intensive (it is equal for systems in thermal equilibrium, regardless of size).
- For most systems, $S$ increases with $E$, giving $T > 0$.
- It agrees with the empirical (ideal gas) definition of temperature.

(ch_mce_ss_law_of_large_numbers)=
### Sharpness of the Energy Distribution

The law of large numbers guarantees that the energy distribution $P(E_1) \propto e^{[S_1(E_1) + S_2(E-E_1)]/k_B}$ is **extremely sharply peaked** around its maximum value $E_1^*$ for large $N$. The width of the peak scales as $\sqrt{E_1^*} \propto \sqrt{N}$, while $E_1^* \propto N$. The relative width is:

$$\frac{\delta E_1}{E_1^*} \sim \frac{\sqrt{N}}{N} = \frac{1}{\sqrt{N}} \to 0 \quad \text{as } N \to \infty$$

For $N \sim 10^{23}$, relative fluctuations are of order $10^{-11\text{.}5}$ — completely unobservable. This is why thermodynamics is so precise.

(ch_mce_s_physical_temperature)=
## Physical Interpretation of Temperature

Temperature has a beautiful physical interpretation. From the definition $1/T = \partial S/\partial E|_{V,N}$:

$$T = \left(\frac{\partial E}{\partial S}\right)_{V,N}$$

This says that temperature measures **how much the energy of a system increases when entropy is added to it at fixed volume and particle number**.

More intuitively: entropy counts accessible microstates. A large $1/T$ means that adding a small amount of energy dramatically increases the number of accessible microstates — the system is "very responsive" to energy. A small $1/T$ (high temperature) means the system is already exploring so many microstates that adding more energy doesn't much increase the count.

**Direction of heat flow:** Suppose system 1 has $\partial S_1/\partial E_1 > \partial S_2/\partial E_2$, i.e., $T_1 < T_2$. Then transferring a small amount of energy $\delta E > 0$ from system 2 to system 1 increases the total entropy:

$$\delta S = \delta S_1 + \delta S_2 = \left(\frac{1}{T_1} - \frac{1}{T_2}\right)\delta E > 0$$

since $T_1 < T_2$. Heat flows spontaneously from the hotter system (higher $T$) to the cooler one (lower $T$) — exactly as we observe. This is the **second law of thermodynamics** emerging from the statistics.

**Kinetic interpretation:** For an ideal gas (derived in detail below), one finds $E = \frac{3}{2}Nk_BT$, so:

$$k_BT = \frac{2E}{3N} = \frac{2}{3} \times \text{average kinetic energy per particle}$$

Temperature is thus proportional to the **average kinetic energy per degree of freedom** — a result of the equipartition theorem. At room temperature ($T = 300$ K), $k_BT \approx 4 \times 10^{-21}$ J $\approx 0.025$ eV, setting the scale of thermal motion.

(ch_mce_s_validity)=
## Validity of the Statistical Description

The statistical mechanics framework rests on several assumptions. It is important to understand when these are valid.

**1. Large system (thermodynamic limit):** The statistical description becomes exact only in the limit $N \to \infty$, $V \to \infty$ with $N/V = $ const. For finite systems, fluctuations are suppressed only as $1/\sqrt{N}$, and corrections to thermodynamics are of order $(\ln N)/N$. In practice, for $N \gtrsim 10^{20}$, the description is excellent.

**2. Ergodicity:** The fundamental postulate (equal a priori probabilities) implicitly assumes that the system explores all accessible microstates equally over long times. This is the **ergodic hypothesis**: time averages equal ensemble averages. This is extremely difficult to prove rigorously for general Hamiltonians, but is believed to hold for most many-body systems and is supported by overwhelming experimental evidence.

**3. Classical vs. quantum treatment:** The classical phase-space description is valid when:

$$\lambda_{\text{dB}} = \frac{h}{\sqrt{2\pi m k_B T}} \ll n^{-1/3} = \left(\frac{V}{N}\right)^{1/3}$$

i.e., when the thermal de Broglie wavelength $\lambda_{\text{dB}}$ is much smaller than the mean inter-particle spacing. For an ideal gas at room temperature and atmospheric pressure, $\lambda_{\text{dB}} \sim 10^{-11}$ m while the mean spacing is $\sim 3 \times 10^{-9}$ m, so the classical treatment is justified.

**4. Equilibrium:** Statistical mechanics as developed here describes only **equilibrium** states. Non-equilibrium phenomena (transport, relaxation) require additional theoretical machinery.

**5. The $N!$ factor:** Even in the classical regime, the factor $N!$ from quantum indistinguishability must be included. Without it, the entropy is not extensive (the Gibbs paradox, discussed below).

+++ { "page-break": true }
+++

(ch_mce_s_pressure)=
## Definition of Pressure and the First Law

(ch_mce_ss_pressure_def)=
### Pressure from Entropy

By an argument entirely analogous to the derivation of temperature, we can define **pressure** through the entropy. Consider two systems separated by a movable, thermally insulating wall. Energy is fixed in each subsystem, but volumes $V_1$ and $V_2$ (with $V_1 + V_2 = V$ fixed) can adjust. Maximising the total entropy:

$$\frac{\partial S_1}{\partial V_1}\bigg|_{E_1, N_1} = \frac{\partial S_2}{\partial V_2}\bigg|_{E_2, N_2}$$

The common value of this derivative at equilibrium defines the pressure:

$$\boxed{\frac{P}{T} \equiv \left(\frac{\partial S}{\partial V}\right)_{E, N}} \tag{8}$$

Similarly, for a permeable wall (allowing particle exchange), the **chemical potential** is defined by:

$$\boxed{-\frac{\mu}{T} \equiv \left(\frac{\partial S}{\partial N}\right)_{E, V}} \tag{9}$$

The sign conventions in (8) and (9) are chosen to match the standard thermodynamic definitions.

(ch_mce_ss_first_law)=
### The First Law of Thermodynamics

From the definitions (7), (8), (9), we can write the total differential of the entropy $S(E, V, N)$:

$$dS = \left(\frac{\partial S}{\partial E}\right)_{V,N} dE + \left(\frac{\partial S}{\partial V}\right)_{E,N} dV + \left(\frac{\partial S}{\partial N}\right)_{E,V} dN$$

$$dS = \frac{1}{T}\,dE + \frac{P}{T}\,dV - \frac{\mu}{T}\,dN$$

Rearranging (identifying the internal energy $U = E$ for an isolated system):

$$\boxed{dU = T\,dS - P\,dV + \mu\,dN} \tag{10}$$

This is the **fundamental thermodynamic relation** — the combined first and second laws. It encodes:

- $T\,dS$: energy added as **heat** (reversibly)
- $-P\,dV$: work done **by** the system
- $\mu\,dN$: energy change due to **particle exchange**

The entire framework of thermodynamics — all thermodynamic potentials, Maxwell relations, equations of state — follows from equation (10) once $S(E, V, N)$ is known.

(ch_mce_ss_calculation_recipe)=
### Recipe for Microcanonical Calculations

Given a Hamiltonian $H(\mathbf{q}, \mathbf{p})$, the microcanonical calculation proceeds as:

1. **Write the Hamiltonian** $H(\mathbf{q}, \mathbf{p})$.
2. **Compute the phase-space volume** $\Phi(E)$ enclosed by the energy surface $H = E$.
3. **Obtain the entropy:** $S(E, V, N) = k_B \ln[\Phi(E) / \Gamma_0(N)]$.
4. **Differentiate** to get temperature, pressure, and chemical potential via (7), (8), (9).
5. **Invert** $1/T = \partial S/\partial E$ to get energy as a function of temperature: $U = U(T, V, N)$ (caloric equation of state).
6. **Use** $P/T = \partial S/\partial V$ to get the equation of state $P = P(T, V, N)$.
7. **Legendre transform** to obtain other thermodynamic potentials: $F = U - TS$, $H = U + PV$, $G = U + PV - TS$.

(ch_mce_s_statistical_entropy)=
## Statistical Interpretation of Entropy

The statistical definition $S = k_B \ln \Gamma$ gives entropy a profound physical meaning. Let us explore several aspects of this.

(ch_mce_ss_extensivity)=
### Extensivity of Entropy

For a system of $N$ weakly interacting particles in volume $V$ with energy $E$, the phase-space volume scales as:

$$\Gamma(E, V, N) \sim V^N \gamma^N(E/N, V/N)$$

The $V^N$ factor comes from the $N$ independent position integrals. With $\Gamma_0(N) = h^{3N}N!$ and Stirling's formula $\ln N! \approx N(\ln N - 1)$:

$$S = k_B \ln\frac{\Gamma}{\Gamma_0} = k_B \ln\frac{V^N \gamma^N}{h^{3N}N!} = k_B N \left[\ln\!\left(\frac{V}{N}\frac{\gamma}{h^3}\right) + 1\right] \equiv N\, s(E/N, V/N)$$

The entropy is **extensive**: $S$ is proportional to $N$ for fixed intensive quantities $E/N$ and $V/N$. This is the correct thermodynamic behaviour.

(ch_mce_ss_second_law)=
### Second Law of Thermodynamics

The statistical entropy satisfies the second law: **the entropy of an isolated system can only increase (or remain constant) for any spontaneous process**.

Consider **free expansion**: a gas initially confined to volume $V_1$ is allowed to expand into volume $V_2 > V_1$ (into vacuum, no heat or work exchange). The total energy $E$ and $N$ are unchanged, but the accessible phase-space volume increases because more positions are now available for each particle:

$$\Gamma(E, V_2, N) > \Gamma(E, V_1, N) \implies S(E, V_2, N) > S(E, V_1, N)$$

More generally:

```{important}
**Equivalence:** The second law of thermodynamics is equivalent to the statement that removing dynamical constraints on a system (opening a partition, allowing a wall to move) increases the number of accessible microstates, and hence increases the entropy.
```

The entropy is maximised at equilibrium: $S(x)$ is maximal with respect to any internal constraint $x$ that is free to adjust.

(ch_mce_ss_additivity)=
### Additivity and Thermal Equilibrium

When two systems with the same intensive properties (same $T$, $P$, $\mu$) are combined, entropy is additive:

$$S(E, V, N) = S(E_1, V_1, N_1) + S(E_2, V_2, N_2)$$

For systems in thermal contact (different initial temperatures), the combined system evolves to the state of maximum entropy, which corresponds to equal temperatures. The correction term $\sim k_B \ln(\sqrt{E_{\text{max}}})$ from the width of the energy distribution is negligible compared to the extensive entropy $\sim Nk_B$, confirming additivity in the thermodynamic limit.

+++ { "page-break": true }
+++

(ch_mce_s_ideal_gas)=
## Classical Ideal Gas in the Microcanonical Ensemble

The ideal gas is the paradigmatic application of the microcanonical ensemble. It consists of $N$ non-interacting, structureless classical particles in a box of volume $V$, described by the Hamiltonian:

$$H(\mathbf{q}, \mathbf{p}) = \sum_{i=1}^{N} \frac{\mathbf{p}_i^2}{2m} \tag{11}$$

Note that $H$ depends only on momenta — the position degrees of freedom are free (confined only to the box volume $V$).

(ch_mce_ss_phase_space_volume)=
### Phase-Space Volume

The total phase-space volume enclosed by the energy surface $H = E$ is:

$$\Phi(E) = \int_{H(\mathbf{q},\mathbf{p}) \leq E} d^{3N}q\, d^{3N}p$$

The position integrals are free (each particle ranges over the volume $V$), giving $V^N$:

$$\Phi(E) = V^N \int_{\sum_i \mathbf{p}_i^2/(2m) \leq E} d^{3N}p = V^N \int_{\sum_i \mathbf{p}_i^2 \leq 2mE} d^{3N}p$$

The remaining integral is the volume of a $3N$-dimensional hypersphere of radius $R = \sqrt{2mE}$:

$$\int_{\sum_i \mathbf{p}_i^2 \leq 2mE} d^{3N}p = V_{3N}(\sqrt{2mE}) = \frac{\pi^{3N/2}}{\Gamma(3N/2 + 1)}\,(2mE)^{3N/2}$$

Therefore:

$$\boxed{\Phi(E, V, N) = V^N \cdot \frac{\pi^{3N/2}}{\Gamma(3N/2 + 1)}\,(2mE)^{3N/2}} \tag{12}$$

(ch_mce_ss_sackur_tetrode)=
### Entropy: The Sackur-Tetrode Equation

Using $S = k_B \ln[\Phi(E)/\Gamma_0(N)]$ with $\Gamma_0(N) = h^{3N}N!$:

$$S = k_B \ln\left[\frac{V^N (2mE)^{3N/2} \pi^{3N/2}}{h^{3N} N! \,\Gamma(3N/2+1)}\right]$$

For large $N$, we use Stirling's formula in two places:

- $\ln(N!) \approx N\ln N - N$
- $\ln\Gamma(3N/2+1) = \ln(3N/2)! \approx \frac{3N}{2}\ln\frac{3N}{2} - \frac{3N}{2}$

A careful expansion gives the **Sackur-Tetrode equation** (1912):

$$\boxed{S(E, V, N) = k_B N \left[\ln\left(\frac{V}{N}\left(\frac{4\pi m E}{3Nh^2}\right)^{3/2}\right) + \frac{5}{2}\right]} \tag{13}$$

This is exact in the thermodynamic limit. It is a remarkable result: from purely mechanical counting, we have derived the full entropy of an ideal gas.

```{note}
The argument of the logarithm in the Sackur-Tetrode equation can also be written using the thermal de Broglie wavelength $\lambda = h/\sqrt{2\pi m k_B T}$:

$$S = Nk_B\left[\ln\frac{1}{\rho\lambda^3} + \frac{5}{2}\right]$$

where $\rho = N/V$ is the number density. The condition for the classical description to be valid is $\rho\lambda^3 \ll 1$, i.e., the de Broglie wavelength is much smaller than the mean inter-particle spacing $\rho^{-1/3}$.
```

(ch_mce_ss_equations_of_state)=
### Equations of State

From the Sackur-Tetrode equation, all thermodynamic quantities follow by differentiation.

**Caloric equation of state** (energy–temperature relation):

$$\frac{1}{T} = \left(\frac{\partial S}{\partial E}\right)_{V,N} = k_B N \cdot \frac{3}{2} \cdot \frac{1}{E} \implies \boxed{U = E = \frac{3}{2}Nk_BT} \tag{14}$$

Each particle contributes $\frac{3}{2}k_BT$ to the energy — $\frac{1}{2}k_BT$ per translational degree of freedom. This is the **equipartition theorem** for the ideal gas.

**Thermal equation of state** (pressure–volume–temperature relation):

$$\frac{P}{T} = \left(\frac{\partial S}{\partial V}\right)_{E,N} = \frac{k_B N}{V} \implies \boxed{PV = Nk_BT} \tag{15}$$

This is the **ideal gas law**, here derived from first principles in statistical mechanics.

**Heat capacity at constant volume:**

$$C_V = \left(\frac{\partial U}{\partial T}\right)_{V,N} = \frac{3}{2}Nk_B \tag{16}$$

**Chemical potential:**

$$\mu = -T\left(\frac{\partial S}{\partial N}\right)_{E,V} = k_BT\ln(\rho\lambda^3) \tag{17}$$

where $\lambda = h/\sqrt{2\pi m k_B T}$ is the thermal de Broglie wavelength. Note that $\mu < 0$ in the classical regime (where $\rho\lambda^3 \ll 1$).

(ch_mce_ss_python_ideal_gas)=
### Numerical Illustration: The Ideal Gas

Let us verify the Sackur-Tetrode equation numerically and explore how entropy, temperature, and pressure depend on the thermodynamic variables.

````{example} Plotting entropy and equations of state

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln

# Constants
kB = 1.381e-23   # J/K  Boltzmann constant
h  = 6.626e-34   # J·s  Planck constant
m  = 4.65e-26    # kg   mass of N2 molecule (~28 amu)
N  = 1e20        # number of particles (macroscopic but tractable)

def sackur_tetrode(E, V, N):
    """Return S/N k_B (dimensionless entropy per particle)."""
    prefactor = (V/N) * (4*np.pi*m*E / (3*N*h**2))**(3/2)
    return np.log(prefactor) + 5/2

# ── Panel 1: S vs E at fixed V, N ──────────────────────────────────────────
E_vals = np.linspace(0.01, 5, 300) * N * kB * 300   # energies around 300K scale

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

V0 = 1e-3  # 1 litre
S_vals = sackur_tetrode(E_vals, V0, N)

axes[0].plot(E_vals / (N * kB), S_vals, 'royalblue', lw=2)
axes[0].set_xlabel(r'$E\;/\;Nk_B$ (K)', fontsize=12)
axes[0].set_ylabel(r'$S\;/\;Nk_B$', fontsize=12)
axes[0].set_title('Entropy vs Energy\n(ideal gas, fixed $V$, $N$)', fontsize=11)
axes[0].grid(True, alpha=0.3)

# ── Panel 2: T vs E (caloric equation of state) ────────────────────────────
T_vals = (2/3) * E_vals / (N * kB)   # E = (3/2) N kB T  →  T = 2E/(3NkB)

axes[1].plot(E_vals / (N * kB), T_vals, 'tomato', lw=2)
axes[1].set_xlabel(r'$E\;/\;Nk_B$ (K)', fontsize=12)
axes[1].set_ylabel(r'Temperature $T$ (K)', fontsize=12)
axes[1].set_title('Caloric Equation of State\n$E = \\frac{3}{2}Nk_BT$', fontsize=11)
axes[1].grid(True, alpha=0.3)

# ── Panel 3: PV/NkBT vs T (should = 1 for ideal gas) ─────────────────────
T_range = np.linspace(100, 1000, 200)
# P/T = kB N / V  →  PV/(NkBT) = 1  (ideal gas law check)
PV_over_NkT = np.ones_like(T_range)  # identically 1 by construction

axes[2].plot(T_range, PV_over_NkT, 'seagreen', lw=2)
axes[2].axhline(1, color='k', lw=0.8, linestyle='--', alpha=0.5)
axes[2].set_ylim(0.8, 1.2)
axes[2].set_xlabel(r'Temperature $T$ (K)', fontsize=12)
axes[2].set_ylabel(r'$PV\;/\;Nk_BT$', fontsize=12)
axes[2].set_title('Ideal Gas Law\n$PV = Nk_BT$', fontsize=11)
axes[2].grid(True, alpha=0.3)

plt.suptitle('Classical Ideal Gas — Microcanonical Ensemble', fontsize=13, y=1.02)
plt.tight_layout()
plt.show()
```
````

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt

# ─── Entropy per particle as a function of V/N (at fixed T) ───────────────
kB = 1.381e-23
h  = 6.626e-34
m  = 4.65e-26   # N2

T = 300   # K
# From E = (3/2)N kB T:  E/N = (3/2) kB T
# Sackur-Tetrode: S/NkB = ln[(V/N)(4πm(3/2 NkBT)/(3Nh²))^{3/2}] + 5/2
#                       = ln[(V/N) / λ³] + 5/2
# thermal de Broglie wavelength
lam = h / np.sqrt(2 * np.pi * m * kB * T)
print(f"Thermal de Broglie wavelength at {T} K: λ = {lam:.3e} m")

rho_range = np.logspace(20, 27, 300)   # particles per m³
V_over_N  = 1.0 / rho_range

S_over_NkB = np.log(V_over_N / lam**3) + 5/2

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.semilogx(rho_range, S_over_NkB, 'darkorchid', lw=2)
ax.axvline(1/lam**3, color='red', linestyle='--', lw=1.2,
           label=r'$\rho\lambda^3 = 1$ (classical breakdown)')
ax.axhline(0, color='gray', linestyle=':', lw=0.8, label='$S = 0$')
ax.fill_betweenx([-5, S_over_NkB.min()], 1/lam**3, rho_range.max(),
                 alpha=0.15, color='red', label='Quantum regime')
ax.set_xlabel(r'Number density $\rho = N/V$ (m$^{-3}$)', fontsize=12)
ax.set_ylabel(r'$S\;/\;Nk_B$', fontsize=12)
ax.set_title(f'Sackur-Tetrode Entropy vs Density at $T = {T}$ K', fontsize=12)
ax.legend(fontsize=10)
ax.set_ylim(-5, 20)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Atmospheric density
rho_atm = 2.7e25   # m^{-3}
S_atm   = np.log(1/rho_atm / lam**3) + 5/2
print(f"At atmospheric density: S/(NkB) = {S_atm:.3f}")
print(f"ρλ³ = {rho_atm * lam**3:.3e}  (≪ 1: classical regime valid)")
```

````

The red dashed line marks $\rho\lambda^3 = 1$, where the classical ideal gas description breaks down and quantum statistics must be used. At atmospheric conditions for N₂, $\rho\lambda^3 \approx 10^{-6} \ll 1$, confirming that classical mechanics is an excellent approximation.

(ch_mce_s_gibbs_paradox)=
## The Gibbs Paradox

(ch_mce_ss_problem)=
### The Problem: Non-Extensive Entropy

If we naively ignore the factor of $N!$ (i.e., treat the particles as distinguishable), the entropy of the ideal gas becomes:

$$S_{\text{classical}} = k_BN\left[\ln\!\left(\frac{(2\pi m E)^{3/2} V}{h^{3N}}\cdot\frac{1}{\Gamma(3N/2+1)}\right)\right] \approx k_BN\left[\ln\!\left(V\left(\frac{4\pi mE}{3Nh^2}\right)^{3/2}\right) + \frac{3}{2}\right]$$

Note that the volume $V$ appears instead of $V/N$. This means $S_{\text{classical}}$ is **not extensive**: if we double both $E$, $V$, and $N$, we get:

$$S_{\text{classical}}(2E, 2V, 2N) = 2Nk_B\left[\ln\!\left(2V\left(\frac{4\pi m \cdot 2E}{3 \cdot 2N h^2}\right)^{3/2}\right) + \frac{3}{2}\right] \neq 2\,S_{\text{classical}}(E, V, N)$$

The entropy contains a spurious additive term $2Nk_B\ln 2$.

(ch_mce_ss_mixing_paradox)=
### Mixing Entropy and the Paradox

The paradox is sharpest when considering the **mixing** of two identical ideal gases. Imagine a container divided by a partition, with the same gas at the same temperature, pressure, and density on both sides. When the partition is removed, intuitively **nothing changes** — the gas is the same on both sides. The entropy should not change.

However, with distinguishable particles:

$$\Delta S_{\text{mix}} = Nk_B\ln 2 + Nk_B\ln 2 = 2Nk_B\ln 2 > 0$$

This is the **Gibbs paradox** (1875): the entropy increases upon mixing two identical gases, which is physically wrong. If the gases were different, a positive $\Delta S_{\text{mix}}$ would be correct (and is the entropy of mixing familiar from chemistry). But for identical gases, there should be no entropy change.

```{note}
The entropy of mixing two *different* ideal gases (Dalton's law situation) is genuinely positive:
$$\Delta S_{\text{mix}} = -Nk_B(x_1 \ln x_1 + x_2 \ln x_2) > 0$$
where $x_i = N_i/N$ are mole fractions. This is physically correct and measurable. The Gibbs paradox arises only when the two gases are *identical*.
```

(ch_mce_ss_resolution)=
### Resolution: Indistinguishability and the $N!$ Factor

The resolution of the Gibbs paradox is the **indistinguishability** of identical particles. If particles are quantum mechanically identical, then permuting two particles does not create a new microstate — it produces the same physical state. The number of **physically distinct** microstates is therefore $\Gamma_{\text{physical}} = \Gamma_{\text{naive}}/N!$.

With the correct $\Gamma_0(N) = h^{3N}N!$, the Sackur-Tetrode entropy (13) is explicitly extensive (it depends on $V/N$, not $V$), and:

$$\Delta S_{\text{mix}} = 0$$

for mixing of identical gases at equal $T$, $P$ — as it should be.

```{important}
The Gibbs paradox reveals that even a classical system cannot be treated without incorporating a quantum mechanical input: the indistinguishability of identical particles. The factor $N!$ is not an arbitrary "fudge" — it is dictated by the deep structure of quantum mechanics and is absolutely necessary for a self-consistent thermodynamics.
```

(ch_mce_ss_python_gibbs)=
### Numerical Illustration: Entropy of Mixing

````{example} Gibbs paradox visualised

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt

def S_sackur(N, V, T, m, kB, h):
    """Sackur-Tetrode entropy (correct, with N!)."""
    lam = h / np.sqrt(2 * np.pi * m * kB * T)
    return N * kB * (np.log((V/N) / lam**3) + 5/2)

def S_classical(N, V, T, m, kB, h):
    """Entropy without N! (wrong for identical particles)."""
    lam = h / np.sqrt(2 * np.pi * m * kB * T)
    return N * kB * (np.log(V / lam**3) + 5/2 - np.log(N))
    # Note: the -ln(N) above comes from Stirling but not the N! denominator;
    # without N! at all: S_wrong = N kB [ln(V/λ³) + 3/2]

kB = 1.381e-23
h  = 6.626e-34
m  = 4.65e-26  # N2
T  = 300.0     # K

N_vals = np.logspace(10, 23, 200)
V      = 1e-3   # 1 litre total (half on each side before mixing)

# Before mixing: two halves, each with N/2 particles in V/2
# After mixing: N particles in V  (same gas, same T and rho)
dS_correct  = np.array([S_sackur(N, V, T, m, kB, h)
                         - 2*S_sackur(N/2, V/2, T, m, kB, h)
                         for N in N_vals])

dS_wrong    = np.array([
    (N * kB * np.log(V) - 2*(N/2)*kB*np.log(V/2))   # simplified S_wrong ~ NkB ln V
    for N in N_vals])

fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogx(N_vals, dS_correct / kB, 'royalblue', lw=2,
            label=r'Correct ($N!$ included): $\Delta S = 0$')
ax.semilogx(N_vals, dS_wrong / kB, 'tomato', lw=2, linestyle='--',
            label=r'Wrong (no $N!$): $\Delta S = N\ln 2$')
ax.axhline(0, color='k', lw=0.8, linestyle=':')
ax.set_xlabel('Number of particles $N$', fontsize=12)
ax.set_ylabel(r'$\Delta S_{\rm mix}\;/\;k_B$', fontsize=12)
ax.set_title('Gibbs Paradox: Entropy of Mixing Identical Ideal Gases', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```
````

The blue curve (correct result with $N!$) gives exactly $\Delta S = 0$ for all $N$, while the red curve (without $N!$) gives a spurious entropy increase growing with system size — the Gibbs paradox in full display.

+++ { "page-break": true }
+++

(ch_mce_s_summary)=
## Summary

In this chapter we constructed the microcanonical ensemble from first principles and derived the fundamental results of equilibrium statistical mechanics. The key results are:

**The microcanonical ensemble** describes an isolated system $(E, V, N)$ fixed. All accessible microstates are equally probable:

$$\rho(\mathbf{q}, \mathbf{p}) = \frac{1}{\Gamma(E,V,N)}\;\text{ for }\; E < H < E+\Delta$$

**Boltzmann entropy:**

$$S(E, V, N) = k_B \ln\!\left[\frac{\Gamma(E,V,N)}{h^{3N}N!}\right]$$

**Thermodynamic quantities from entropy** (with $U \equiv E$):

$$\frac{1}{T} = \left(\frac{\partial S}{\partial U}\right)_{V,N}, \quad \frac{P}{T} = \left(\frac{\partial S}{\partial V}\right)_{U,N}, \quad -\frac{\mu}{T} = \left(\frac{\partial S}{\partial N}\right)_{U,V}$$

**First law (fundamental relation):**

$$dU = T\,dS - P\,dV + \mu\,dN$$

**Classical ideal gas — Sackur-Tetrode equation:**

$$S(E,V,N) = k_BN\left[\ln\!\left(\frac{V}{N}\!\left(\frac{4\pi mE}{3Nh^2}\right)^{3/2}\right) + \frac{5}{2}\right]$$

**Equations of state for the ideal gas:**

$$U = \frac{3}{2}Nk_BT, \qquad PV = Nk_BT, \qquad C_V = \frac{3}{2}Nk_B$$

**Gibbs paradox:** Treating identical particles as distinguishable leads to a non-extensive entropy and a spurious entropy of mixing. The resolution requires the factor $N!$ from quantum-mechanical indistinguishability, even within an otherwise classical treatment.

```{note} What comes next
The microcanonical ensemble is conceptually fundamental but often computationally cumbersome — fixing $E$ exactly makes calculations difficult for realistic interacting systems. In the next chapters, we introduce:

- The **canonical ensemble** (fixed $T$, $V$, $N$): the system is in contact with a heat bath at temperature $T$. The constraint on energy is relaxed in favour of fixing $T$. This is usually much easier to work with and gives equivalent results in the thermodynamic limit.
- The **grand-canonical ensemble** (fixed $T$, $V$, $\mu$): both energy and particle number can fluctuate. This is essential for quantum gases (Bose-Einstein condensation, Fermi gases).

All three ensembles are equivalent in the thermodynamic limit $N \to \infty$, but differ in computational convenience and in how they handle fluctuations.
```

```{important} Suggested Books

- R. K. Pathria and P. D. Beale, *Statistical Mechanics*, 4th Ed. (Elsevier)
- L. D. Landau and E. M. Lifshitz, *Statistical Physics*, Part 1 (Butterworth-Heinemann)
- M. Kardar, *Statistical Physics of Particles* (Cambridge University Press)
- K. Huang, *Statistical Mechanics*, 2nd Ed. (Wiley)
- H. B. Callen, *Thermodynamics and an Introduction to Thermostatistics*, 2nd Ed. (Wiley)
- F. Reif, *Fundamentals of Statistical and Thermal Physics* (Waveland Press)
- C. Kittel and H. Kroemer, *Thermal Physics*, 2nd Ed. (W. H. Freeman)
- L. Bocquet, *Lecture Notes on Statistical Physics* (ENS-PSL)
```

(ch_mce_s_exercises)=
## Exercises

```{exercise} Energy shell and density of states
:label: ex_mce_dos

For a single classical particle of mass $m$ in a 3D box of volume $V$:

1. Write down the Hamiltonian $H(\mathbf{q}, \mathbf{p})$.
2. Compute the cumulative phase-space volume $\Phi(E)$ enclosed by the energy surface.
3. From $\Phi(E)$, compute the density of states $\Omega(E) = \partial \Phi/\partial E$.
4. Sketch $\Omega(E)$ vs $E$. What power law does it follow?
```

```{exercise} Entropy and temperature of a classical ideal gas
:label: ex_mce_entropy

Starting from the Sackur-Tetrode equation (13):

1. Verify that $1/T = (\partial S/\partial E)_{V,N}$ gives $E = \frac{3}{2}Nk_BT$.
2. Verify that $P/T = (\partial S/\partial V)_{E,N}$ gives $PV = Nk_BT$.
3. Compute the chemical potential $\mu = -T(\partial S/\partial N)_{E,V}$ and express it in terms of the de Broglie wavelength $\lambda$.
4. Show that $S$ as given by (13) is explicitly extensive: $S(\lambda E, \lambda V, \lambda N) = \lambda S(E, V, N)$ for any $\lambda > 0$.
```

```{exercise} Thickness of the energy shell
:label: ex_mce_shell

Show explicitly that the entropy $S = k_B \ln \Gamma$ does not depend on the shell width $\Delta$ in the thermodynamic limit:

1. Write $\Gamma(E) = \Omega(E)\,\Delta$ where $\Omega(E) \propto E^{3N/2 - 1}$.
2. Compute $S_\Delta = k_B \ln(\Omega(E)\Delta)$ and $S_\Phi = k_B \ln \Phi(E)$ where $\Phi(E) \propto E^{3N/2}$.
3. Show that $(S_\Delta - S_\Phi)/S_\Phi \to 0$ as $N \to \infty$ for any finite $\Delta$.
```

```{exercise} Physical temperature interpretation
:label: ex_mce_temp_interp

An ideal gas is in equilibrium at temperature $T$.

1. Show that the average kinetic energy per particle is $\langle \epsilon \rangle = \frac{3}{2}k_BT$.
2. If the gas is heated from $T_1 = 300$ K to $T_2 = 600$ K at constant volume, by what factor does the number of accessible microstates change? (*Hint*: use $\Gamma \propto E^{3N/2}$ and $E = \frac{3}{2}Nk_BT$.)
3. Explain qualitatively why heat flows from hot to cold in terms of the total entropy.
```

```{exercise} Gibbs paradox
:label: ex_mce_gibbs

Consider a container of volume $V$ divided into two equal halves, each containing $N/2$ identical molecules at temperature $T$ and pressure $P$.

1. Compute the entropy change $\Delta S$ when the partition is removed, using (a) the Sackur-Tetrode formula and (b) the formula without the $N!$ factor.
2. Show that with the correct formula, $\Delta S = 0$.
3. Now repeat the calculation for the case where the two halves contain *different* gases (A and B). Show that $\Delta S = -Nk_B(x_A \ln x_A + x_B \ln x_B)$ where $x_A = x_B = 1/2$.
4. Explain physically why the two cases give different results.
```

```{exercise} Hypersphere volumes
:label: ex_mce_hypersphere

Using the gamma-function result $V_n(R) = \pi^{n/2}R^n/\Gamma(n/2+1)$:

1. Verify the formula for $n = 1, 2, 3$.
2. Show that for large $n$, the surface area $A_n(R) = dV_n/dR$ satisfies $A_n(R)/V_n(R) = n/R$.
3. Show that for large $n$, most of the volume of a hypersphere lies within a thin shell of relative thickness $\sim 1/n$ near the surface. (*Hint*: compute the fraction of volume in the shell $[R-\delta, R]$ for small $\delta$.)
4. Explain the physical significance of this result for statistical mechanics.
```

```{exercise} Heat capacities and adiabatic processes
:label: ex_mce_adiabatic

For an ideal gas with $N$ particles in $d$ dimensions, the energy is $E = \frac{d}{2}Nk_BT$.

1. Derive $C_V = \frac{d}{2}Nk_B$ and $C_P = \frac{d+2}{2}Nk_B$.
2. Define $\gamma = C_P/C_V$ and show that for an adiabatic process ($dS = 0$), $PV^\gamma = \text{const}$.
3. For $d = 3$ (monatomic ideal gas), $\gamma = 5/3$. Evaluate the ratio for a diatomic gas treated classically (5 quadratic degrees of freedom per molecule).
```
```

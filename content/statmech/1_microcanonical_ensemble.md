---
title: Introduction to Statistical Physics

numbering:
  title:
    enabled: false
    offset: 0

kernelspec:
  name: python3
  display_name: 'Python 3'
---

```{card}
In this part we introduce the foundations of **Statistical Physics**. We begin with a simple but deep question: how can the behaviour of a gas, a magnet, or a liquid — each made of an astronomically large number of atoms — be understood from the microscopic rules governing each individual particle?

The central idea is that we do not need to track every atom. Instead, we use *probability* and *statistics* to extract the macroscopic behaviour that emerges when enormously many particles interact together. The key concepts developed here are the *microstate*, the *macrostate*, the *fundamental postulate of equal probabilities*, and the celebrated *Boltzmann entropy* $S = k_B \ln \Omega$.

From these foundations, the familiar quantities of thermodynamics — *temperature*, *pressure*, *energy*, and *entropy* — are given deep microscopic meaning for the first time. Special attention is given to the *classical ideal gas* as the paradigmatic application, which lets us derive the ideal gas law and the equipartition theorem from scratch.

This part aims to build physical intuition first, introducing mathematical tools only when they genuinely illuminate the physics, and preparing the reader for the study of the canonical and grand-canonical ensembles.
```

(ch-statphys)=
# Introduction to Statistical Physics

Imagine you have a bottle of air sitting on your desk. Inside that bottle are roughly $10^{22}$ molecules — nitrogen, oxygen, a little argon — all flying around, bouncing off each other and off the walls, changing direction billions of times per second. If you tried to write down the position and velocity of every single molecule and track how each one evolves according to Newton's laws, you would need more computer memory than exists on Earth, and more time than the age of the universe.

And yet, with just three numbers — **temperature**, **pressure**, and **volume** — you can describe that bottle of air with extraordinary precision. How is that possible?

The answer is **statistical physics**: the art of extracting exact macroscopic predictions from microscopic chaos, not *despite* the enormous number of particles, but *because* of it. When $10^{22}$ particles are involved, statistical fluctuations become negligibly small, and average behaviour becomes essentially certain. This is the deep reason thermodynamics works so well.

This chapter introduces the ideas, vocabulary, and first results of statistical physics. We will:

1. Understand the difference between **microstates** and **macrostates**.
2. State and explore the **fundamental postulate** of statistical mechanics.
3. Introduce **phase space** as the arena where statistical mechanics lives.
4. Define the **Boltzmann entropy** and understand what it really means.
5. Derive **temperature**, **pressure**, and **the first law** from statistics alone.
6. Apply everything to the **classical ideal gas**.
7. Resolve the **Gibbs paradox** and understand why particles must be treated as indistinguishable.

(ch_sp_s_micro_macro)=
## Microstates and Macrostates

(ch_sp_ss_the_gap)=
### The Gap Between Micro and Macro

The most fundamental distinction in statistical physics is the difference between a **microstate** and a **macrostate**.

A **microstate** is the complete, precise description of a system at the atomic level — the exact position and velocity of every single particle, at every instant. For a gas of $N$ particles in three dimensions, this is a list of $6N$ numbers: three position coordinates and three momentum components per particle. The microstate changes billions of times per second as particles collide and move.

A **macrostate** is what we actually measure in the lab — quantities like temperature $T$, pressure $P$, volume $V$, and total energy $E$. A macrostate is specified by just a handful of numbers, regardless of whether the system contains $10$ or $10^{23}$ particles.

```{important}
The central observation of statistical mechanics is that a single macrostate corresponds to an **astronomically large number** of microstates. All of those microstates look identical to any macroscopic measurement, yet they are genuinely different at the microscopic level.

The job of statistical mechanics is to figure out which macrostate a system will adopt, given what we know about the microstates.
```

(ch_sp_ss_coins)=
### A Warm-Up: Coins and Counting

Before diving into atoms and gases, let us build intuition with the simplest possible example: flipping coins.

Imagine you flip $N = 4$ coins. Each coin can be Heads (H) or Tails (T). The **microstate** is the full outcome of every coin — something like HHTH. There are $2^4 = 16$ equally probable microstates in total.

Now suppose you only care about the **total number of heads** — this is the macrostate. There are only 5 possible macrostates: 0, 1, 2, 3, or 4 heads.

How many microstates correspond to each macrostate?

| Macrostate (# heads) | Number of microstates | Example microstates |
|:---:|:---:|:---|
| 0 | 1 | TTTT |
| 1 | 4 | HTTT, THTT, TTHT, TTTH |
| 2 | 6 | HHTT, HTHT, HTTH, THHT, THTH, TTHH |
| 3 | 4 | HHHT, HHTH, HTHH, THHH |
| 4 | 1 | HHHH |

The macrostate with 2 heads is **6 times more likely** than the macrostate with 0 heads, simply because it has 6 times as many corresponding microstates. If we observed the system at a random moment, we would almost certainly find it in the most probable macrostate.

Now scale this up: for $N = 100$ coins, the number of microstates for 50 heads is roughly $10^{29}$, while 100 heads has just 1. For $N = 10^{23}$ coins (the scale of real physics), the most probable macrostate is so overwhelmingly more likely than any other that it is, for all practical purposes, **the only macrostate we ever observe**. This is why thermodynamics is so precise.

````{example} Coin flipping and the emergence of a macrostate

The code below simulates flipping $N$ coins many times and shows how the distribution of outcomes sharpens dramatically as $N$ grows.

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

for ax, N in zip(axes, [10, 100, 1000]):
    k = np.arange(0, N + 1)
    # Exact binomial probabilities  p = 1/2
    prob = comb(N, k, exact=False) / 2**N

    ax.bar(k / N, prob * N, width=1/N, color='steelblue',
           alpha=0.75, align='center')
    ax.set_xlabel('Fraction of Heads', fontsize=12)
    ax.set_ylabel('Probability', fontsize=11)
    ax.set_title(f'$N = {N}$ coins', fontsize=13)
    ax.set_xlim(0, 1)
    ax.axvline(0.5, color='crimson', linestyle='--', lw=1.5,
               label='Most probable')
    ax.legend(fontsize=9)

plt.suptitle('Distribution of Heads — the Peak Sharpens with N',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.show()
```

Notice how the distribution becomes sharper and sharper as $N$ increases. For large $N$, the fraction of heads is essentially always exactly $1/2$ — the fluctuations become vanishingly small. This is the statistical origin of the certainty we find in thermodynamics.

````

+++ { "page-break": true }
+++

(ch_sp_s_fundamental_postulate)=
## The Fundamental Postulate

(ch_sp_ss_equal_probability)=
### All Microstates Are Equally Likely

The entire edifice of statistical mechanics rests on one foundational assumption:

```{important}
**The Fundamental Postulate of Statistical Mechanics:**

An isolated system in equilibrium is equally likely to be found in any one of its accessible microstates.

No accessible microstate is favoured over any other.
```

This might seem like a very bold claim. Why should we believe it? There are several reasons:

**It is the most unbiased assumption.** If we have no information that distinguishes one microstate from another, it is logically most honest to assign them equal probability. Any other choice would require us to explain why some microstates are preferred.

**It is consistent with the equations of motion.** Liouville's theorem in classical mechanics — and its quantum analogue — shows that the density of phase-space trajectories is conserved as the system evolves. A uniform distribution over accessible microstates is therefore a stationary (time-independent) distribution, consistent with equilibrium.

**It works.** Every prediction derived from this postulate agrees with experiment to extraordinary precision. That is ultimately its strongest justification.

```{note}
The word "accessible" is important. A microstate is accessible if it is consistent with the known macroscopic constraints — in particular, if it has the right total energy $E$, is confined to the right volume $V$, and contains the right number of particles $N$. Microstates that violate these constraints are not accessible and have zero probability.
```

(ch_sp_ss_why_most_probable)=
### Why Systems Settle in the Most Probable Macrostate

Given the fundamental postulate, why do we always observe the *most probable* macrostate?

The answer is pure arithmetic. In a system of $10^{23}$ particles, the most probable macrostate has so many more corresponding microstates than any other macrostate that the probability of ever observing a different macrostate is not just small — it is essentially zero.

To put it concretely: if you released a gas into one half of a box, the probability that all molecules spontaneously return to the left half is $2^{-10^{23}}$ — a number so tiny that even if you waited for a trillion times the age of the universe, it would almost certainly never happen.

This is the statistical origin of the **second law of thermodynamics**: systems evolve towards states of higher probability (higher entropy) and never spontaneously return to states of lower probability.

+++ { "page-break": true }
+++

(ch_sp_s_phase_space)=
## Phase Space: The Arena of Statistical Mechanics

(ch_sp_ss_what_is_phase_space)=
### What is Phase Space?

To count microstates for a continuous system (like a gas of real atoms), we need a mathematical arena. That arena is **phase space**.

For a single particle moving in three dimensions, we need 6 numbers to specify its state completely: three position coordinates $(x, y, z)$ and three momentum components $(p_x, p_y, p_z)$. Together, these define a single **point** in a 6-dimensional space — the phase space of one particle.

For a system of $N$ particles, the phase space is $6N$-dimensional. A **single point** in this high-dimensional space specifies the complete microstate of the entire system — where every particle is and how fast it is moving.

```{note}
Phase space is not a physical space — you cannot point to it. It is a mathematical construction that lets us visualise and count microstates. As the system evolves in time, its representative point traces a **trajectory** through phase space, never crossing itself (because the equations of motion have unique solutions).
```

(ch_sp_ss_energy_surface)=
### The Energy Surface and the Energy Shell

For an isolated system with fixed total energy $E$, all accessible microstates must satisfy:

$$H(\mathbf{q}, \mathbf{p}) = E$$

where $H$ is the Hamiltonian (total energy function). This condition defines a surface in phase space — the **energy surface** — and the system's trajectory is forever confined to it.

In practice, we always allow a tiny energy uncertainty $\Delta$ (much smaller than $E$), so we work with a thin shell between energies $E$ and $E + \Delta$ rather than an infinitely thin surface. The **volume** of this energy shell is the key quantity we need:

$$\Gamma(E, V, N) = \text{volume of phase space with } E < H(\mathbf{q},\mathbf{p}) < E + \Delta$$

This volume $\Gamma$ counts (in a precise sense) **how many microstates** the system has access to. The larger $\Gamma$, the more microstates are accessible, and — as we will see — the larger the entropy.

::::{tab-set}

:::{tab-item} One Particle in 1D
For a single particle of mass $m$ moving in one dimension with energy $E$, the Hamiltonian is $H = p^2/2m$. The "energy surface" is just two points: $p = \pm\sqrt{2mE}$. The energy shell has width $\Delta p = m\Delta E / \sqrt{2mE}$ at each point.
:::

:::{tab-item} One Particle in 3D
For a single particle in 3D, the energy surface in momentum space is a **sphere** of radius $R = \sqrt{2mE}$. The energy shell is a thin spherical shell of volume $4\pi R^2 \cdot \Delta p \propto E^{1/2}\,\Delta$.
:::

:::{tab-item} N Particles in 3D
For $N$ particles, each contributes 3 momentum components. The energy surface in the $3N$-dimensional momentum space is a hypersphere of radius $R = \sqrt{2mE}$. Its "surface area" grows as $E^{(3N/2)-1}$ — an astronomically large number for macroscopic $N$.
:::

::::

+++ { "page-break": true }
+++

(ch_sp_s_entropy)=
## Entropy: Counting Disorder

(ch_sp_ss_boltzmann)=
### The Boltzmann Definition

We have seen that some macrostates are more probable than others simply because they correspond to more microstates. Boltzmann had the brilliant idea of capturing this in a single number. In 1877, he proposed:

$$\boxed{S = k_B \ln \Omega}$$

where:
- $S$ is the **entropy** of the macrostate,
- $\Omega$ is the **number of microstates** compatible with that macrostate (proportional to the phase-space volume $\Gamma$),
- $k_B = 1.381 \times 10^{-23}$ J/K is **Boltzmann's constant**, which sets the scale.

This formula is engraved on Boltzmann's tombstone in Vienna. It is one of the most important equations in all of physics.

```{important}
Entropy is not a mysterious or abstract quantity. It is simply a measure of **how many microstates** are consistent with what we observe macroscopically.

- **High entropy** = many accessible microstates = the system could be in any of a large number of arrangements.
- **Low entropy** = few accessible microstates = the system is in one of only a small number of highly specific arrangements.
```

(ch_sp_ss_why_log)=
### Why the Logarithm?

The logarithm in $S = k_B \ln \Omega$ is not arbitrary. It is chosen for a crucial physical reason: **entropy must be additive**.

If you have two independent systems A and B, their combined entropy should be $S_{A+B} = S_A + S_B$ (total entropy is the sum of parts). But the total number of microstates for two independent systems is the *product*: $\Omega_{A+B} = \Omega_A \times \Omega_B$.

The logarithm converts this product into a sum:

$$S_{A+B} = k_B \ln(\Omega_A \times \Omega_B) = k_B \ln \Omega_A + k_B \ln \Omega_B = S_A + S_B \;\checkmark$$

Any other function of $\Omega$ would fail this requirement.

(ch_sp_ss_entropy_everyday)=
### Entropy in Everyday Terms

Let us build some intuition with a few everyday examples before turning to physics.

::::{tab-set}

:::{tab-item} A Messy Room
A tidy room has a very small number of arrangements (everything in exactly the right place). A messy room has an enormous number of arrangements (socks could be on the floor, the desk, the chair...). Entropy is higher for the messy room. The second law says the room will spontaneously become messier over time, never spontaneously tidier — exactly as anyone who has lived in student accommodation can confirm.
:::

:::{tab-item} Ink in Water
When you drop a spot of ink into a glass of water, it spreads out and mixes uniformly. The mixed state has vastly more microstates (each ink molecule could be anywhere in the full volume) than the unmixed state (each ink molecule confined to a small spot). The system moves from low entropy to high entropy. It never spontaneously unmixes.
:::

:::{tab-item} A Gas Expanding
When a gas expands from a small volume $V_1$ to a larger volume $V_2$, each molecule has more positions available to it. With $N$ molecules, the number of accessible microstates increases by $(V_2/V_1)^N$ — an unimaginably large factor for any macroscopic gas. The entropy increase is $\Delta S = Nk_B \ln(V_2/V_1)$.
:::

:::{tab-item} A Perfect Crystal at 0 K
A perfect crystal at absolute zero has just one accessible microstate: every atom is sitting at its precise lattice site, not vibrating. $\Omega = 1$, so $S = k_B \ln 1 = 0$. This is the **third law of thermodynamics**: the entropy of a perfect crystal at absolute zero is zero.
:::

::::

````{example} Entropy of mixing — visualised

The code below computes and plots how the entropy of a gas increases when its volume expands — a direct illustration of $S = k_B \ln \Omega$.

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt

kB = 1.381e-23   # J/K
N  = 1e20        # number of molecules

V1 = 1.0   # initial volume (arbitrary units)
V_ratios = np.linspace(1, 10, 300)

# Entropy change when volume increases from V1 to V2 = ratio * V1
delta_S = N * kB * np.log(V_ratios)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# ── Left: ΔS vs volume ratio ──────────────────────────────────────────────
axes[0].plot(V_ratios, delta_S / kB, color='steelblue', lw=2.5)
axes[0].fill_between(V_ratios, delta_S / kB, alpha=0.15, color='steelblue')
axes[0].set_xlabel(r'Volume ratio $V_2 / V_1$', fontsize=12)
axes[0].set_ylabel(r'$\Delta S \;/\; k_B$', fontsize=12)
axes[0].set_title('Entropy Increase on Expansion\n'
                  r'$\Delta S = N k_B \ln(V_2/V_1)$', fontsize=12)
axes[0].grid(True, alpha=0.3)

# ── Right: Cartoon — accessible phase space grows ─────────────────────────
ax = axes[1]
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Accessible microstates grow\nwhen volume increases',
             fontsize=12)

# Small box (before expansion)
small = plt.Rectangle((0.5, 1.5), 2, 2, linewidth=2,
                       edgecolor='steelblue', facecolor='#dbeafe')
ax.add_patch(small)
ax.text(1.5, 2.5, r'$V_1$', ha='center', va='center',
        fontsize=14, color='steelblue')
ax.text(1.5, 1.2, r'$\Omega_1$ microstates', ha='center',
        fontsize=9, color='steelblue')

# Arrow
ax.annotate('', xy=(4.5, 2.5), xytext=(3.0, 2.5),
            arrowprops=dict(arrowstyle='->', color='gray', lw=2))
ax.text(3.75, 2.8, 'expand', ha='center', fontsize=10, color='gray')

# Large box (after expansion)
large = plt.Rectangle((4.5, 0.5), 4, 4, linewidth=2,
                       edgecolor='tomato', facecolor='#fee2e2')
ax.add_patch(large)
ax.text(6.5, 2.5, r'$V_2 > V_1$', ha='center', va='center',
        fontsize=14, color='tomato')
ax.text(6.5, 0.2, r'$\Omega_2 \gg \Omega_1$ microstates', ha='center',
        fontsize=9, color='tomato')

plt.tight_layout()
plt.show()
```
````

+++ { "page-break": true }
+++

(ch_sp_s_temperature)=
## Temperature: What Does It Really Mean?

(ch_sp_ss_two_systems)=
### Heat Flow Between Two Systems

We now ask: if two systems are placed in thermal contact and allowed to exchange energy, what determines how they share that energy at equilibrium?

Consider two systems, 1 and 2, with energies $E_1$ and $E_2$ and entropies $S_1$ and $S_2$. Their total energy is fixed: $E_1 + E_2 = E = \text{const}$. Energy can flow back and forth, but the total is conserved.

The combined system will settle into the **most probable state**, which means the state with the **maximum total number of microstates**, i.e., the state that maximises the total entropy $S_1 + S_2$.

When we work out the mathematics of this maximisation, we find that at equilibrium:

$$\left(\frac{\partial S_1}{\partial E_1}\right)_{V,N} = \left(\frac{\partial S_2}{\partial E_2}\right)_{V,N}$$

In other words, equilibrium requires that both systems have **the same rate of change of entropy with energy**. This common value is what we define as the inverse temperature:

$$\boxed{\frac{1}{T} \equiv \left(\frac{\partial S}{\partial E}\right)_{V,N}}$$

(ch_sp_ss_temperature_meaning)=
### The Physical Meaning of Temperature

This definition might look abstract, but it has a beautifully simple physical interpretation.

$\partial S / \partial E$ measures how rapidly the number of accessible microstates *increases* when you add a little energy to the system. Think of it as the system's **appetite for energy**:

- A system with a **large** $\partial S/\partial E$ (low $T$) has a strong appetite — adding energy dramatically increases the number of accessible microstates.
- A system with a **small** $\partial S/\partial E$ (high $T$) has a weak appetite — it already has so many accessible microstates that adding more energy barely changes the count.

**Heat flows from high $T$ to low $T** because it increases the total entropy:**

When a hot system (small $\partial S/\partial E$) transfers energy to a cold system (large $\partial S/\partial E$), the cold system gains more entropy than the hot system loses. The total entropy increases. This process is spontaneous. The reverse — heat flowing from cold to hot — would *decrease* total entropy and never happens spontaneously.

```{note}
For an ideal gas, this definition gives $T = \frac{2E}{3Nk_B}$, which means $k_BT \approx \frac{2}{3}$ of the average kinetic energy per particle. Temperature is therefore a measure of the **average kinetic energy of the particles** — the familiar result from kinetic theory, but now derived from first principles.
```

(ch_sp_ss_zeroth_law)=
### The Zeroth Law of Thermodynamics

The definition of temperature as the quantity equalised at thermal equilibrium immediately implies the **zeroth law of thermodynamics**:

```{important}
**Zeroth Law:** If system A is in thermal equilibrium with system B, and system B is in thermal equilibrium with system C, then A and C are also in thermal equilibrium with each other.

This makes temperature a well-defined, transitive property — the basis for the very concept of a thermometer.
```

+++ { "page-break": true }
+++

(ch_sp_s_pressure_first_law)=
## Pressure, Chemical Potential, and the First Law

(ch_sp_ss_pressure)=
### Where Does Pressure Come From?

Just as temperature was defined by asking what is equalised when two systems can exchange *energy*, we can define **pressure** by asking what is equalised when two systems can exchange *volume*.

Imagine two chambers separated by a movable, thermally insulating piston. The total volume $V_1 + V_2 = V$ is fixed, but the piston can slide, exchanging volume between the two sides. The system maximises total entropy, which gives the equilibrium condition:

$$\left(\frac{\partial S_1}{\partial V_1}\right)_{E,N} = \left(\frac{\partial S_2}{\partial V_2}\right)_{E,N}$$

The common value of this derivative defines pressure through $P/T = (\partial S/\partial V)_{E,N}$. Microscopically, pressure arises because moving the wall gives the particles more space to occupy — more accessible microstates — and the system pushes in the direction that increases its entropy.

(ch_sp_ss_chemical_potential)=
### Chemical Potential

Similarly, if two systems can exchange **particles**, the equilibrium condition involves the **chemical potential** $\mu$, defined through:

$$-\frac{\mu}{T} \equiv \left(\frac{\partial S}{\partial N}\right)_{E,V}$$

Particles flow from regions of high chemical potential to low chemical potential, just as heat flows from high temperature to low temperature.

(ch_sp_ss_first_law)=
### The First Law of Thermodynamics

Putting together the definitions of temperature, pressure, and chemical potential, we can write a single equation that expresses how the internal energy $U$ changes when the macroscopic state of the system changes:

$$\boxed{dU = T\,dS - P\,dV + \mu\,dN}$$

This is the **first law of thermodynamics** — the conservation of energy — but written in a form that makes explicit the three ways energy can change:

::::{tab-set}

:::{tab-item} Heat: $T\,dS$
Adding heat $dQ = T\,dS$ to the system increases its internal energy by increasing the entropy (increasing the number of accessible microstates). This is energy transferred *randomly* at the microscopic level.
:::

:::{tab-item} Work: $-P\,dV$
Compressing the system (decreasing $V$) does work $dW = -P\,dV$ on it, increasing its energy. This is energy transferred *systematically* — by pushing the wall.
:::

:::{tab-item} Particles: $\mu\,dN$
Adding particles changes the energy by $\mu\,dN$. In many problems $N$ is fixed ($dN = 0$) and this term drops out, giving the simpler form $dU = T\,dS - P\,dV$.
:::

::::

```{important}
The first law is a statement of **energy conservation**. Energy cannot be created or destroyed — it can only be converted between heat, work, and the energy carried by particles. Statistical mechanics gives this law a microscopic foundation: it emerges from the definitions of $T$, $P$, and $\mu$ through the entropy function $S(E, V, N)$.
```

+++ { "page-break": true }
+++

(ch_sp_s_second_law)=
## The Second Law: Entropy Always Increases

(ch_sp_ss_statement)=
### The Statement

The second law of thermodynamics is one of the most profound statements in all of science:

```{important}
**Second Law of Thermodynamics:** The total entropy of an isolated system can never decrease. For any spontaneous process:

$$\Delta S_{\text{total}} \geq 0$$

Equality holds only for perfectly reversible (ideal, quasi-static) processes.
```

From the statistical viewpoint, this is not mysterious at all. An isolated system evolves randomly through its accessible microstates. There are enormously many more microstates corresponding to the high-entropy macrostate than to any low-entropy macrostate. So the system almost always moves towards higher entropy — not because of any driving force, but simply because of the overwhelming weight of numbers.

(ch_sp_ss_arrow_of_time)=
### The Arrow of Time

The second law gives time its **direction**. The microscopic laws of physics (Newton's laws, quantum mechanics) are perfectly symmetric in time — if you filmed a single particle bouncing, you could not tell if the film was playing forwards or backwards. But if you filmed a gas expanding into a vacuum and played it backwards, the absurdity would be immediately obvious.

This asymmetry is entirely statistical. The expanding gas visits ever more probable macrostates. The reversed film would show a journey towards an astronomically improbable macrostate — possible in principle, but never observed in practice.

(ch_sp_ss_free_expansion)=
### Free Expansion as an Example

A gas occupies the left half of a box. The right half is empty (vacuum). The partition is suddenly removed.

- **Before:** all $N$ molecules confined to volume $V/2$. Few accessible microstates.
- **After:** all $N$ molecules free to move in volume $V$. Each molecule has twice as many positions available. The number of accessible microstates increases by $(2)^N = 2^{10^{23}}$ — an unimaginably enormous factor.

The entropy increase is:

$$\Delta S = Nk_B \ln 2$$

For one mole of gas, this is about $5.8$ J/K — a macroscopically measurable entropy increase. The gas never spontaneously contracts back to the left half.

````{example} Free expansion — entropy change with volume

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt

kB   = 1.381e-23    # J/K
N_A  = 6.022e23     # Avogadro's number
n_moles = np.array([0.001, 0.01, 0.1, 1.0, 10.0])  # moles of gas
N_particles = n_moles * N_A

V_ratios = np.linspace(1.0, 4.0, 300)

fig, ax = plt.subplots(figsize=(8, 5))

colors = ['#93c5fd', '#60a5fa', '#3b82f6', '#1d4ed8', '#1e3a8a']
for n_mol, N, color in zip(n_moles, N_particles, colors):
    dS = N * kB * np.log(V_ratios)   # in J/K
    ax.plot(V_ratios, dS, color=color, lw=2,
            label=f'{n_mol} mol  ($N={n_mol:.3g}N_A$)')

ax.axvline(2.0, color='gray', linestyle='--', alpha=0.6, lw=1.2)
ax.text(2.05, ax.get_ylim()[0] * 0.9 if ax.get_ylim()[0] > 0 else 0.1,
        'Doubles\n(free expansion\ninto equal vacuum)',
        fontsize=8.5, color='gray')

ax.set_xlabel(r'Volume ratio $V_{\rm final} / V_{\rm initial}$', fontsize=12)
ax.set_ylabel(r'Entropy change $\Delta S$  (J/K)', fontsize=12)
ax.set_title('Entropy Increase During Expansion', fontsize=13)
ax.legend(fontsize=9, loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```
````

+++ { "page-break": true }
+++

(ch_sp_s_ideal_gas)=
## The Classical Ideal Gas

Now we apply the full statistical mechanics framework to the most important example: the **classical ideal gas**. This is a gas of $N$ identical, non-interacting particles in a box of volume $V$, with total energy $E$.

"Non-interacting" means the particles do not push or pull on each other — they only feel the walls of the box. This is an excellent approximation for low-density gases at high temperature (like ordinary air at room conditions).

(ch_sp_ss_counting)=
### Counting Microstates

For the ideal gas, the total energy is just the sum of kinetic energies of all particles:

$$E = \sum_{i=1}^{N} \frac{\mathbf{p}_i^2}{2m}$$

where $\mathbf{p}_i$ is the momentum of particle $i$ and $m$ is its mass.

The constraint $H = E$ defines a hypersphere in the $3N$-dimensional momentum space with radius $R = \sqrt{2mE}$. Multiplying by the position space volume $V^N$ (each of the $N$ particles can be anywhere in the box), the total accessible phase-space volume is:

$$\Phi(E, V, N) = V^N \times \left(\text{volume of a }3N\text{-dimensional sphere of radius }\sqrt{2mE}\right)$$

The volume of a hypersphere in $n$ dimensions with radius $R$ is $\frac{\pi^{n/2}}{\Gamma(n/2+1)} R^n$, giving:

$$\Phi(E, V, N) = V^N \cdot \frac{\pi^{3N/2}}{\Gamma\!\left(\frac{3N}{2}+1\right)} \cdot (2mE)^{3N/2}$$

(ch_sp_ss_sackur_tetrode)=
### The Sackur-Tetrode Equation

Taking the logarithm and using the correct normalisation $\Gamma_0 = h^{3N} N!$ (more on the $N!$ shortly), the entropy of the ideal gas works out to be:

$$\boxed{S(E, V, N) = Nk_B \left[\,\ln\!\left(\frac{V}{N}\left(\frac{4\pi m E}{3Nh^2}\right)^{3/2}\right) + \frac{5}{2}\right]}$$

This is the **Sackur-Tetrode equation** (1912), independently derived by Hugo Sackur and Otto Tetrode. It is the complete entropy of a monatomic classical ideal gas.

From this single equation, all thermodynamic properties of the ideal gas follow by differentiation.

(ch_sp_ss_equations_of_state)=
### Everything Follows from Entropy

::::{tab-set}

:::{tab-item} Temperature
Differentiating $S$ with respect to $E$ at fixed $V$ and $N$:

$$\frac{1}{T} = \frac{\partial S}{\partial E}\bigg|_{V,N} \implies \boxed{U = E = \frac{3}{2}Nk_BT}$$

Each particle contributes $\frac{3}{2}k_BT$ to the total energy — $\frac{1}{2}k_BT$ for each of its three translational degrees of freedom. This is the **equipartition theorem**.
:::

:::{tab-item} Pressure
Differentiating $S$ with respect to $V$ at fixed $E$ and $N$:

$$\frac{P}{T} = \frac{\partial S}{\partial V}\bigg|_{E,N} \implies \boxed{PV = Nk_BT}$$

This is the **ideal gas law** — derived entirely from counting microstates, with no empirical input!
:::

:::{tab-item} Heat Capacity
Since $U = \frac{3}{2}Nk_BT$, the heat capacity at constant volume is:

$$\boxed{C_V = \frac{\partial U}{\partial T}\bigg|_{V,N} = \frac{3}{2}Nk_B}$$

This means it takes $\frac{3}{2}k_B \approx 2.07 \times 10^{-23}$ J of heat to raise the temperature of each molecule by 1 K.
:::

:::{tab-item} Speed of Sound
Using $C_P = C_V + Nk_B = \frac{5}{2}Nk_B$ and $\gamma = C_P/C_V = 5/3$, the speed of sound in the gas is:

$$v_s = \sqrt{\frac{\gamma k_BT}{m}} = \sqrt{\frac{5k_BT}{3m}}$$

For nitrogen at room temperature ($T = 293$ K, $m = 4.65\times10^{-26}$ kg), this gives $v_s \approx 350$ m/s — in good agreement with the measured value of 343 m/s.
:::

::::

````{example} The ideal gas law from statistical mechanics

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt

kB  = 1.381e-23   # J/K
N   = 6.022e23    # 1 mole of molecules
m   = 4.65e-26    # kg (N2)

T_vals = np.linspace(100, 1000, 300)   # K
V      = 22.4e-3                        # m³ (molar volume at STP)

# Ideal gas law: P = N kB T / V
P_vals = N * kB * T_vals / V

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# ── Left: P vs T (ideal gas law) ──────────────────────────────────────────
axes[0].plot(T_vals, P_vals / 1e5, color='royalblue', lw=2.5)
axes[0].axvline(273, color='gray', linestyle='--', lw=1.2, alpha=0.7,
                label='$T = 273$ K (STP)')
axes[0].axhline(1.0, color='tomato', linestyle='--', lw=1.2, alpha=0.7,
                label='$P = 1$ atm')
axes[0].set_xlabel('Temperature $T$ (K)', fontsize=12)
axes[0].set_ylabel('Pressure $P$ (bar)', fontsize=12)
axes[0].set_title('Ideal Gas Law: $PV = Nk_BT$\n'
                  '(1 mole, $V = 22.4$ L fixed)', fontsize=12)
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.3)

# ── Right: Average kinetic energy per particle ────────────────────────────
KE_per_particle = 1.5 * kB * T_vals   # J
axes[1].plot(T_vals, KE_per_particle / kB, color='seagreen', lw=2.5,
             label=r'$\langle \epsilon \rangle = \frac{3}{2}k_BT$')
axes[1].set_xlabel('Temperature $T$ (K)', fontsize=12)
axes[1].set_ylabel(r'$\langle \epsilon \rangle / k_B$  (K)', fontsize=12)
axes[1].set_title('Equipartition Theorem\n'
                  r'Average KE per particle $= \frac{3}{2}k_BT$', fontsize=12)
axes[1].legend(fontsize=11)
axes[1].grid(True, alpha=0.3)

plt.suptitle('Classical Ideal Gas — Statistical Mechanics Results',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.show()

# Print numerical checks
print(f"At T = 273 K, P = {N*kB*273/V/1e5:.4f} bar  (should be ~1.00)")
print(f"Average KE per N2 molecule at 300 K: "
      f"{1.5*kB*300:.3e} J  = {1.5*kB*300/(1.602e-19)*1000:.2f} meV")
```
````

+++ { "page-break": true }
+++

(ch_sp_s_gibbs_paradox)=
## The Gibbs Paradox and Indistinguishability

(ch_sp_ss_the_puzzle)=
### A Surprising Puzzle

Here is a thought experiment that reveals something deep about the nature of identical particles.

Take a box divided by a partition. On the left: $N/2$ molecules of gas A at temperature $T$ and pressure $P$. On the right: $N/2$ molecules of gas B, *also* at temperature $T$ and pressure $P$. Now remove the partition.

If A and B are **different gases** (say, nitrogen and oxygen), the gases mix, and the entropy increases. This is the **entropy of mixing** — a real, measurable effect. You can verify it by noting that you now need to do work to separate them again.

But what if A and B are **the same gas** — say, nitrogen on both sides? When you remove the partition, nothing observable happens. The gas on both sides was already at the same temperature, pressure, and density. The system is in equilibrium both before and after.

**The entropy should not change.** And yet if you naively count microstates by treating the molecules as distinguishable labelled particles (as classical mechanics suggests), you find a spurious entropy increase of $\Delta S = Nk_B \ln 2 > 0$ — as if mixing two identical gases were somehow different from having one big container of gas. This is the **Gibbs paradox**, noticed by Josiah Willard Gibbs in 1875.

(ch_sp_ss_resolution)=
### The Resolution: Particles Are Indistinguishable

The resolution is both simple and profound: **identical particles are not distinguishable**.

In quantum mechanics, two electrons, two nitrogen molecules, or any two identical particles are not merely similar — they are in principle completely indistinguishable. There is no experiment, even in principle, that can tell "molecule number 347" from "molecule number 8,193" if both are nitrogen molecules in the same quantum state.

This means that if you swap two identical particles, you have not created a new microstate — you have the **same** microstate. To correctly count microstates, we must divide by the number of ways we can permute $N$ identical particles, which is $N!$ (N factorial):

$$\Omega_{\text{correct}} = \frac{\Omega_{\text{naive}}}{N!}$$

With this correction, the entropy becomes:

$$S_{\text{correct}}(E, V, N) = Nk_B \left[\ln\left(\frac{V}{N}\cdot \left(\frac{4\pi mE}{3Nh^2}\right)^{3/2}\right) + \frac{5}{2}\right]$$

Notice that the volume appears as $V/N$ (volume per particle), not $V$. This makes the entropy **extensive** — proportional to $N$ when $E/N$ and $V/N$ are held fixed — and the mixing entropy for two identical gases is exactly zero.

```{important}
The Gibbs paradox reveals that statistical mechanics cannot be done purely classically. Even when the particles themselves move according to classical mechanics (which is fine at high temperature and low density), the **counting** of microstates must respect quantum indistinguishability. The factor $N!$ is a quantum input into an otherwise classical calculation.

Without it, the entropy is not extensive, and thermodynamics breaks down.
```

(ch_sp_ss_mixing_summary)=
### What Happens for Different Gases?

When genuinely different gases mix, the entropy increase is real and is given by:

$$\Delta S_{\text{mix}} = -Nk_B\bigl(x_1 \ln x_1 + x_2 \ln x_2\bigr)$$

where $x_1 = N_1/N$ and $x_2 = N_2/N$ are the mole fractions. For equal amounts ($x_1 = x_2 = 1/2$), this gives $\Delta S = Nk_B \ln 2 > 0$. This entropy of mixing is the thermodynamic reason why spontaneous mixing of different gases is irreversible.

````{example} Gibbs paradox visualised

```{code-cell} Python
:tag: hide-input

import numpy as np
import matplotlib.pyplot as plt

kB = 1.381e-23
N  = 1e20   # particles

x1_vals = np.linspace(0.001, 0.999, 300)
x2_vals = 1.0 - x1_vals

# Entropy of mixing for DIFFERENT gases (real, positive)
dS_mix_different = -N * kB * (x1_vals * np.log(x1_vals) +
                               x2_vals * np.log(x2_vals))

# Spurious entropy for SAME gas WITHOUT N! correction
# ΔS_wrong = N kB [x1 ln(1/x1) + x2 ln(1/x2)]  same formula!
dS_wrong_same = -N * kB * (x1_vals * np.log(x1_vals) +
                             x2_vals * np.log(x2_vals))

# Correct entropy for SAME gas WITH N! correction = 0
dS_correct_same = np.zeros_like(x1_vals)

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x1_vals, dS_mix_different / (N * kB), lw=2.5, color='royalblue',
        label='Different gases (correct: $\\Delta S > 0$)')
ax.plot(x1_vals, dS_wrong_same / (N * kB), lw=2, color='tomato',
        linestyle='--', label='Same gas without $N!$ (wrong: $\\Delta S > 0$)')
ax.plot(x1_vals, dS_correct_same, lw=2.5, color='seagreen',
        label='Same gas with $N!$ (correct: $\\Delta S = 0$)')

ax.set_xlabel(r'Mole fraction $x_1$', fontsize=12)
ax.set_ylabel(r'$\Delta S_{\rm mix}\;/\;Nk_B$', fontsize=12)
ax.set_title('Gibbs Paradox: Entropy of Mixing', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1)
ax.set_ylim(-0.05, 0.75)

ax.annotate('Gibbs paradox:\nspurious entropy\nfor identical gases',
            xy=(0.5, 0.693), xytext=(0.65, 0.45),
            fontsize=9, color='tomato',
            arrowprops=dict(arrowstyle='->', color='tomato', lw=1.5))

plt.tight_layout()
plt.show()
```
````

+++ { "page-break": true }
+++

(ch_sp_s_ensembles_preview)=
## A Preview: The Three Ensembles

The microcanonical ensemble — the framework we have built throughout this chapter — describes a system that is **completely isolated**: fixed $E$, $V$, and $N$. It is the most fundamental ensemble, but it is often the most difficult to work with in practice, because real experiments rarely fix the energy exactly.

Statistical mechanics has developed two other ensembles, each suited to a different experimental situation. All three give the same thermodynamic results in the limit of large $N$.

```{list-table} The Three Ensembles of Statistical Mechanics
:header-rows: 1
:widths: 20 20 20 20 20

* - Ensemble
  - Fixed quantities
  - Fluctuating quantity
  - Thermodynamic potential
  - Best suited for
* - **Microcanonical**
  - $E$, $V$, $N$
  - (nothing)
  - Entropy $S$
  - Isolated systems; conceptual foundations
* - **Canonical**
  - $T$, $V$, $N$
  - Energy $E$
  - Free energy $F = U - TS$
  - Systems in contact with a heat bath; most practical calculations
* - **Grand Canonical**
  - $T$, $V$, $\mu$
  - Energy $E$ and particle number $N$
  - Grand potential $\Omega = F - \mu N$
  - Open systems; quantum gases; chemical reactions
```

The microcanonical ensemble answers the question: *given fixed energy, what are the equilibrium properties?*

The canonical ensemble answers: *given fixed temperature (heat bath), what are the equilibrium properties?*

The grand canonical ensemble answers: *given fixed temperature and chemical potential (heat and particle reservoir), what are the equilibrium properties?*

```{note}
In the **thermodynamic limit** ($N \to \infty$), all three ensembles give identical predictions for macroscopic quantities. They differ only in how they handle fluctuations: the canonical ensemble allows energy to fluctuate (but with tiny relative fluctuations $\sim 1/\sqrt{N}$), while the grand canonical ensemble allows both energy and particle number to fluctuate.
```

+++ { "page-break": true }
+++

(ch_sp_s_summary)=
## Summary

In this chapter we have laid the foundations of statistical physics. The central ideas to take away are:

::::{tab-set}

:::{tab-item} The Big Picture
Statistical physics bridges the microscopic world of atoms and the macroscopic world of thermodynamics. It does so through probability: by counting microstates, we identify the most probable macrostate, which is overwhelmingly more probable than any other for large systems. This is why thermodynamics is so precise.
:::

:::{tab-item} Key Definitions
- A **microstate** is the complete atomic-level description. A **macrostate** is the small set of macroscopic variables we measure.
- The **fundamental postulate** says all accessible microstates are equally likely.
- **Entropy** is $S = k_B \ln \Omega$ — a measure of how many microstates correspond to a given macrostate.
- **Temperature** is $1/T = \partial S / \partial E$ — the quantity equalised at thermal equilibrium.
- **Pressure** is $P/T = \partial S/\partial V$ — the quantity equalised when volumes can adjust.
:::

:::{tab-item} Key Results
- The **first law** $dU = T\,dS - P\,dV + \mu\,dN$ follows directly from these definitions.
- The **second law** $\Delta S \geq 0$ is a statement of overwhelming probability, not a fundamental law of nature.
- The **ideal gas law** $PV = Nk_BT$ and the **equipartition theorem** $U = \frac{3}{2}Nk_BT$ are derived entirely from counting microstates.
- The **Gibbs paradox** is resolved by recognising that identical particles are indistinguishable, requiring the $N!$ correction.
:::

:::{tab-item} What Comes Next
The microcanonical ensemble is conceptually fundamental but computationally awkward. In practice, we usually work with the **canonical ensemble** (fixed $T$, $V$, $N$), which is mathematically more tractable and directly applicable to systems in contact with a heat bath — which describes almost all real laboratory situations.
:::

::::

```{important} Recommended Reading

- F. Reif, *Fundamentals of Statistical and Thermal Physics* — the most accessible and physically intuitive introduction.
- C. Kittel & H. Kroemer, *Thermal Physics* — clear, concept-first approach, excellent for undergraduates.
- R. K. Pathria & P. D. Beale, *Statistical Mechanics* — comprehensive graduate-level reference.
- M. Kardar, *Statistical Physics of Particles* (Cambridge) — modern, elegant, with excellent physical discussions.
- L. D. Landau & E. M. Lifshitz, *Statistical Physics, Part 1* — terse but deep; rewarding for advanced students.
- L. Bocquet, *Lecture Notes on Statistical Physics* (ENS-PSL) — very clear modern notes, freely available.
```

(ch_sp_s_exercises)=
## Exercises

```{exercise} Counting microstates
:label: ex_sp_counting

You have a system of 3 distinguishable particles. Each particle can be in one of 4 equally probable energy levels: $\epsilon_0, \epsilon_1, \epsilon_2, \epsilon_3$.

1. How many total microstates are there?
2. List all the microstates in which the total energy is $E = 2\epsilon_1$ (assume $\epsilon_n = n\epsilon_1$).
3. If the fundamental postulate holds, what is the probability of each macrostate (defined by total energy)?
4. Which macrostate is most probable?
```

```{exercise} Entropy and information
:label: ex_sp_entropy_info

A system can be in one of $\Omega$ equally probable microstates.

1. If $\Omega = 1$, what is the entropy? What does this tell you physically?
2. If $\Omega$ doubles (e.g., the volume doubles for an ideal gas), by how much does $S$ increase?
3. Two independent systems have $\Omega_A = 10^{20}$ and $\Omega_B = 10^{30}$. What is the total entropy of the combined system in units of $k_B$?
4. Why is it the logarithm of $\Omega$ that appears in the entropy, rather than $\Omega$ itself?
```

```{exercise} Temperature from entropy
:label: ex_sp_temp

A simple model system has entropy $S = aN^{2/3} E^{1/3}$ where $a$ is a positive constant, $N$ is the number of particles, and $E$ is the total energy.

1. Find the temperature $T$ as a function of $E$ and $N$.
2. Invert this to find the energy $E$ as a function of $T$ and $N$.
3. Find the heat capacity $C_V = \partial E / \partial T$ at fixed $N$.
4. Is the heat capacity extensive (proportional to $N$)? Should it be?
```

```{exercise} Heat flow and entropy
:label: ex_sp_heat_flow

System A has $N_A = 10^{23}$ particles at temperature $T_A = 400$ K. System B has $N_B = 10^{23}$ particles at temperature $T_B = 200$ K. Both are monatomic ideal gases.

1. What is the total energy $E = E_A + E_B$?
2. When the two systems are brought into thermal contact, what is the final equilibrium temperature $T_f$? (*Hint*: total energy is conserved; assume equal heat capacities.)
3. Calculate the entropy change of system A, system B, and the total entropy change when equilibrium is reached. Use $\Delta S = Nk_B \ln(T_f/T_i)$ for each system (derive this if you can).
4. Confirm that $\Delta S_{\text{total}} > 0$.
```

```{exercise} Free expansion
:label: ex_sp_free_expansion

One mole of an ideal gas (monatomic, $C_V = \frac{3}{2}R$) is confined to the left half of an insulated box of total volume $V = 10$ L. The right half is vacuum. The partition is removed.

1. What is the change in temperature? (*Hint*: no work is done and no heat is exchanged with the surroundings.)
2. What is the entropy change? Express your answer in J/K.
3. Is this process reversible? Explain using the concept of microstates.
4. What would you need to do to restore the gas to its original state?
```

```{exercise} The Gibbs paradox
:label: ex_sp_gibbs

A box of volume $V$ is divided into two equal halves. The left half contains $N/2$ molecules of gas A and the right half contains $N/2$ molecules of gas B, both at the same temperature $T$ and pressure $P$.

1. If A and B are **different** gases, calculate the entropy of mixing when the partition is removed. Express your answer in terms of $N$ and $k_B$.
2. If A and B are **identical** gases and you treat the molecules as distinguishable, what entropy change do you get? Why is this result wrong?
3. Explain in words why the $N!$ correction resolves the paradox.
4. What is the physical content of the statement that "identical particles are indistinguishable"?
```

```{exercise} The ideal gas from scratch
:label: ex_sp_ideal_gas

Starting from the Sackur-Tetrode entropy:

$$S(E, V, N) = Nk_B \left[\ln\!\left(\frac{V}{N}\left(\frac{4\pi m E}{3Nh^2}\right)^{3/2}\right) + \frac{5}{2}\right]$$

1. Derive the ideal gas law $PV = Nk_BT$ by computing $P/T = (\partial S/\partial V)_{E,N}$.
2. Derive the caloric equation of state $E = \frac{3}{2}Nk_BT$ by computing $1/T = (\partial S/\partial E)_{V,N}$.
3. Show that $S$ is extensive: verify that $S(\lambda E, \lambda V, \lambda N) = \lambda S(E, V, N)$ for any positive $\lambda$.
4. Compute the entropy change when the temperature doubles at constant volume and particle number.
```

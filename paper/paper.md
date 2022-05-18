---
title: 'pypbomb: A Python package with tools for the design of detonation tubes'
tags:
  - detonation
authors:
  - name: Mick Carter
    affiliation: 1
  - name: David Blunck
    affiliation: 1
affiliations:
 - name: School of Mechanical, Industrial, and Manufacturing Engineering, Oregon State University, Corvallis, OR, USA
   index: 1

date: 20 January, 2021
bibliography: paper.bib
---

# Statement of need 

A detonation is supersonic combustion in which a reaction front is coupled with a pressure shock front [@Lee2008]. Adiabatic heating from the shock front helps to sustain the combustion, which in turn accelerates the shock front and allows it to propagate supersonically. Because the products of a detonation are at a higher pressure than the reactants, thermodynamic cycles using detonations, such as the Humphrey cycle, have the potential for higher thermodynamic efficiency than deflagration-based cycles, such as the Brayton cycle [@Coleman2001]. However, knowledge of a reactant mixture's characteristic detonation structure is needed in order to make practical use of it with detonation-based cycles.

The behavior of a detonation can be investigated either via simulation or physical experimentation. Given that detonation is an inherently three-dimensional process [@Lee2008; @ciccarelli], accurate simulations can be computationally expensive [@Kessler2010], requiring runs on the order of days to months for a single 2-D simulation depending on the hardware involved [@Kessler2011; @Radulescu2007]. Alternatively, a detonation tube is often used to experimentally study the behavior of detonations. In a detonation tube, data for a given set of initial conditions can be collected anywhere from multiple times per second (e.g., in a pulse detonation engine) to multiple times per hour (e.g., in a closed-end detonation tube).

The design of a detonation tube requires many considerations, including estimation of the required length for deflagration-to-detonation transition (DDT), tube material and size selection (including the effects of transient pressures), fastener failure calculations (including bolt pull-out), flange class selection, viewing window sizing (if optical access is required), and prediction of safe operating conditions (including accounting for detonation reflection in the case of a closed-end tube). All of this is specific to the mixture being detonated, therefore it is important to be able to quickly re-perform the analysis for new reactant mixtures. The tools within ``pypbomb`` provide a first-order analysis meant to serve as the basis for the previously detailed analysis.

# Summary

``pypbomb`` contains a series of tools that can be used to iterate on initial design parameters and obtain an estimate for the operational envelope of a closed-end detonation tube. This package is not meant to replace the design process entirely, but rather to provide an initial design of a detonation tube through a series of simplified analyses. A conservative safety factor is recommended given the dynamic nature of detonations, as well as a more in-depth analysis of the tube's individual components; the analysis in this package assumes steady state propagation at equilibrium conditions and does not account for transient effects other than through a dynamic load factor. The first iteration of this package was written during the design of a detonation tube that has been used to measure cell sizes of gaseous detonations, and will be used in the near future for the study of detonations in two-phase mixtures.

`pypbomb.Tube` allows the user to quickly iterate on the design of the piping portion of a detonation tube and determine its safe operational limits. Nominal Pipe Size (NPS) lookups are included, allowing the user to quickly assess different tube geometries. Once the tube geometry is determined, the stress limits of the pipe material are used to evaluate the maximum allowable static pressure within the tube. For convenience, maximum allowable stress values of stainless steels are available as a function of temperature [@asmeb311]. Alternative allowable stress values, such as those calculated using the ASME Boiler and Pressure Vessel code, may also be manually supplied by the user. However, it is important to keep in mind that no code currently exists for the design of detonation tubes (because of the dynamic transient pressures they experience as well as with their experimental nature), and any estimate using existing code should be accompanied by a more in-depth evaluation and/or a large safety factor. Once the user has determined the allowable stress, it is used to calculate the maximum operating pressure of the tube [@megyesy]. Accounting for the tube geometry and predicted mixture detonation behavior, a dynamic load factor is calculated and applied [@Shepherd2009]. The dynamic load factor is used to adjust the static pressure limit to account for the tube's response to the transient pressure caused by the detonation wave. The maximum initial reactant pressure is then iteratively determined for a given tube geometry, reactant mixture, and initial temperature using the maximum operating pressure. Detonation wave speeds and reflection properties for the desired reactant mixture are calculated using selected functions adapted from the shock and detonation toolbox (SDToolbox) [@sdt]. The curve fitting portion of the wave speed estimation in SDToolbox has been parallelized in order to speed up computation.

Once the operational limits of a tube are determined, flanges can be sized. `pypbomb.Flange` looks up the minimum necessary flange class based on the maximum tube pressure and temperature based on the standards set forth in ASME B16.5-2003 [@asmeb165]. Although they are used here to provide an estimate of minimum flange class, ASME codes do not account for impulsive loads such as those caused by detonations. Therefore further analysis should be conducted on a per-flange basis, using the recommended flange size as an initial design.

A successful detonation tube design must account for the deflagration-to-detonation transition (DDT). DDT is usually achieved using a series of blockages, which causes the combustion wave to undergo local accelerations, thereby aiding in the DDT process [@ciccarelli]. The blockages must be properly sized, and must continue for a minimum (mixture specific) distance in order to maximize the probability of a successful transition to detonation. To this end, `pypbomb.DDT` contains tools for Shchelkin spiral blockage ratio and diameter calculations, and allows the user to estimate the necessary DDT run-up length for a desired mixture using Cantera [@ciccarelli; @cantera]. For blockage ratios $BR \leq 0.1$ the run-up length, $X_{S}$, inside a tube of diameter $D$ is estimated to be

$$X_{S} = \frac{D \gamma}{C} \left[ \frac{1}{\kappa} \ln \left( \gamma \frac{D}{h} \right) + K \right]$$

where $\kappa=0.4$, $K=5.5$, $C=0.2$,

$$ \frac{D}{h} = \frac{2}{1 - \sqrt{1-BR}} $$

and

$$ \gamma = \left[ \frac{a_{p}}{\eta (\sigma - 1)^{2}S_{L}} \left( \frac{\delta}{D} \right)^{\frac{1}{3}} \right]^{\frac{1}{2m + 7/3}} $$

where $S_{L}$ is the mixture's laminar flame speed, $\delta = \nu / S_{L}$ is the mixture's laminar flame thickness, $\nu$ is the kinematic viscosity, $\sigma$ is the product/reactant density ratio, $\eta=2.1$, and $m=-0.18$ [@ciccarelli]. For blockage ratios $0.3 \leq BR \leq 0.75$, the run-up length is estimated to be

$$ X_{S} \approx a \frac{D a_{p}(1 - BR)}{20 S_{L} (1 + b BR)(\sigma - 1)} $$

where $a=2$, $b=1.5$, and $a_{p}$ is the speed of sound within the products [@ciccarelli]. Run-up length for blockage ratios $0.1 < BR < 0.3$ are estimated by linearly interpolating between the two formulas, each evaluated at the relevant endpoint. Blockage ratios greater than 0.75 are not considered.

Finally, `pypbomb` provides some tools to facilitate the inclusion of optical access in the detonation tube. Historically, the structure of detonations have typically been studied using soot covered foils inserted along the wall or end-cap of detonation tubes [@Lee2008]. More recently, however, researchers have begun using high speed photography to study detonation waves, including PLIF and focusing schlieren methods [@Pintgen2003; @Mevel2015; @Rankin2016; @Radulescu2007; @Stevens2015]. In some cases, soot foil and schlieren techniques have been used simultaneously [@Kellenberger2017]. If optical access is desired, window thickness and factor of safety calculations can be quickly performed for clamped rectangular windows using `pypbomb.Window` [@crystran]. These calculations do not account for loads applied to the window due to contact with the detonation tube or window retainers; it is critical that windows be isolated from contact with any hard surfaces. In our tube this was accomplished using rubber gaskets on the faces of the windows as well as around the periphery. In addition to window calculations, `pypbomb.Bolt` allows the user to estimate bolt stress areas and safety factors in order to keep the windows intact and prevent bolts from pulling out of the tube [@machinery].

# Acknowledgements

This work is supported by the Office of Naval Research, contract N000141612429. The authors would like to thank the Detonation Group at the Air Force Research Laboratory for providing access to a detonation tube, as well as for their invaluable input and technical advice. Additionally, the authors would like to thank Kyle Niemeyer for his instruction and encouragement in the area of open-source software development for engineering research. Finally, this package and the research that it precedes are enabled by the resources made available to the public by Joseph Shepherd's [Explosion Dynamics Laboratory](https://shepherd.caltech.edu/EDL/) at Caltech; a debt of gratitude is owed to everyone who has contributed to that group, and insightful comments about this manuscript by Dr. Shepherd are gratefully acknowledged.

# References
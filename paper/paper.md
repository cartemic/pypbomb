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

A detonation is supersonic combustion in which a reaction front is coupled with a pressure shock front [@Lee2008]. Adiabatic heating from the shocks helps sustain the combustion front, which in turn accelerates the shock front and allows it to propagate supersonically. Because the products of a detonation are at a higher pressure than the reactants, thermodynamic cycles using detonations, such as the Humphrey cycle, have the potential for higher thermodynamic efficiency than deflagration-based cycles, such as the Brayton cycle [@Coleman2001]. However, knowledge of the characteristic structural properties of a given reactant mixture is needed in order to make practical use of detonation cycles.

The behavior of a detonation can be investigated either via simulation or physical experimentation. Given that detonation is an inherently three-dimensional process [@Lee2008; @ciccarelli], accurate simulations can be computationally expensive [@Kessler2010], requiring runs on the order of days to months for a single 2-D simulation depending on the hardware involved [@Kessler2011; @Radulescu2007]. Alternatively, a detonation tube is often used to experimentally study the behavior of detonations. In a detonation tube, data for a given set of initial conditions can be collected anywhere from multiple times per second (e.g. in a pulse detonation engine) to multiple times per hour (e.g. in a closed-end detonation tube).

The design of a detonation tube requires many considerations, including estimation of the required length for deflagration-to-detonation transition (DDT), tube material and size selection (including the effects of transient pressures), fastener failure calculations (including bolt pull-out), flange class selection, viewing window sizing (if optical access is required), and prediction of safe operating conditions (including accounting for detonation reflection in the case of a closed-end tube). All of this is specific to the mixture being detonated, therefore it is important to be able to quickly re-perform the analysis for new reactant mixtures. ``pypbomb`` provides tools to address the previously detailed analysis and combines them into one easy-to-use package.

# Summary

``pypbomb`` contains a series of tools that can be used to design or evaluate the operational envelope of closed-end detonation tubes. The first iteration of this package was written during the design of just such a tube, which is currently being used to measure cell sizes of gaseous detonations, and will be used in the near future for the study of detonations in two-phase mixtures.

`pypbomb.Tube` allows the user to quickly iterate on the design of the piping portion of a detonation tube and determine its safe operational limits. Nominal Pipe Size (NPS) lookups are included, allowing the user to quickly assess different tube geometries. Once the tube geometry is determined, the stress limits of the pipe material are used to evaluate the maximum allowable static pressure within the tube. For convenience, maximum allowable stress values of stainless steels are available as a function of temperature [@asmeb311]; other values may also be manually supplied by the user, and the maximum operating pressure of the tube is calculated [@megyesy]. Accounting for the tube geometry and predicted mixture detonation behavior, Shepherd's dynamic load factor is calculated [@Shepherd2009]. The dynamic load factor is used to adjust the static pressure limit to account for the tube's response to the transient pressure caused by the detonation wave. The maximum initial reactant pressure is then iteratively determined for a given tube geometry, reactant mixture, and initial temperature using the reduced pressure. Detonation wave speeds and reflection properties for the desired reactant mixture are calculated using selected functions from the shock and detonation toolbox (SDToolbox) [@sdt]. The curve fitting portion of the wave speed estimation in SDToolbox has been parallelized in order to speed up computation.

Once the operational limits of a tube are determined, flanges can be sized. `pypbomb.Flange` looks up the minimum necessary flange class based on the maximum tube pressure and temperature based on the standards set forth in ASME B16.5-2003 [@asmeb165].

A successful detonation tube design must account for the deflagration-to-detonation transition (DDT). DDT is usually achieved using a series of blockages, which causes the combustion wave to undergo local accelerations, thereby aiding in the DDT process [@ciccarelli]. The blockages must be properly sized, and must continue for a minimum (mixture specific) distance in order to maximize the probability of a successful transition to detonation. To this end, `pypbomb.DDT` contains tools for Shchelkin spiral blockage ratio and diameter calculations, and allows the user to estimate the necessary DDT run-up length for a desired mixture using Cantera [@ciccarelli; @cantera].

One feature that researchers may desire in a detonation tube is optical access. Historically, the structure of detonations have typically been studied using soot covered foils inserted along the wall or end-cap of detonation tubes [@Lee2008]. More recently, however, researchers have begun using high speed photography to study detonation waves, including PLIF and focusing schlieren methods [@Mevel2016; @Rankin2016; @Radulescu2007; @Stevens2015]. In some cases, soot foil and schlieren techniques have been used simultaneously [@Kellenberger2017]. If optical access is desired, window thickness and safety factor calculations can be quickly performed for clamped rectangular windows using `pypbomb.Window` [@crystran]. Additionally, `pypbomb.Bolt` allows the user to calculate bolt stress areas and safety factors in order to keep the windows intact  and prevent bolts from pulling out of the tube [@machinery].

# Acknowledgements

This work is supported by the Office of Naval Research, contract N000141612429. The authors would like to thank the Detonation Group at the Air Force Research Laboratory for providing access to a detonation tube, as well as for their invaluable input and technical advice. Additionally, the authors would like to thank Kyle Niemeyer for his instruction and encouragement in the area of open-source software development for engineering research. Finally, this package and the research that it precedes rely heavily on the resources made available to the public by Joe Shepherd's [Explosion Dynamics Laboratory](https://shepherd.caltech.edu/EDL/) at Caltech, and an enormous debt of gratitude is owed to everyone who has contributed to that group.

# References
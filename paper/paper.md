---
title: 'HBAT2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures'
abstract: |
  Hydrogen bonds and other non-covalent interactions play a crucial role in maintaining the structural integrity and functionality of biological macromolecules such as proteins and nucleic acids. Accurate identification and analysis of hydrogen bonds and other non-covalent interactions are essential for understanding molecular interactions, protein folding, and drug design. HBAT (Hydrogen Bond Analysis Tool) is a widely used software for analysing hydrogen bonds and other weak interactions in macromolecular structures. In this paper, we present HBAT2, an updated Python reimplementation of the original HBAT tool.
tags:
  - Python
  - structural biology
  - hydrogen bonds
  - molecular interactions
  - protein structures
  - bioinformatics
authors:
  - name: Abhishek Tiwari
    orcid: 0000-0003-2222-2395
    affiliation: 1
affiliations:
 - name: Independent Researcher
   index: 1
date: 3 August 2025
bibliography: paper.bib
---

# Summary
HBAT2 is a Python package for automated analysis of hydrogen bonds and other non-covalent interactions in macromolecular structures, available in Protein Data Bank (PDB) file format. Originally developed in Perl/Tk and published in 2007 [@tiwari2007hbat], HBAT2 has been completely rewritten in Python. The software identifies and analyses traditional hydrogen bonds, weak hydrogen bonds, halogen bonds, X-H$\cdots$$\pi$, and $\pi$-$\pi$ stacking, and n$\rightarrow$$\pi$* interactions using geometric criteria. It also detects cooperativity and anticooperativity chains and renders them as 2D visualisations. The latest version of HBAT offers improved cross-platform `tkinter`-based graphical user interface (GUI), a simple command-line interface (CLI), and a developer-friendly API making it accessible to users with different computational backgrounds and needs.

# Statement of need
Hydrogen bonds and other non-covalent interactions are fundamental to protein structure, stability, and function. With over 200,000 structures in the Protein Data Bank [@berman2000protein], there is an increasing need for automated tools to analyse these interactions systematically.

The landscape of hydrogen bond analysis tools is diverse but fragmented. Classic tools like HBPLUS [@mcdonald1994satisfying] and HBexplore [@lindauer1996hbexplore] pioneered automated H-bond detection but lack modern interfaces and support for a diverse range of interactions. More recent tools serve specialized niches: PLIP [@salentin_plip_2015] and Arpeggio [@jubb_arpeggio_2017] excel at protein-ligand interactions but are web-based without standalone GUI options; HBonanza [@durrant_hbonanza_2011], HBCalculator [@wang_hbcalculator_2024], and BRIDGE2 [@siemers_interactive_2021] focus on molecular dynamics trajectories rather than static structures; MDAnalysis [@noauthor_431_nodate], GROMACS [@noauthor_gmx_nodate], and AMBER [@noauthor_hbond_2020] provide H-bond analysis within larger MD suites. Tools like VMD [@noauthor_vmd_nodate] and ChimeraX [@noauthor_tool_nodate] offer interactive hydrogen bond visualisation but limited statistical analysis. ProteinTools [@ferruz_proteintools_2021] provides web-based network analysis but lacks cross-platform desktop capabilities.

Despite this ecosystem, there remains a gap for a comprehensive and reliable cross-platform desktop tool that: (1) analyses diverse interaction types beyond canonical hydrogen bonds, (2) provides both graphical and command-line interfaces for different user workflows, (3) identifies potential cooperativity, (4) integrates seamlessly with the existing scientific ecosystem, and (5) supports flexible parameter customization with domain-specific presets. 

HBAT2 addresses these limitations by providing a modern, Python implementation that integrates seamlessly with contemporary structural biology workflows. The software is particularly valuable for researchers in structural biology, computational chemistry, and drug design who need detailed analysis of molecular interactions.

The original HBAT [@tiwari2007hbat] was developed in Perl/Tk with a Windows-only GUI, limiting its adoption in modern computational environments. This update addresses this limitation by making HBAT2 cross-platform with support for Windows, Linux, and Mac.

![The latest update to HBAT2 uses tkinter to provide a cross-platform graphical user interface (GUI)](https://static.abhishek-tiwari.com/hbat/hbat-window-v3.png)

# Key Enhancements

HBAT2 introduces several key improvements over the original 2007 version:

## Structure Preparation
HBAT2 uses PDBFixer [@noauthor_openmmpdbfixer_2025; @eastman_openmm_2013] and OpenBabel [@oboyle_pybel_2008] to automatically enhance macromolecular structures by adding missing atoms, converting residues, and cleaning up structural issues. These capabilities are particularly valuable when working with crystal structures missing hydrogen atoms, low-resolution structures with incomplete side chains, structures containing non-standard amino acid residues, and structures with unwanted ligands or contaminants. This integration improves the quality of analysis.

## Interaction Coverage

HBAT2 analyses a broad spectrum of interactions including hydrogen bonds (O-H$\cdots$O, N-H$\cdots$O, N-H$\cdots$N, C-H$\cdots$O), halogen bonds (C-X$\cdots$Y where X=F,Cl,Br,I), X-H$\cdots$$\pi$ interactions with aromatic systems, $\pi$-$\pi$ stacking [@mcgaughey_pi-stacking_1998; @vernon_pi-pi_nodate], carbonyl-carbonyl n$\rightarrow$$\pi$* interactions [@rahim_reciprocal_2017; @newberry_n_2017], and n$\rightarrow$$\pi$* interactions [@choudhary_nature_2009]. This comprehensive approach addresses the growing recognition of weak interactions' importance in protein structure and stability [@desiraju_weak_2001; @cavallo_halogen_2016; @brandl_c-h-interactions_2001].

## Cooperativity Chains
HBAT2 offers two ways to visualise hydrogen bond networks: NetworkX [@hagberg2008networkx]/Matplotlib [@hunter2007matplotlib] and GraphViz [@graphviz2024]. Unlike tools that provide only basic visualisation (VMD, ChimeraX) or focus solely on MD trajectory dynamics (BRIDGE2, HBonanza), HBAT2 emphasises cooperativity chains and network topology in static structures with customizable layouts and high-resolution export (PNG, SVG, PDF).

## Parameter Presets
Built-in parameter sets optimised for different experimental conditions (high-resolution X-ray, NMR, membrane proteins, drug design) address a common challenge in hydrogen bond analysis. While other tools require manual parameter specification, HBAT's preset system makes it accessible to experimental structural biologists while maintaining flexibility for computational experts.

## Flexible Output
Multiple export formats (text, CSV, JSON) enable integration with downstream analysis pipelines and statistical software. Combined with optional PDBFixer and OpenBabel integration for automated hydrogen addition, HBAT2 provides a complete workflow from raw PDB files to publication-ready analyses.

![An example visualisation of potential cooperativity chain generated by HBAT2 software for Protein Data Bank (PDB) entry 6RSA](https://static.abhishek-tiwari.com/hbat/6rsa-pdb-chain-6.png)

# Implementation

HBAT2 employs a modular architecture with separate components for PDB parsing, geometric analysis, statistical computation, and visualisation. The core analysis engine uses efficient nearest-neighbor searching with configurable distance cutoffs, followed by geometric filtering based on distance and angular criteria.

The software implements the same fundamental geometric approach as the original version [@tiwari2007hbat] but with optimised algorithms and improved handling.

# Impact and Applications

Since its original publication, HBAT has been [cited](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC) in numerous studies of protein structure and molecular recognition [@tiwari2007hbat]. The latest update, however, extends this impact by providing modern tools for structure-based drug design, protein engineering, molecular dynamics analysis, crystallographic studies, and comparative structural analysis.

The software's preset configurations and flexible parameter system make it accessible to both computational experts and experimental structural biologists, broadening its potential user base compared to the original version.

# Availability
HBAT2 is freely available to download from [GitHub](https://github.com/abhishektiwari/hbat) and [PyPI](https://pypi.org/project/hbat) under the MIT license with detailed [user and API documentation](https://hbat.abhishek-tiwari.com/). The software can be installed via PyPI (`pip install hbat`) or via Conda (`conda install -c hbat hbat`), with optional GraphViz integration for advanced visualization features.

# Acknowledgements

The author thanks the original co-developer Sunil K. Panigrahi and acknowledges the structural biology community for feedback that guided the modernization of HBAT.

# References
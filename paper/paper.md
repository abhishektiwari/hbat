---
title: 'HBAT 2: A Python Package to Analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures'
authors:
  - name: Abhishek Tiwari
    orcid: 0000-0003-2222-2395
    affiliation: 1
affiliations:
 - name: Independent Researcher, New York, United States
   index: 1
date: 3 August 2025
bibliography: paper.bib
tags:
  - python
  - structural biology
  - hydrogen bonds
  - molecular interactions
  - protein structures
  - bioinformatics
  - non-covalent interactions
  - PDB analysis
---

# Summary

HBAT 2 is a Python package for automated analysis of hydrogen bonds and other non-covalent interactions in macromolecular structures available in Protein Data Bank (PDB) format. The software identifies and analyzes traditional hydrogen bonds, weak hydrogen bonds, halogen bonds, X-H$\cdots$$\pi$ interactions, $\pi$-$\pi$ stacking, and n$\rightarrow$$\pi$* interactions using geometric criteria. HBAT 2 also detects networks of linked hydrogen bonds and renders them as interactive 2D network visualizations. Originally developed in Perl/Tk and published in 2007 [@tiwari2007hbat], HBAT 2 has been completely rewritten in Python with a modern cross-platform `tkinter`-based graphical user interface (GUI), web-based interface, command-line interface (CLI), and developer-friendly APIs enabling integration with Jupyter notebooks and custom workflows. The software also provides JavaScript-based 3D visualization available both on the web server and as a widget in Jupyter notebooks, enabling researchers to interactively explore detected interactions without requiring separate visualization software.

![The latest update to HBAT 2 uses tkinter to provide a cross-platform graphical user interface (GUI)](./images/hbat-gui.png)

# Statement of Need

Hydrogen bonds [@desiraju_weak_2001] and other non-covalent interactions, including halogen bonds [@cavallo_halogen_2016], are fundamental to the structure, stability, and function of macromolecular systems including proteins, nucleic acids, and their complexes. These interactions determine how biological molecules fold, recognize their targets, and bind ligands. With over 200,000 macromolecular structures archived in the Protein Data Bank [@berman2000protein], there is an increasing need for automated tools to systematically analyze these critical interactions across diverse structural types.

The target audience includes structural biologists studying molecular mechanisms, computational chemists designing new drugs, protein engineers improving enzyme properties, and researchers analyzing molecular dynamics simulations. These researchers need to quickly identify and quantify weak interactions in diverse macromolecular structures to answer biological questions about molecular function and evolution.

The original HBAT [@tiwari2007hbat] was developed in Perl/Tk with a Windows-only GUI, limiting its adoption in modern computational environments where researchers use diverse operating systems (Windows, Linux, macOS) and integrate analyses into Python-based scientific workflows. This update addresses these limitations while expanding the types of interactions analyzed.

# State of the Field

The landscape of hydrogen bond analysis tools is diverse but fragmented, with most solutions designed for specific use cases or as secondary features within larger systems. Classic tools like HBPLUS [@mcdonald1994satisfying] and HBexplore [@lindauer1996hbexplore] pioneered automated H-bond detection but lack modern interfaces and support for diverse interaction types.

Current tools are typically constrained: PLIP [@salentin_plip_2015] and Arpeggio [@jubb_arpeggio_2017] focus on protein-ligand interactions; HBonanza [@durrant_hbonanza_2011], HBCalculator [@wang_hbcalculator_2024], and BRIDGE2 [@siemers_interactive_2021] are optimized for molecular dynamics trajectories; MDAnalysis [@noauthor_431_nodate], GROMACS [@noauthor_gmx_nodate], and AMBER [@noauthor_hbond_2020] provide H-bond analysis as secondary features within larger simulation suites. Visualization tools like VMD [@noauthor_vmd_nodate] and ChimeraX [@noauthor_tool_nodate] excel at interactive exploration but offer limited quantitative analysis. ProteinTools [@ferruz_proteintools_2021] provides web-based network analysis but lacks cross-platform accessibility.

HBAT 2 uniquely provides a **general-purpose solution** with: (1) detection of diverse interaction types including hydrogen bonds, halogen bonds, π-interactions, and n→π* interactions; (2) multiple access modes (desktop GUI, command-line, web server, Python API) supporting diverse research workflows; (3) identification of hydrogen bond networks - both potential cooperativity and anticooperativity chains and water-mediated hydrogen bond networks - in static structures; (4) seamless Jupyter integration; and (5) flexible parameters with domain-specific presets. The modern Python implementation with cross-platform support and web accessibility makes comprehensive interaction analysis available to researchers regardless of computational background or infrastructure.

# Software Design

HBAT 2 employs a modular architecture with distinct components for PDB parsing, geometric analysis, statistical computation, and network visualization. The core analysis engine uses efficient nearest-neighbor searching with configurable distance cutoffs, followed by geometric filtering based on distance and angular criteria. This fundamental approach matches the original HBAT [@tiwari2007hbat] but with optimized algorithms for analyzing large macromolecular structures.

![HBAT 2 analysis workflow showing the complete pipeline from PDB input to multiple output formats. \label{fig:workflow}](./images/overall-workflow-graphviz.pdf){width=100%}

Key design decisions include: (1) Python for accessibility and integration with scientific ecosystems; (2) integration with PDBFixer [@noauthor_openmmpdbfixer_2025; @eastman_openmm_2013] and OpenBabel [@oboyle_pybel_2008] for automated structure preparation, addressing real-world challenges of missing atoms and incomplete side chains; (3) support for multiple interaction types reflecting advances in weak interaction chemistry [@desiraju_weak_2001; @cavallo_halogen_2016; @brandl_c-h-interactions_2001]; (4) preset parameter systems for common experimental conditions (high-resolution X-ray, NMR, membrane proteins, drug design) that balance accessibility for experimentalists with flexibility for computational experts; and (5) multiple output formats (CSV, JSON, network graphics in PNG/SVG/PDF) enabling integration into downstream analysis pipelines.

The software analyzes hydrogen bonds (O-H$\cdots$O, N-H$\cdots$O, N-H$\cdots$N, C-H$\cdots$O), halogen bonds (C-X$\cdots$Y where X=F,Cl,Br,I), X-H$\cdots$$\pi$ interactions, $\pi$-$\pi$ stacking [@mcgaughey_pi-stacking_1998; @vernon_pi-pi_nodate], carbonyl-carbonyl n$\rightarrow$$\pi$* interactions [@rahim_reciprocal_2017; @newberry_n_2017], and n$\rightarrow$$\pi$* interactions [@choudhary_nature_2009]. Network visualization uses either NetworkX/Matplotlib [@hagberg2008networkx; @hunter2007matplotlib] or GraphViz [@graphviz2024], emphasizing cooperativity chains and network topology with customizable layouts and high-resolution export.

![C-H$\cdots$$\pi$ interaction between A:ALA:20 → A:TYR:25 of PDB Structure 6RSA visualized using 3Dmol.js integration in HBAT 2](./images/6rsa_A_ALA_20_to_A_TYR_25_pi.png)

## 3D Visualization with 3Dmol.js Integration

HBAT 2 integrates 3Dmol.js [@rego_3dmol_2014] for interactive 3D molecular visualization in both web and Jupyter notebook environments. The JavaScript-based viewer enables researchers to interactively explore detected interactions with customizable color schemes, allowing rotation, zoom, and inspection of specific interaction geometries. Within Jupyter notebooks, the 3D visualization widget is embedded directly alongside analysis code, bridging automated analysis and manual structural inspection.

# Research Impact Statement

Since its original publication, HBAT has been [cited more than 80 times](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC) in numerous studies spanning multiple research domains [@tiwari2007hbat]. HBAT has been used for:

**Structure-based drug design**: Researchers have used HBAT to characterize protein-ligand interactions for anticancer compounds targeting ABL kinase [@diaz-cervantes_molecular_2020], analyze cardiovascular disease receptors with metabolite interactions [@ravindran_interaction_2015], predict drug resistance based on hydrogen bond patterns [@zhou_prediction_2013], and elucidate anticoagulant interactions [@russo_krauss_thrombinaptamer_2011].

**Protein engineering**: HBAT enabled identification of stabilizing C-H$\cdots$$\pi$ and C-H$\cdots$O interactions in thermostable enzyme variants [@wang_simultaneously_2020] and characterization of hydrogen bonding networks for functional specificity [@fournier_vanadium_2014].

**Molecular dynamics analysis**: The tool has been used to characterize interactions in viral proteins important for understanding RNA encapsulation [@huang_kernel_2009], glycoprotein interactions for cellular quality control [@jayaprakash_atomic_2025], and DNA structural dynamics [@barzegar_stabilization_2025].

**Comparative structural analysis**: HBAT has enabled systematic mutation analysis identifying disease-associated variants in cancer genes [@kavitha_insights_2024; @khan_prediction_2017; @khan_identification_2018; @abdulazeez_rs61742690_2019; @abdulazeez_-silico_2016].

**Recent citations**: Recent applications include biosorption and molecular docking studies utilizing HBAT 2 for hydrogen bond analysis in environmental remediation [@mehmet_karadayı_removal_2026].

![An example visualization of hydrogen-bonds network detected by HBAT 2 software for Protein Data Bank (PDB) entry 6RSA](./images/6rsa-network-hbonds.pdf)

HBAT 2 extends this demonstrated impact by providing modern cross-platform access, expanded interaction coverage, and integration with contemporary Python workflows, supporting growing adoption across diverse research applications.

# Availability

HBAT 2 is freely available to download from [GitHub](https://github.com/abhishektiwari/hbat) and [PyPI](https://pypi.org/project/hbat) under the MIT license with detailed [user and API documentation](https://hbat.abhishek-tiwari.com/). The software can be installed via [PyPI](https://pypi.org/project/hbat) (`pip install hbat`) or via [Conda](https://anaconda.org/hbat/hbat) (`conda install -c hbat hbat`), with optional GraphViz integration for advanced visualization features. A [web-based interface](https://hbat-web.abhishek-tiwari.com) is available for immediate browser-based analysis without installation. Interactive [Jupyter notebooks](https://github.com/abhishektiwari/hbat/tree/main/notebooks) with Google Colab integration are provided for computational workflows.

# AI Usage Disclosure

AI tools were used for porting original HBAT from Perl to Python, as well as for refactoring, test scaffolding, and documentation generation for portions of the codebase. AI tools were also used to assist with editing portions of this manuscript. All AI-assisted code and documentation outputs were reviewed, edited, validated, and tested by the human author, who takes full responsibility for the final software and paper. No figures or data were generated by AI.

# Acknowledgements

The author thanks the original co-developer Sunil K. Panigrahi and acknowledges the structural biology community for feedback that guided the modernization of HBAT.

# References
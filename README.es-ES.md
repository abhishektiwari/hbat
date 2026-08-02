

![HBAT](https://github.com/abhishektiwari/hbat/raw/main/hbat.svg)

# HBAT 2 (Hydrogen Bond Analysis Tool 2) 

Un paquete de Python para automatizar el análisis de enlaces de hidrógeno potenciales y otros tipos de interacciones débiles en estructuras macromoleculares del Banco de Datos de Proteínas (PDB). HBAT 2 admite los formatos de archivo `.pdb` y `.cif` (mmCIF) y utiliza un enfoque geométrico para identificar interacciones moleculares analizando criterios de distancia y ángulo.

**Tipos de interacciones admitidos:**

- **Enlaces de hidrógeno**: Interacciones clásicas `N-H···O`, `O-H···O` y débiles `C-H···O`
- **Enlaces de halógeno**: Interacciones `C-X···A` (`X = Cl, Br, I`)
- **Interacciones π**: Interacciones X-H···π y `C-X···π` con anillos aromáticos (`Phe`, `Tyr`, `Trp`, `His`, etc.)
- **Apilamiento π-π**: Interacciones entre anillos aromáticos (paralelo, en forma de T, desplazado)
- **Interacciones de carbonilo**: Interacciones `n→π*` entre grupos carbonilo
- **Interacciones n-π**: Interacciones de pares de electrones solitarios con sistemas `π` aromáticos
- **Puentes de agua**: Redes de enlaces de hidrógeno mediadas por agua que conectan residuos de proteínas/ligandos
- **Interacciones de ligandos**: Detección completa de todos los tipos de interacciones entre ligandos y residuos de proteínas/ácidos nucleicos
 
> **¡La interfaz web de HBAT 2 ya está disponible!** Pruébala en [hbat-web.abhishek-tiwari.com](https://hbat-web.abhishek-tiwari.com)


![GitHub Release](https://img.shields.io/github/v/release/abhishektiwari/hbat)
![GitHub Actions Test Workflow Status](https://img.shields.io/github/actions/workflow/status/abhishektiwari/hbat/test.yml?label=tests)
![PyPI - Version](https://img.shields.io/pypi/v/hbat)
![Python Wheels](https://img.shields.io/pypi/wheel/hbat)
![Python Versions](https://img.shields.io/pypi/pyversions/hbat?logo=python&logoColor=white)
![GitHub last commit](https://img.shields.io/github/last-commit/abhishektiwari/hbat)
![PyPI - Status](https://img.shields.io/pypi/status/hbat)
![Conda Version](https://img.shields.io/conda/v/hbat/hbat)
![License](https://img.shields.io/github/license/abhishektiwari/hbat)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/abhishektiwari/hbat/total?label=GitHub%20Downloads)
![SourceForge Downloads](https://img.shields.io/sourceforge/dt/hbat?label=SourceForge%20Downloads)
![PyPI Downloads](https://img.shields.io/pepy/dt/hbat?label=PyPI%20Downloads)
[![codecov](https://codecov.io/gh/abhishektiwari/hbat/graph/badge.svg?token=QSKYLB3M1V)](https://codecov.io/gh/abhishektiwari/hbat)
[![Socket](https://socket.dev/api/badge/pypi/package/hbat/2.2.11?artifact_id=py3-none-any-whl)](https://socket.dev/pypi/package/hbat/overview/2.2.11/py3-none-any-whl)
[![CodeFactor](https://www.codefactor.io/repository/github/abhishektiwari/hbat/badge/main)](https://www.codefactor.io/repository/github/abhishektiwari/hbat/overview/main)
[![DOI HBAT](https://img.shields.io/badge/10.3233%2FISI-2007-00337?logo=doi&label=10.3233%2FISI-2007-00337&link=https%3A%2F%2Fdoi.org%2F10.3233%2FISI-2007-00337)](https://doi.org/10.3233/ISI-2007-00337)
[![Google Scholar Citation](https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3233%2FISI-2007-00337&link=https%3A%2F%2Fscholar.google.com%2Fcitations%3Fview_op%3Dview_citation%26hl%3Den%26user%3DMb7eYKYAAAAJ%26citation_for_view%3DMb7eYKYAAAAJ%3Au-x6o8ySG0sC)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:u-x6o8ySG0sC)
[![arXiv](https://img.shields.io/badge/arXiv-2602.17712-b31b1b.svg)](https://doi.org/10.48550/arXiv.2602.17712)
[![Chemrxiv](https://img.shields.io/badge/chemrxiv-15000141--v1-green)](https://doi.org/10.26434/chemrxiv.15000141/v1)
[![Google Scholar Citation](https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.26434%2Fchemrxiv.15000141%2Fv1&link=https%3A%2F%2Fscholar.google.com%2Fcitations%3Fview_op%3Dview_citation%26hl%3Den%26user%3DMb7eYKYAAAAJ%26citation_for_view%3DMb7eYKYAAAAJ%3A3bvyWxjaHKcC)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Mb7eYKYAAAAJ&citation_for_view=Mb7eYKYAAAAJ:3bvyWxjaHKcC)

**HBAT Desktop (Mac, Windows, Linux):**
![HBAT Desktop](https://static.abhishek-tiwari.com/hbat/hbat-window-v3.png)

**HBAT Web:** Pruébala en [hbat-web.abhishek-tiwari.com](https://hbat-web.abhishek-tiwari.com)

![HBAT Web](https://static.abhishek-tiwari.com/hbat/hbat-2-web-v2.png)

**Visualización de puente de agua con HBAT Web usando PyMOL (Entrada PDB 6RSA):**
![Water Bridge in PDB Entry 6RSA](https://static.abhishek-tiwari.com/hbat/6rsa_A_ARG_10_to_A_ASP_38_pymol.png)

**Visualización de interacciones con HBAT Web usando D3MOl (Entrada PDB 6RSA):**

![Pi Interaction in PDB Entry 6RSA](https://static.abhishek-tiwari.com/hbat/6rsa_A_MET_29_to_A_PHE_46_pi.png)

**Visualización de interacciones de ligandos en 2D usando LigPlots de HBAT Web (Entrada PDB 2IZF)** 
![Ligand Interaction in 2IZF](https://static.abhishek-tiwari.com/hbat/2izf_B_BTN_300_ligplot_static.svg)

**Detección y visualización de cadenas de cooperatividad (Entrada PDB 6RSA):**

![Cooperativity chain detection and visualization](https://static.abhishek-tiwari.com/hbat/6rsa_chain_H_bond_chain_10.png)

## Antecedentes
HBAT 2 es una reimplementación moderna en Python de la herramienta original basada en Perl desarrollada por [Abhishek Tiwari](https://www.abhishek-tiwari.com) y Sunil Kumar Panigrahi. La versión HBAT v1 aún se puede descargar desde [SourceForge](https://sourceforge.net/projects/hbat/files/HBAT/), sin embargo, la versión de Perl ya no se mantiene. 


## Características principales de HBAT 2

- Detección y análisis de enlaces de hidrógeno potenciales, enlaces de halógeno, interacciones π, apilamiento π-π, interacciones de carbonilo, interacciones n-π, puentes de agua e interacciones de ligandos
- Corrección automatizada de archivos PDB con integración de OpenBabel y PDBFixer
- Admite interfaces gráficas (tkinter), de línea de comandos y API de programación
- Utilice interfaces gráficas para análisis interactivo, CLI/API para procesamiento por lotes y automatización
- Análisis de interacciones de ligandos con visualización y filtrado específicos por residuo
- Detección y análisis de puentes de agua con visualización de la ruta del puente
- Visualización de redes de enlaces de hidrógeno (cadenas potenciales de cooperatividad/anticooperatividad y redes de enlaces de hidrógeno mediadas por agua) utilizando NetworkX/matplotlib y GraphViz
- Exportación de visualizaciones de redes de enlaces de hidrógeno a formatos PNG, SVG y PDF
- Visualización 3D de interacciones usando 3Dmol.js en cuadernos Jupyter y la interfaz web de HBAT
- Exportación y visualización de interacciones en PyMOL desde la interfaz web de HBAT
- Presets integrados para diferentes tipos de estructuras (alta resolución, RMN, proteínas de membrana, etc.)
- Límites de distancia personalizables, umbrales de ángulo y modos de análisis.
- Múltiples formatos de salida: opciones de exportación a texto, CSV y JSON
- Algoritmos optimizados para un análisis eficiente de estructuras grandes
- Multiplataforma: Funciona en Windows, macOS y Linux.

Consulte la [documentación de HBAT](https://hbat.abhishek-tiwari.com/) para obtener más detalles.

## Instalación

### Opción 1: Instalar desde PyPI (Recomendado)

```bash
pip install hbat
```

Ejecute la interfaz de línea de comandos (CLI) de HBAT usando `hbat` o inicie la GUI de HBAT usando `hbat-gui`.

### Opción 2: Instalar desde el código fuente

```bash
git clone https://github.com/abhishektiwari/hbat.git
cd hbat
pip install -e .
```

Alternativamente,  

```bash
pip install git+https://github.com/abhishektiwari/hbat.git
```

Ejecute la interfaz de línea de comandos (CLI) de HBAT usando `hbat` o inicie la GUI de HBAT usando `hbat-gui`.

### Opción 3: Instalar desde Conda

```
conda install -c hbat hbat
```

### Requisitos

#### Requisitos del sistema
- Python: 3.9 o superior
- tkinter: tkinter se incluye en la biblioteca estándar de Python en la mayoría de los sistemas. Sin embargo, en Mac instale Python y tkinter usando `brew`. 

```
brew install python python3-tk
```

- GraphViz (Opcional): Requerido para la visualización avanzada de cadenas de cooperatividad con renderizado de gráficos de alta calidad. HBAT volverá automáticamente a la visualización de NetworkX/matplotlib si GraphViz no está disponible.

Instalar GraphViz:

En Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install graphviz
```

En macOS (usando Homebrew):
```bash
brew install graphviz
```

En Windows:
- Descargue e instale desde el [sitio web oficial de GraphViz](https://graphviz.org/download/)
- O usando Chocolatey: `choco install graphviz`
- O usando conda: `conda install -c conda-forge graphviz`

> Nota: Después de instalar GraphViz, reinicie su terminal/consola de comandos antes de ejecutar HBAT para asegurarse de que los ejecutables de GraphViz estén disponibles en su PATH.

## Uso

### Interfaz gráfica

Inicie la aplicación GUI:

```bash
hbat-gui
```

La GUI proporciona:
- Explorador de archivos para cargar archivos PDB
- Paneles de configuración de parámetros
- Visualización de resultados con pestañas
- Opciones de exportación y visualización

### Interfaz de línea de comandos

Uso básico:

```bash
hbat input.pdb
hbat input.cif
```

#### Opciones de formato de salida

HBAT admite múltiples formatos de salida:

```bash
# Sin bandera de salida: muestra los resultados en la consola
hbat input.pdb
                    # Muestra los resultados en la consola

# Salidas de archivo único (formato detectado automáticamente desde la extensión)
hbat input.pdb -o results.txt     # Formato de texto (resumen legible por humanos + detalles)
hbat input.pdb -o results.json    # Formato JSON (archivo único con todas las interacciones)

# Salidas de múltiples archivos (archivos separados por tipo de interacción)
hbat input.pdb --csv results      # Crea results_h_bonds.csv, results_x_bonds.csv, etc.
hbat input.pdb --json results     # Crea results_h_bonds.json, results_x_bonds.json, etc.
```

Con parámetros personalizados:

```bash
hbat input.pdb -o results.txt --hb-distance 3.0 --mode inter
```

Modos de inclusión de interacciones:

| Modo | Inter-residuos | Intra-residuos |
|---|---:|---:|
| `inter` | Sí | No |
| `all` | Sí | Sí |

> **Cambio importante:** Los valores de modo anteriores `local` y `complete` ya no son válidos. Reemplace `local` por `inter` y `complete` por `all`.

#### Listar presets disponibles

```bash
hbat --list-presets
```

#### Usar un preset específico

```bash
hbat protein.pdb --preset high_resolution
hbat membrane_protein.pdb --preset membrane_proteins
```

#### Usar preset con anulaciones personalizadas

```bash
hbat protein.pdb --preset drug_design_strict --hb-distance 3.0 --verbose
```

#### Opciones de CLI

```
argumentos posicionales:
  input                 Archivo PDB de entrada

argumentos opcionales:
  -h, --help            mostrar este mensaje de ayuda y salir
  -o OUTPUT, --output OUTPUT
                        Archivo de salida (formato detectado automáticamente desde la extensión: .txt, .json)
  --json JSON           Exportar a múltiples archivos JSON (nombre base para los archivos)
  --csv CSV             Exportar a múltiples archivos CSV (nombre base para los archivos)

Opciones de preset:
  --preset PRESET       Cargar parámetros desde archivo preset (.hbat o .json)
  --list-presets        Listar presets de ejemplo disponibles y salir

Parámetros de análisis:
  Parámetros de enlaces de hidrógeno:
  --hb-distance HB_DISTANCE
                        Límite de distancia H...A del enlace de hidrógeno en Å (predeterminado: 2.5)
  --hb-angle HB_ANGLE   Límite de ángulo D-H...A del enlace de hidrógeno en grados (predeterminado: 120)
  --da-distance DA_DISTANCE
                        Límite de distancia donador-aceptor en Å (predeterminado: 3.5)

  Parámetros de enlaces de halógeno:
  --xb-distance XB_DISTANCE
                        Límite de distancia X...A del enlace de halógeno en Å (predeterminado: 3.9)
  --xb-angle XB_ANGLE   Límite de ángulo C-X...A del enlace de halógeno en grados (predeterminado: 150)

  Parámetros de interacciones π:
  --pi-distance PI_DISTANCE
                        Límite de distancia H...π de la interacción π en Å (predeterminado: 3.5)
  --pi-angle PI_ANGLE   Límite de ángulo D-H...π de la interacción π en grados (predeterminado: 110)

  Parámetros de apilamiento π-π:
  --pi-pi-distance PI_PI_DISTANCE
                        Límite de distancia centroide-a-centroide del apilamiento π-π en Å (predeterminado: 3.8)
  --pi-pi-parallel-angle PI_PI_PARALLEL_ANGLE
                        Ángulo máximo para apilamiento π-π paralelo en grados (predeterminado: 30.0)
  --pi-pi-tshaped-angle-min PI_PI_TSHAPED_ANGLE_MIN
                        Ángulo mínimo para apilamiento π-π en forma de T en grados (predeterminado: 60.0)
  --pi-pi-tshaped-angle-max PI_PI_TSHAPED_ANGLE_MAX
                        Ángulo máximo para apilamiento π-π en forma de T en grados (predeterminado: 90.0)
  --pi-pi-offset PI_PI_OFFSET
                        Desplazamiento lateral máximo para apilamiento π-π paralelo en Å (predeterminado: 2.0)

  Parámetros de interacciones de carbonilo (n→π*):
  --carbonyl-distance CARBONYL_DISTANCE
                        Límite de distancia O···C del carbonilo en Å (predeterminado: 3.2)
  --carbonyl-angle-min CARBONYL_ANGLE_MIN
                        Ángulo mínimo O···C=O para interacciones de carbonilo en grados (predeterminado: 95.0)
  --carbonyl-angle-max CARBONYL_ANGLE_MAX
                        Ángulo máximo O···C=O para interacciones de carbonilo en grados (predeterminado: 125.0)

  Parámetros de interacciones n→π*:
  --n-pi-distance N_PI_DISTANCE
                        Límite de distancia par solitario al centro π en Å (predeterminado: 3.6)
  --n-pi-sulfur-distance N_PI_SULFUR_DISTANCE
                        Límite de distancia específico para azufre en Å (predeterminado: 4.0)
  --n-pi-angle-min N_PI_ANGLE_MIN
                        Ángulo mínimo al plano π en grados (predeterminado: 0.0)
  --n-pi-angle-max N_PI_ANGLE_MAX
                        Ángulo máximo al plano π en grados (predeterminado: 45.0)

  Parámetros generales:
  --covalent-factor COVALENT_FACTOR
                        Factor de detección de enlaces covalentes (predeterminado: 0.85)
  --mode {inter,all}
                        Modo de inclusión de interacciones: inter incluye interacciones
                        solo entre diferentes residuos; all también incluye
                        interacciones intra-residuo

Control de salida:
  --verbose, -v         Salida detallada con progreso exhaustivo
  --quiet, -q           Modo silencioso con salida mínima
  --summary-only        Mostrar solo estadísticas resumidas

Filtros de análisis:
  --no-hydrogen-bonds   Omitir análisis de enlaces de hidrógeno
  --no-halogen-bonds    Omitir análisis de enlaces de halógeno
  --no-pi-interactions  Omitir análisis de interacciones π
  --no-pi-pi-stacking   Omitir análisis de apilamiento π-π
  --no-carbonyl-interactions
                        Omitir análisis de interacciones n→π* de carbonilo
  --no-n-pi-interactions
                        Omitir análisis de interacciones n→π*
```

## Cuadernos de ejemplo

Cuadernos Jupyter interactivos que demuestran el uso de HBAT con visualizaciones 3D usando Py3DMol.

| Cuaderno | Descripción | Colab |
|----------|-------------|-------|
| [01_analyze_6rsa_with_visualization.ipynb](notebooks/01_analyze_6rsa_with_visualization.ipynb) | Análisis de enlaces de hidrógeno de 6RSA (Ribonucleasa A) con visualización en py3Dmol | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/abhishektiwari/hbat/blob/main/notebooks/01_analyze_6rsa_with_visualization.ipynb) |
| [02_halogen_bonds_4x21.ipynb](notebooks/02_halogen_bonds_4x21.ipynb) | Detección y visualización de enlaces de halógeno en la estructura 4X21 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/abhishektiwari/hbat/blob/main/notebooks/02_halogen_bonds_4x21.ipynb) |
| [03_pdbfixer_vs_openbabel_comparison.ipynb](notebooks/03_pdbfixer_vs_openbabel_comparison.ipynb) | Comparación de PDBFixer frente a OpenBabel para la adición de hidrógenos | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/abhishektiwari/hbat/blob/main/notebooks/03_pdbfixer_vs_openbabel_comparison.ipynb) |


Consulte el [directorio de cuadernos](notebooks/) para más detalles.

## Licencia

Este proyecto está licenciado bajo la Licencia MIT - consulte el archivo [LICENSE](LICENSE) para más detalles.

## Citar HBAT y HBAT 2

Si utiliza HBAT 2 en su investigación, por favor cite:

[![arXiv](https://img.shields.io/badge/arXiv-2602.17712-b31b1b.svg)](https://doi.org/10.48550/arXiv.2602.17712)
```
@article{tiwari_2026_hbat_arxiv,
  author       = {Tiwari, Abhishek},
  title        = {HBAT 2: A Python Package to analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
  year         = 2026,
  publisher    = {arXiv},
  doi          = {10.48550/arXiv.2602.17712},
  url          = {https://arxiv.org/abs/2602.17712}, 
}
```

o 

[![Chemrxiv](https://img.shields.io/badge/chemrxiv-15000141--v1-green)](https://doi.org/10.26434/chemrxiv.15000141/v1)
```
@article{tiwari_2026_hbat_chemrxiv,
  author = {Abhishek Tiwari },
  title = {HBAT 2: A Python Package to Analyse Hydrogen Bonds and Other Non-covalent Interactions in Macromolecular Structures},
  publisher = {ChemRxiv},
  year = {2026},
  doi = {10.26434/chemrxiv.15000141/v1},
  URL = {https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15000141/v1},
  eprint = {https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.15000141/v1},
}
```

Si utiliza HBAT 1.0 o 1.1 en su investigación, por favor cite:

[![DOI HBAT](https://img.shields.io/badge/10.3233%2FISI-2007-00337?logo=doi&label=10.3233%2FISI-2007-00337&link=https%3A%2F%2Fdoi.org%2F10.3233%2FISI-2007-00337)](https://doi.org/10.3233/ISI-2007-00337)

```
@article{tiwari2007hbat,
author = {Tiwari, Abhishek and Panigrahi, Sunil Kumar},
doi = {10.3233/ISI-2007-00337},
journal = {In Silico Biology},
month = dec,
number = {6},
title = {{HBAT: A Complete Package for Analysing Strong and Weak Hydrogen Bonds in Macromolecular Crystal Structures}},
volume = {7},
year = {2007}
}
```

## Contribuciones 

Consulte nuestra [guía de contribuciones](CONTRIBUTING.md) y la [guía de desarrollo](https://hbat.abhishek-tiwari.com/development). A grandes rasgos,

1. Bifurque (fork) el repositorio
2. Cree una rama de características
3. Realice sus cambios
4. Agregue pruebas si es aplicable
5. Envíe una solicitud de extracción (pull request)

"""
Jupyter notebook utilities for molecular interaction visualization.

This module provides convenient functions for displaying 3D molecular
visualizations in Jupyter notebooks using py3Dmol.
"""

from typing import TYPE_CHECKING

from .pymol3d import (
    generate_carbonyl_interaction_viewer_js,
    generate_halogen_bond_viewer_js,
    generate_hydrogen_bond_viewer_js,
    generate_n_pi_interaction_viewer_js,
    generate_pi_interaction_viewer_js,
    generate_pi_pi_stacking_viewer_js,
)

if TYPE_CHECKING:
    from ..core.interactions import (
        CarbonylInteraction,
        HalogenBond,
        HydrogenBond,
        NPiInteraction,
        PiInteraction,
        PiPiStacking,
    )


def display_hydrogen_bond(
    hb: "HydrogenBond",
    pdb_content: str,
    viewer_id: str = "hb_viewer",
    width: int = 800,
    height: int = 600,
):
    """Display hydrogen bond 3D visualization in Jupyter notebook.

    :param hb: Hydrogen bond interaction
    :type hb: HydrogenBond
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div (default: "hb_viewer")
    :type viewer_id: str
    :param width: Viewer width in pixels (default: 800)
    :type width: int
    :param height: Viewer height in pixels (default: 600)
    :type height: int

    Example::

        >>> from hbat.visualization.notebook import display_hydrogen_bond
        >>> # After running analysis
        >>> for hb in analyzer.hydrogen_bonds:
        ...     display_hydrogen_bond(hb, pdb_content)
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    js_code = generate_hydrogen_bond_viewer_js(hb, pdb_content, viewer_id)

    # Display viewer container
    display(
        HTML(
            f'<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>'
        )
    )

    # Load 3Dmol library if not already loaded
    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )

    # Display viewer JavaScript
    display(HTML(f"<script>{js_code}</script>"))


def display_halogen_bond(
    xb: "HalogenBond",
    pdb_content: str,
    viewer_id: str = "xb_viewer",
    width: int = 800,
    height: int = 600,
):
    """Display halogen bond 3D visualization in Jupyter notebook.

    :param xb: Halogen bond interaction
    :type xb: HalogenBond
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div (default: "xb_viewer")
    :type viewer_id: str
    :param width: Viewer width in pixels (default: 800)
    :type width: int
    :param height: Viewer height in pixels (default: 600)
    :type height: int

    Example::

        >>> from hbat.visualization.notebook import display_halogen_bond
        >>> for xb in analyzer.halogen_bonds:
        ...     display_halogen_bond(xb, pdb_content)
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    js_code = generate_halogen_bond_viewer_js(xb, pdb_content, viewer_id)

    display(
        HTML(
            f'<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>'
        )
    )
    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )
    display(HTML(f"<script>{js_code}</script>"))


def display_pi_interaction(
    pi: "PiInteraction",
    pdb_content: str,
    viewer_id: str = "pi_viewer",
    width: int = 800,
    height: int = 600,
):
    """Display π interaction 3D visualization in Jupyter notebook.

    :param pi: π interaction
    :type pi: PiInteraction
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div (default: "pi_viewer")
    :type viewer_id: str
    :param width: Viewer width in pixels (default: 800)
    :type width: int
    :param height: Viewer height in pixels (default: 600)
    :type height: int

    Example::

        >>> from hbat.visualization.notebook import display_pi_interaction
        >>> for pi in analyzer.pi_interactions:
        ...     display_pi_interaction(pi, pdb_content)
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    js_code = generate_pi_interaction_viewer_js(pi, pdb_content, viewer_id)

    display(
        HTML(
            f'<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>'
        )
    )
    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )
    display(HTML(f"<script>{js_code}</script>"))


def display_pi_pi_stacking(
    pi_pi: "PiPiStacking",
    pdb_content: str,
    viewer_id: str = "pipi_viewer",
    width: int = 800,
    height: int = 600,
):
    """Display π-π stacking 3D visualization in Jupyter notebook.

    :param pi_pi: π-π stacking interaction
    :type pi_pi: PiPiStacking
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div (default: "pipi_viewer")
    :type viewer_id: str
    :param width: Viewer width in pixels (default: 800)
    :type width: int
    :param height: Viewer height in pixels (default: 600)
    :type height: int

    Example::

        >>> from hbat.visualization.notebook import display_pi_pi_stacking
        >>> for pi_pi in analyzer.pi_pi_interactions:
        ...     display_pi_pi_stacking(pi_pi, pdb_content)
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    js_code = generate_pi_pi_stacking_viewer_js(pi_pi, pdb_content, viewer_id)

    display(
        HTML(
            f'<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>'
        )
    )
    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )
    display(HTML(f"<script>{js_code}</script>"))


def display_carbonyl_interaction(
    carbonyl: "CarbonylInteraction",
    pdb_content: str,
    viewer_id: str = "carbonyl_viewer",
    width: int = 800,
    height: int = 600,
):
    """Display carbonyl n→π* interaction 3D visualization in Jupyter notebook.

    :param carbonyl: Carbonyl interaction
    :type carbonyl: CarbonylInteraction
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div (default: "carbonyl_viewer")
    :type viewer_id: str
    :param width: Viewer width in pixels (default: 800)
    :type width: int
    :param height: Viewer height in pixels (default: 600)
    :type height: int

    Example::

        >>> from hbat.visualization.notebook import display_carbonyl_interaction
        >>> for carbonyl in analyzer.carbonyl_interactions:
        ...     display_carbonyl_interaction(carbonyl, pdb_content)
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    js_code = generate_carbonyl_interaction_viewer_js(carbonyl, pdb_content, viewer_id)

    display(
        HTML(
            f'<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>'
        )
    )
    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )
    display(HTML(f"<script>{js_code}</script>"))


def display_n_pi_interaction(
    n_pi: "NPiInteraction",
    pdb_content: str,
    viewer_id: str = "npi_viewer",
    width: int = 800,
    height: int = 600,
):
    """Display n→π* interaction 3D visualization in Jupyter notebook.

    :param n_pi: n→π* interaction
    :type n_pi: NPiInteraction
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div (default: "npi_viewer")
    :type viewer_id: str
    :param width: Viewer width in pixels (default: 800)
    :type width: int
    :param height: Viewer height in pixels (default: 600)
    :type height: int

    Example::

        >>> from hbat.visualization.notebook import display_n_pi_interaction
        >>> for n_pi in analyzer.n_pi_interactions:
        ...     display_n_pi_interaction(n_pi, pdb_content)
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    js_code = generate_n_pi_interaction_viewer_js(n_pi, pdb_content, viewer_id)

    display(
        HTML(
            f'<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>'
        )
    )
    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )
    display(HTML(f"<script>{js_code}</script>"))


def load_3dmol_library():
    """Load the 3Dmol.js library in a Jupyter notebook.

    Call this once at the beginning of your notebook to ensure the 3Dmol library is loaded.

    Example::

        >>> from hbat.visualization.notebook import load_3dmol_library
        >>> load_3dmol_library()
    """
    try:
        from IPython.display import HTML, display
    except ImportError:
        raise ImportError(
            "IPython is required for notebook visualization. Install with: pip install ipython"
        )

    display(
        HTML('<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    )

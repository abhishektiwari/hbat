"""
3D molecular visualization utilities using py3Dmol.

This module provides standalone functions for visualizing protein structures
and molecular interactions using py3Dmol, extracted from notebook examples
for reuse in Jupyter notebooks and scripts.
"""

from typing import List, Optional

from hbat.constants import WATER_MOLECULES

try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError:
    PY3DMOL_AVAILABLE = False


def visualize_structure_with_interactions(
    pdb_file: str,
    analyzer,
    width: int = 800,
    height: int = 600,
    max_hbonds: int = 20,
    max_pi_interactions: int = 10,
    show_all_halogen_bonds: bool = True,
    show_water: bool = True
) -> Optional['py3Dmol.view']:
    """Visualize protein structure with molecular interactions.

    Creates a 3D visualization showing the protein structure with hydrogen bonds,
    halogen bonds, and π interactions highlighted as colored cylinders.

    :param pdb_file: Path to PDB file
    :type pdb_file: str
    :param analyzer: HBAT analyzer object with analysis results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param width: Viewer width in pixels
    :type width: int
    :param height: Viewer height in pixels
    :type height: int
    :param max_hbonds: Maximum number of hydrogen bonds to display
    :type max_hbonds: int
    :param max_pi_interactions: Maximum number of π interactions to display
    :type max_pi_interactions: int
    :param show_all_halogen_bonds: Whether to show all halogen bonds
    :type show_all_halogen_bonds: bool
    :param show_water: Whether to show water molecules as red spheres. Uses WATER_MOLECULES constant (HOH, WAT, DOD, TIP3, TIP4, TIP5, W)
    :type show_water: bool
    :returns: py3Dmol viewer object if available, None otherwise
    :rtype: Optional[py3Dmol.view]

    Example::

        >>> from hbat.visualization import visualize_structure_with_interactions
        >>> viewer = visualize_structure_with_interactions('protein.pdb', analyzer)
        >>> viewer.show()  # In Jupyter notebook
    """
    if not PY3DMOL_AVAILABLE:
        raise ImportError(
            "py3Dmol is required. Install with: pip install py3Dmol"
        )

    # Read PDB file
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()

    # Create viewer
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(pdb_data, 'pdb')

    # Style the protein
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})

    # Show water molecules if requested
    if show_water:
        viewer.addStyle(
            {'resn': WATER_MOLECULES},
            {'sphere': {'color': 'red', 'radius': 0.3}}
        )

    # Add hydrogen bonds as dashed yellow lines
    for hb in analyzer.hydrogen_bonds[:max_hbonds]:
        viewer.addCylinder({
            'start': {
                'x': hb.donor.coords.x,
                'y': hb.donor.coords.y,
                'z': hb.donor.coords.z
            },
            'end': {
                'x': hb.acceptor.coords.x,
                'y': hb.acceptor.coords.y,
                'z': hb.acceptor.coords.z
            },
            'radius': 0.1,
            'color': 'yellow',
            'dashed': True
        })

    # Add halogen bonds as orange lines
    halogen_bonds = analyzer.halogen_bonds if show_all_halogen_bonds else analyzer.halogen_bonds[:max_hbonds]
    for xb in halogen_bonds:
        viewer.addCylinder({
            'start': {
                'x': xb.halogen.coords.x,
                'y': xb.halogen.coords.y,
                'z': xb.halogen.coords.z
            },
            'end': {
                'x': xb.acceptor.coords.x,
                'y': xb.acceptor.coords.y,
                'z': xb.acceptor.coords.z
            },
            'radius': 0.15,
            'color': 'orange',
            'dashed': True
        })

    # Add π interactions as cyan lines
    for pi in analyzer.pi_interactions[:max_pi_interactions]:
        if hasattr(pi, 'donor') and hasattr(pi, 'aromatic_center'):
            viewer.addCylinder({
                'start': {
                    'x': pi.donor.coords.x,
                    'y': pi.donor.coords.y,
                    'z': pi.donor.coords.z
                },
                'end': {
                    'x': pi.aromatic_center[0],
                    'y': pi.aromatic_center[1],
                    'z': pi.aromatic_center[2]
                },
                'radius': 0.12,
                'color': 'cyan',
                'dashed': True
            })

    viewer.zoomTo()
    return viewer


def visualize_residue_interactions(
    pdb_file: str,
    analyzer,
    chain: str,
    residue_number: int,
    residue_name: str,
    width: int = 800,
    height: int = 600,
    show_labels: bool = True,
    show_water: bool = True
) -> Optional['py3Dmol.view']:
    """Visualize interactions involving a specific residue.

    Creates a focused 3D visualization showing hydrogen bonds involving
    a specific residue, with the target residue and interacting partners
    highlighted.

    :param pdb_file: Path to PDB file
    :type pdb_file: str
    :param analyzer: HBAT analyzer object with analysis results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param chain: Chain identifier (e.g., 'A', 'B')
    :type chain: str
    :param residue_number: Residue sequence number
    :type residue_number: int
    :param residue_name: Residue name (e.g., 'GLY', 'ALA')
    :type residue_name: str
    :param width: Viewer width in pixels
    :type width: int
    :param height: Viewer height in pixels
    :type height: int
    :param show_labels: Whether to show residue labels
    :type show_labels: bool
    :param show_water: Whether to show water molecules as red spheres. Uses WATER_MOLECULES constant (HOH, WAT, DOD, TIP3, TIP4, TIP5, W)
    :type show_water: bool
    :returns: py3Dmol viewer object if available, None otherwise
    :rtype: Optional[py3Dmol.view]

    Example::

        >>> from hbat.visualization import visualize_residue_interactions
        >>> viewer = visualize_residue_interactions('protein.pdb', analyzer, 'A', 10, 'GLY')
        >>> viewer.show()  # In Jupyter notebook
    """
    if not PY3DMOL_AVAILABLE:
        raise ImportError(
            "py3Dmol is required. Install with: pip install py3Dmol"
        )

    # Construct target residue identifier (format: chain_id + res_seq + res_name)
    target_residue_id = f"{chain}{residue_number}{residue_name}"

    # Filter hydrogen bonds involving the target residue using get_donor_residue/get_acceptor_residue
    residue_hbonds = [
        hb for hb in analyzer.hydrogen_bonds
        if hb.get_donor_residue() == target_residue_id or hb.get_acceptor_residue() == target_residue_id
    ]

    # Read PDB file
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()

    # Create viewer
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(pdb_data, 'pdb')

    # Show entire protein as semi-transparent cartoon
    viewer.setStyle({'cartoon': {'color': 'lightgray', 'opacity': 0.3}})

    # Show water molecules if requested
    if show_water:
        viewer.addStyle(
            {'resn': WATER_MOLECULES},
            {'sphere': {'color': 'red', 'radius': 0.3}}
        )

    # Highlight target residue in green
    viewer.addStyle(
        {'chain': chain, 'resi': residue_number},
        {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}}
    )

    # Add hydrogen bonds and highlight interacting residues
    for hb in residue_hbonds:
        # Highlight interacting residues in cyan
        viewer.addStyle(
            {'chain': hb.donor.chain_id, 'resi': hb.donor.res_seq},
            {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15}}
        )
        viewer.addStyle(
            {'chain': hb.acceptor.chain_id, 'resi': hb.acceptor.res_seq},
            {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15}}
        )

        # Draw hydrogen bond
        viewer.addCylinder({
            'start': {
                'x': hb.donor.coords.x,
                'y': hb.donor.coords.y,
                'z': hb.donor.coords.z
            },
            'end': {
                'x': hb.acceptor.coords.x,
                'y': hb.acceptor.coords.y,
                'z': hb.acceptor.coords.z
            },
            'radius': 0.15,
            'color': 'yellow',
            'dashed': True
        })

        # Add labels if requested
        if show_labels:
            viewer.addLabel(
                f"{hb.donor.res_name}{hb.donor.res_seq}",
                {
                    'position': {
                        'x': hb.donor.coords.x,
                        'y': hb.donor.coords.y,
                        'z': hb.donor.coords.z
                    },
                    'backgroundColor': 'white',
                    'fontColor': 'black',
                    'fontSize': 10
                }
            )

    # Zoom to target residue
    viewer.zoomTo({'chain': chain, 'resi': residue_number})
    return viewer


def visualize_residue_halogen_bonds(
    pdb_file: str,
    analyzer,
    chain: str,
    residue_number: int,
    residue_name: str,
    width: int = 800,
    height: int = 600,
    show_labels: bool = True,
    show_water: bool = True
) -> Optional['py3Dmol.view']:
    """Visualize halogen bonds involving a specific residue.

    Creates a focused 3D visualization showing halogen bonds involving
    a specific residue (typically a ligand with halogen atoms), with the
    target residue and interacting partners highlighted.

    :param pdb_file: Path to PDB file
    :type pdb_file: str
    :param analyzer: HBAT analyzer object with analysis results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param chain: Chain identifier (e.g., 'A', 'B')
    :type chain: str
    :param residue_number: Residue sequence number
    :type residue_number: int
    :param residue_name: Residue name (e.g., '3WH', 'CLU')
    :type residue_name: str
    :param width: Viewer width in pixels
    :type width: int
    :param height: Viewer height in pixels
    :type height: int
    :param show_labels: Whether to show residue labels
    :type show_labels: bool
    :param show_water: Whether to show water molecules as red spheres. Uses WATER_MOLECULES constant (HOH, WAT, DOD, TIP3, TIP4, TIP5, W)
    :type show_water: bool
    :returns: py3Dmol viewer object if available, None otherwise
    :rtype: Optional[py3Dmol.view]

    Example::

        >>> from hbat.visualization import visualize_residue_halogen_bonds
        >>> # Visualize halogen bonds from ligand residue 501 in chain A
        >>> viewer = visualize_residue_halogen_bonds('protein.pdb', analyzer, 'A', 501, '3WH')
        >>> viewer.show()  # In Jupyter notebook
    """
    if not PY3DMOL_AVAILABLE:
        raise ImportError(
            "py3Dmol is required. Install with: pip install py3Dmol"
        )

    # Construct target residue identifier (format: chain_id + res_seq + res_name)
    target_residue_id = f"{chain}{residue_number}{residue_name}"

    # Filter halogen bonds involving the target residue using get_donor_residue/get_acceptor_residue
    residue_xbonds = [
        xb for xb in analyzer.halogen_bonds
        if xb.get_donor_residue() == target_residue_id or xb.get_acceptor_residue() == target_residue_id
    ]

    # Read PDB file
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()

    # Create viewer
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(pdb_data, 'pdb')

    # Show entire protein as semi-transparent cartoon
    viewer.setStyle({'cartoon': {'color': 'lightgray', 'opacity': 0.3}})

    # Show water molecules if requested
    if show_water:
        viewer.addStyle(
            {'resn': WATER_MOLECULES},
            {'sphere': {'color': 'red', 'radius': 0.3}}
        )

    # Highlight target residue in purple (halogen donor)
    viewer.addStyle(
        {'chain': chain, 'resi': residue_number},
        {'stick': {'colorscheme': 'purpleCarbon', 'radius': 0.25}}
    )

    # Add halogen bonds and highlight interacting residues
    for xb in residue_xbonds:
        # Highlight halogen donor residue in purple
        viewer.addStyle(
            {'chain': xb.halogen.chain_id, 'resi': xb.halogen.res_seq},
            {'stick': {'colorscheme': 'purpleCarbon', 'radius': 0.25}}
        )
        # Highlight acceptor residue in cyan
        viewer.addStyle(
            {'chain': xb.acceptor.chain_id, 'resi': xb.acceptor.res_seq},
            {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.2}}
        )

        # Draw halogen bond
        viewer.addCylinder({
            'start': {
                'x': xb.halogen.coords.x,
                'y': xb.halogen.coords.y,
                'z': xb.halogen.coords.z
            },
            'end': {
                'x': xb.acceptor.coords.x,
                'y': xb.acceptor.coords.y,
                'z': xb.acceptor.coords.z
            },
            'radius': 0.2,
            'color': 'orange',
            'dashed': True
        })

        # Add labels if requested
        if show_labels:
            viewer.addLabel(
                f"{xb.halogen.res_name}{xb.halogen.res_seq}\n{xb.halogen.name}",
                {
                    'position': {
                        'x': xb.halogen.coords.x,
                        'y': xb.halogen.coords.y,
                        'z': xb.halogen.coords.z
                    },
                    'backgroundColor': 'purple',
                    'fontColor': 'white',
                    'fontSize': 10
                }
            )
            viewer.addLabel(
                f"{xb.acceptor.res_name}{xb.acceptor.res_seq}",
                {
                    'position': {
                        'x': xb.acceptor.coords.x,
                        'y': xb.acceptor.coords.y,
                        'z': xb.acceptor.coords.z
                    },
                    'backgroundColor': 'cyan',
                    'fontColor': 'black',
                    'fontSize': 10
                }
            )

    # Zoom to target residue
    viewer.zoomTo({'chain': chain, 'resi': residue_number})
    return viewer


def visualize_pi_pi_stacking(
    pdb_file: str,
    analyzer,
    width: int = 800,
    height: int = 600,
    max_interactions: int = 5,
    show_water: bool = True
) -> Optional['py3Dmol.view']:
    """Visualize π-π stacking interactions.

    Creates a 3D visualization showing π-π stacking interactions between
    aromatic residues, with aromatic rings highlighted and connections shown.

    :param pdb_file: Path to PDB file
    :type pdb_file: str
    :param analyzer: HBAT analyzer object with analysis results
    :type analyzer: NPMolecularInteractionAnalyzer
    :param width: Viewer width in pixels
    :type width: int
    :param height: Viewer height in pixels
    :type height: int
    :param max_interactions: Maximum number of π-π interactions to display
    :type max_interactions: int
    :param show_water: Whether to show water molecules as red spheres. Uses WATER_MOLECULES constant (HOH, WAT, DOD, TIP3, TIP4, TIP5, W)
    :type show_water: bool
    :returns: py3Dmol viewer object if available, None otherwise
    :rtype: Optional[py3Dmol.view]

    Example::

        >>> from hbat.visualization import visualize_pi_pi_stacking
        >>> viewer = visualize_pi_pi_stacking('protein.pdb', analyzer)
        >>> if viewer:
        ...     viewer.show()  # In Jupyter notebook
    """
    if not PY3DMOL_AVAILABLE:
        raise ImportError(
            "py3Dmol is required. Install with: pip install py3Dmol"
        )

    # Check if π-π stacking interactions exist
    if not hasattr(analyzer, 'pi_pi_stacking') or len(analyzer.pi_pi_stacking) == 0:
        return None

    # Read PDB file
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()

    # Create viewer
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(pdb_data, 'pdb')

    # Show protein as semi-transparent cartoon
    viewer.setStyle({'cartoon': {'color': 'lightgray', 'opacity': 0.2}})

    # Show water molecules if requested
    if show_water:
        viewer.addStyle(
            {'resn': WATER_MOLECULES},
            {'sphere': {'color': 'red', 'radius': 0.3}}
        )

    # Visualize π-π interactions
    for pi_pi in analyzer.pi_pi_stacking[:max_interactions]:
        # Highlight aromatic residues in magenta
        if hasattr(pi_pi, 'residue1_number'):
            viewer.addStyle(
                {'resi': pi_pi.residue1_number},
                {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.2}}
            )
        if hasattr(pi_pi, 'residue2_number'):
            viewer.addStyle(
                {'resi': pi_pi.residue2_number},
                {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.2}}
            )

        # Draw connection between aromatic centers
        if hasattr(pi_pi, 'centroid1') and hasattr(pi_pi, 'centroid2'):
            viewer.addCylinder({
                'start': {
                    'x': pi_pi.centroid1[0],
                    'y': pi_pi.centroid1[1],
                    'z': pi_pi.centroid1[2]
                },
                'end': {
                    'x': pi_pi.centroid2[0],
                    'y': pi_pi.centroid2[1],
                    'z': pi_pi.centroid2[2]
                },
                'radius': 0.2,
                'color': 'magenta',
                'dashed': False
            })

    viewer.zoomTo()
    return viewer

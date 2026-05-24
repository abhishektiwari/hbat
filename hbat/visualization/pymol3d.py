"""
py3Dmol visualization utilities for molecular interactions.

This module provides functions to generate py3Dmol viewer JavaScript code
for visualizing different types of molecular interactions. These functions
can be used in web applications, Jupyter notebooks, and other contexts.
"""

from typing import TYPE_CHECKING

from ..constants.pdb_constants import WATER_MOLECULES

if TYPE_CHECKING:
    from ..core.interactions import (
        CarbonylInteraction,
        HalogenBond,
        HydrogenBond,
        NPiInteraction,
        PiInteraction,
        PiPiStacking,
    )


def _escape_pdb_content(pdb_content: str) -> str:
    """Escape PDB content for JavaScript.

    :param pdb_content: Raw PDB file content
    :type pdb_content: str
    :returns: Escaped PDB content safe for JavaScript template literals
    :rtype: str
    """
    return pdb_content.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")


def _create_viewer_init_wrapper(viewer_id: str, pdb_escaped: str, viewer_init_code: str) -> str:
    """Create standard viewer initialization wrapper JavaScript.

    Common pattern for all interaction viewers: initialize 3Dmol, add PDB model,
    execute viewer-specific code, render, and resize.

    :param viewer_id: Unique ID for the viewer instance
    :param pdb_escaped: Escaped PDB content
    :param viewer_init_code: The viewer-specific JavaScript code to execute inside try block
    :returns: Complete wrapped JavaScript with init, render, and resize logic
    :rtype: str
    """
    return f"""
    (function() {{
        function init3Dmol() {{
            if (typeof $3Dmol === 'undefined') {{
                setTimeout(init3Dmol, 100);
                return;
            }}

            let element = document.getElementById("{viewer_id}");
            if (!element) {{
                setTimeout(init3Dmol, 100);
                return;
            }}

            try {{
                let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});
                window.{viewer_id}_instance = viewer;
                let pdbData = `{pdb_escaped}`;
                viewer.addModel(pdbData, "pdb", {{keepH: true}});

                // Viewer-specific code
                {viewer_init_code}

                viewer.zoomTo();
                viewer.render();

                // Force multiple resizes to ensure canvas is sized correctly
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 100);
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 300);
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 500);
            }} catch (error) {{
                console.error("Error creating viewer:", error);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """


def _create_cartoon_style() -> str:
    """Generate JavaScript for protein/RNA cartoon style (gray, low opacity).

    :returns: JavaScript style code
    :rtype: str
    """
    return "viewer.setStyle({}, {cartoon: {color: 'lightgray', opacity: 0.3}});"


def _create_residue_stick_style(chain: str, resi: int, color: str = "cyanCarbon") -> str:
    """Generate JavaScript for residue stick style.

    :param chain: Chain identifier
    :param resi: Residue sequence number
    :param color: Color scheme (e.g., 'cyanCarbon', 'orangeCarbon')
    :returns: JavaScript style code
    :rtype: str
    """
    return f"viewer.addStyle({{chain: '{chain}', resi: {resi}}}, {{stick: {{colorscheme: '{color}'}}}}); "


def _create_dashed_line(start_x: float, start_y: float, start_z: float,
                       end_x: float, end_y: float, end_z: float) -> str:
    """Generate JavaScript for dashed cylinder line between atoms.

    Shortens line by 0.4 Å on each end to avoid overlapping with atom spheres.

    :param start_x, start_y, start_z: Starting atom coordinates
    :param end_x, end_y, end_z: Ending atom coordinates
    :returns: JavaScript code to draw dashed line
    :rtype: str
    """
    return f"""
                const dx = {end_x} - {start_x}, dy = {end_y} - {start_y}, dz = {end_z} - {start_z};
                const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                const offset = 0.4;
                const ratio_start = offset / dist;
                const ratio_end = (dist - offset) / dist;

                viewer.addCylinder({{
                    start: {{x: {start_x} + dx*ratio_start, y: {start_y} + dy*ratio_start, z: {start_z} + dz*ratio_start}},
                    end: {{x: {start_x} + dx*ratio_end, y: {start_y} + dy*ratio_end, z: {start_z} + dz*ratio_end}},
                    radius: 0.15,
                    color: 'yellow',
                    dashed: true
                }});
    """


def _create_label(text: str, x: float, y: float, z: float,
                  bg_color: str = "cyan", font_color: str = "black", font_size: int = 14) -> str:
    """Generate JavaScript for residue label.

    :param text: Label text
    :param x, y, z: Label position coordinates
    :param bg_color: Background color
    :param font_color: Font color
    :param font_size: Font size in pixels
    :returns: JavaScript code to add label
    :rtype: str
    """
    return f"""
                viewer.addLabel('{text}',
                               {{position: {{x: {x}, y: {y}, z: {z}}},
                                backgroundColor: '{bg_color}', fontColor: '{font_color}', fontSize: {font_size}}});
    """


def _create_pi_center_sphere(x: float, y: float, z: float, color: str = "green") -> str:
    """Generate JavaScript for aromatic ring center sphere.

    Used for π-interactions and π-π stacking to visualize ring centers.

    :param x, y, z: Sphere center coordinates
    :param color: Sphere color
    :returns: JavaScript code to add sphere
    :rtype: str
    """
    return f"""
                viewer.addSphere({{
                    center: {{x: {x}, y: {y}, z: {z}}},
                    radius: 0.3,
                    color: '{color}',
                    alpha: 0.7
                }});
    """


def _create_dashed_line_to_center(start_x: float, start_y: float, start_z: float,
                                  center_x: float, center_y: float, center_z: float) -> str:
    """Generate JavaScript for dashed cylinder from atom to ring center.

    Special handling for π-interactions: shortens line to avoid overlapping with spheres.

    :param start_x, start_y, start_z: Starting atom coordinates
    :param center_x, center_y, center_z: Ring center coordinates
    :returns: JavaScript code to draw dashed line to center
    :rtype: str
    """
    return f"""
                const pix = {start_x}, piy = {start_y}, piz = {start_z};
                const pcx = {center_x}, pcy = {center_y}, pcz = {center_z};
                const dpix = pcx - pix, dpiy = pcy - piy, dpiz = pcz - piz;
                const dist_pi = Math.sqrt(dpix*dpix + dpiy*dpiy + dpiz*dpiz);
                const offset_pi = 0.4;
                const ratio_start_pi = offset_pi / dist_pi;
                const ratio_end_pi = (dist_pi - offset_pi) / dist_pi;

                viewer.addCylinder({{
                    start: {{x: pix + dpix*ratio_start_pi, y: piy + dpiy*ratio_start_pi, z: piz + dpiz*ratio_start_pi}},
                    end: {{x: pix + dpix*ratio_end_pi, y: piy + dpiy*ratio_end_pi, z: piz + dpiz*ratio_end_pi}},
                    radius: 0.1,
                    color: 'yellow',
                    dashed: true
                }});
    """


def _create_distance_label(text: str, x1: float, y1: float, z1: float,
                          x2: float, y2: float, z2: float) -> str:
    """Generate JavaScript for distance label at midpoint between two positions.

    Used for π-interactions and π-π stacking to show distance at midpoint.

    :param text: Label text (e.g., "3.50 Å")
    :param x1, y1, z1: First position coordinates
    :param x2, y2, z2: Second position coordinates
    :returns: JavaScript code to add distance label at midpoint
    :rtype: str
    """
    return f"""
                viewer.addLabel('{text}',
                              {{position: {{x: ({x1} + {x2}) / 2,
                                          y: ({y1} + {y2}) / 2,
                                          z: ({z1} + {z2}) / 2}},
                               backgroundColor: 'black', fontColor: 'white', fontSize: 12}});
    """


def generate_png_export_js(viewer_id: str, filename: str) -> str:
    """Generate JavaScript code to export viewer as PNG.

    :param viewer_id: Unique ID of the viewer instance
    :type viewer_id: str
    :param filename: Filename for the downloaded PNG
    :type filename: str
    :returns: JavaScript code to export viewer as PNG
    :rtype: str

    Example::

        >>> js_code = generate_png_export_js("viewer1", "interaction.png")
        >>> # In web app with nicegui:
        >>> ui.run_javascript(js_code)
    """
    return f"""
        (function() {{
            const viewer = window.{viewer_id}_instance;
            if (viewer) {{
                viewer.render();
                const imgData = viewer.pngURI();
                const link = document.createElement('a');
                link.href = imgData;
                link.download = '{filename}';
                link.click();
            }} else {{
                console.error('Viewer not found');
            }}
        }})();
    """


def generate_hydrogen_bond_viewer_js(
    hb: "HydrogenBond", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for hydrogen bond 3D visualization.

    :param hb: Hydrogen bond interaction
    :type hb: HydrogenBond
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str

    Example::

        >>> js_code = generate_hydrogen_bond_viewer_js(hb, pdb_content, "viewer1")
        >>> # In Jupyter notebook:
        >>> from IPython.display import HTML, display
        >>> display(HTML(f'<div id="viewer1" style="width:800px;height:600px;"></div>'))
        >>> display(HTML(f'<script>{js_code}</script>'))
    """
    donor_chain = hb.donor.chain_id
    donor_resi = hb.donor.res_seq
    acceptor_chain = hb.acceptor.chain_id
    acceptor_resi = hb.acceptor.res_seq

    pdb_escaped = _escape_pdb_content(pdb_content)

    # Build viewer-specific code using helpers
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(donor_chain, donor_resi, "cyanCarbon") +
        _create_residue_stick_style(acceptor_chain, acceptor_resi, "orangeCarbon") +
        _create_dashed_line(hb.hydrogen.coords.x, hb.hydrogen.coords.y, hb.hydrogen.coords.z,
                           hb.acceptor.coords.x, hb.acceptor.coords.y, hb.acceptor.coords.z) +
        _create_label(hb.get_donor_residue(), hb.donor.coords.x, hb.donor.coords.y, hb.donor.coords.z,
                     "cyan", "black", 14) +
        _create_label(hb.get_acceptor_residue(), hb.acceptor.coords.x, hb.acceptor.coords.y, hb.acceptor.coords.z,
                     "orange", "white", 14) +
        f"viewer.zoomTo({{chain: ['{donor_chain}', '{acceptor_chain}'], resi: [{donor_resi}, {acceptor_resi}]}});"
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)


def generate_halogen_bond_viewer_js(
    xb: "HalogenBond", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for halogen bond 3D visualization.

    :param xb: Halogen bond interaction
    :type xb: HalogenBond
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    donor_chain = xb.donor_atom.chain_id
    donor_resi = xb.donor_atom.res_seq
    acceptor_chain = xb.acceptor.chain_id
    acceptor_resi = xb.acceptor.res_seq

    pdb_escaped = _escape_pdb_content(pdb_content)

    # Build viewer-specific code using helpers
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(donor_chain, donor_resi, "purpleCarbon") +
        _create_residue_stick_style(acceptor_chain, acceptor_resi, "orangeCarbon") +
        _create_dashed_line(xb.halogen.coords.x, xb.halogen.coords.y, xb.halogen.coords.z,
                           xb.acceptor.coords.x, xb.acceptor.coords.y, xb.acceptor.coords.z) +
        _create_label(xb.get_donor_residue(), xb.donor_atom.coords.x, xb.donor_atom.coords.y, xb.donor_atom.coords.z,
                     "purple", "white", 14) +
        _create_label(xb.get_acceptor_residue(), xb.acceptor.coords.x, xb.acceptor.coords.y, xb.acceptor.coords.z,
                     "orange", "white", 14) +
        f"viewer.zoomTo({{chain: ['{donor_chain}', '{acceptor_chain}'], resi: [{donor_resi}, {acceptor_resi}]}});"
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)


def generate_pi_interaction_viewer_js(
    pi: "PiInteraction", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for π interaction 3D visualization.

    :param pi: π interaction
    :type pi: PiInteraction
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    donor_chain = pi.donor.chain_id
    donor_resi = pi.donor.res_seq
    pi_res_seq = pi.pi_atoms[0].res_seq if pi.pi_atoms else ""
    pi_chain_id = pi.pi_atoms[0].chain_id if pi.pi_atoms else ""

    pdb_escaped = _escape_pdb_content(pdb_content)

    # Build viewer-specific code using helpers
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(donor_chain, donor_resi, "cyanCarbon") +
        _create_residue_stick_style(pi_chain_id, pi_res_seq, "greenCarbon") +
        _create_pi_center_sphere(pi.pi_center.x, pi.pi_center.y, pi.pi_center.z, "green") +
        _create_dashed_line_to_center(pi.hydrogen.coords.x, pi.hydrogen.coords.y, pi.hydrogen.coords.z,
                                      pi.pi_center.x, pi.pi_center.y, pi.pi_center.z) +
        _create_label(pi.get_donor_residue(), pi.donor.coords.x, pi.donor.coords.y, pi.donor.coords.z,
                     "cyan", "black", 14) +
        _create_label(pi.get_acceptor_residue(), pi.pi_center.x, pi.pi_center.y, pi.pi_center.z,
                     "green", "white", 14) +
        _create_distance_label(f"{pi.distance:.2f} Å", pi.hydrogen.coords.x, pi.hydrogen.coords.y, pi.hydrogen.coords.z,
                              pi.pi_center.x, pi.pi_center.y, pi.pi_center.z) +
        f"viewer.zoomTo({{chain: ['{donor_chain}', '{pi_chain_id}'], resi: [{donor_resi}, {pi_res_seq}]}});"
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)


def generate_pi_pi_stacking_viewer_js(
    pi_pi: "PiPiStacking", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for π-π stacking 3D visualization.

    :param pi_pi: π-π stacking interaction
    :type pi_pi: PiPiStacking
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    ring1_res = pi_pi.ring1_atoms[0].res_seq if pi_pi.ring1_atoms else ""
    ring1_chain = pi_pi.ring1_atoms[0].chain_id if pi_pi.ring1_atoms else ""
    ring2_res = pi_pi.ring2_atoms[0].res_seq if pi_pi.ring2_atoms else ""
    ring2_chain = pi_pi.ring2_atoms[0].chain_id if pi_pi.ring2_atoms else ""

    pdb_escaped = _escape_pdb_content(pdb_content)

    # Build viewer-specific code using helpers
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(ring1_chain, ring1_res, "cyanCarbon") +
        _create_residue_stick_style(ring2_chain, ring2_res, "magentaCarbon") +
        _create_pi_center_sphere(pi_pi.ring1_center.x, pi_pi.ring1_center.y, pi_pi.ring1_center.z, "cyan") +
        _create_pi_center_sphere(pi_pi.ring2_center.x, pi_pi.ring2_center.y, pi_pi.ring2_center.z, "magenta") +
        _create_dashed_line(pi_pi.ring1_center.x, pi_pi.ring1_center.y, pi_pi.ring1_center.z,
                           pi_pi.ring2_center.x, pi_pi.ring2_center.y, pi_pi.ring2_center.z) +
        _create_label(pi_pi.ring1_residue, pi_pi.ring1_center.x, pi_pi.ring1_center.y, pi_pi.ring1_center.z,
                     "cyan", "black", 14) +
        _create_label(pi_pi.ring2_residue, pi_pi.ring2_center.x, pi_pi.ring2_center.y, pi_pi.ring2_center.z,
                     "magenta", "white", 14) +
        _create_distance_label(f"{pi_pi._distance:.2f} Å ({pi_pi.stacking_type})",
                              pi_pi.ring1_center.x, pi_pi.ring1_center.y, pi_pi.ring1_center.z,
                              pi_pi.ring2_center.x, pi_pi.ring2_center.y, pi_pi.ring2_center.z) +
        f"viewer.zoomTo({{chain: ['{ring1_chain}', '{ring2_chain}'], resi: [{ring1_res}, {ring2_res}]}});"
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)


def generate_carbonyl_interaction_viewer_js(
    carbonyl: "CarbonylInteraction", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for carbonyl n→π* interaction 3D visualization.

    :param carbonyl: Carbonyl interaction
    :type carbonyl: CarbonylInteraction
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    pdb_escaped = _escape_pdb_content(pdb_content)

    # Build viewer-specific code using helpers
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(carbonyl.donor_oxygen.chain_id, carbonyl.donor_oxygen.res_seq, "redCarbon") +
        _create_residue_stick_style(carbonyl.acceptor_carbon.chain_id, carbonyl.acceptor_carbon.res_seq, "blueCarbon") +
        _create_dashed_line(carbonyl.donor_oxygen.coords.x, carbonyl.donor_oxygen.coords.y, carbonyl.donor_oxygen.coords.z,
                           carbonyl.acceptor_carbon.coords.x, carbonyl.acceptor_carbon.coords.y, carbonyl.acceptor_carbon.coords.z) +
        _create_label(carbonyl.get_donor_residue(), carbonyl.donor_oxygen.coords.x, carbonyl.donor_oxygen.coords.y, carbonyl.donor_oxygen.coords.z,
                     "red", "white", 14) +
        _create_label(carbonyl.get_acceptor_residue(), carbonyl.acceptor_carbon.coords.x, carbonyl.acceptor_carbon.coords.y, carbonyl.acceptor_carbon.coords.z,
                     "blue", "white", 14) +
        _create_distance_label(f"{carbonyl.distance:.2f} Å",
                              carbonyl.donor_oxygen.coords.x, carbonyl.donor_oxygen.coords.y, carbonyl.donor_oxygen.coords.z,
                              carbonyl.acceptor_carbon.coords.x, carbonyl.acceptor_carbon.coords.y, carbonyl.acceptor_carbon.coords.z) +
        f"viewer.zoomTo({{chain: ['{carbonyl.donor_oxygen.chain_id}', '{carbonyl.acceptor_carbon.chain_id}'], resi: [{carbonyl.donor_oxygen.res_seq}, {carbonyl.acceptor_carbon.res_seq}]}});"
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)


def generate_n_pi_interaction_viewer_js(
    n_pi: "NPiInteraction", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for n→π* interaction 3D visualization.

    :param n_pi: n→π* interaction
    :type n_pi: NPiInteraction
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    pi_res_seq = n_pi.pi_atoms[0].res_seq if n_pi.pi_atoms else ""
    pi_chain_id = n_pi.pi_atoms[0].chain_id if n_pi.pi_atoms else ""

    pdb_escaped = _escape_pdb_content(pdb_content)

    # Build viewer-specific code using helpers
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(n_pi.lone_pair_atom.chain_id, n_pi.lone_pair_atom.res_seq, "orangeCarbon") +
        _create_residue_stick_style(pi_chain_id, pi_res_seq, "tealCarbon") +
        _create_pi_center_sphere(n_pi.pi_center.x, n_pi.pi_center.y, n_pi.pi_center.z, "teal") +
        _create_dashed_line_to_center(n_pi.lone_pair_atom.coords.x, n_pi.lone_pair_atom.coords.y, n_pi.lone_pair_atom.coords.z,
                                      n_pi.pi_center.x, n_pi.pi_center.y, n_pi.pi_center.z) +
        _create_label(n_pi.get_donor_residue(), n_pi.lone_pair_atom.coords.x, n_pi.lone_pair_atom.coords.y, n_pi.lone_pair_atom.coords.z,
                     "orange", "white", 14) +
        _create_label(n_pi.get_acceptor_residue(), n_pi.pi_center.x, n_pi.pi_center.y, n_pi.pi_center.z,
                     "teal", "white", 14) +
        _create_distance_label(f"{n_pi.distance:.2f} Å",
                              n_pi.lone_pair_atom.coords.x, n_pi.lone_pair_atom.coords.y, n_pi.lone_pair_atom.coords.z,
                              n_pi.pi_center.x, n_pi.pi_center.y, n_pi.pi_center.z) +
        f"viewer.zoomTo({{chain: ['{n_pi.lone_pair_atom.chain_id}', '{pi_chain_id}'], resi: [{n_pi.lone_pair_atom.res_seq}, {pi_res_seq}]}});"
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)

def generate_water_bridge_viewer_js(
    water_bridge: "WaterBridge", pdb_content: str, viewer_id: str
) -> str:
    """Generate JavaScript code for water bridge 3D visualization.

    Visualizes water-mediated hydrogen bond networks with donor and acceptor
    residues highlighted and water molecules shown as stick models with all atoms.

    :param water_bridge: Water bridge interaction
    :type water_bridge: WaterBridge
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    donor = water_bridge.get_donor()
    acceptor = water_bridge.get_acceptor()
    donor_chain = donor.chain_id
    donor_resi = donor.res_seq
    acceptor_chain = acceptor.chain_id
    acceptor_resi = acceptor.res_seq

    # Get donor atom coordinates
    donor_x = donor.coords.x
    donor_y = donor.coords.y
    donor_z = donor.coords.z

    # Get acceptor atom coordinates
    acceptor_x = acceptor.coords.x
    acceptor_y = acceptor.coords.y
    acceptor_z = acceptor.coords.z

    # Extract water residue data (resi, oxygen coordinates)
    water_resi_list = []
    water_coords = []

    # Get bridge path to find water oxygen atoms
    for hbond in water_bridge.bridge_path:
        donor_hb = hbond.get_donor()
        acceptor_hb = hbond.get_acceptor()

        # Check if either side is water
        if hasattr(donor_hb, 'res_name') and donor_hb.res_name in WATER_MOLECULES:
            if donor_hb.res_seq not in water_resi_list:
                water_resi_list.append(donor_hb.res_seq)
                water_coords.append({
                    'resi': donor_hb.res_seq,
                    'x': donor_hb.coords.x,
                    'y': donor_hb.coords.y,
                    'z': donor_hb.coords.z
                })
        if hasattr(acceptor_hb, 'res_name') and acceptor_hb.res_name in WATER_MOLECULES:
            if acceptor_hb.res_seq not in water_resi_list:
                water_resi_list.append(acceptor_hb.res_seq)
                water_coords.append({
                    'resi': acceptor_hb.res_seq,
                    'x': acceptor_hb.coords.x,
                    'y': acceptor_hb.coords.y,
                    'z': acceptor_hb.coords.z
                })

    pdb_escaped = _escape_pdb_content(pdb_content)
    water_resi_json = "[" + ", ".join(str(r) for r in water_resi_list) + "]"

    # Create JavaScript water coordinate data and labels
    water_coords_js = "{"
    water_labels_js = "{"
    for i, wc in enumerate(water_coords):
        water_coords_js += f"{wc['resi']}: {{x: {wc['x']}, y: {wc['y']}, z: {wc['z']}}}, "
        # Get water residue label from water_residues list
        if i < len(water_bridge.water_residues):
            water_label = water_bridge.water_residues[i]
        else:
            water_label = f"W{wc['resi']}"
        water_labels_js += f"{wc['resi']}: '{water_label}', "
    water_coords_js = water_coords_js.rstrip(", ") + "}"
    water_labels_js = water_labels_js.rstrip(", ") + "}"

    # Build viewer-specific code (inline water-specific logic with helpers)
    viewer_code = (
        _create_cartoon_style() +
        _create_residue_stick_style(donor_chain, donor_resi, "cyanCarbon") +
        _create_residue_stick_style(acceptor_chain, acceptor_resi, "orangeCarbon") +
        f"""
                // Show water molecules as sticks
                let waterResis = {water_resi_json};
                waterResis.forEach(function(resi) {{
                    viewer.addStyle({{resi: resi}},
                                   {{stick: {{colorscheme: 'lightblueCarbon', radius: 0.25, opacity: 0.5}}}});
                }});

                // Draw dashed lines connecting donor → water → acceptor
                const dx = {donor_x}, dy = {donor_y}, dz = {donor_z};
                const ax = {acceptor_x}, ay = {acceptor_y}, az = {acceptor_z};
                let waterCoords = {water_coords_js};

                // Draw connections from donor to each water
                for (let waterResi in waterCoords) {{
                    let waterO = waterCoords[waterResi];
                    const wdx = waterO.x - dx, wdy = waterO.y - dy, wdz = waterO.z - dz;
                    const wdist = Math.sqrt(wdx*wdx + wdy*wdy + wdz*wdz);
                    if (wdist > 0) {{
                        const offset = 0.4;
                        const ratio_start = offset / wdist;
                        const ratio_end = (wdist - offset) / wdist;
                        viewer.addCylinder({{
                            start: {{x: dx + wdx*ratio_start, y: dy + wdy*ratio_start, z: dz + wdz*ratio_start}},
                            end: {{x: dx + wdx*ratio_end, y: dy + wdy*ratio_end, z: dz + wdz*ratio_end}},
                            radius: 0.15, color: 'yellow', dashed: true
                        }});
                    }}
                }}

                // Draw connections from each water to acceptor
                for (let waterResi in waterCoords) {{
                    let waterO = waterCoords[waterResi];
                    const wadx = ax - waterO.x, wady = ay - waterO.y, wadz = az - waterO.z;
                    const wadist = Math.sqrt(wadx*wadx + wady*wady + wadz*wadz);
                    if (wadist > 0) {{
                        const offset = 0.4;
                        const ratio_start = offset / wadist;
                        const ratio_end = (wadist - offset) / wadist;
                        viewer.addCylinder({{
                            start: {{x: waterO.x + wadx*ratio_start, y: waterO.y + wady*ratio_start, z: waterO.z + wadz*ratio_start}},
                            end: {{x: waterO.x + wadx*ratio_end, y: waterO.y + wady*ratio_end, z: waterO.z + wadz*ratio_end}},
                            radius: 0.15, color: 'yellow', dashed: true
                        }});
                    }}
                }}
        """ +
        _create_label(water_bridge.get_donor_residue(), donor_x, donor_y, donor_z,
                     "cyan", "black", 12) +
        _create_label(water_bridge.get_acceptor_residue(), acceptor_x, acceptor_y, acceptor_z,
                     "orange", "white", 12) +
        f"""
                // Add labels for water molecules
                let waterLabels = {water_labels_js};
                for (let waterResi in waterCoords) {{
                    let waterO = waterCoords[waterResi];
                    let waterLabel = waterLabels[waterResi] || ('W' + waterResi);
                    viewer.addLabel(waterLabel,
                                   {{position: {{x: waterO.x, y: waterO.y, z: waterO.z}},
                                     backgroundColor: 'lightblue', fontColor: 'black', fontSize: 11}});
                }}
                viewer.zoomTo();
        """
    )

    return _create_viewer_init_wrapper(viewer_id, pdb_escaped, viewer_code)


def generate_ligand_interactions_viewer_js(
    interactions_data: list, pdb_content: str, viewer_id: str, ligand_res: str = None
) -> str:
    """Generate JavaScript code for unified 3D visualization of all ligand interactions.

    :param interactions_data: List of interaction dicts with donor/acceptor atom info
    :type interactions_data: list
    :param pdb_content: PDB file content
    :type pdb_content: str
    :param viewer_id: Unique ID for the viewer div
    :type viewer_id: str
    :param ligand_res: Ligand residue name (e.g., "GTP") for styling
    :type ligand_res: str
    :returns: JavaScript code to initialize the viewer
    :rtype: str
    """
    import json
    pdb_escaped = _escape_pdb_content(pdb_content)
    ligand_name = ligand_res or ""
    interactions_json = json.dumps(interactions_data or [])

    # Build viewer-specific code for dynamic interaction drawing
    viewer_code = f"""
                // Show protein/RNA as cartoon
                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.4}}}});

                // Highlight ligand atoms
                viewer.addStyle({{resn: '{ligand_name}'}},
                               {{stick: {{colorscheme: 'cyanCarbon', radius: 0.20, showNonBonded: true}}}},
                               {{cartoon: {{color: 'yellow', thickness: 0.8, opacity: 0.8}}}});

                // Show interacting residues as sticks
                viewer.addStyle({{not: {{resn: '{ligand_name}'}}}},
                               {{stick: {{colorscheme: 'whiteCarbon', radius: 0.15}}}});

                // Draw dashed interaction lines using actual interaction data
                let interactions = {interactions_json};
                let atoms = model.atoms;

                // Helper function to find atom by name and residue
                function findAtom(atomName, residueId) {{
                    // Parse residue ID (e.g., "A:GTP:300" or "B:ASP:42")
                    let parts = residueId.split(':');
                    let chain = parts[0];
                    let resName = parts[1];
                    let resSeq = parts[2];

                    return atoms.find(a =>
                        a.atom === atomName &&
                        a.chain === chain &&
                        a.resn === resName &&
                        a.resi == resSeq
                    );
                }}

                // Track labeled residues to avoid duplicates
                let labeledResidues = new Set();

                // Draw dashed lines for each interaction and collect residues
                interactions.forEach(interaction => {{
                    let donorAtom = findAtom(interaction.donor_atom, interaction.donor_res);
                    let acceptorAtom = findAtom(interaction.acceptor_atom, interaction.acceptor_res);

                    if (donorAtom) {{
                        try {{
                            // Determine endpoint: π-center for π-interactions, acceptor atom otherwise
                            let endPoint;
                            let endLabel = interaction.acceptor_res;
                            let endColor = 'orange';

                            if (interaction.type && (interaction.type.includes('π') || interaction.type.includes('pi')) && interaction.pi_center) {{
                                // π-interaction: draw to ring center
                                endPoint = interaction.pi_center;
                                endLabel = interaction.acceptor_res + ' (π)';
                                endColor = 'green';

                                // Add green sphere at π center
                                viewer.addSphere({{
                                    center: endPoint,
                                    radius: 0.3,
                                    color: 'green',
                                    alpha: 0.7
                                }});
                            }} else if (acceptorAtom) {{
                                // Regular interaction: draw to acceptor atom
                                endPoint = {{x: acceptorAtom.x, y: acceptorAtom.y, z: acceptorAtom.z}};
                            }} else {{
                                // No endpoint found, skip
                                return;
                            }}

                            // Draw dashed cylinder
                            viewer.addCylinder({{
                                start: {{x: donorAtom.x, y: donorAtom.y, z: donorAtom.z}},
                                end: endPoint,
                                radius: 0.15,
                                color: 'yellow',
                                dashed: true
                            }});

                            // Add label for donor residue if not already labeled
                            if (!labeledResidues.has(interaction.donor_res)) {{
                                viewer.addLabel(interaction.donor_res,
                                               {{position: {{x: donorAtom.x, y: donorAtom.y, z: donorAtom.z}},
                                                 backgroundColor: 'cyan', fontColor: 'black', fontSize: 12}});
                                labeledResidues.add(interaction.donor_res);
                            }}

                            // Add label for acceptor residue if not already labeled
                            if (!labeledResidues.has(endLabel)) {{
                                viewer.addLabel(endLabel,
                                               {{position: endPoint,
                                                 backgroundColor: endColor, fontColor: 'white', fontSize: 12}});
                                labeledResidues.add(endLabel);
                            }}
                        }} catch (e) {{
                            console.warn('Could not draw interaction line:', interaction, e);
                        }}
                    }}
                }});

                viewer.zoomTo();
    """

    # Create custom wrapper that uses model for atom access
    javascript = f"""
    (function() {{
        function init3Dmol() {{
            if (typeof $3Dmol === 'undefined') {{
                setTimeout(init3Dmol, 100);
                return;
            }}

            let element = document.getElementById("{viewer_id}");
            if (!element) {{
                setTimeout(init3Dmol, 100);
                return;
            }}

            try {{
                let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: 'white'}});
                window.{viewer_id}_instance = viewer;
                let pdbData = `{pdb_escaped}`;

                let model = viewer.addModel(pdbData, "pdb", {{keepH: true}});

                // Viewer-specific code
                {viewer_code}

                viewer.render();

                // Force multiple resizes to ensure proper rendering
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 100);
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 300);
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 500);
            }} catch(error) {{
                console.error('Error initializing ligand interactions viewer:', error);
                console.error(error.stack);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript

"""
py3Dmol visualization utilities for molecular interactions.

This module provides functions to generate py3Dmol viewer JavaScript code
for visualizing different types of molecular interactions. These functions
can be used in web applications, Jupyter notebooks, and other contexts.
"""

from typing import TYPE_CHECKING

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

    javascript = f"""
    (function() {{
        // Wait for both 3Dmol to load and div to be ready
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
                window.{viewer_id}_instance = viewer;  // Store for PNG export
                let pdbData = `{pdb_escaped}`;

                viewer.addModel(pdbData, "pdb", {{keepH: true}});
                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                // Show donor residue with hydrogens
                viewer.addStyle({{chain: '{donor_chain}', resi: {donor_resi}}},
                               {{stick: {{colorscheme: 'cyanCarbon', radius: 0.20, showNonBonded: false}}}});

                // Show acceptor residue with hydrogens
                viewer.addStyle({{chain: '{acceptor_chain}', resi: {acceptor_resi}}},
                               {{stick: {{colorscheme: 'orangeCarbon', radius: 0.20, showNonBonded: false}}}});

                // Add dashed line for hydrogen bond (shortened to avoid overlapping atoms)
                const hx = {hb.hydrogen.coords.x}, hy = {hb.hydrogen.coords.y}, hz = {hb.hydrogen.coords.z};
                const ax = {hb.acceptor.coords.x}, ay = {hb.acceptor.coords.y}, az = {hb.acceptor.coords.z};
                const dx = ax - hx, dy = ay - hy, dz = az - hz;
                const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                const offset = 0.4; // Shorten by 0.4 Å on each end
                const ratio_start = offset / dist;
                const ratio_end = (dist - offset) / dist;

                viewer.addCylinder({{
                    start: {{x: hx + dx*ratio_start, y: hy + dy*ratio_start, z: hz + dz*ratio_start}},
                    end: {{x: hx + dx*ratio_end, y: hy + dy*ratio_end, z: hz + dz*ratio_end}},
                    radius: 0.15,
                    color: 'yellow',
                    dashed: true
                }});

                // Add labels
                viewer.addLabel("{hb.get_donor_residue()}",
                               {{position: {{x: {hb.donor.coords.x}, y: {hb.donor.coords.y}, z: {hb.donor.coords.z}}},
                                backgroundColor: 'cyan', fontColor: 'black', fontSize: 14}});

                viewer.addLabel("{hb.get_acceptor_residue()}",
                               {{position: {{x: {hb.acceptor.coords.x}, y: {hb.acceptor.coords.y}, z: {hb.acceptor.coords.z}}},
                                backgroundColor: 'orange', fontColor: 'white', fontSize: 14}});

                viewer.zoomTo({{chain: ['{donor_chain}', '{acceptor_chain}'], resi: [{donor_resi}, {acceptor_resi}]}});
                viewer.render();
                viewer.zoom(1.2);

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

        // Start with a delay to ensure div is in DOM and sized
        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript


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

    javascript = f"""
    (function() {{
        // Wait for both 3Dmol to load and div to be ready
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
                window.{viewer_id}_instance = viewer;  // Store for PNG export
                let pdbData = `{pdb_escaped}`;

                viewer.addModel(pdbData, "pdb", {{keepH: true}});
                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                viewer.addStyle({{chain: '{donor_chain}', resi: {donor_resi}}},
                               {{stick: {{colorscheme: 'purpleCarbon', radius: 0.25}}}});

                viewer.addStyle({{chain: '{acceptor_chain}', resi: {acceptor_resi}}},
                               {{stick: {{colorscheme: 'orangeCarbon', radius: 0.25}}}});

                // Add dashed line for halogen bond (shortened to avoid overlapping atoms)
                const xx = {xb.halogen.coords.x}, xy = {xb.halogen.coords.y}, xz = {xb.halogen.coords.z};
                const axb = {xb.acceptor.coords.x}, ayb = {xb.acceptor.coords.y}, azb = {xb.acceptor.coords.z};
                const dxb = axb - xx, dyb = ayb - xy, dzb = azb - xz;
                const dist_xb = Math.sqrt(dxb*dxb + dyb*dyb + dzb*dzb);
                const offset_xb = 0.4;
                const ratio_start_xb = offset_xb / dist_xb;
                const ratio_end_xb = (dist_xb - offset_xb) / dist_xb;

                viewer.addCylinder({{
                    start: {{x: xx + dxb*ratio_start_xb, y: xy + dyb*ratio_start_xb, z: xz + dzb*ratio_start_xb}},
                    end: {{x: xx + dxb*ratio_end_xb, y: xy + dyb*ratio_end_xb, z: xz + dzb*ratio_end_xb}},
                    radius: 0.15,
                    color: 'orange',
                    dashed: true
                }});

                // Add labels
                viewer.addLabel("{xb.get_donor_residue()}",
                               {{position: {{x: {xb.donor_atom.coords.x}, y: {xb.donor_atom.coords.y}, z: {xb.donor_atom.coords.z}}},
                                backgroundColor: 'purple', fontColor: 'white', fontSize: 14}});

                viewer.addLabel("{xb.get_acceptor_residue()}",
                               {{position: {{x: {xb.acceptor.coords.x}, y: {xb.acceptor.coords.y}, z: {xb.acceptor.coords.z}}},
                                backgroundColor: 'orange', fontColor: 'white', fontSize: 14}});

                viewer.zoomTo({{chain: ['{donor_chain}', '{acceptor_chain}'], resi: [{donor_resi}, {acceptor_resi}]}});
                viewer.render();
                viewer.zoom(1.2);

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

        // Start with a delay to ensure div is in DOM and sized
        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript


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

    # Get pi system residue
    pi_res_seq = pi.pi_atoms[0].res_seq if pi.pi_atoms else ""
    pi_chain_id = pi.pi_atoms[0].chain_id if pi.pi_atoms else ""

    pdb_escaped = _escape_pdb_content(pdb_content)

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
                window.{viewer_id}_instance = viewer;  // Store for PNG export
                let pdb_data = `{pdb_escaped}`;
                viewer.addModel(pdb_data, "pdb", {{keepH: true}});

                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                // Highlight donor
                viewer.setStyle({{resi: '{donor_resi}', chain: '{donor_chain}'}}, {{stick: {{colorscheme: 'cyanCarbon'}}}});

                // Highlight π system
                viewer.setStyle({{resi: '{pi_res_seq}', chain: '{pi_chain_id}'}}, {{stick: {{colorscheme: 'greenCarbon'}}}});

                // Add sphere at π center to show the interaction point
                viewer.addSphere({{
                    center: {{x: {pi.pi_center.x}, y: {pi.pi_center.y}, z: {pi.pi_center.z}}},
                    radius: 0.3,
                    color: 'green',
                    alpha: 0.7
                }});

                // Highlight the hydrogen/X-atom as a sphere
                viewer.addStyle({{serial: {pi.hydrogen.serial}}},
                               {{sphere: {{color: 'yellow', radius: 0.5}}}});

                // Add dashed line from X-atom to π center (shortened to avoid overlapping)
                const pix = {pi.hydrogen.coords.x}, piy = {pi.hydrogen.coords.y}, piz = {pi.hydrogen.coords.z};
                const pcx = {pi.pi_center.x}, pcy = {pi.pi_center.y}, pcz = {pi.pi_center.z};
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

                // Add labels
                viewer.addLabel('{pi.get_donor_residue()}',
                              {{position: {{x: {pi.donor.coords.x}, y: {pi.donor.coords.y}, z: {pi.donor.coords.z}}},
                               backgroundColor: 'cyan', fontColor: 'black', fontSize: 14}});

                viewer.addLabel('{pi.get_acceptor_residue()}',
                              {{position: {{x: {pi.pi_center.x}, y: {pi.pi_center.y}, z: {pi.pi_center.z}}},
                               backgroundColor: 'green', fontColor: 'white', fontSize: 14}});

                // Distance label (at midpoint between hydrogen and pi center)
                viewer.addLabel('{pi.distance:.2f} Å',
                              {{position: {{x: ({pi.hydrogen.coords.x} + {pi.pi_center.x}) / 2,
                                          y: ({pi.hydrogen.coords.y} + {pi.pi_center.y}) / 2,
                                          z: ({pi.hydrogen.coords.z} + {pi.pi_center.z}) / 2}},
                               backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                viewer.zoomTo({{chain: ['{donor_chain}', '{pi_chain_id}'], resi: ['{donor_resi}', '{pi_res_seq}']}});
                viewer.render();

                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
            }} catch (error) {{
                console.error("Error creating viewer:", error);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript


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
    # Get residue info from ring atoms
    ring1_res = pi_pi.ring1_atoms[0].res_seq if pi_pi.ring1_atoms else ""
    ring1_chain = pi_pi.ring1_atoms[0].chain_id if pi_pi.ring1_atoms else ""
    ring2_res = pi_pi.ring2_atoms[0].res_seq if pi_pi.ring2_atoms else ""
    ring2_chain = pi_pi.ring2_atoms[0].chain_id if pi_pi.ring2_atoms else ""

    pdb_escaped = _escape_pdb_content(pdb_content)

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
                window.{viewer_id}_instance = viewer;  // Store for PNG export

                let pdb_data = `{pdb_escaped}`;
                viewer.addModel(pdb_data, "pdb", {{keepH: true}});

                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                // Highlight both aromatic rings
                viewer.setStyle({{resi: '{ring1_res}', chain: '{ring1_chain}'}}, {{stick: {{colorscheme: 'cyanCarbon'}}}});
                viewer.setStyle({{resi: '{ring2_res}', chain: '{ring2_chain}'}}, {{stick: {{colorscheme: 'magentaCarbon'}}}});

                // Add residue labels
                viewer.addLabel('{pi_pi.ring1_residue}',
                              {{position: {{x: {pi_pi.ring1_center.x}, y: {pi_pi.ring1_center.y}, z: {pi_pi.ring1_center.z}}},
                               backgroundColor: 'cyan', fontColor: 'black', fontSize: 14}});

                viewer.addLabel('{pi_pi.ring2_residue}',
                              {{position: {{x: {pi_pi.ring2_center.x}, y: {pi_pi.ring2_center.y}, z: {pi_pi.ring2_center.z}}},
                               backgroundColor: 'magenta', fontColor: 'white', fontSize: 14}});

                // Add distance label
                viewer.addLabel('{pi_pi._distance:.2f} Å ({pi_pi.stacking_type})',
                              {{position: {{x: {pi_pi.midpoint.x}, y: {pi_pi.midpoint.y}, z: {pi_pi.midpoint.z}}},
                               backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                // Add line between ring centers (shortened to avoid overlapping)
                const r1x = {pi_pi.ring1_center.x}, r1y = {pi_pi.ring1_center.y}, r1z = {pi_pi.ring1_center.z};
                const r2x = {pi_pi.ring2_center.x}, r2y = {pi_pi.ring2_center.y}, r2z = {pi_pi.ring2_center.z};
                const drx = r2x - r1x, dry = r2y - r1y, drz = r2z - r1z;
                const dist_pipi = Math.sqrt(drx*drx + dry*dry + drz*drz);
                const offset_pipi = 0.3; // Smaller offset for ring centers
                const ratio_start_pipi = offset_pipi / dist_pipi;
                const ratio_end_pipi = (dist_pipi - offset_pipi) / dist_pipi;

                viewer.addCylinder({{
                    start: {{x: r1x + drx*ratio_start_pipi, y: r1y + dry*ratio_start_pipi, z: r1z + drz*ratio_start_pipi}},
                    end: {{x: r1x + drx*ratio_end_pipi, y: r1y + dry*ratio_end_pipi, z: r1z + drz*ratio_end_pipi}},
                    radius: 0.1,
                    color: 'purple',
                    dashed: true
                }});

                viewer.zoomTo({{chain: ['{ring1_chain}', '{ring2_chain}'], resi: ['{ring1_res}', '{ring2_res}']}});
                viewer.render();

                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
            }} catch (error) {{
                console.error("Error creating viewer:", error);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript


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
                window.{viewer_id}_instance = viewer;  // Store for PNG export

                let pdb_data = `{pdb_escaped}`;
                viewer.addModel(pdb_data, "pdb", {{keepH: true}});

                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                // Highlight donor and acceptor residues
                viewer.setStyle({{resi: '{carbonyl.donor_oxygen.res_seq}', chain: '{carbonyl.donor_oxygen.chain_id}'}},
                              {{stick: {{colorscheme: 'redCarbon'}}}});
                viewer.setStyle({{resi: '{carbonyl.acceptor_carbon.res_seq}', chain: '{carbonyl.acceptor_carbon.chain_id}'}},
                              {{stick: {{colorscheme: 'blueCarbon'}}}});

                // Add residue labels
                viewer.addLabel('{carbonyl.get_donor_residue()}',
                              {{position: {{x: {carbonyl.donor_oxygen.coords.x}, y: {carbonyl.donor_oxygen.coords.y}, z: {carbonyl.donor_oxygen.coords.z}}},
                               backgroundColor: 'red', fontColor: 'white', fontSize: 14}});

                viewer.addLabel('{carbonyl.get_acceptor_residue()}',
                              {{position: {{x: {carbonyl.acceptor_carbon.coords.x}, y: {carbonyl.acceptor_carbon.coords.y}, z: {carbonyl.acceptor_carbon.coords.z}}},
                               backgroundColor: 'blue', fontColor: 'white', fontSize: 14}});

                // Add distance label (at midpoint)
                let carbonyl_midpoint = {{
                    x: ({carbonyl.donor_oxygen.coords.x} + {carbonyl.acceptor_carbon.coords.x}) / 2,
                    y: ({carbonyl.donor_oxygen.coords.y} + {carbonyl.acceptor_carbon.coords.y}) / 2,
                    z: ({carbonyl.donor_oxygen.coords.z} + {carbonyl.acceptor_carbon.coords.z}) / 2
                }};
                viewer.addLabel('{carbonyl.distance:.2f} Å',
                              {{position: carbonyl_midpoint,
                               backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                // Add line from donor O to acceptor C (shortened to avoid overlapping)
                const cox = {carbonyl.donor_oxygen.coords.x}, coy = {carbonyl.donor_oxygen.coords.y}, coz = {carbonyl.donor_oxygen.coords.z};
                const cax = {carbonyl.acceptor_carbon.coords.x}, cay = {carbonyl.acceptor_carbon.coords.y}, caz = {carbonyl.acceptor_carbon.coords.z};
                const dcx = cax - cox, dcy = cay - coy, dcz = caz - coz;
                const dist_co = Math.sqrt(dcx*dcx + dcy*dcy + dcz*dcz);
                const offset_co = 0.4;
                const ratio_start_co = offset_co / dist_co;
                const ratio_end_co = (dist_co - offset_co) / dist_co;

                viewer.addCylinder({{
                    start: {{x: cox + dcx*ratio_start_co, y: coy + dcy*ratio_start_co, z: coz + dcz*ratio_start_co}},
                    end: {{x: cox + dcx*ratio_end_co, y: coy + dcy*ratio_end_co, z: coz + dcz*ratio_end_co}},
                    radius: 0.1,
                    color: 'orange',
                    dashed: true
                }});

                viewer.zoomTo({{chain: ['{carbonyl.donor_oxygen.chain_id}', '{carbonyl.acceptor_carbon.chain_id}'], resi: ['{carbonyl.donor_oxygen.res_seq}', '{carbonyl.acceptor_carbon.res_seq}']}});
                viewer.render();

                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
            }} catch (error) {{
                console.error("Error creating viewer:", error);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript


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
    # Get the π system residue sequence number from the first pi atom
    pi_res_seq = n_pi.pi_atoms[0].res_seq if n_pi.pi_atoms else ""
    pi_chain_id = n_pi.pi_atoms[0].chain_id if n_pi.pi_atoms else ""

    pdb_escaped = _escape_pdb_content(pdb_content)

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
                window.{viewer_id}_instance = viewer;  // Store for PNG export

                let pdb_data = `{pdb_escaped}`;
                viewer.addModel(pdb_data, "pdb", {{keepH: true}});

                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.3}}}});

                // Highlight lone pair atom residue
                viewer.setStyle({{resi: '{n_pi.lone_pair_atom.res_seq}', chain: '{n_pi.lone_pair_atom.chain_id}'}},
                              {{stick: {{colorscheme: 'orangeCarbon'}}}});

                // Highlight π system residue
                viewer.setStyle({{resi: '{pi_res_seq}', chain: '{pi_chain_id}'}},
                              {{stick: {{colorscheme: 'tealCarbon'}}}});

                // Add residue labels
                viewer.addLabel('{n_pi.get_donor_residue()}',
                              {{position: {{x: {n_pi.lone_pair_atom.coords.x}, y: {n_pi.lone_pair_atom.coords.y}, z: {n_pi.lone_pair_atom.coords.z}}},
                               backgroundColor: 'orange', fontColor: 'white', fontSize: 14}});

                viewer.addLabel('{n_pi.get_acceptor_residue()}',
                              {{position: {{x: {n_pi.pi_center.x}, y: {n_pi.pi_center.y}, z: {n_pi.pi_center.z}}},
                               backgroundColor: 'teal', fontColor: 'white', fontSize: 14}});

                // Add distance label (at midpoint)
                let npi_midpoint = {{
                    x: ({n_pi.lone_pair_atom.coords.x} + {n_pi.pi_center.x}) / 2,
                    y: ({n_pi.lone_pair_atom.coords.y} + {n_pi.pi_center.y}) / 2,
                    z: ({n_pi.lone_pair_atom.coords.z} + {n_pi.pi_center.z}) / 2
                }};
                viewer.addLabel('{n_pi.distance:.2f} Å',
                              {{position: npi_midpoint,
                               backgroundColor: 'black', fontColor: 'white', fontSize: 12}});

                // Add line from lone pair atom to π center (shortened to avoid overlapping)
                const npx = {n_pi.lone_pair_atom.coords.x}, npy = {n_pi.lone_pair_atom.coords.y}, npz = {n_pi.lone_pair_atom.coords.z};
                const npcx = {n_pi.pi_center.x}, npcy = {n_pi.pi_center.y}, npcz = {n_pi.pi_center.z};
                const dnpx = npcx - npx, dnpy = npcy - npy, dnpz = npcz - npz;
                const dist_npi = Math.sqrt(dnpx*dnpx + dnpy*dnpy + dnpz*dnpz);
                const offset_npi = 0.4;
                const ratio_start_npi = offset_npi / dist_npi;
                const ratio_end_npi = (dist_npi - offset_npi) / dist_npi;

                viewer.addCylinder({{
                    start: {{x: npx + dnpx*ratio_start_npi, y: npy + dnpy*ratio_start_npi, z: npz + dnpz*ratio_start_npi}},
                    end: {{x: npx + dnpx*ratio_end_npi, y: npy + dnpy*ratio_end_npi, z: npz + dnpz*ratio_end_npi}},
                    radius: 0.1,
                    color: 'green',
                    dashed: true
                }});

                viewer.zoomTo({{chain: ['{n_pi.lone_pair_atom.chain_id}', '{pi_chain_id}'], resi: ['{n_pi.lone_pair_atom.res_seq}', '{pi_res_seq}']}});
                viewer.render();

                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 100);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 300);
                setTimeout(function() {{ viewer.resize(); viewer.render(); }}, 500);
            }} catch (error) {{
                console.error("Error creating viewer:", error);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript

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
        if hasattr(donor_hb, 'res_name') and donor_hb.res_name in ['HOH', 'WAT', 'DOD', 'TIP3', 'TIP4', 'TIP5', 'W']:
            if donor_hb.res_seq not in water_resi_list:
                water_resi_list.append(donor_hb.res_seq)
                water_coords.append({
                    'resi': donor_hb.res_seq,
                    'x': donor_hb.coords.x,
                    'y': donor_hb.coords.y,
                    'z': donor_hb.coords.z
                })
        if hasattr(acceptor_hb, 'res_name') and acceptor_hb.res_name in ['HOH', 'WAT', 'DOD', 'TIP3', 'TIP4', 'TIP5', 'W']:
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

    javascript = f"""
    (function() {{
        // Wait for both 3Dmol to load and div to be ready
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
                window.{viewer_id}_instance = viewer;  // Store for PNG export
                let pdbData = `{pdb_escaped}`;

                viewer.addModel(pdbData, "pdb", {{keepH: true}});
                viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.2}}}});

                // Show donor residue in cyan with stick representation
                viewer.addStyle({{chain: '{donor_chain}', resi: {donor_resi}}},
                               {{stick: {{colorscheme: 'cyanCarbon', radius: 0.25}},
                                 cartoon: {{color: 'cyan', opacity: 0.7}}}});

                // Show acceptor residue in orange with stick representation
                viewer.addStyle({{chain: '{acceptor_chain}', resi: {acceptor_resi}}},
                               {{stick: {{colorscheme: 'orangeCarbon', radius: 0.25}},
                                 cartoon: {{color: 'orange', opacity: 0.7}}}});

                // Show water molecules as sticks with all atoms visible
                // Use opacity 0.5 to match PDB occupancy visualization
                let waterResis = {water_resi_json};
                waterResis.forEach(function(resi) {{
                    viewer.addStyle({{resi: resi}},
                                   {{stick: {{colorscheme: 'lightblueCarbon', radius: 0.25, opacity: 0.5}}}});
                }});

                // Draw dashed lines connecting donor → water → acceptor
                // Use embedded coordinates (same approach as hydrogen bond viewer)
                const dx = {donor_x}, dy = {donor_y}, dz = {donor_z};
                const ax = {acceptor_x}, ay = {acceptor_y}, az = {acceptor_z};

                let waterCoords = {water_coords_js};

                // Draw connections from donor to each water
                for (let waterResi in waterCoords) {{
                    let waterO = waterCoords[waterResi];
                    const wdx = waterO.x - dx;
                    const wdy = waterO.y - dy;
                    const wdz = waterO.z - dz;
                    const wdist = Math.sqrt(wdx*wdx + wdy*wdy + wdz*wdz);
                    if (wdist > 0) {{
                        const offset = 0.4;
                        const ratio_start = offset / wdist;
                        const ratio_end = (wdist - offset) / wdist;
                        viewer.addCylinder({{
                            start: {{x: dx + wdx*ratio_start, y: dy + wdy*ratio_start, z: dz + wdz*ratio_start}},
                            end: {{x: dx + wdx*ratio_end, y: dy + wdy*ratio_end, z: dz + wdz*ratio_end}},
                            radius: 0.15,
                            color: 'yellow',
                            dashed: true
                        }});
                    }}
                }}

                // Draw connections from each water to acceptor
                for (let waterResi in waterCoords) {{
                    let waterO = waterCoords[waterResi];
                    const wadx = ax - waterO.x;
                    const wady = ay - waterO.y;
                    const wadz = az - waterO.z;
                    const wadist = Math.sqrt(wadx*wadx + wady*wady + wadz*wadz);
                    if (wadist > 0) {{
                        const offset = 0.4;
                        const ratio_start = offset / wadist;
                        const ratio_end = (wadist - offset) / wadist;
                        viewer.addCylinder({{
                            start: {{x: waterO.x + wadx*ratio_start, y: waterO.y + wady*ratio_start, z: waterO.z + wadz*ratio_start}},
                            end: {{x: waterO.x + wadx*ratio_end, y: waterO.y + wady*ratio_end, z: waterO.z + wadz*ratio_end}},
                            radius: 0.15,
                            color: 'yellow',
                            dashed: true
                        }});
                    }}
                }}

                // Add labels for donor, water molecules, and acceptor
                viewer.addLabel("{water_bridge.get_donor_residue()}",
                               {{position: {{x: {donor_x}, y: {donor_y}, z: {donor_z}}},
                                 backgroundColor: 'cyan', fontColor: 'black', fontSize: 12}});
                viewer.addLabel("{water_bridge.get_acceptor_residue()}",
                               {{position: {{x: {acceptor_x}, y: {acceptor_y}, z: {acceptor_z}}},
                                 backgroundColor: 'orange', fontColor: 'white', fontSize: 12}});

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
                viewer.render();
                viewer.zoom(1.2);

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
                console.error('Error initializing water bridge viewer:', error);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript


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
    ligand_name = ligand_res or ""  # Use residue name if provided

    # Convert interactions_data to JSON for embedding in JavaScript
    interactions_json = json.dumps(interactions_data or [])

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

                            if (interaction.type && interaction.type.includes('pi') && interaction.pi_center) {{
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
                viewer.render();

                // Force resize to ensure proper rendering
                setTimeout(function() {{
                    viewer.resize();
                    viewer.render();
                }}, 100);

            }} catch(error) {{
                console.error('Error initializing ligand interactions viewer:', error);
                console.error(error.stack);
            }}
        }}

        setTimeout(init3Dmol, 400);
    }})();
    """
    return javascript

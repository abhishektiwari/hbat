API Usage Examples
------------------

While the HBAT command-line interface (CLI) and graphical user interface (GUI) are sufficient for most standard analyses, the HBAT's analysis API provides additional capabilities for advanced users who need to:

- **Integrate HBAT into automated pipelines** - Process hundreds of structures programmatically without manual intervention
- **Customize analysis workflows** - Apply complex filtering criteria or combine multiple analysis steps
- **Extract specific data subsets** - Focus on particular residues, regions, or interaction types
- **Perform statistical analyses** - Calculate distributions, correlations, or other metrics across multiple structures
- **Create custom visualizations** - Generate publication-quality figures or interactive plots
- **Combine with other tools** - Integrate HBAT results with molecular dynamics trajectories, docking scores, or structural databases

This section provides practical examples demonstrating how to use the analysis API for various analysis tasks, from simple single-structure analysis to complex comparative studies and high-throughput screening workflows.

Basic Analysis
--------------

Simple Hydrogen Bond Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from hbat.core.analysis import HBondAnalyzer
   
   # Initialize analyzer
   analyzer = HBondAnalyzer()
   
   # Analyze PDB file
   analyzer.analyze_file("1abc.pdb")
   
   # Print summary
   print(analyzer.get_results_summary())
   
   # Access individual results
   for hb in analyzer.hydrogen_bonds:
       print(f"H-bond: {hb.donor_residue} -> {hb.acceptor_residue}")
       print(f"Distance: {hb.distance:.2f} Å")
       print(f"Angle: {hb.angle * 180 / 3.14159:.1f}°")

Advanced Parameter Customization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from hbat.core.analysis import HBondAnalyzer, AnalysisParameters
   
   # Create custom parameters for drug design
   drug_params = AnalysisParameters(
       hb_distance_cutoff=2.8,        # Stricter distance
       hb_angle_cutoff=140.0,         # Stricter angle
       hb_donor_acceptor_cutoff=3.5,  # Tighter D...A distance
       analysis_mode="global"          # Include all interactions
   )
   
   analyzer = HBondAnalyzer(parameters=drug_params)
   analyzer.analyze_file("drug_target_complex.pdb")
   
   # Focus on strong interactions only
   strong_hbonds = [hb for hb in analyzer.hydrogen_bonds 
                    if hb.distance < 2.5]
   
   print(f"Found {len(strong_hbonds)} strong hydrogen bonds")

Multi-Type Interaction Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze all interaction types
   analyzer = HBondAnalyzer()
   analyzer.analyze_file("membrane_protein.pdb")
   
   print("=== Interaction Summary ===")
   print(f"Hydrogen bonds: {len(analyzer.hydrogen_bonds)}")
   print(f"Halogen bonds: {len(analyzer.halogen_bonds)}")
   print(f"π interactions: {len(analyzer.pi_interactions)}")
   print(f"Cooperative chains: {len(analyzer.cooperativity_chains)}")
   
   # Analyze cooperativity
   if analyzer.cooperativity_chains:
       print("\\n=== Cooperative Chains ===")
       for chain in analyzer.cooperativity_chains:
           print(f"Chain length: {chain.chain_length}")
           print(f"Chain type: {chain.chain_type}")
           print(chain)

Batch Processing
----------------

Processing Multiple Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import glob
   import csv
   from hbat.core.analysis import HBondAnalyzer
   
   # Process all PDB files in directory
   pdb_files = glob.glob("structures/*.pdb")
   results = []
   
   for pdb_file in pdb_files:
       analyzer = HBondAnalyzer()
       if analyzer.analyze_file(pdb_file):
           stats = analyzer.get_statistics()
           results.append({
               'file': pdb_file,
               'hydrogen_bonds': stats['hydrogen_bonds'],
               'halogen_bonds': stats['halogen_bonds'],
               'pi_interactions': stats['pi_interactions'],
               'total_interactions': stats['total_interactions']
           })
   
   # Save results to CSV
   with open('batch_results.csv', 'w', newline='') as csvfile:
       fieldnames = ['file', 'hydrogen_bonds', 'halogen_bonds', 
                     'pi_interactions', 'total_interactions']
       writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
       writer.writeheader()
       writer.writerows(results)

Comparative Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Compare native vs mutant structures
   structures = {
       'native': 'wild_type.pdb',
       'mutant': 'mutant_Y123F.pdb'
   }
   
   results = {}
   
   for name, pdb_file in structures.items():
       analyzer = HBondAnalyzer()
       analyzer.analyze_file(pdb_file)
       
       results[name] = {
           'hydrogen_bonds': len(analyzer.hydrogen_bonds),
           'cooperativity_chains': len(analyzer.cooperativity_chains),
           'avg_hb_distance': sum(hb.distance for hb in analyzer.hydrogen_bonds) / 
                             len(analyzer.hydrogen_bonds) if analyzer.hydrogen_bonds else 0
       }
   
   # Compare results
   print("Structure Comparison:")
   for metric in results['native'].keys():
       native_val = results['native'][metric]
       mutant_val = results['mutant'][metric]
       change = mutant_val - native_val
       print(f"{metric}: Native={native_val:.2f}, Mutant={mutant_val:.2f}, Change={change:+.2f}")

Data Export and Visualization
------------------------------

Detailed CSV Export
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import csv
   from hbat.core.analysis import HBondAnalyzer
   
   analyzer = HBondAnalyzer()
   analyzer.analyze_file("protein.pdb")
   
   # Export hydrogen bonds to CSV
   with open('hydrogen_bonds.csv', 'w', newline='') as csvfile:
       fieldnames = ['donor_residue', 'donor_atom', 'acceptor_residue', 
                     'acceptor_atom', 'distance', 'angle_degrees', 
                     'da_distance', 'bond_type']
       writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
       writer.writeheader()
       
       for hb in analyzer.hydrogen_bonds:
           writer.writerow({
               'donor_residue': hb.donor_residue,
               'donor_atom': hb.donor.name,
               'acceptor_residue': hb.acceptor_residue,
               'acceptor_atom': hb.acceptor.name,
               'distance': round(hb.distance, 3),
               'angle_degrees': round(hb.angle * 180 / 3.14159, 1),
               'da_distance': round(hb.donor_acceptor_distance, 3),
               'bond_type': hb.bond_type
           })

JSON Export with Metadata
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import json
   from datetime import datetime
   from hbat.core.analysis import HBondAnalyzer
   
   analyzer = HBondAnalyzer()
   analyzer.analyze_file("complex.pdb")
   
   # Create comprehensive results dictionary
   results = {
       'metadata': {
           'analysis_date': datetime.now().isoformat(),
           'pdb_file': 'complex.pdb',
           'parameters': {
               'hb_distance_cutoff': analyzer.parameters.hb_distance_cutoff,
               'hb_angle_cutoff': analyzer.parameters.hb_angle_cutoff,
               'analysis_mode': analyzer.parameters.analysis_mode
           }
       },
       'statistics': analyzer.get_statistics(),
       'interactions': {
           'hydrogen_bonds': [
               {
                   'donor': hb.donor_residue,
                   'acceptor': hb.acceptor_residue,
                   'distance': hb.distance,
                   'angle': hb.angle,
                   'type': hb.bond_type
               }
               for hb in analyzer.hydrogen_bonds
           ],
           'cooperativity_chains': [
               {
                   'length': chain.chain_length,
                   'type': chain.chain_type,
                   'description': str(chain)
               }
               for chain in analyzer.cooperativity_chains
           ]
       }
   }
   
   # Save to JSON with pretty formatting
   with open('analysis_results.json', 'w') as f:
       json.dump(results, f, indent=2, default=str)

Specialized Analysis Tasks
--------------------------

Drug-Target Interaction Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from hbat.core.analysis import HBondAnalyzer, AnalysisParameters
   
   # Custom parameters for drug analysis
   drug_params = AnalysisParameters(
       hb_distance_cutoff=3.2,
       hb_angle_cutoff=120.0,
       analysis_mode="global"
   )
   
   analyzer = HBondAnalyzer(parameters=drug_params)
   analyzer.analyze_file("drug_target.pdb")
   
   # Filter interactions involving the drug (assuming it's a HET residue)
   drug_interactions = []
   
   for hb in analyzer.hydrogen_bonds:
       # Check if either donor or acceptor is from drug
       if ('HET' in hb.donor_residue or 'HET' in hb.acceptor_residue or
           'LIG' in hb.donor_residue or 'LIG' in hb.acceptor_residue):
           drug_interactions.append(hb)
   
   print(f"Drug-target hydrogen bonds: {len(drug_interactions)}")
   for interaction in drug_interactions:
       print(f"  {interaction}")

Membrane Protein Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze interactions in membrane proteins
   analyzer = HBondAnalyzer()
   analyzer.analyze_file("membrane_protein.pdb")
   
   # Categorize interactions by region (transmembrane vs extracellular)
   # This assumes Z-coordinate indicates membrane position
   
   tm_interactions = []  # Transmembrane region
   ec_interactions = []  # Extracellular region
   
   for hb in analyzer.hydrogen_bonds:
       # Simple Z-coordinate based classification
       donor_z = hb.donor.coords.z
       acceptor_z = hb.acceptor.coords.z
       avg_z = (donor_z + acceptor_z) / 2
       
       if -20 < avg_z < 20:  # Transmembrane region
           tm_interactions.append(hb)
       elif avg_z > 20:      # Extracellular region
           ec_interactions.append(hb)
   
   print(f"Transmembrane H-bonds: {len(tm_interactions)}")
   print(f"Extracellular H-bonds: {len(ec_interactions)}")

Integration with Other Tools
----------------------------

Using with Pandas for Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   from hbat.core.analysis import HBondAnalyzer
   
   analyzer = HBondAnalyzer()
   analyzer.analyze_file("protein.pdb")
   
   # Convert results to pandas DataFrame
   hb_data = []
   for hb in analyzer.hydrogen_bonds:
       hb_data.append({
           'donor_res': hb.donor_residue,
           'acceptor_res': hb.acceptor_residue,
           'distance': hb.distance,
           'angle': hb.angle * 180 / 3.14159,
           'bond_type': hb.bond_type
       })
   
   df = pd.DataFrame(hb_data)
   
   # Perform statistical analysis
   print("Distance Statistics:")
   print(df['distance'].describe())
   
   print("\\nBond Type Distribution:")
   print(df['bond_type'].value_counts())
   
   # Find strongest interactions
   strongest = df.nsmallest(5, 'distance')
   print("\\nStrongest hydrogen bonds:")
   print(strongest)
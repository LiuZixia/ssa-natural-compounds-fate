| Variable | Full Name | Unit | Process Relevance | Remarks | References |
|:---|:---|:---|:---|:---|:---|
| `molecule_id` | Molecule ID | — | Not relevant | Identifier only | — |
| `total_atom_count` | Total Atom Count | Count | Both | Size proxy; can correlate with hydrophobicity/diffusivity and interfacial packing | [10.1016/j.envres.2021.112555](https://doi.org/10.1016/j.envres.2021.112555), [10.5194/acp-22-5223-2022](https://doi.org/10.5194/acp-22-5223-2022) |
| `heavy_atom_count` | Heavy Atom Count | Count | Both | Size proxy; similar to total atom count but excludes H | [10.1016/j.envres.2021.112555](https://doi.org/10.1016/j.envres.2021.112555), [10.5194/acp-22-5223-2022](https://doi.org/10.5194/acp-22-5223-2022) |
| `molecular_weight` | Molecular Weight | g/mol | Both | Size proxy; affects volatility & aerosol partitioning | [10.1016/j.envres.2021.112555](https://doi.org/10.1016/j.envres.2021.112555), [10.5194/acp-22-5223-2022](https://doi.org/10.5194/acp-22-5223-2022) |
| `exact_molecular_weight` | Exact Molecular Weight | g/mol | Both | More precise MW; same process links as MW | [10.1016/j.envres.2021.112555](https://doi.org/10.1016/j.envres.2021.112555), [10.5194/acp-22-5223-2022](https://doi.org/10.5194/acp-22-5223-2022) |
| `molecular_formula` | Molecular Formula | — | Not relevant | Mostly identifier (unless you engineer elemental ratios) | — |
| `alogp` | ALogP | Log Unit | Both | Hydrophobicity/amphiphilicity proxy → surface enrichment and organic film partitioning | [10.1021/acsearthspacechem.1c00080](https://doi.org/10.1021/acsearthspacechem.1c00080), [10.5194/acp-22-5223-2022](https://doi.org/10.5194/acp-22-5223-2022) |
| `topological_polar_surface_area` | Topological Polar Surface Area | Å² | Enrichment-related | Polarity proxy → influences solubility, interfacial adsorption, enrichment | [10.1021/acsearthspacechem.1c00080](https://doi.org/10.1021/acsearthspacechem.1c00080), [10.5194/acp-23-6571-2023](https://doi.org/10.5194/acp-23-6571-2023) |
| `rotatable_bond_count` | Rotatable Bond Count | Count | Enrichment-related | Flexibility; can affect packing at interfaces & phase state | [10.1021/acsearthspacechem.1c00080](https://doi.org/10.1021/acsearthspacechem.1c00080), [10.1039/C7CS00008A](https://doi.org/10.1039/C7CS00008A) |
| `hydrogen_bond_acceptors` | Hydrogen Bond Acceptors | Count | Enrichment-related | Polarity & water affinity; affects SML/SSA partitioning and hygroscopicity | [10.1021/acsearthspacechem.1c00080](https://doi.org/10.1021/acsearthspacechem.1c00080), [10.1039/C7CS00008A](https://doi.org/10.1039/C7CS00008A) |
| `hydrogen_bond_donors` | Hydrogen Bond Donors | Count | Enrichment-related | As above; H-bonding impacts solubility & surface activity | [10.1021/acsearthspacechem.1c00080](https://doi.org/10.1021/acsearthspacechem.1c00080), [10.1039/C7CS00008A](https://doi.org/10.1039/C7CS00008A) |
| `hydrogen_bond_acceptors_lipinski` | Hydrogen Bond Acceptors (Lipinski) | Count | Not relevant | Redundant with HBA; Lipinski framing not marine-transfer-specific | — |
| `hydrogen_bond_donors_lipinski` | Hydrogen Bond Donors (Lipinski) | Count | Not relevant | Redundant with HBD; Lipinski framing not marine-transfer-specific | — |
| `lipinski_rule_of_five_violations` | Lipinski Rule of Five Violations | Count | Not relevant | Drug-likeness heuristic; not mechanistic for SML/SSA transfer | — |
| `aromatic_rings_count` | Aromatic Rings Count | Count | Enrichment-related | Structure proxy; influences hydrophobicity/surface activity and chemical reactivity | [10.1039/C7CS00008A](https://doi.org/10.1039/C7CS00008A), [10.5194/acp-23-6571-2023](https://doi.org/10.5194/acp-23-6571-2023) |
| `qed_drug_likeliness` | QED Drug Likeliness | Score | Not relevant | Drug-likeness score; not mechanistic | — |
| `formal_charge` | Formal Charge | e | Both | Ionization strongly affects solubility & interaction with salts/ions in films | [10.1002/2016GL069070](https://doi.org/10.1002/2016GL069070), [10.1021/acs.est.6b02988](https://doi.org/10.1021/acs.est.6b02988) |
| `fractioncsp3` | Fraction Csp3 | Ratio | Enrichment-related | Structural saturation proxy; relates to lipid-like character | [10.5194/acp-23-6571-2023](https://doi.org/10.5194/acp-23-6571-2023), [10.1039/C7CS00008A](https://doi.org/10.1039/C7CS00008A) |
| `number_of_minimal_rings` | Number of Minimal Rings | Count | Both | Structure proxy; affects hydrophobicity & film behavior | [10.1039/C7CS00008A](https://doi.org/10.1039/C7CS00008A), [10.1021/acsomega.8b01157](https://doi.org/10.1021/acsomega.8b01157) |
| `van_der_walls_volume` | Van der Waals Volume | — | Enrichment-related | Size/packing proxy; can influence interfacial adsorption and enrichment | [10.1016/j.envres.2021.112555](https://doi.org/10.1016/j.envres.2021.112555), [10.1021/acsomega.8b01157](https://doi.org/10.1021/acsomega.8b01157) |
| `contains_sugar` | Contains Sugar | — | Both | Sugary/glyco- motif → water affinity; saccharides can be transferred to SSA | [10.1021/acs.est.6b02988](https://doi.org/10.1021/acs.est.6b02988), [10.1073/pnas.0908905107](https://doi.org/10.1073/pnas.0908905107), [10.1002/2016GL069070](https://doi.org/10.1002/2016GL069070) |
| `contains_ring_sugars` | Contains Ring Sugars | — | Both | As above; cyclic saccharides are common monosaccharide motifs | [10.1021/acs.est.6b02988](https://doi.org/10.1021/acs.est.6b02988), [10.1073/pnas.0908905107](https://doi.org/10.1073/pnas.0908905107), [10.1002/2016GL069070](https://doi.org/10.1002/2016GL069070) |
| `contains_linear_sugars` | Contains Linear Sugars | — | Both | As above; linear saccharide motifs in (poly)saccharides | [10.1021/acs.est.6b02988](https://doi.org/10.1021/acs.est.6b02988), [10.1073/pnas.0908905107](https://doi.org/10.1073/pnas.0908905107), [10.1002/2016GL069070](https://doi.org/10.1002/2016GL069070) |
| `fragments` | Fragments | Count | Not relevant | Substructure feature; not a direct physchem driver (still useful for clustering) | — |
| `fragments_with_sugar` | Fragments with Sugar | Count | Not relevant | Substructure feature; not a direct physchem driver | — |
| `murcko_framework` | Murcko Framework | — | Not relevant | Scaffold identifier; not a direct physchem driver | — |
| `np_likeness` | Natural Product Likeliness | Score | Not relevant | Commonly tied mechanistically to SML/SSA transfer in reviews | — |
| `chemical_class` | ClassyFire Chemical Class | — | Not relevant | Class labels used for grouping/interpretation; not a direct physchem parameter | [10.1186/s13321-016-0174-y](https://doi.org/10.1186/s13321-016-0174-y) |
| `chemical_sub_class` | ClassyFire Chemical Sub-Class | — | Not relevant | Class labels used for grouping/interpretation; not a direct physchem parameter | [10.1186/s13321-016-0174-y](https://doi.org/10.1186/s13321-016-0174-y) |
| `chemical_super_class` | ClassyFire Chemical Super-Class | — | Not relevant | Class labels used for grouping/interpretation; not a direct physchem parameter | [10.1186/s13321-016-0174-y](https://doi.org/10.1186/s13321-016-0174-y) |
| `direct_parent_classification` | Direct Parent Classification | — | Not relevant | Class labels used for grouping/interpretation; not a direct physchem parameter | [10.1186/s13321-016-0174-y](https://doi.org/10.1186/s13321-016-0174-y) |
| `np_classifier_pathway` | NPClassifier Pathway | — | Not relevant | Pathway label; interpretive grouping not a direct physchem driver | [10.1021/acs.jnatprod.1c00399](https://doi.org/10.1021/acs.jnatprod.1c00399) |
| `np_classifier_superclass` | NPClassifier Superclass | — | Not relevant | Pathway label; interpretive grouping not a direct physchem driver | [10.1021/acs.jnatprod.1c00399](https://doi.org/10.1021/acs.jnatprod.1c00399) |
| `np_classifier_class` | NPClassifier Class | — | Not relevant | Pathway label; interpretive grouping not a direct physchem driver | [10.1021/acs.jnatprod.1c00399](https://doi.org/10.1021/acs.jnatprod.1c00399) |
| `np_classifier_is_glycoside` | NPClassifier Is Glycoside | — | Not relevant | Glycoside tag; relevant to sugar transfer but still a label | [10.1021/acs.jnatprod.1c00399](https://doi.org/10.1021/acs.jnatprod.1c00399) |
| `standard_inchi` | Standard InChI | — | Not relevant | Identifier | — |
| `standard_inchi_key` | Standard InChI Key | — | Not relevant | Identifier | — |
| `canonical_smiles` | Canonical SMILES | — | Not relevant | Identifier/structure encoding | — |
| `sugar_free_smiles` | Sugar Free SMILES | — | Not relevant | Identifier/structure encoding | — |
| `identifier` | Identifier | — | Not relevant | Identifier | — |
| `name` | Chemical Name | — | Not relevant | Identifier | — |
| `cas` | CAS Registry Number | — | Not relevant | Identifier | — |
| `synonyms` | Synonyms | — | Not relevant | Identifier list | — |
| `iupac_name` | IUPAC Name | — | Not relevant | Identifier | — |
| `murko_framework` | Murko Framework | — | Not relevant | Scaffold identifier; not a direct physchem driver | — |
| `structural_comments` | Structural Comments | — | Not relevant | Free text | — |
| `name_trust_level` | Name Trust Level | Score | Not relevant | Metadata | — |
| `annotation_level` | Annotation Level | Score | Not relevant | Metadata | — |
| `parent_id` | Parent ID | — | Not relevant | Metadata | — |
| `variants_count` | Variants Count | Count | Not relevant | Metadata | — |
| `ticker` | Ticker | — | Not relevant | Metadata | — |
| `status` | Status | — | Not relevant | Metadata | — |
| `active` | Active | — | Not relevant | Metadata | — |
| `has_variants` | Has Variants | — | Not relevant | Metadata | — |
| `has_stereo` | Has Stereo | — | Not relevant | Metadata | — |
| `is_tautomer` | Is Tautomer | — | Not relevant | Metadata | — |
| `is_parent` | Is Parent | — | Not relevant | Metadata | — |
| `is_placeholder` | Is Placeholder | — | Not relevant | Metadata | — |
| `comment` | Comment | — | Not relevant | Metadata | — |
| `organism_count` | Organism Count | Count | Not relevant | Sampling/DB metadata | — |
| `geo_count` | Geographic Count | Count | Not relevant | Sampling/DB metadata | — |
| `citation_count` | Citation Count | Count | Not relevant | Sampling/DB metadata | — |
| `collection_count` | Collection Count | Count | Not relevant | Sampling/DB metadata | — |
| `synonym_count` | Synonym Count | Count | Not relevant | Sampling/DB metadata | — |
| `is_duplicate` | Is Duplicate | — | Not relevant | QC metadata | — |
| `is_marine` | Is Marine | — | Not relevant | Metadata flag (not a physchem driver) | — |
| `aop_oh_rate_constant` | AOP OH Rate Constant | cm³/molecule | Not relevant | Atmospheric oxidation proxy (gas-phase), not marine transfer | [10.1021/cr0206420](https://doi.org/10.1021/cr0206420), [US EPA EPI Suite (AOPWIN)](https://www.epa.gov/tsca-screening-tools/epi-suitetm-estimation-program-interface) |
| `aop_oh_half_life_hours` | AOP OH Half Life (Hours) | h | Not relevant | Derived from OH rate constant; atmospheric persistence | [10.1021/cr0206420](https://doi.org/10.1021/cr0206420), [US EPA EPI Suite (AOPWIN)](https://www.epa.gov/tsca-screening-tools/epi-suitetm-estimation-program-interface) |
| `aop_ozone_rate_constant` | AOP Ozone Rate Constant | cm³/molecule | Not relevant | Atmospheric oxidation proxy (gas-phase) | [10.1021/cr0206420](https://doi.org/10.1021/cr0206420) |
| `aop_ozone_half_life_hours` | AOP Ozone Half Life (Hours) | h | Not relevant | Derived from ozone rate constant; atmospheric persistence | [10.1021/cr0206420](https://doi.org/10.1021/cr0206420) |
| `bp_est` | Boiling Point (Estimated) | °C | Aerosolization-related | Volatility proxy; impacts partitioning & phase | [US EPA EPI Suite (MPBPWIN)](https://www.epa.gov/tsca-screening-tools/epi-suitetm-estimation-program-interface), [10.1080/00986448708960487](https://doi.org/10.1080/00986448708960487) |
| `mp_est` | Melting Point (Estimated) | °C | Aerosolization-related | Phase state proxy (solid vs liquid) | [US EPA EPI Suite (MPBPWIN)](https://www.epa.gov/tsca-screening-tools/epi-suitetm-estimation-program-interface), [10.1039/C1CP22617G](https://doi.org/10.1039/C1CP22617G) |
| `vp_est` | Vapor Pressure (Estimated) | mm Hg | Aerosolization-related | Volatility proxy; key driver of gas/particle partitioning | [10.1016/j.fluid.2008.04.020](https://doi.org/10.1016/j.fluid.2008.04.020), [10.1016/1352-2310(94)90093-0](https://doi.org/10.1016/1352-2310(94)90093-0) |

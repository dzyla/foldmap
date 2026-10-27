# Reading a structure file

## Scope

- `fetch.fetch`: a PDB id, a UniProt or AlphaFold DB id, a local path or uploaded bytes, to a file on disk (cached).
- `io.read_model`: the file read by gemmi, first model only, one conformer per residue, and the assembly chosen
  (`asu`, an assembly id, or `auto`).
- `io.load_backbone`: the protein backbone (N, CA, C, O) of each chain, as `model.Backbone`.
- `io.load_nucleic`: the nucleic-acid strands, their base pairs and the protein residues touching them
  (see [nucleic-acids.md](nucleic-acids.md)).
- `scope.scope_note`: what an automatic choice left out of the figure, said in the notes.

## Assumes

- The file is PDB or mmCIF that gemmi reads. Only the first model is used.
- A residue in alternative locations is read in one conformer: the model's first altloc letter (`io._one_conformer`).
- A protein residue needs N, CA, C and O. A missing carbonyl O is rebuilt in the peptide plane
  (`io._rebuild_oxygens`). A residue missing N, CA or C is skipped and breaks the chain.
- Two residues are bonded when C to N is at most `io._PEPTIDE_BOND_MAX` (2.5 Å), as DSSP takes it.
- A CA-only or P-only model cannot be drawn. foldmap refuses it with a one-line message.
- Under `--assembly auto`, the first assembly is built when it adds copies and stays within `io.AUTO_MAX` (24)
  chains, or when it is one of several copies of a complex in the asymmetric unit (`io._one_complex_of_several`).
  An assembly that would leave out the DNA or RNA is not chosen.
- Two nucleic-acid chains held by base pairs (at least `io.PAIRED_MIN` pairs within `io._PAIR_DIST`) are one object
  (`io._paired`).

## Guarantees

- First model only: `tests/test_io.py::test_only_first_model_used`.
- One conformer per residue: `tests/test_io.py::test_a_residue_is_read_in_one_conformer`.
- PDB and mmCIF give the same backbone: `tests/test_io.py::test_pdb_format_matches_mmcif`.
- A residue missing a backbone atom breaks the chain:
  `tests/test_io.py::test_residue_missing_backbone_atom_is_skipped_and_breaks_chain`.
- Missing carbonyl oxygens are rebuilt: `tests/test_robustness.py::test_missing_carbonyl_oxygens_are_rebuilt`.
- A CA-only model is refused with its reason: `tests/test_robustness.py::test_ca_only_model_is_explained`.
- An assembly that leaves out the DNA keeps the asymmetric unit:
  `tests/test_assembly.py::test_dna_outside_assembly_one_keeps_the_asymmetric_unit`.
- A lattice contact of two base pairs does not pair two chains:
  `tests/test_assembly.py::test_a_lattice_contact_of_two_base_pairs_is_no_pairing`.
- A rigid motion, a chain renaming, shuffled atoms and the PDB writer leave the structure facts exactly equal on five
  fixtures (1UBQ, 1TIM, 1C3W, 1LMB, 1EHZ): `tests/test_metamorphic.py::test_structure_facts_are_exact`, and the
  transforms table of `docs/validation/metamorphic/report.md`.
- Coordinate noise of 0.05 and 0.1 Å changes the facts by a bounded rate:
  `tests/test_metamorphic.py::test_noise_changes_are_bounded`.
- What `auto` leaves out is said: `tests/test_scope.py::test_scope_notes_name_the_settings_that_draw_everything`.

## Fails

- A CA-only or P-only entry is refused, never drawn. The campaigns meet such entries and replace them in the sample
  (`docs/validation/fourth-sample-protocol.md`).
- 1C3W written as PDB: the BIOMT operators have 6 decimals, so the assembly moves by up to 0.0024 page units
  (`tests/test_metamorphic.py::test_pdb_assembly_operators_are_rounded`).
- 5K7Z: DNA bound in two orientations. Before `io._one_conformer`, atoms were taken from both conformers and the strand
  was drawn across the page. Fixed; held by the one-conformer test above.
- Which conformer is drawn is a rule (the first altloc), not a choice by occupancy. Whether another conformer would
  read better is not measured.
- `io.PAIRED_MIN` and `io._PAIR_DIST` were set for 4PLX and 3OWW. No view of the third random sample reads them
  (`docs/validation/sensitivity/report.md`, "not exercised"). Only tests hold them (`docs/validation/prune/report.md`).
- Protein-nucleic acid hydrogen bonds (`hbonds.find_hbonds`) class an atom at a mid-chain gap as terminal, so it may
  donate. How often this happens is not measured.

## Seen by

- `A1`: a crash, a traceback or a Python warning in any view (campaign only).
- `A1c`: a CA-only or P-only entry must exit non-zero with foldmap's message (campaign only).
- `A9`: chains or nucleotides left out without a note.
- `A10`: the secondary structure against mkdssp (campaign only). A misread backbone would show there; that use is not
  measured.
- No check sees a plausible but wrong assembly choice: the notes say what `auto` drew, and nothing judges whether it
  was the biological one (not measured).
- No check sees which altloc was drawn.

## Evidence

- [metamorphic/report.md](../validation/metamorphic/report.md): transforms, noise and the 1C3W PDB gap.
- [sensitivity/report.md](../validation/sensitivity/report.md): `io:PAIRED_MIN` and `io:_PAIR_DIST` not exercised.
- [prune/report.md](../validation/prune/report.md): both kept by the test suite.
- [inventory/inventory.md](../validation/inventory/inventory.md): the reading constants and their origins.

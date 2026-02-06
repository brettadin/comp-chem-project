# HybridQSE - Claude Build Instructions (Consolidated)
Version 1.0 - February 05, 2026

**Source tags used inline**
- **[P-Docs]**: repo docs (START_HERE / PROJECT_STATUS / DEVELOPMENT_GUIDE / ACHIEVEMENT_SUMMARY / PROJECT_TREE)
- **[BP]**: *Blueprint for a Hybrid Multi-Scale Modeling Engine* (PDF)
- **[GQ]**: *Designing a Generalized Quantum Simulation Engine for Atoms, Molecules, and Subatomic Particles* (PDF)
- **[ExtDocs]**: official external library docs (PySCF/Psi4/etc.)

---

## 0) What exists now (do not break) [P-Docs]
- Working `Molecule/Atom` model + XYZ IO + nuclear repulsion
- Basis-set + integral interface (PySCF-backed)
- Restricted HF (RHF) SCF solver returning energy + orbitals + dipole + HOMO/LUMO gap
- Validation framework with reference DB and auto comparison
- Example: `examples/01_h2_hartree_fock.py`
- Unit tests for core components

Read order:
1. `START_HERE.md` [P-Docs]
2. `ACHIEVEMENT_SUMMARY.md` [P-Docs]
3. `PROJECT_STATUS.md` [P-Docs]
4. `DEVELOPMENT_GUIDE.md` [P-Docs]
5. `PROJECT_TREE.txt` [P-Docs]

## 1) Non-negotiables [BP, P-Docs]
- **Validation-first:** every new method must match an external reference on at least one small system before expanding features. [BP]
- **Reproducibility:** deterministic with fixed seed; log versions + inputs + tolerances. [BP]
- **Modularity:** add methods via plugin registry / entry points; avoid core edits. [BP]
- **Units + provenance everywhere:** geometry units, basis/functional IDs, tolerances, backend versions embedded in results. [BP]

## 2) Scope boundaries (v1.x) [P-Docs, GQ]
Primary: molecular electronic structure (Gaussian basis) via PySCF backend. [P-Docs]  
Secondary: BO-MD (classical nuclei on quantum PES). [P-Docs]  
Tertiary: vibrations + IR spectrum simulation. [GQ]  
Deferred: full lattice QCD (separate module; start with toy validation). [GQ]

### IR -> geometry (truthful framing) [GQ]
IR encodes **normal modes** (mass-weighted Hessian) around an equilibrium geometry. [GQ]  
But "IR -> unique 3D geometry" is generally ill-posed. Implement **model-based inversion**: propose geometries -> compute spectra -> fit/optimize candidates with priors + uncertainty.

## 3) Milestones (do in this order)

### M1 - DFT wrapper (RKS then UKS) [P-Docs, ExtDocs]
Deliver: `hybrid_qse/methods/dft.py` (wrapper around PySCF DFT). [P-Docs, ExtDocs: PySCF]
- Support common functionals (B3LYP, PBE, PBE0, etc.) by name/alias. [P-Docs]
- Return energy, converged, mo_coeff, mo_energy, (optional density handle). [P-Docs]
Acceptance:
- DFT(H2O, B3LYP/6-31G) converges.
- Energies match PySCF within tight tolerance.
- Unit tests + add reference points in BenchmarkDatabase.

### M2 - Geometry optimization [P-Docs, BP]
Deliver: `hybrid_qse/methods/optimizer.py`
- Start with SciPy BFGS using numerical gradients; later move to analytic gradients. [P-Docs, ExtDocs: SciPy/PySCF]
- Add constraints: freeze atoms, fix bonds/angles. [BP]
Acceptance:
- Perturbed water converges; energy decreases.
- Frequency analysis after optimization shows no imaginary modes. [GQ]

### M3 - Vibrations + IR spectrum [GQ, P-Docs]
Deliver: `hybrid_qse/properties/vibrations.py`
- Numerical Hessian -> mass-weighted Hessian -> normal modes/frequencies (cm^-1). [GQ, P-Docs]
- IR intensities from dipole derivatives; broaden peaks to form simulated spectrum. [GQ]
Acceptance:
- H2O frequencies physically reasonable; peaks qualitatively align with references.
- Tests: Hessian symmetry, unit conversions, stable diagonalization.

### M4 - Born-Oppenheimer MD (toy) [P-Docs, BP]
Deliver: `hybrid_qse/dynamics/md.py`
- Velocity Verlet; forces via numerical gradient of energy. [P-Docs]
- Track energy drift; only add thermostats after NVE is stable. [BP]
Acceptance:
- Short NVE run on H2 shows bounded energy drift; reproducible with fixed seed.

## 4) Core representations (keep extensible) [BP, GQ]
Define minimal interfaces:
- `System` (Molecule now; later PeriodicSystem / FieldSystem) [GQ]
- `Method.run(system, settings) -> Result` [BP]
- `Result` includes values + convergence + provenance dict [BP]
- `Task` objects for energy/gradient/optimize/vib/md workflows [BP]

Implement a plugin registry / entry points. [BP]

## 5) "Shortcuts" that preserve truth (and quantify error) [BP]
Physics-first:
- symmetry, RI/density fitting, Cholesky ERI, pseudopotentials (later), adaptive basis refinement. [BP]

ML (only with uncertainty):
- delta-learning corrections, ML potentials with uncertainty-triggered fallback/active learning, optional neural-VMC later. [BP]

## 6) Validation playbook [BP, P-Docs, GQ]
- Validate vs PySCF at identical settings first. [P-Docs]
- Add literature/experimental refs next. [BP]
- For spectra: validate normal modes + intensities vs curated reference spectra; document scaling. [GQ]
- Single command runs full suite and outputs pass/fail. [BP]

Minimum structured log schema (JSON):
```json
{
  "engine_version": "...",
  "backend": {"pyscf": "...", "numpy": "...", "scipy": "..."},
  "system": {"formula":"H2O","charge":0,"multiplicity":1,"units":"angstrom"},
  "method": {"name":"DFT","functional":"B3LYP","basis":"6-31G"},
  "settings": {"scf_tol":1e-9,"max_iter":50},
  "results": {"energy_hartree":-76.4,"converged":true},
  "provenance": {"input_hash":"...","seed":1234}
}
```

## 7) UX + docs (community-ready) [BP]
- "5-minute path": pip install -> run example -> validated output. [P-Docs]
- How-to pages per milestone with expected outputs. [BP]
- Actionable errors (SCF convergence help, etc.). [BP]
- Ship notebooks + CLI (`hybridqse run config.yaml`). [BP]

## 8) Official docs (use these as ground truth) [ExtDocs]
- PySCF docs: https://pyscf.org
- Psi4 docs: https://psicode.org/psi4manual/master/
- OpenFermion: https://quantumai.google/openfermion
- ASE: https://wiki.fysik.dtu.dk/ase/
- LAMMPS: https://www.lammps.org
- Quantum ESPRESSO: https://www.quantum-espresso.org
- JAX: https://jax.readthedocs.io
- QUDA (lattice QCD): https://lattice.github.io/quda/

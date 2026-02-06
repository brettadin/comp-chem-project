# Quick Start Development Guide

## 🚀 **RIGHT NOW - What You Can Do**

### Option 1: Run Existing Code (Recommended First Step)

Install dependencies and run the example:

```bash
# Install required packages
pip install numpy scipy pyscf matplotlib pytest

# Run the H₂ Hartree-Fock example
cd /mnt/project
python examples/01_h2_hartree_fock.py
```

**Expected output:**
- H₂ geometry display
- SCF convergence in ~5 iterations
- Final energy: ~-1.117 Hartree
- Validation against PySCF and literature
- Molecular properties (dipole, HOMO-LUMO gap)

### Option 2: Run Tests

```bash
cd /mnt/project
pytest tests/ -v
```

### Option 3: Interactive Exploration

```python
from hybrid_qse import Molecule
from hybrid_qse.methods import HartreeFock
from hybrid_qse.validation import BenchmarkDatabase

# Try different molecules
he = Molecule.helium()
h2 = Molecule.h2(bond_length=0.75)  # Slightly stretched
water = Molecule.water()

# Run calculations
hf = HartreeFock(water, basis_name='sto-3g')
results = hf.run()

# Check against references
db = BenchmarkDatabase()
refs = db.list_references_for("H2O")
print(refs)
```

---

## 📝 **NEXT CODING SESSION - Add DFT** (2-3 hours)

Create `hybrid_qse/methods/dft.py`:

```python
"""
Density Functional Theory implementation.

Wrapper around PySCF's DFT module with consistent interface.
"""

import numpy as np
try:
    from pyscf import dft, gto, scf
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False

from ..core.molecule import Molecule


class DFT:
    """
    Density Functional Theory wrapper.
    
    Supports common functionals:
    - LDA: SVWN
    - GGA: PBE, BLYP
    - Hybrid: B3LYP, PBE0
    """
    
    def __init__(
        self,
        molecule: Molecule,
        functional: str = 'B3LYP',
        basis_name: str = 'sto-3g'
    ):
        self.molecule = molecule
        self.functional = functional
        self.basis_name = basis_name
        
        # Results
        self.energy = None
        self.converged = False
    
    def run(self, verbose=True):
        """Run DFT calculation."""
        if not PYSCF_AVAILABLE:
            raise RuntimeError("PySCF required for DFT")
        
        # Build PySCF molecule
        atom_string = []
        for atom in self.molecule.atoms:
            x, y, z = atom.position
            atom_string.append([atom.symbol, (x, y, z)])
        
        mol = gto.M(
            atom=atom_string,
            basis=self.basis_name,
            charge=self.molecule.charge,
            spin=self.molecule.multiplicity - 1,
            unit='angstrom'
        )
        
        # Run DFT
        mf = dft.RKS(mol)
        mf.xc = self.functional
        mf.verbose = 4 if verbose else 0
        
        self.energy = mf.kernel()
        self.converged = mf.converged
        
        return {
            'energy_total': self.energy,
            'energy_nuclear': self.molecule.nuclear_repulsion_energy(),
            'converged': self.converged,
            'mo_coefficients': mf.mo_coeff,
            'orbital_energies': mf.mo_energy
        }
```

**Test it:**
```python
from hybrid_qse import Molecule
from hybrid_qse.methods import DFT

water = Molecule.water()

# Compare functionals
for func in ['B3LYP', 'PBE', 'PBE0']:
    dft = DFT(water, functional=func, basis_name='6-31g')
    results = dft.run(verbose=False)
    print(f"{func:10s}: {results['energy_total']:12.6f} Hartree")
```

---

## 🎯 **SESSION 2 - Geometry Optimization** (3-4 hours)

Create `hybrid_qse/methods/optimizer.py`:

```python
"""Geometry optimization."""

from scipy.optimize import minimize
import numpy as np

class GeometryOptimizer:
    """Optimize molecular geometry to find minimum energy."""
    
    def __init__(self, molecule, method, basis='sto-3g'):
        self.molecule = molecule
        self.method = method  # 'HF' or 'DFT'
        self.basis = basis
    
    def _energy_function(self, coords_flat):
        """Compute energy for given geometry."""
        # Reshape coords and update molecule
        coords = coords_flat.reshape(-1, 3)
        for i, atom in enumerate(self.molecule.atoms):
            atom.position = coords[i]
        
        # Calculate energy
        if self.method == 'HF':
            from .hf import HartreeFock
            calc = HartreeFock(self.molecule, self.basis)
        elif self.method.startswith('DFT'):
            from .dft import DFT
            func = self.method.split(':')[1]  # e.g., 'DFT:B3LYP'
            calc = DFT(self.molecule, func, self.basis)
        
        results = calc.run(verbose=False)
        return results['energy_total']
    
    def optimize(self):
        """Run optimization."""
        # Initial coordinates
        coords_init = np.array([atom.position for atom in self.molecule.atoms])
        
        # Optimize
        result = minimize(
            self._energy_function,
            coords_init.flatten(),
            method='BFGS',
            options={'disp': True}
        )
        
        return result
```

**Test on water:**
```python
# Start with approximate geometry
water = Molecule.water()
# Perturb slightly
water.atoms[1].position += np.array([0.1, 0.0, 0.0])

opt = GeometryOptimizer(water, method='HF', basis='sto-3g')
result = opt.optimize()

print(f"Optimized geometry:")
print(water)
```

---

## 🔬 **SESSION 3 - Spectral Calculations** (4-5 hours)

### Vibrational Frequencies

Create `hybrid_qse/properties/vibrations.py`:

```python
"""Vibrational frequency calculations."""

import numpy as np

def compute_hessian_numerical(molecule, method, basis, delta=0.001):
    """
    Compute Hessian by finite differences.
    
    H_ij = ∂²E/∂x_i∂x_j
    """
    n_atoms = len(molecule.atoms)
    n_coords = 3 * n_atoms
    
    hessian = np.zeros((n_coords, n_coords))
    
    # Get energy function
    def energy_func(coords):
        for i, atom in enumerate(molecule.atoms):
            atom.position = coords[i]
        
        if method == 'HF':
            from ..methods.hf import HartreeFock
            calc = HartreeFock(molecule, basis)
        
        return calc.run(verbose=False)['energy_total']
    
    # Central difference for each element
    coords_0 = np.array([atom.position for atom in molecule.atoms]).flatten()
    
    for i in range(n_coords):
        for j in range(i, n_coords):
            # Four-point formula for ∂²E/∂x_i∂x_j
            coords = coords_0.copy()
            
            # E(x_i+δ, x_j+δ)
            coords[i] += delta
            coords[j] += delta
            E_pp = energy_func(coords.reshape(-1, 3))
            
            # ... (similar for other combinations)
            
            hessian[i, j] = (E_pp - E_pm - E_mp + E_mm) / (4 * delta**2)
            hessian[j, i] = hessian[i, j]
    
    return hessian

def compute_frequencies(hessian, masses):
    """
    Compute vibrational frequencies from mass-weighted Hessian.
    
    Returns frequencies in cm⁻¹
    """
    # Mass-weight the Hessian
    mass_matrix = np.repeat(masses, 3)
    M_sqrt = np.sqrt(np.outer(mass_matrix, mass_matrix))
    H_mw = hessian / M_sqrt
    
    # Diagonalize
    eigvals, eigvecs = np.linalg.eigh(H_mw)
    
    # Convert to frequencies (cm⁻¹)
    # ω = sqrt(λ) * conversion_factor
    freq_cm = np.sqrt(np.abs(eigvals)) * 5140.48  # a.u. to cm⁻¹
    
    return freq_cm, eigvecs
```

---

## 🧪 **SESSION 4 - Real-Time Dynamics** (5-6 hours)

### Born-Oppenheimer Molecular Dynamics

Create `hybrid_qse/dynamics/md.py`:

```python
"""Molecular dynamics on quantum potential energy surface."""

class BornOppenheimerMD:
    """
    Classical nuclear dynamics on quantum PES.
    
    Uses velocity Verlet integration.
    """
    
    def __init__(self, molecule, method, basis, dt=1.0):
        self.molecule = molecule
        self.method = method  # Method for energy/forces
        self.basis = basis
        self.dt = dt  # Time step in femtoseconds
        
        self.positions = []
        self.velocities = []
        self.energies = []
    
    def _compute_forces(self):
        """Compute forces by numerical gradient."""
        forces = np.zeros((len(self.molecule.atoms), 3))
        delta = 0.001
        
        for i, atom in enumerate(self.molecule.atoms):
            for j in range(3):
                # Forward
                atom.position[j] += delta
                E_plus = self._get_energy()
                
                # Backward
                atom.position[j] -= 2*delta
                E_minus = self._get_energy()
                
                # Restore
                atom.position[j] += delta
                
                # Force = -dE/dx
                forces[i, j] = -(E_plus - E_minus) / (2 * delta)
        
        return forces
    
    def run(self, n_steps, temperature=300):
        """Run MD simulation."""
        # Initialize velocities (Maxwell-Boltzmann)
        self._initialize_velocities(temperature)
        
        for step in range(n_steps):
            # Velocity Verlet
            forces = self._compute_forces()
            
            # Update positions
            for i, atom in enumerate(self.molecule.atoms):
                atom.position += self.v[i] * self.dt + 0.5 * forces[i] * self.dt**2
            
            # Update velocities
            forces_new = self._compute_forces()
            self.v += 0.5 * (forces + forces_new) * self.dt
            
            # Store trajectory
            self.positions.append([a.position.copy() for a in self.molecule.atoms])
            self.energies.append(self._get_energy())
```

---

## 🤖 **SESSION 5+ - Machine Learning** (2-3 weeks)

Following Blueprint Section 2 (Physics-Informed ML):

1. **Neural Network Potentials**
   - Train on CCSD(T) data
   - SchNet or ANI architecture
   - Symmetry-preserving layers

2. **Active Learning**
   - Uncertainty estimation
   - On-the-fly training
   - Adaptive sampling

3. **Neural Wavefunctions**
   - FermiNet-style architecture
   - Variational Monte Carlo
   - AutodiffPyTorch integration

---

## 📚 **REFERENCE WHILE CODING**

### Key Equations to Implement:

**HF Fock Matrix:**
```
F_μν = H_μν + Σ_λσ P_λσ [(μν|λσ) - 0.5(μλ|νσ)]
```

**DFT Exchange-Correlation:**
```
E_xc = ∫ ε_xc(ρ(r)) ρ(r) dr
```

**Vibrational Frequencies:**
```
ω = sqrt(eigenvalues of mass-weighted Hessian)
```

**MD Force:**
```
F_i = -∇_i E_quantum
```

### Code Style:
- Docstrings for everything
- Type hints where helpful
- Validation against references
- Tests for each new feature

---

## 🎯 **SUCCESS METRICS**

After each session, validate:

✅ **Code runs without errors**
✅ **Results match PySCF reference**
✅ **Literature values reproduced**
✅ **Tests pass**
✅ **Documentation updated**

---

## 🔥 **LET'S BUILD THIS!**

You have:
- ✅ Solid foundation
- ✅ Clear roadmap
- ✅ Working examples
- ✅ Validation framework

Each session builds on the last. Each feature is **immediately usable** and **scientifically validated**.

**Next command to run:**
```bash
python examples/01_h2_hartree_fock.py
```

Then start building the next feature! 🚀⚛️

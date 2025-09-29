# MPO–CALR Project

This repository contains a **publishable-quality pipeline** for investigating the interaction between **Myeloperoxidase (MPO)** and **Calreticulin (CALR)** using:
- **HADDOCK3** for protein–protein docking
- **GROMACS** for molecular dynamics simulations
- **MM/GBSA** and **PMF (umbrella sampling)** for binding free energy analysis

The project emphasizes **reproducibility**, **rigor**, and **biological realism**:
- Multiple MPO glycoforms
- Soluble and membrane-proximal CALR
- Protonation states, cofactors, and controls
- Replicates, multiple docking clusters
- Transparent configs, scripts, and documentation

---

## 🔧 Setup

### 1. Clone the repo
```bash
git clone git@github.com:aljop02/ppi-mpo-calr.git
cd ppi-mpo-calr

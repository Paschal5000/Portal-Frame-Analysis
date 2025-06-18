🏗 2D Portal Frame Structural Analysis – Stiffness Matrix Method

This Python program performs structural analysis of 2D portal frames using the stiffness matrix method. It’s built to support engineers, researchers, and students in analyzing statically indeterminate structures efficiently and accurately.

---

✨ Features
- Stiffness matrix assembly for 2D frames
- Supports:
  - Point loads
  - Uniformly distributed loads (UDLs)
  - Applied moments
- User-defined:
  - Node coordinates and element connectivity
  - Material properties (E, A)
  - Boundary conditions
- Global stiffness matrix assembly with transformation
- Solves for nodal displacements and member forces
- Outputs key matrices and results in a clean, report-style format using tabulate

---

🛠 Applications
- Educational demonstrations
- Preliminary structural design
- Learning and teaching matrix-based analysis

---

🐍 Requirements
- Python 3.x
- numpy
- tabulate

Install with:
bash
pip install numpy tabulate


---

📁 How to Use
1. Clone or download the repo.
2. Define frame geometry, material properties, loads, and supports in the input section.
3. Run the script.
4. View detailed output including: 
  - Global stiffness matrix
   - Displacements
   - Reactions and member forces.
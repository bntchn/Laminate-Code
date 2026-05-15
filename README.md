# Classical Lamination Theory MATLAB Tool

A MATLAB class-based implementation of Classical Lamination Theory (CLT) for analyzing laminated composite materials. This program computes laminate stiffness properties, transformed ply stiffness matrices, ABD matrices, and effective laminate engineering constants for arbitrary stacking sequences.

This tool was developed to perform a sensitivity analysis of composite laminates to variations in fiber orientation. The analysis consists of a Monte-Carlo simulation using the CLT tool.

---

# Features

- Define multiple composite material systems
- Arbitrary laminate stacking sequences
- Computes:
  - Ply stiffness matrices \(Q\)
  - Transformed stiffness matrices \(\bar{Q}\)
  - Transformation matrices \(T\)
  - Laminate ABD matrix
  - Effective laminate properties:
    - \(E_1\)
    - \(E_2\)
    - \(G_{12}\)
    - \(\nu_{12}\)
- Automatic ply z-coordinate generation
- Generates formatted output reports (`Results.txt`)
- Supports multiple material types in a single laminate

---

# File Structure

```text
laminateCLT_700.m    % Main CLT MATLAB class
RunMe.m              % Script used to generate Monte-Carlo data and post-processing
```

---

# Inputs

<details>
<summary>Material List</summary>

`materialList` is an \(n \times 5\) matrix:

```matlab
[E1, E2, G12, nu12, thickness]
```

Each row corresponds to a unique fiber material.

Example:

```matlab
materialList = [
    19e6  2e6  1e6  0.3  0.005
];
```

</details>

---

<details>
<summary>Ply Orientations</summary>

Ply angles in degrees:

```matlab
orientations = [0 45 -45 90];
```

</details>

---

<details>
<summary>Fiber Assignments</summary>

Defines which material row is used for each ply:

```matlab
fibers = [1 1 1 1];
```

</details>

---

# Example Usage

<details>
<summary>Click to expand example</summary>

```matlab
materialList = [
    19e6  2e6  1e6  0.3  0.005
];

orientations = [0 45 -45 90];
fibers = [1 1 1 1];

% Create laminate object
laminate = laminateCLT_700(materialList,orientations,fibers);

% Compute laminate properties
laminate = laminate.properties();

% Generate results file
laminate = laminate.ResultsFile();
```

</details>

---

# Outputs

The program computes:

- Ply stiffness matrices
- Rotated stiffness matrices
- Transformation matrices
- Ply z-locations
- ABD matrix
- Effective laminate properties

A formatted report is automatically generated:

```text
Results.txt
```

---

# Current Capabilities

- Symmetric and asymmetric laminates
- Multiple material systems
- Arbitrary stacking sequences
- Automated laminate stiffness calculations

---

# Planned Improvements

- Failure criteria implementation
  - Maximum stress
  - Maximum strain
- Thermal loading effects
- Progressive ply failure

---

# References

- Jones, R. M. *Mechanics of Composite Materials*
- Classical Lamination Theory (CLT)

---

# Acknowledgments

Portions of this project documentation and development were completed with assistance from ChatGPT, an AI language model developed by OpenAI.

---

# License

This project is open-source and intended for educational and research purposes.
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

---

#  Acknowledgement

```text
Portions of this project were developed with assistance from ChatGPT, an AI language model developed by OpenAI.
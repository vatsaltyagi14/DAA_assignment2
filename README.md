# CS F364 - Design and Analysis of Algorithms

## Assignment 2: Densest Subgraph Discovery

---


## Dataset Preparation for Algorithm 1 (Exact Algorithm)

- Uncompress the datasets if not already done and place them in the same folder as the code files.
- Remove the first few comment lines before giving them as input to the code files.
- Make sure the dataset filenames exactly match the ones specified in the code.
- The datasets are :
  - as-caida.txt 
  - as20000102.txt 
  - CA-HepTh.txt 
  - netscience.txt 

---


## Dataset Preparation for Algorithm 4 (CoreExact Algorithm)

- Uncompress the datasets if not already done and place them in the same folder as the code files.
- Remove the first few comment lines before giving them as input to the code files.
- Make sure the dataset filenames exactly match the ones specified in the code.
- In the first line write <H,Number of vertices,Number of edges>
- The datasets are :
  - as-caida.txt 
  - as20000102.txt 
  - CA-HepTh.txt 
  - netscience.txt 

---

## Execution Instructions

- Before executing the code, change the dataset path if necessary.
- For compilation g++ -O3 -std=c++17 exact.cpp -o output
- For execution ./output (for Algorithm 1)
- For execution ./output input.txt (for Algorithm 4)

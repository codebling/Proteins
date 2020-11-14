# Proteins
This is the working implementation of the Branch and Bound HP lattice protein folding algorithm developed by [Chen and Huang](https://www.brown.edu/Research/Istrail_Lab/_proFolding/papers/2005/bran-06.pdf). The general gist of all HP protein folding algorithms is:
1. Given a linear sequence of amino acids, first convert each acid to either a Hydrophobic (H) or Hydrophilic (P) type acid. Since these proteins are in aqueous solutions, the protein chain will want to fold such that the H acids are as far away from the exterior environment as possible (ie the H acids will tend to clump towards the center of the folded protein).
2. This can be equivalently thought of as maximizing the number of H acids that are adjacent to each other *excluding* consecutive H acids along the chain.
3. The negative number of these adjacent pairs can be thought of as the total energy state of the protein.
4. The allowable locations for each of these acids are bound to the vertices of a desired grid, typically square or triangular. 
5. No two acids can be on the same vertex at the same time.

The algorithm proposed by Chen and Huang implements a Branch-and-Bound growing search technique. Each acid in the chain is placed consecutively, and the energy state is compared to the historic average and historic minimum of all potential chain conformations of the same length. Their paper shows a 2D square lattice implementation, and uses a depth-first search technique. A big drawback of the depth-first implementation, is that not all conformations of length *k* are evaluated for their low-energy potential at the same time--the potential for a given partial conformation can not be fully evaluated until all possible branches from this partial conformation is completed. 

This implementation uses a *breadth*-first search technique, which evaluates *layers* of the protein chain: in other words, it tests all possible conformations of length k, evaluates each energy potential, and then begins a search for branches of length k+1 only using high-potential partial chains. This has sped up the implementation significantly. Another implementation is the 2D triangular lattice, as well as 3-Dimensional implementations of the square and triangular lattice shapes. Here are some examples on a 200-acid long chain.

![All possible lattice + dimension types on a 200-acid long chain.](Figures/chain200.png)

Another key implementation is an additional threshold culling routine. 
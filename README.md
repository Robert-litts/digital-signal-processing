# Digital Signal Processing Coursework
Coursework based on *Digital Signal Processing, 4th edition* by John G. Proakis and Dimitris G. Manolakis

## Term Project: *Distributed Compressed Sensing of Correlated Sparse Signals*
- Compressed sensing (CS) deals with reconstruction of a “sparse signal” using only a few samples, much lower than required by the Nyquist rate
- When operating in a network model, previous work focused on multiple nodes gathering data, but signal reconstruction took place at central node
- However, distributed spectrum estimation, multiple-sensor image or sound capturing require the ability to reconstruct signals at individual nodes
- Led to need for development of Distributed Compressed Sensing (DCS), network of users each individually reconstruct sparse signals
- Algorithm based off the Greedy Pursuit principle, balance between complexity and quality
- Goal is to build a common, or joint support set by merging the individual support sets of correlated sparse signals at individual nodes
- Support set information shared over network
- DIPP Algorithm: distributed parallel pursuit with side information
- Data fusion from individual nodes using democratic voting
- Improves signal reconstruction at each node as network size increases
- Perfect signal reconstruction in theory [1]

# References
1. D. Sundman, S. Chatterjee, and M. Skoglund, “Design and Analysis of a Greedy Pursuit for Distributed Compressed Sensing,” IEEE Transactions on Signal Processing, vol. 64, no. 11, pp. 2803–2818, 2016.
2. l1-Magic. [Online]. Available: https://statweb.stanford.edu/~candes/software/l1magic/index.html#links. [Accessed: 25-Apr-2020].

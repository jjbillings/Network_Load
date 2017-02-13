# Network_Load
C++ app to simulate the load on an optical network.
###To-Do:
- [x] 1. Rall Symposium Abstract 
- [ ] 2. Phase 0.5
	- [ ] * Finish CPU implementation
		- [x] * Compute all simple paths of length N between all pairs of vertices
		- [x] * Store results in memory
		- [ ] * Copy each SimplePath into a Path to represent a primary path
		- [ ] * For each primary path, compute cost of backup path using stored simple paths/channels/load
		- [ ] * Select cheapest primary/backup combo
		- [ ] * Increase Network Load
	- [ ] * Measure CPU running time
	- [ ] * **This needs to be done for Rall.**
- [ ] 3. Phase 1
	- [ ] * Develop technique for splitting up work accross cores
	- [ ] * Implement algorithm for GPU
	- [ ] * Measure GPU running time
- [ ] 4. Phase 2
	- [ ] * Compare CPU/GPU running time
	- [ ] * **If we have this for Rall, that would be dope.**

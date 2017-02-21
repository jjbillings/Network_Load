# Network_Load
C++ app to simulate the load on an optical network.
###To-Do:
- [x] 1. Rall Symposium Abstract 
- [ ] 2. Phase 0.5
	- [x] * Finish CPU implementation
		- [x] * Compute all simple paths of length N between all pairs of vertices
		- [x] * Store results in memory
		- [x] * Copy each SimplePath into a Path to represent a primary path
		- [x] * For each primary path, compute cost of backup path using stored simple paths/channels/load
		- [x] * Select cheapest primary/backup combo
		- [x] * Increase Network Load
	- [ ] * Measure CPU running time
	- [ ] * **This needs to be done for Rall.**
- [ ] 3. Phase 1
	- [x] * Develop technique for splitting up work accross cores
	- [ ] * Implement algorithm for GPU
		- [x] * Generate connection dataset for testing
		- [ ] * Memory Management
		- [x] * Copy SimplePaths array to GPU
		- [ ] * Write Kernel for determineCompatibleBackups()
			- [ ] * Combine pathCosts array and potPathInd array into one
	- [ ] * Measure GPU running time
- [ ] 4. Phase 2
	- [ ] * Compare CPU/GPU running time
	- [ ] * **If we have this for Rall, that would be dope.**

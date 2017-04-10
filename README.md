# Network_Load
C++ app to simulate the load on an optical network.
###To-Do:
- [x] 1. Rall Symposium Abstract 
- [x] 2. Phase 0.5
	- [x] * Finish CPU implementation
		- [x] * Compute all simple paths of length N between all pairs of vertices
		- [x] * Store results in memory
		- [x] * Copy each SimplePath into a Path to represent a primary path
		- [x] * For each primary path, compute cost of backup path using stored simple paths/channels/load
		- [x] * Select cheapest primary/backup combo
		- [x] * Increase Network Load
	- [x] * Measure CPU running time
	- [ ] * **This needs to be done for Rall.**
- [ ] 3. Phase 1
	- [x] * Develop technique for splitting up work accross cores
	- [ ] * Implement algorithm for GPU
		- [x] * Generate connection dataset for testing
		- [ ] * Memory Management - not like it's ever actually done...
		- [x] * Copy SimplePaths array to GPU
		- [x] * Write Kernel for determineCompatibleBackups()
			- [x] * Figure out how to index into the arrays via threads/blocks
			- [x] * Combine pathCosts array and potPathInd array into one
			- [x] * Copy potPathCosts array back to Host
		- [x] * Write Kernel for computeCost()
			- [x] * Copy channels array to the GPU
			- [x] * Determine cost per primary/backup combo
			- [x] * Copy results back to Host
	- [x] * Measure GPU running time
- [ ] 4. Phase 2
	- [x] * Compare CPU/GPU running time
	- [ ] * **If we have this for Rall, that would be dope.**
- [ ] 5. Try a larger network (COSTS Network was recommended by Kim)
- [ ] 6. Use a much larger set of test data.
- [ ] 7. Drastically change the parameters - Allocate like 100 random connections and then measure the CPU vs. GPU performance over a large number of samples (1 sample being an allocation of 100 random connections). Then increase or decrease the parameters one at a time. MAX_CHANNELS,NUM_CONNECTIONS I think would be appropriate ones to start with.


Data Collection: Keep MAX_CHANNELS and NUM_CONNECTIONS constant, allocate varying numbers of random connections on the network, measure cpu vs gpu performance. Take a large number of samples for each quantity of random connections.

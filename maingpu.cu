/*
 * File:   maingpu.cu
 * Author: jjbillings
 *
 * Created on October 16, 2016, 9:09 PM
 */

#include<cstdlib>
#include<stdio.h>
#include<queue>
#include<stack>
#include<iostream>
#include<fstream>
#include<ctime>
#include"nets.h"

using namespace std;

#define NUM_CONNECTIONS 500
#define MAX_CHANNELS 30
#define SAMPLES 1

struct SimplePath;
struct Path;
struct Edge;
struct Connection;
struct Connection;
struct Channel;

struct Channel{
    bool primary; //is this channel used for a primary path?
    int numBackups; //total protected;
    Connection *backupsOnChannel[NUM_CONNECTIONS];//Realistically, there will be far fewer than NUM_CONNECTIONS
    Connection *d_backupsOnChannel[NUM_CONNECTIONS];
};

struct Edge {
    int edgeNum;
    int v1;
    int v2;
    int load; //load <= MAX_CHANNELS. Also, load is the sum of the primary AND backups paths using it.
    int totalProtected;
};

struct SimplePath {
    int sourceNode;
    int destNode;

    int hops;
    int index;
    int edgeNums[N_NODES];

    Edge *edges[N_NODES];

  SimplePath() {
    for(int i = 0; i < N_NODES; ++i) {
      edgeNums[i] = -1;
      edges[i] = 0;
    }
    sourceNode = -1;
    destNode = -1;
    hops = -1;
    index = -1;
  };
};

struct Path {
    int sourceNode;
    int destNode;

    int hops;
    int index;
    int cost;

    //Every path that uses a particular edge just has a reference to it (not a copy), so they can each manipulate it.
    Edge *edges[N_NODES];
    bool freeEdges[N_NODES]; //whether or not that edge has a cost of 0
    int channelNum[N_NODES]; //Channel number for each edge that it uses
    int edgeNums[N_NODES];
    bool primary;
    bool active;
};

struct Connection {
    int sourceNode;
    int destNode;
    int combinedCost;
    bool validBackup;
    bool validPrimary;
    Path *backupPath;
    Path *primaryPath;
};

void readGraphReorderEdgeList(int vertexList[],Edge compactEdgeList[2*N_EDGES],Edge reorderedEdgeList[2*N_NODES]);

int computeAllSimplePathsN(SimplePath **ps, int *vertexList, Edge *edgeList, int sourceNode, int destNode, int hops);
void simulate(int *vertexList, Edge *edgeList);
void simulate_GPU(int *vertexList, Edge *edgeList);
void computeCostForBackupsWithGPU(SimplePath *p, int *potPathCosts, int primaryInd, Channel cs[2*N_EDGES][MAX_CHANNELS]);
int determineCompatibleBackups(SimplePath *p, int *potPathInd, int numPossiblePaths, int pInd);
void computeCostForBackups(SimplePath *p, int *potPathInd, int numPotPaths, int backupIndex, int *pathCosts,Channel cs[2*N_EDGES][MAX_CHANNELS]);
void selectChannels(Connection *c, Channel chan[2*N_EDGES][MAX_CHANNELS]);
void increaseLoad(Connection *connection, Channel channels[2*N_EDGES][MAX_CHANNELS], Connection *d_con);

int vertexList[N_NODES+1];
Edge edgeList[2*N_EDGES];
Edge reorderedEdgeList[2*N_EDGES];
Connection cons[NUM_CONNECTIONS];
Channel channels[2*N_EDGES][MAX_CHANNELS];



//-----------Kernel for Determining which Backups are compatible with which Primaries. WORKING---------//
__global__ void determineCompatibleBackups(SimplePath *ps, int *potPathCosts,int conInd){
 
  int p_ind = (conInd * NUM_CONNECTIONS) +  blockIdx.x;
  int b_ind = (conInd * NUM_CONNECTIONS) +  threadIdx.x;
  int output_ind = (blockIdx.x * NUM_CONNECTIONS) + threadIdx.x;

  int primIndex = ps[p_ind].index;
  int backIndex = ps[b_ind].index;

  int primHops = ps[p_ind].hops;
  int backHops = ps[b_ind].hops;
  
  if(primHops > 0 && backHops > 0) {
    bool disjoint = true;

    for(int e1 = 0; disjoint && e1 <= primIndex; ++e1) {
      for(int e2 = 0; disjoint && e2 <= backIndex; ++e2){
	if(ps[p_ind].edgeNums[e1] == ps[b_ind].edgeNums[e2]) {
	  disjoint = false;
	}
      }
    }
    if(disjoint) {
      potPathCosts[output_ind] = 1;
    }else {
      potPathCosts[output_ind] = -1;
    }
  }else {
    potPathCosts[output_ind] = -1;
  }
}

//-----------TEST_KERNEL_FOR_WARPS---------//
__global__ void determineCompatibleBackups2(SimplePath *ps, int *potPathCosts,int conInd){

  int warp_id = threadIdx.x / 32;
  int warp_offset = threadIdx.x % 32;
  
  int p_ind = (conInd * NUM_CONNECTIONS) +  blockIdx.x;
  int b_ind = (conInd * NUM_CONNECTIONS) +  threadIdx.x;
  int output_ind = (blockIdx.x * NUM_CONNECTIONS) + threadIdx.x;

  int primIndex = ps[p_ind].index;
  int backIndex = ps[b_ind].index;

  int primHops = ps[p_ind].hops;
  int backHops = ps[b_ind].hops;
  
  if(primHops > 0 && backHops > 0) {
    bool disjoint = true;

    for(int e1 = 0; disjoint && e1 <= primIndex; ++e1) {
      for(int e2 = 0; disjoint && e2 <= backIndex; ++e2){
	if(ps[p_ind].edgeNums[e1] == ps[b_ind].edgeNums[e2]) {
	  disjoint = false;
	}
      }
    }
    if(disjoint) {
      potPathCosts[output_ind] = 1;
    }else {
      potPathCosts[output_ind] = -1;
    }
  }else {
    potPathCosts[output_ind] = -1;
  }
}

__global__ void costsKernel(SimplePath *p, int *potPathCosts, int conInd , Channel *cs) {

  int p_ind = (conInd * NUM_CONNECTIONS) + blockIdx.x;
  int b_ind = (conInd * NUM_CONNECTIONS) + threadIdx.x;
  int index = (blockIdx.x * NUM_CONNECTIONS) + threadIdx.x;


  //If we already know that this combo is unusable, just quit.
        if(potPathCosts[index] == -1) {
	    return;
        }
        
        int cost = 0;

        for(int e = 0; e <= p[b_ind].index; ++e) {
            bool free = false;
            int edgeNum = p[b_ind].edgeNums[e];
	    
            int firstOpenChannel = MAX_CHANNELS+1;

            for(int c = 0; !free && c < MAX_CHANNELS; ++c) {

	        int channelIndex = (edgeNum * MAX_CHANNELS)+c;
	        if(cs[channelIndex].primary == true) {
                    continue;
                }

                //At this point, we know that there are no primary paths on this channel
                //Thus we must check and see if it is "free".

                //we COULD use this channel, but there may be a "free" one further down.
                if(cs[channelIndex].numBackups == 0) {
                    if(c < firstOpenChannel) {
                        firstOpenChannel = c;
                    }
                    continue;
                }

                bool disjoint = true;

                //Check every connection currently on protected on the channel
                for(int bup = 0; disjoint && bup < cs[channelIndex].numBackups; ++bup) {

                    //At this point, we know that there is at least one path protected on this channel.
                    //Technically, we should also know that it's not a primary path.

                    //for each edge of the protected connection's primary path
                    for(int e2 = 0; disjoint && e2 <= (*(*cs[channelIndex].d_backupsOnChannel[bup]).primaryPath).index; ++e2) {

                        //see if its the same edge as used by our primary path.
                        for(int e3 = 0; disjoint && e3 <= p[p_ind].index; ++e3 ) {

                            if((*(*cs[channelIndex].d_backupsOnChannel[bup]).primaryPath).edgeNums[e2] == p[p_ind].edgeNums[e3]) {
                                //There is a non-disjoint primary path on this channel, so it is unusable.

                                disjoint = false;
                            }
                        }
                    }
                }

                if(disjoint) {
                    //This channel is free
                    free = true;
                }
            }

            if(!free) {
                if(firstOpenChannel < MAX_CHANNELS) {
                    cost++;
                }else {
                    cost = 1000000;
                    break;
                }

            }

        }

        potPathCosts[index] = cost;
    
}

/*
 *TODO: I totally thought I made the algorithm be based on BFS, but it is in fact based on DFS.
 *So REVERSE the order of the edge list. Currently, the neighbor with the lowest degree gets pushed
 *to the "bottom" of the stack, so we end up computing the path with high-degree nodes in it...
 */
int main(int argc, char** argv) {
    cout <<"Welcome to main\n";

    cudaFree(0);
    for(int f = 0; f < (2*N_EDGES); ++f){
        for(int g = 0; g < MAX_CHANNELS; ++g) {
            channels[f][g].numBackups = 0;
            channels[f][g].primary = false;
        }
    }
    readGraphReorderEdgeList(vertexList,edgeList,reorderedEdgeList);

    srand(time(NULL));

    simulate_GPU(vertexList,edgeList);
    //simulate(vertexList,edgeList);
    return 0;
}

void simulate_GPU(int *vertexList, Edge *edgeList){

    clock_t cpu_startTime, cpu_endTime;
    double cpu_elapsedTime = 0;
    float gpu_totalTime = 0;
    cpu_startTime = clock();
    
    int connectionNum = 0;
    const size_t sp_size = sizeof(SimplePath);
    const size_t potPathCosts_size = (NUM_CONNECTIONS * NUM_CONNECTIONS) * sizeof(int);
    const size_t ps_size = ((N_NODES*N_NODES)*NUM_CONNECTIONS)*sp_size; //Size of the entire 2D array
    const size_t row_size = NUM_CONNECTIONS*sp_size; //Size of a SINGLE row in the array of SimplePaths

    const size_t channels_size = ((2*N_EDGES)*MAX_CHANNELS)*sizeof(Channel);
    
    //Test Data
    int v1[40] = {9, 5, 6, 1, 3, 5, 4, 9, 9, 9, 7, 8, 2, 10, 3, 5, 9, 3, 2, 3, 5, 2, 3, 3, 10, 9, 10, 2, 1, 1, 3, 2, 9, 5, 4, 6, 10, 5, 0, 1};
    int v2[40] = {3, 8, 4, 3, 8, 3, 7, 1, 5, 6, 0, 6, 10, 5, 8, 2, 3, 6, 5, 4, 2, 3, 9, 7, 9, 5, 6, 5, 0, 2, 5, 5, 10, 3, 9, 3, 4, 1, 10, 2};
    
    SimplePath **ps = new SimplePath*[N_NODES * N_NODES]; //Storage for paths

    SimplePath *d_ps; //Device pointer for the array of SimplePaths
    int *d_potPathCosts; //Device pointer for the array of Potential Path Costs
    int *h_potPathCosts; //Host pointer for the array of potential path costs.

    Connection *d_cons; //Device pointer to the array of connections.
    Channel *d_channels; //Device pointer for the array of channels.
    
    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
        ps[i] = new SimplePath[NUM_CONNECTIONS];
    }


    if(cudaSuccess != cudaMalloc((void **)&d_ps,ps_size)) {
    	cout << "Malloc Error\n";
    }else {
      cout << "allocated SimplePaths array on Device\n";
    }

    
    if(cudaSuccess != cudaMalloc((void **)&d_channels,channels_size)) {
	cout << "Error Allocating channels on GPU\n";
    }else {
	cout << "Allocated Channels array on GPU\n";
    }

    cudaMalloc((void **)&d_cons,sizeof(Connection)*NUM_CONNECTIONS);
    

    cudaMalloc((void **)&d_potPathCosts,potPathCosts_size);
    cout << "Allocated potential Path Costs array on device\n";

    cudaMemcpy(d_channels,&channels,channels_size,cudaMemcpyHostToDevice);

    
    h_potPathCosts = (int *)malloc(potPathCosts_size);

    //We COULD parallelize this by giving a thread a source/dest combo to compute the paths of. potentially beneficial for large graphs
    for(int src = 0; src < N_NODES; ++src) {
        for(int dest = 0; dest < N_NODES; ++dest) {
            if(src != dest) {
                int index = (src*N_NODES)+dest;
                computeAllSimplePathsN(ps,vertexList,edgeList,src,dest,N_NODES);
                //cout <<"All simple paths computed and stored! " << npaths[index] << " paths between " << src << " and " << dest << "\n";
            }
        }
    }


    //Copy Simple paths to the GPU
    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
      cudaMemcpy(d_ps + (i*(NUM_CONNECTIONS)),ps[i],row_size,cudaMemcpyHostToDevice);
    }

    //Setup components for GPU benchmarking.
    cudaEvent_t start, stop;
    

    //cpu_startTime = clock();
    for(int c = 0; c < 5; ++c) {
      
      //Attempt to allocate SOME connection onto the network
      int s = v1[connectionNum];
      int d = v2[connectionNum];

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //Allocate storage for the potential primary/backup path combos
    int index = (s*N_NODES) + d;

    
    //BENCHMARKING
    cudaEventRecord(start);
    
    //-----------Launch the Kernel-------------//
    determineCompatibleBackups<<<NUM_CONNECTIONS,NUM_CONNECTIONS>>>(d_ps, d_potPathCosts,index);
    //    cudaDeviceSynchronize();

    //BENCHMARKING
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milli = 0;
    cudaEventElapsedTime(&milli,start,stop);
    gpu_totalTime += milli;
    //cout << "Kernel Execution took: " << milli << " milliseconds\n";

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    
    if(cudaSuccess != cudaGetLastError()) {
      cout << "CUDA ERROR IN KERNEL: " << cudaGetLastError() << "\n";
    }

    //---------Launch the Kernel----------//
    costsKernel<<<NUM_CONNECTIONS,NUM_CONNECTIONS>>>(d_ps, d_potPathCosts, index,d_channels);
    
    //---------Copy the Results back to the host ---//
    cudaMemcpy(h_potPathCosts,d_potPathCosts,potPathCosts_size,cudaMemcpyDeviceToHost);    

    //for(int i = 0; i < NUM_CONNECTIONS; ++i) {
    //  computeCostForBackupsWithGPU(ps[index],h_potPathCosts,i,channels);
    //}

    //-----------Select the cheapest combo using GPU Results-----------//
    int minCostGPU = 100000000;
    int minPrimIndGPU = -1;
    int minBackIndGPU = -1;

    for(int p = 0; p < NUM_CONNECTIONS; ++p) {
        int primaryCostGPU = ps[index][p].hops;

        for(int b = 0; b < NUM_CONNECTIONS; ++b) {
	  int f = (p*NUM_CONNECTIONS)+b;
	  if(h_potPathCosts[(p*NUM_CONNECTIONS)+b] < 0) {
	      continue;
	  }
          if((h_potPathCosts[(p*NUM_CONNECTIONS)+b] + primaryCostGPU) < minCostGPU) {
                minCostGPU = (h_potPathCosts[(p*NUM_CONNECTIONS)+b] + primaryCostGPU);
                minPrimIndGPU = p;
                minBackIndGPU = b;
          }
        }
    }

    cout << "Min cost on GPU is: " << minCostGPU << "\n";
    cout << "PRIM: "<<minPrimIndGPU << "\n";
    for(int i = 0; i <= ps[index][minPrimIndGPU].index; ++i) {
      cout << (*ps[index][minPrimIndGPU].edges[i]).v1 << " -> " << (*ps[index][minPrimIndGPU].edges[i]).v2 << "\n";
    }
    cout << "BACK: " << minBackIndGPU << "\n";
    for(int i = 0; i <= ps[index][minBackIndGPU].index; ++i) {
      cout << (*ps[index][minBackIndGPU].edges[i]).v1 << " -> " << (*ps[index][minBackIndGPU].edges[i]).v2 << "\n";
    }


    //--------------Store the connection--------------//
    cons[connectionNum].sourceNode = s;
    cons[connectionNum].destNode = d;
    cons[connectionNum].combinedCost = minCostGPU;
    cons[connectionNum].validBackup = true;
    cons[connectionNum].validPrimary = true;
    cons[connectionNum].backupPath = new Path();
    cons[connectionNum].primaryPath = new Path();
    (*cons[connectionNum].primaryPath).hops = ps[index][minPrimIndGPU].hops;
    (*cons[connectionNum].primaryPath).index = ps[index][minPrimIndGPU].index;
    (*cons[connectionNum].primaryPath).primary = true;
    (*cons[connectionNum].backupPath).hops = ps[index][minBackIndGPU].hops;
    (*cons[connectionNum].backupPath).index = ps[index][minBackIndGPU].index;

    for(int p = 0; p <= ps[index][minPrimIndGPU].index; ++p) {
        (*cons[connectionNum].primaryPath).edges[p] = ps[index][minPrimIndGPU].edges[p];
        (*cons[connectionNum].primaryPath).freeEdges[p] = false;
	(*cons[connectionNum].primaryPath).edgeNums[p] = ps[index][minPrimIndGPU].edgeNums[p];
    }
    for(int p = 0; p <= ps[index][minBackIndGPU].index; ++p) {
        (*cons[connectionNum].backupPath).edges[p] = ps[index][minBackIndGPU].edges[p];
	(*cons[connectionNum].backupPath).edgeNums[p] = ps[index][minBackIndGPU].edgeNums[p];
    }

    

    //Select Channels
    selectChannels(&cons[connectionNum],channels);

    //Increase the network load
    increaseLoad(&cons[connectionNum],channels,&d_cons[connectionNum]);

    //NOTE: We can 100% only copy individual channels to the GPU. i.e. if only channels 3 and 41 were updated, we can copy ONLY those channels if we want to
    cudaMemcpy(d_channels,&channels,channels_size,cudaMemcpyHostToDevice);

    //TESTING FOR COSTS KERNEL
    cudaMemcpy(d_cons,&cons,sizeof(Connection)*NUM_CONNECTIONS,cudaMemcpyHostToDevice);
    
    //--------------Print Network Load--------------//
    for(int m = 0; m < 2*N_EDGES; ++m) {
        cout << "LOAD: " << edgeList[m].v1 << " -> " << edgeList[m].v2 << ": " << edgeList[m].load << " | TP: " << edgeList[m].totalProtected << " | ";
        if(edgeList[m].load > 0) {
            for(int c = 0; c < edgeList[m].load; ++c) {
                cout << "C" << c << ": " << channels[m][c].numBackups << " ";
                if(channels[m][c].primary == true) {
                    cout << "P ";
                }
            }
        }
        cout << "\n";

    }
    

    connectionNum++;
    }//ENDFOR

    cudaDeviceSynchronize();
    cpu_endTime = clock();


    //--------------Clean up memory--------------//
    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
        delete[] ps[i];
    }
    delete[] ps;
    
    cudaFree(d_ps);
    cudaFree(d_potPathCosts);
    cudaFree(d_channels);
    cudaFree(d_cons);
    
    free(h_potPathCosts);
    //cpu_endTime = clock();
    cpu_elapsedTime = ((double) (cpu_endTime - cpu_startTime)/CLOCKS_PER_SEC) * 1000;

        cout << "Kernel Execution took: " << gpu_totalTime << " milliseconds\n";
	cout << "Total time: " << cpu_elapsedTime << " milliseconds\n";
	//cout << "CPU Start: " << cpu_startTime << " CPU End: " << cpu_endTime << "\n";
}


void computeCostForBackupsWithGPU(SimplePath *p, int *potPathCosts, int primaryInd, Channel cs[2*N_EDGES][MAX_CHANNELS]) {

    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
        int pid = (primaryInd * NUM_CONNECTIONS) + i;
        if(potPathCosts[pid] == -1) {
            continue;
        }
        
        int cost = 0;

        for(int e = 0; e <= p[i].index; ++e) {
            bool free = false;
            int edgeNum = (*p[i].edges[e]).edgeNum;
            int firstOpenChannel = MAX_CHANNELS+1;

            for(int c = 0; !free && c < MAX_CHANNELS; ++c) {

                if(cs[edgeNum][c].primary == true) {
                    continue;
                }

                //At this point, we know that there are no primary paths on this channel
                //Thus we must check and see if it is "free".

                //we COULD use this channel, but there may be a "free" one further down.
                if(cs[edgeNum][c].numBackups == 0) {
                    if(c < firstOpenChannel) {
                        firstOpenChannel = c;
                    }
                    continue;
                }

                bool disjoint = true;

                //Check every connection currently on protected on the channel
                for(int bup = 0; disjoint && bup < channels[edgeNum][c].numBackups; ++bup) {

                    //At this point, we know that there is at least one path protected on this channel.
                    //Technically, we should also know that it's not a primary path.

                    //for each edge of the protected connection's primary path
                    for(int e2 = 0; disjoint && e2 <= (*(*channels[edgeNum][c].backupsOnChannel[bup]).primaryPath).index; ++e2) {

                        //see if its the same edge as used by our primary path.
                        for(int e3 = 0; disjoint && e3 <= p[primaryInd].index; ++e3 ) {

                            if((*(*channels[edgeNum][c].backupsOnChannel[bup]).primaryPath).edges[e2] == p[primaryInd].edges[e3]) {
                                //There is a non-disjoint primary path on this channel, so it is unusable.

                                disjoint = false;
                            }
                        }
                    }
                }

                if(disjoint) {
                    //This channel is free
                    free = true;
                }
            }

            if(!free) {
                if(firstOpenChannel < MAX_CHANNELS) {
                    cost++;
                }else {
                    cost = 1000000;
                    break;
                }

            }

        }

        potPathCosts[pid] = cost;
    }
}








void simulate(int *vertexList, Edge *edgeList){
    clock_t cpu_startTime, cpu_endTime;
    double cpu_elapsedTime = 0;
    cpu_startTime = clock();

    //Test Data
    int v1[40] = {9, 5, 6, 1, 3, 5, 4, 9, 9, 9, 7, 8, 2, 10, 3, 5, 9, 3, 2, 3, 5, 2, 3, 3, 10, 9, 10, 2, 1, 1, 3, 2, 9, 5, 4, 6, 10, 5, 0, 1};
    int v2[40] = {3, 8, 4, 3, 8, 3, 7, 1, 5, 6, 0, 6, 10, 5, 8, 2, 3, 6, 5, 4, 2, 3, 9, 7, 9, 5, 6, 5, 0, 2, 5, 5, 10, 3, 9, 3, 4, 1, 10, 2};
    
    int connectionNum = 0;
    //We want to compute and store all possible paths between our source and desitination.
    SimplePath **ps = new SimplePath*[N_NODES * N_NODES]; //Storage for paths
    int *npaths = new int[N_NODES*N_NODES];

    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
        ps[i] = new SimplePath[NUM_CONNECTIONS];
    }


    //We COULD parallelize this by giving a thread a source/dest combo to compute the paths of. potentially beneficial for large graphs
    for(int src = 0; src < N_NODES; ++src) {
        for(int dest = 0; dest < N_NODES; ++dest) {
            if(src != dest) {
                int index = (src*N_NODES)+dest;
                npaths[index] = computeAllSimplePathsN(ps,vertexList,edgeList,src,dest,N_NODES);
                //cout <<"All simple paths computed and stored! " << npaths[index] << " paths between " << src << " and " << dest << "\n";
            }
        }
    }
    //At this point, we COULD delete[] any paths in the array that we didn't use.


    for(int num = 0; num < 40; ++num) {
    //Attempt to allocate SOME connection onto the network
    int s = v1[connectionNum];
    int d = v2[connectionNum];
    //int s = rand() % N_NODES;
    //int d = rand() % N_NODES;
    //while(s == d) {
    //    s = rand()%N_NODES;
    //    d = rand()%N_NODES;
    //}

    //Allocate storage for the potential primary/backup path combos
    int index = (s*N_NODES) + d;
    int numPossiblePaths = npaths[index];

    //Stores indices into the ps[index][] array for each disjoint backup path.
    //potPathInd[i][j] = k where ps[index][k] is a path that is edge-disjoint from ps[index][i].
    int ** potPathInd = new int*[NUM_CONNECTIONS];
    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
        potPathInd[i] = new int[NUM_CONNECTIONS];
    }


    //--------------Find all paths which are edge-disjoint from this primary--------------//
    int k = -1;
    //On the GPU, instead of iterating i..numPossiblePaths, we would give thread_i backup_i
    for(int i = 0; i < numPossiblePaths; ++i) {
        k = determineCompatibleBackups(ps[index],potPathInd[i],numPossiblePaths,i);
        //cout << "Number of paths which are disjoint from this primary path: " << k << "\n";
    }



    //--------------Compute Cost for each backup path--------------//
    int ** pathCosts = new int*[numPossiblePaths];
    for(int i = 0; i < numPossiblePaths; ++i) {
        pathCosts[i] = new int[numPossiblePaths];
    }

    for(int i = 0; i < numPossiblePaths; ++i) {
        computeCostForBackups(ps[index],potPathInd[i],numPossiblePaths,i,pathCosts[i],channels);
    }



    //--------------Select cheapest connection--------------//
    int minCost = 100000000;
    int minPrimInd = -1;
    int minBackInd = -1;

    for(int p = 0; p < numPossiblePaths; ++p) {
        int backInd = 0;
        int primaryCost = ps[index][p].hops;

        while(pathCosts[p][backInd] != -1) {
            if((pathCosts[p][backInd] + primaryCost) < minCost) {
                minCost = (pathCosts[p][backInd] + primaryCost);
                minPrimInd = p;
                minBackInd = backInd;
            }
            backInd++;
        }
    }
    cout << "Min cost is: " << minCost << "\n";



    //--------------Store the connection--------------//
    cons[connectionNum].sourceNode = s;
    cons[connectionNum].destNode = d;
    cons[connectionNum].combinedCost = minCost;
    cons[connectionNum].validBackup = true;
    cons[connectionNum].validPrimary = true;
    cons[connectionNum].backupPath = new Path();
    cons[connectionNum].primaryPath = new Path();
    (*cons[connectionNum].primaryPath).hops = ps[index][minPrimInd].hops;
    (*cons[connectionNum].primaryPath).index = ps[index][minPrimInd].index;
    (*cons[connectionNum].primaryPath).primary = true;
    (*cons[connectionNum].backupPath).hops = ps[index][potPathInd[minPrimInd][minBackInd]].hops;
    (*cons[connectionNum].backupPath).index = ps[index][potPathInd[minPrimInd][minBackInd]].index;

    for(int p = 0; p <= ps[index][minPrimInd].index; ++p) {
        (*cons[connectionNum].primaryPath).edges[p] = ps[index][minPrimInd].edges[p];
        (*cons[connectionNum].primaryPath).freeEdges[p] = false;
    }
    for(int p = 0; p <= ps[index][potPathInd[minPrimInd][minBackInd]].index; ++p) {
        (*cons[connectionNum].backupPath).edges[p] = ps[index][potPathInd[minPrimInd][minBackInd]].edges[p];
    }

    //Select Channels
    selectChannels(&cons[connectionNum],channels);

    //Increase the network load
    //    increaseLoad(&cons[connectionNum],channels); TODO: REIMPLEMENT WITH THE CORRECT NUM OF PARAMETERS FOR NON-GPU version.


    //--------------Print Network Load--------------//
    for(int m = 0; m < 2*N_EDGES; ++m) {
        cout << "LOAD: " << edgeList[m].v1 << " -> " << edgeList[m].v2 << ": " << edgeList[m].load << " | TP: " << edgeList[m].totalProtected << " | ";
        if(edgeList[m].load > 0) {
            for(int c = 0; c < edgeList[m].load; ++c) {
                cout << "C" << c << ": " << channels[m][c].numBackups << " ";
                if(channels[m][c].primary == true) {
                    cout << "P ";
                }
            }
        }
        cout << "\n";

    }


    //--------------Clean up memory--------------//
    for(int i = 0; i < numPossiblePaths; ++i) {
        delete[] potPathInd[i];
    }
    delete[] potPathInd;

    for(int i = 0; i < numPossiblePaths; ++i) {
        delete[] pathCosts[i];
    }
    delete[] pathCosts;
    connectionNum++;
}//end loop

    cpu_endTime = clock();
    
    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
        delete[] ps[i];
    }
    delete[] ps;
    delete[] npaths;
    cout << "ps and npaths deleted\n";

    //cpu_endTime = clock();
    cpu_elapsedTime = ((double)(cpu_endTime-cpu_startTime)/CLOCKS_PER_SEC)*1000;
    cout << "CPU Total Elapsed Time: " << cpu_elapsedTime << "\n";
}

void increaseLoad(Connection *connection, Channel channels[2*N_EDGES][MAX_CHANNELS], Connection *d_con) {
    if((*(*connection).primaryPath).index < 0) {
        cout << "Primary Path DNE?\n";
        return;
    }
    //Increment the network load; put the backup on the channels

    //Here we are incrementing the network load for the PRIMARY PATH
    for(int i = 0; i <= (*(*connection).primaryPath).index; ++i) {

        //Every edge in the primary path gets its load increased
        channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*connection).primaryPath).channelNum[i]].primary = true;
        channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*connection).primaryPath).channelNum[i]].backupsOnChannel[0] = connection;
	channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*connection).primaryPath).channelNum[i]].d_backupsOnChannel[0] = d_con;
        channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*connection).primaryPath).channelNum[i]].numBackups += 1;
        (*(*(*connection).primaryPath).edges[i]).load += 1;
        (*(*(*connection).primaryPath).edges[i]).totalProtected += 1;
    }

    //Here we are increasing the network load for the BACKUP PATH
    for(int i = 0; i <= (*(*connection).backupPath).index; ++i) {
        //Temp
        Edge *e = (*(*connection).backupPath).edges[i];
        int cNum = (*(*connection).backupPath).channelNum[i];

        //first path to use this channel, or this is not a free edge for the backup path.
        //if(channels[(*e).edgeNum][cNum].numBackups == 0 || (*(*connection).backupPath).freeEdges[i] == false) {
        if((*(*connection).backupPath).freeEdges[i] == false) {
            (*e).load += 1;
        }

        //Marks that the connection is protected on this channel.
        int en = (*e).edgeNum;
        int numbs = channels[en][cNum].numBackups;
        channels[en][cNum].primary = false;
        channels[en][cNum].backupsOnChannel[numbs] = connection;
	channels[en][cNum].d_backupsOnChannel[numbs] = d_con;
        channels[en][cNum].numBackups += 1;
        (*e).totalProtected +=1;
    }

}

//TODO: This method contains a lot of redundant code that is also in computeCostForBackups. Consider combining.
//I wanted to modularize the code as much as possible this time around, which is why there's so much redundancy in this method.
void selectChannels(Connection *c, Channel chan[2*N_EDGES][MAX_CHANNELS]) {

  /*
    cout << "prim\n";
    for(int i = 0; i <= (*(*c).primaryPath).index; ++i) {
        cout << (*(*(*c).primaryPath).edges[i]).v1 << " -> " << (*(*(*c).primaryPath).edges[i]).v2 << "\n";
    }
    cout << "back\n";
    for(int i = 0; i <= (*(*c).backupPath).index; ++i) {
        cout << (*(*(*c).backupPath).edges[i]).v1 << " -> " << (*(*(*c).backupPath).edges[i]).v2 << "\n";
    }

  */
    int edgeNum = -1;
    //Select Primary path channels;
    for(int p = 0; p <= (*(*c).primaryPath).index; ++p){
        edgeNum = (*(*(*c).primaryPath).edges[p]).edgeNum;
        bool allSet = false;
        for(int ch = 0; !allSet && ch < MAX_CHANNELS; ++ch) {
            if(chan[edgeNum][ch].numBackups == 0) {
                allSet = true;
                (*(*c).primaryPath).channelNum[p] = ch;
            }
        }
    }

    for(int e = 0; e <= (*(*c).backupPath).index; ++e) {
        bool free = false;
        edgeNum = (*(*(*c).backupPath).edges[e]).edgeNum;
        int firstOpenChannel = MAX_CHANNELS+1;

        for(int ch = 0; !free && ch < MAX_CHANNELS; ++ch) {

            if(chan[edgeNum][ch].primary == true) {
                continue;
            }

            //At this point, we know that there are no primary paths on this channel
            //Thus we must check and see if it is "free".

            //we COULD use this channel, but there may be a "free" one further down.
            if(chan[edgeNum][ch].numBackups == 0) {
                if(ch < firstOpenChannel) {
                    firstOpenChannel = ch;
                }
                continue;
            }

            bool disjoint = true;

            //Check every connection currently on protected on the channel
            for(int bup = 0; disjoint && bup < chan[edgeNum][ch].numBackups; ++bup) {

                //At this point, we know that there is at least one path protected on this channel.
                //Technically, we should also know that it's not a primary path.

                //for each edge of the protected connection's primary path
                for(int e2 = 0; disjoint && e2 <= (*(*chan[edgeNum][ch].backupsOnChannel[bup]).primaryPath).index; ++e2) {

                    //see if its the same edge as used by our primary path.
                    for(int e3 = 0; disjoint && e3 <= (*(*c).primaryPath).index; ++e3 ) {

                        if((*(*chan[edgeNum][ch].backupsOnChannel[bup]).primaryPath).edges[e2] == (*(*c).primaryPath).edges[e3]) {
                            //There is a non-disjoint primary path on this channel, so it is unusable.
                            //goto CHANNEL_LOOP_END;
                            disjoint = false;
                        }
                    }
                }
            }

            if(disjoint) {
                //This channel is free
                free = true;
                (*(*c).backupPath).channelNum[e] = ch;
                (*(*c).backupPath).freeEdges[e] = true;
            }
        }

        if((*(*c).backupPath).freeEdges[e] == false) {
            (*(*c).backupPath).channelNum[e] = firstOpenChannel;
        }
    }
    //cout << "all set?\n";
}

void computeCostForBackups(SimplePath *p, int *potPathInd, int numPossiblePaths, int primaryInd, int *pathCosts, Channel cs[2*N_EDGES][MAX_CHANNELS]) {

    for(int i = 0; i < numPossiblePaths; ++i) {
        if(potPathInd[i] == -1) {
            pathCosts[i] = -1;
            break;
        }
        int pid = potPathInd[i];
        int cost = 0;

        for(int e = 0; e <= p[pid].index; ++e) {
            bool free = false;
            int edgeNum = (*p[pid].edges[e]).edgeNum;
            int firstOpenChannel = MAX_CHANNELS+1;

            for(int c = 0; !free && c < MAX_CHANNELS; ++c) {

                if(cs[edgeNum][c].primary == true) {
                    continue;
                }

                //At this point, we know that there are no primary paths on this channel
                //Thus we must check and see if it is "free".

                //we COULD use this channel, but there may be a "free" one further down.
                if(cs[edgeNum][c].numBackups == 0) {
                    if(c < firstOpenChannel) {
                        firstOpenChannel = c;
                    }
                    continue;
                }

                bool disjoint = true;

                //Check every connection currently on protected on the channel
                for(int bup = 0; disjoint && bup < channels[edgeNum][c].numBackups; ++bup) {

                    //At this point, we know that there is at least one path protected on this channel.
                    //Technically, we should also know that it's not a primary path.

                    //for each edge of the protected connection's primary path
                    for(int e2 = 0; disjoint && e2 <= (*(*channels[edgeNum][c].backupsOnChannel[bup]).primaryPath).index; ++e2) {

                        //see if its the same edge as used by our primary path.
                        for(int e3 = 0; disjoint && e3 <= p[primaryInd].index; ++e3 ) {

                            if((*(*channels[edgeNum][c].backupsOnChannel[bup]).primaryPath).edges[e2] == p[primaryInd].edges[e3]) {
                                //There is a non-disjoint primary path on this channel, so it is unusable.

                                disjoint = false;
                            }
                        }
                    }
                }

                if(disjoint) {
                    //This channel is free
                    free = true;
                }
            }

            if(!free) {
                if(firstOpenChannel < MAX_CHANNELS) {
                    cost++;
                }else {
                    cost = 1000000;
                    break;
                }

            }

        }

        pathCosts[i] = cost;
    }
}

//TODO: There's some sketchiness going on with numPossiblePaths vs NUM_CONNECTIONS.
int determineCompatibleBackups(SimplePath *p, int *potPathInd, int numPossiblePaths, int pInd) {
    int numDisjoint = 0;
    int numConf = 0;
    //cout << "SRC: " << p[pInd].sourceNode << " DEST: " << p[pInd].destNode << "\n";
    //First pass checks to see which simple paths are disjoint from the primary path.
    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
      if(p[i].hops <= 0 || p[i].index < 0|| p[pInd].hops <= 0 || p[pInd].index < 0){numConf++; continue;}
        bool disjoint = true;
        //Check each edge to make sure they're disjoint
        for(int e1 = 0; disjoint && e1 <= p[pInd].index; ++e1) {
            for(int e2 = 0; disjoint && e2 <= p[i].index; ++e2) {
                if(p[i].edgeNums[e2] == p[pInd].edgeNums[e1]) {
                    disjoint = false;
		    numConf++;
                }
            }
        }
        if(disjoint) {
            potPathInd[numDisjoint] = i;
            numDisjoint++;
        }

    }
    //Mark the end of the array
    potPathInd[numDisjoint] = -1;
    //cout << "disjoint: " << numDisjoint << " out of " << numPossiblePaths <<"\n";
    //cout << "conflicts: " << numConf << "\n";
    return numDisjoint;
}

int computeAllSimplePathsN(SimplePath **ps, int *vertexList, Edge *edgeList, int sourceNode, int destNode, int hops) {
    int index = (sourceNode * N_NODES) + destNode;

    //initialize arrays
    int visited[N_NODES]; //visited[i] is 1 if node i has been visited on this path, 0 otherwise.
    int currentPath = 0;

    //edgeListIndex[i] contains the index into edgeList[] (aka the compact adj list) for node i.
    int edgeListIndex[N_NODES];

    ps[index][currentPath].index = 0;

    //Initialize our search components
    for(int i = 0; i < N_NODES; ++i) {
        visited[i] = 0;
        edgeListIndex[i] = vertexList[i];
    }

    stack <int> st;
    int currentNode;
    int neighbor;
    int currentHop = 1;

    st.push(sourceNode);
    visited[sourceNode] = 1;

    while(st.size() > 0) {
        //use loopCond to get to the beginning of the while loop from inside the for loop.
        bool loopCond = true;
        currentNode = st.top();
        //for each neighbor of currentNode
        for(; loopCond == true && edgeListIndex[currentNode] < vertexList[currentNode+1]; ++edgeListIndex[currentNode]) {
            neighbor = edgeList[edgeListIndex[currentNode]].v2;

            //If we're too far away from our source node, backtrack.
            if(currentHop >= hops) {
                break;
            }

            if(edgeList[edgeListIndex[currentNode]].load == MAX_CHANNELS) {
                continue;
            }

            //If our neighbor is the desired node, AND we're at the correct path length, save this path!
            if(neighbor == destNode && currentHop < hops) {

                ps[index][currentPath].edges[ps[index][currentPath].index] = &edgeList[edgeListIndex[currentNode]];
		ps[index][currentPath].edgeNums[ps[index][currentPath].index] = edgeList[edgeListIndex[currentNode]].edgeNum;
		
                ps[index][currentPath].sourceNode = sourceNode;
                ps[index][currentPath].destNode = destNode;
                ps[index][currentPath].hops = currentHop;

                //Copy the whole path up until the dest node to the next path in the array.
                //Note that we don't copy the COST from the current primary path, as the cost is computed
                //independently for each primary path.
                ps[index][currentPath+1].sourceNode = sourceNode;
                ps[index][currentPath+1].destNode = destNode;
                ps[index][currentPath+1].hops = currentHop;
                ps[index][currentPath+1].index = ps[index][currentPath].index-1;
                for(int i = 0; i < ps[index][currentPath].index; ++i) {
                    ps[index][currentPath+1].edges[i] = ps[index][currentPath].edges[i];
		    ps[index][currentPath+1].edgeNums[i] = ps[index][currentPath].edgeNums[i];
                }

                currentPath += 1;

                ps[index][currentPath].index += 1;
                ++edgeListIndex[currentNode];

                //
                loopCond = false;
                break;
            }

            if(!visited[neighbor]) {

                ps[index][currentPath].edges[ps[index][currentPath].index] = &edgeList[edgeListIndex[currentNode]];
		ps[index][currentPath].edgeNums[ps[index][currentPath].index] = edgeList[edgeListIndex[currentNode]].edgeNum;
                ps[index][currentPath].index += 1;

                st.push(neighbor);
                visited[neighbor] = 1;
                currentHop++;

                //continue the while loop, but increment the ELI first.
                ++edgeListIndex[currentNode];
                loopCond = false;
                break;
            }
        }

        if(loopCond) {
            currentHop--;

            //Once we've visited all of this node's neighbors, we reset it so that a
            //different path involving this node can be explored.
            visited[currentNode] = 0;
            ps[index][currentPath].index -= 1;

            edgeListIndex[currentNode] = vertexList[currentNode];
            st.pop();
        }

    }
    //Last path is invalid
    ps[index][currentPath].hops = 0;

    for(int j = currentPath; j < NUM_CONNECTIONS; ++j) {
      ps[index][j].hops = 0;
    }
    return currentPath;
}



void readGraphReorderEdgeList(int vertexList[],Edge compactEdgeList[2*N_EDGES],Edge reorderedEdgeList[2*N_NODES]) {
    //cout << "Beginning read\n";

    //TODO: We def don't need this extra array... please revise.
    int edgeList[N_NODES][N_NODES];
    for(int i = 0; i < N_NODES; ++i) {
        for(int j = 0; j < N_NODES; ++j) {
            edgeList[i][j] = 0;
        }
    }
    for(int i = 0; i < N_EDGES; ++i) {
        edgeList[base_edges[i][0]][base_edges[i][1]] = 1;
        edgeList[base_edges[i][1]][base_edges[i][0]] = 1;
    }

    int vDegree[N_NODES];

    int counter = 0;
    for(int i = 0; i < N_NODES; ++i) {
        vertexList[i] = counter;
        for(int j = 0; j < N_NODES; ++j) {
            if(edgeList[i][j] != 0) {
                compactEdgeList[counter].v1 = i;
                compactEdgeList[counter].v2 = j;
                compactEdgeList[counter].load = 0;
                compactEdgeList[counter].totalProtected = 0;
                compactEdgeList[counter].edgeNum = counter;

                //for(int x = 0; x < MAX_CHANNELS; ++x) {
                //    compactEdgeList[counter].channels[x].numBackups = 0;
                //}


                counter++;
            }
        }

        vDegree[i] = counter - vertexList[i];

        //cout << i << ": " << vDegree[i] << "\n";
    }
    vertexList[N_NODES] = 2*N_EDGES;

    //THis successfully reorders the edgelist based on the degree of the neighbor.
    //TODO: make this sorting algorithm faster... like WAY faster.
    for(int i = 0; i < N_NODES; ++i) {

        int startInd = vertexList[i];
        int endInd = vertexList[i+1];
        //[startInd,endInd)

        int reorderedInd = startInd;

        while(reorderedInd < endInd) {
            int min = startInd;
            int minVal = 66666; //min degree of the neighbor

            //Find the "smallest" neighbor of this node.
            for(int j = startInd; j < endInd; ++j) {

                bool isReordered = false;

                //Check to see if this node is already in our reordered list.
                for(int k = startInd; k < reorderedInd; ++k) {
                    if(reorderedEdgeList[k].v2 == compactEdgeList[j].v2) {
                        isReordered = true;
                        break;
                    }
                }

                //if its not in our reordered list and it qualifies as the minimum neighbor.
                if(isReordered == false && vDegree[compactEdgeList[j].v2] <= minVal) {
                    min = j;
                    minVal = vDegree[compactEdgeList[j].v2];
                }

            }

            reorderedEdgeList[reorderedInd].v1 = compactEdgeList[min].v1;
            reorderedEdgeList[reorderedInd].v2 = compactEdgeList[min].v2;
            reorderedEdgeList[reorderedInd].load = 0;
            reorderedInd++;
        }
    }
}

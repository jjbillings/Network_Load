
/*
 * File:   main.cpp
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
struct Connection2;
struct Channel;

struct Channel{
    bool primary; //is this channel used for a primary path?
    int numBackups; //total protected;
    Connection2 *backupsOnChannel[NUM_CONNECTIONS];//Realistically, there will be far fewer than NUM_CONNECTIONS
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

    Edge *edges[N_NODES];
};

struct Path {
    //TODO: Re-evaluate if we need source/dest nodes since we have them in the Connection Struct.
    int sourceNode;
    int destNode;

    int hops;
    int index;
    int cost;

    //Every path that uses a particular edge just has a reference to it (not a copy), so they can each manipulate it.
    Edge *edges[N_NODES];
    bool freeEdges[N_NODES]; //whether or not that edge has a cost of 0
    int channelNum[N_NODES]; //Channel number for each edge that it uses
    bool primary;
    bool active;
};

struct Connection2 {
    int sourceNode;
    int destNode;
    int combinedCost;
    bool validBackup;
    bool validPrimary;
    Path *backupPath;
    Path *primaryPath;
};

struct Connection {
	int sourceNode;
	int destNode;
    bool validBackup;
    bool validPrimary;
	Path backupPath;
	Path primaryPath;
};



void randomConnections(int vertexList[],Edge edgeList[2*N_EDGES],int sampleNum);
Connection2 allPrimaryBackupCombos(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode, Connection2 selectedConnection,Channel channels[2*N_EDGES][MAX_CHANNELS]);

bool single_connection_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int hops, int connectionNum, Connection conns[NUM_CONNECTIONS]);

bool computeBackupPath(int vertexList[], Edge edgeList[2*N_EDGES], Connection conns[NUM_CONNECTIONS], int connectionNum, int hops);
bool computePrimaryPath(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, Path *p);

int computeAllPrimaryPaths(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, Path *p[MAX_PATHS], Channel channels[2*N_EDGES][MAX_CHANNELS]);
int computeAllBackupPaths(int vertexList[], Edge edgeList[2*N_EDGES], Path *primaryPath, int hops, Connection2 *p[MAX_PATHS], Channel channels[2*N_EDGES][MAX_CHANNELS]);

int computeAllSimplePathsN(SimplePath **ps, int *vertexList, Edge *edgeList, int sourceNode, int destNode, int hops);
void simulate(int *vertexList, Edge *edgeList);
int determineCompatibleBackups(SimplePath *p, int *potPathInd, int numPossiblePaths, int pInd);
void computeCostForBackups(SimplePath *p, int *potPathInd, int numPotPaths, int backupIndex, int *pathCosts,Channel cs[2*N_EDGES][MAX_CHANNELS]);
void selectChannels(Connection2 *c, Channel chan[2*N_EDGES][MAX_CHANNELS]);
void increaseLoad(Connection2 *connection, Channel channels[2*N_EDGES][MAX_CHANNELS]);

void readGraph(int vertexList[],Edge compactEdgeList[2*N_EDGES]);
void readGraphReorderEdgeList(int vertexList[],Edge compactEdgeList[2*N_EDGES],Edge reorderedEdgeList[2*N_NODES]);
void increaseNetworkLoad(Connection2 *connection,Channel channels[2*N_EDGES][MAX_CHANNELS]);
void printConnection(Connection *connection);
void printPath(Path *path);
void exportNetworkLoad(Connection conns[NUM_CONNECTIONS],Edge edgeList[2*N_EDGES],int sampleNum, int numIncompleteConnections);

bool comparePath(const Path& p1, const Path& p2);

int vertexList[N_NODES+1];
Edge edgeList[2*N_EDGES];
Edge reorderedEdgeList[2*N_EDGES];
Connection2 cons[NUM_CONNECTIONS];
Channel channels[2*N_EDGES][MAX_CHANNELS];

/*
 *TODO: I totally thought I made the algorithm be based on BFS, but it is in fact based on DFS.
 *So REVERSE the order of the edge list. Currently, the neighbor with the lowest degree gets pushed
 *to the "bottom" of the stack, so we end up computing the path with high-degree nodes in it...
 */
int main(int argc, char** argv) {
    cout <<"Welcome to main\n";

    for(int f = 0; f < (2*N_EDGES); ++f){
        for(int g = 0; g < MAX_CHANNELS; ++g) {
            channels[f][g].numBackups = 0;
            channels[f][g].primary = false;
        }
    }
    readGraphReorderEdgeList(vertexList,edgeList,reorderedEdgeList);

    srand(time(NULL));

    simulate(vertexList,edgeList);

    /*Use this for computing all primary/backup combos
    for(int x = 0; x < 100; ++x) {
        int s = rand() % N_NODES;
        int d = rand() % N_NODES;
        //int s = 0;
        //int d = 9;
        while(s == d) {
            s = rand() % N_NODES;
            d = rand() % N_NODES;
        }

        cons[x] = allPrimaryBackupCombos(vertexList,edgeList,s,d,cons[x],channels);
        if(cons[x].validPrimary == false || cons[x].validBackup == false) {
            cout << "Unable to find a valid combo\n";
        }else {
            printPath(cons[x].primaryPath);
            printPath(cons[x].backupPath);
            increaseNetworkLoad(&cons[x],channels);
        }
        cout <<"";
    }


    //print out each edge, with its load and the number of backups protected on each channel.
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
    */
    return 0;
}

void simulate(int *vertexList, Edge *edgeList){
    int connectionNum = 0;
    //We want to compute and store all possible paths between our source and desitination.
    SimplePath **ps = new SimplePath*[N_NODES * N_NODES]; //Storage for paths
    int *npaths = new int[N_NODES*N_NODES];

    int v1[40] = {9, 5, 6, 1, 3, 5, 4, 9, 9, 9, 7, 8, 2, 10, 3, 5, 9, 3, 2, 3, 5, 2, 3, 3, 10, 9, 10, 2, 1, 1, 3, 2, 9, 5, 4, 6, 10, 5, 0, 1};
    int v2[40] = {3, 8, 4, 3, 8, 3, 7, 1, 5, 6, 0, 6, 10, 5, 8, 2, 3, 6, 5, 4, 2, 3, 9, 7, 9, 5, 6, 5, 0, 2, 5, 5, 10, 3, 9, 3, 4, 1, 10, 2};

    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
        ps[i] = new SimplePath[NUM_CONNECTIONS];
    }

    //cout <<"ps created\n";

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
    //cout << "all simple paths computed!\n";


    for(int num = 0; num < 40; ++num) {
    //Attempt to allocate SOME connection onto the network
    //int s = 0;
    //int d = 9;
    //int s = rand() % N_NODES;
    //int d = rand() % N_NODES;
    //while(s == d) {
    //    s = rand()%N_NODES;
    //    d = rand()%N_NODES;
    //}
    int s = v1[num];
    int d = v2[num];

    //Allocate storage for the potential primary/backup path combos
    int index = (s*N_NODES) + d;
    int numPossiblePaths = npaths[index];

    //Stores indices into the ps[index][] array for each disjoint backup path.
    //potPathInd[i][j] = k where ps[index][k] is a path that is edge-disjoint from ps[index][i].
    int ** potPathInd = new int*[numPossiblePaths];
    for(int i = 0; i < numPossiblePaths; ++i) {
        potPathInd[i] = new int[numPossiblePaths];
    }


    //Find all paths which are edge-disjoint from this primary.
    int k = -1;
    //On the GPU, instead of iterating i..numPossiblePaths, we would give thread_i backup_i
    for(int i = 0; i < numPossiblePaths; ++i) {
        k = determineCompatibleBackups(ps[index],potPathInd[i],numPossiblePaths,i);
        //cout << "Number of paths which are disjoint from this primary path: " << k << "\n";
    }



    //Compute the Cost for each backup path.
    int ** pathCosts = new int*[numPossiblePaths];
    for(int i = 0; i < numPossiblePaths; ++i) {
        pathCosts[i] = new int[numPossiblePaths];
    }

    for(int i = 0; i < numPossiblePaths; ++i) {
        computeCostForBackups(ps[index],potPathInd[i],numPossiblePaths,i,pathCosts[i],channels);
    }



    //Select cheapest connection
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
    increaseLoad(&cons[connectionNum],channels);


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

    for(int i = 0; i < (N_NODES*N_NODES); ++i) {
        delete[] ps[i];
    }
    delete[] ps;
    delete[] npaths;
    //cout << "ps and npaths deleted\n";
}

void increaseLoad(Connection2 *connection, Channel channels[2*N_EDGES][MAX_CHANNELS]) {
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
        channels[en][cNum].numBackups += 1;
        (*e).totalProtected +=1;
    }

}

//TODO: This method contains a lot of redundant code that is also in computeCostForBackups. Consider combining.
//I wanted to modularize the code as much as possible this time around, which is why there's so much redundancy in this method.
void selectChannels(Connection2 *c, Channel chan[2*N_EDGES][MAX_CHANNELS]) {

    cout << "prim\n";
    for(int i = 0; i <= (*(*c).primaryPath).index; ++i) {
        cout << (*(*(*c).primaryPath).edges[i]).v1 << " -> " << (*(*(*c).primaryPath).edges[i]).v2 << "\n";
    }
    cout << "back\n";
    for(int i = 0; i <= (*(*c).backupPath).index; ++i) {
        cout << (*(*(*c).backupPath).edges[i]).v1 << " -> " << (*(*(*c).backupPath).edges[i]).v2 << "\n";
    }

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
}

//TODO: Need to test once we actually start loading the network.
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

//TODO: Give each thread an index into the array of simple paths, and have them check to see if "their" path is compatible.
int determineCompatibleBackups(SimplePath *p, int *potPathInd, int numPossiblePaths, int pInd) {
    int numDisjoint = 0;
    //First pass checks to see which simple paths are disjoint from the primary path.
    for(int i = 0; i < numPossiblePaths; ++i) {
        if(i == pInd) {//TODO: Shouldn't even need this, they will clearly be edge-disjoint
            continue;
        }

        bool disjoint = true;
        //Check each edge to make sure they're disjoint
        for(int e1 = 0; disjoint && e1 <= p[pInd].index; ++e1) {
            for(int e2 = 0; disjoint && e2 <= p[i].index; ++e2) {
                if(p[i].edges[e2] == p[pInd].edges[e1]) {
                    disjoint = false;
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
    return currentPath;
}

void increaseNetworkLoad(Connection2 *connection, Channel channels[2*N_EDGES][MAX_CHANNELS]) {

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
        channels[en][cNum].numBackups += 1;
        (*e).totalProtected +=1;
    }

}

Connection2 allPrimaryBackupCombos(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode, Connection2 selectedConnection,Channel channels[2*N_EDGES][MAX_CHANNELS]) {
    Path *paths[MAX_PATHS];
    //Connection2 *cons = new Connection2[MAX_PATHS];
    Connection2 *conns[MAX_PATHS];

    for(int i = 0; i < MAX_PATHS; ++i) {
        //paths[i] = (struct Path*) malloc(sizeof(struct Path));
        paths[i] = new Path();
        (*paths[i]).index = 0;
        (*paths[i]).cost = 0;
    }

    //TODO: have computeAllPrimaryPaths store paths directly in cons[] array, eliminating the need for the paths[] array.
    int k = computeAllPrimaryPaths(vertexList,edgeList,sourceNode,destNode,N_NODES,paths,channels);
    if(k == 0) {
        cout << "NO PRIMARY PATHS POSSIBLE";
        return *(new Connection2());
    }

    for(int i = 0; i < k; ++i) {
        //cons[i] = (struct Connection2*) malloc(sizeof(struct Connection2));
        conns[i] = new Connection2();
        (*conns[i]).sourceNode = (*paths[i]).sourceNode;
        (*conns[i]).destNode = (*paths[i]).destNode;
        (*conns[i]).validBackup = false;
        (*conns[i]).validPrimary = true;
        (*conns[i]).primaryPath = paths[i];
        (*conns[i]).combinedCost = (*(*conns[i]).primaryPath).cost; //TODO: probs remove
    }


    //for each primary path, we need to compute all possible backup paths.
    //bps will store all possible backup paths for THIS Primary Path.
    Connection2 *bps[MAX_PATHS];
    for(int j = 0; j < k; ++j) {

        //Allocate storage for potential backup paths.
        for(int i = 0; i < MAX_PATHS; ++i) {
            //bps[i] = (struct Connection2*) malloc(sizeof(struct Connection2));
            bps[i] = new Connection2();
            //THey all have the same primary path
            (*bps[i]).sourceNode = (*paths[j]).sourceNode;
            (*bps[i]).destNode = (*paths[j]).destNode;
            (*bps[i]).validBackup = false;
            (*bps[i]).validPrimary = false;
            (*bps[i]).primaryPath = paths[j];
            //(*bps[i]).backupPath = (struct Path*) malloc(sizeof(struct Path));
            (*bps[i]).backupPath = new Path();
            (*(*bps[i]).backupPath).index = 0;
        }

        int bp = computeAllBackupPaths(vertexList, edgeList, paths[j], N_NODES, bps,channels);
        //cout << "NUM_BACKUP_PATHS: " << bp <<"\n";

        //Compute the cost of each backup path.
        for(int backup = 0; backup < bp; ++backup) {
            for(int i = 0; i <= (*(*bps[backup]).backupPath).index; ++i) {
                if((*(*bps[backup]).backupPath).freeEdges[i] == false) {
                    (*(*bps[backup]).backupPath).cost += 1;
                }
            }
            //cout << "BACKUP #" << backup << ": COST = " << (*(*bps[backup]).backupPath).cost << "\n";
            cout << "";
        }

        //Find the "Cheapest" backup path for THIS primary path
        int minCost = 100000;
        int minInd = -1;
        Path *cheapest;
        for(int backup = 0; backup < bp; ++backup) {
            if((*(*bps[backup]).backupPath).index > 0 && (*(*bps[backup]).backupPath).cost < minCost) {
                cheapest = (*bps[backup]).backupPath;
                minCost = (*cheapest).cost;
                minInd = backup;
            }
        }

        (*conns[j]).backupPath = cheapest;
        (*conns[j]).validBackup = true;

        (*conns[j]).combinedCost = (*(*conns[j]).backupPath).cost + (*(*conns[j]).primaryPath).cost;

        for(int back = 0; back < MAX_PATHS; ++back) {
            if(back != minInd) {
                delete bps[back];
            }
        }
    }

    //Select the cheapest primary/backup combo.
    int minCost = 100000;
    Connection2 cheapestCon;
    int ind = -1;
    cout << "PREPARING TO FIND THE CHEAPEST\n";
    for(int c = 0; c < k; ++c) {
        if((*conns[c]).combinedCost < minCost) {
            minCost = (*conns[c]).combinedCost;
            cheapestCon = *conns[c];
            ind = c;
        }
    }

    //TODO:THIS IS TESTING
    //int randomInt = rand() % k;
    //cheapestCon = *conns[randomInt];
    //ind = randomInt;
    //TODO:END TEST

    Connection2 ret = cheapestCon;
    for(int i = 0; i < k; ++i) {
        if(i != ind) {
            delete conns[i];
            delete paths[i];
        }
    }
    return ret;
}



int computeAllPrimaryPaths(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, Path *p[MAX_PATHS], Channel channels[2*N_EDGES][MAX_CHANNELS]) {
    //cout << "Preparing to compute all primary paths\n";
    //initialize arrays
    int visited[N_NODES]; //visited[i] is 1 if node i has been visited on this path, 0 otherwise.
    int currentPath = 0;

    //edgeListIndex[i] contains the index into edgeList[] (aka the compact adj list) for node i.
    int edgeListIndex[N_NODES];

    (*p[currentPath]).index = 0;

    //Initialize our search components
    for(int i = 0; i < N_NODES; ++i) {
        visited[i] = 0;
        edgeListIndex[i] = vertexList[i];
    }

    stack <int> st;
    int currentNode;
    int neighbor;
    int currentHop = 0;

    st.push(sourceNode);
    visited[sourceNode] = 1;

    LOOP:
    while(st.size() > 0) {
        currentNode = st.top();
        //for each neighbor of currentNode
        for(; edgeListIndex[currentNode] < vertexList[currentNode+1]; ++edgeListIndex[currentNode]) {
            neighbor = edgeList[edgeListIndex[currentNode]].v2;

            //If we're too far away from our source node, backtrack.
            if(currentHop >= hops) {
                goto NEXT_NODE;
            }

            if(edgeList[edgeListIndex[currentNode]].load == MAX_CHANNELS) {
                //TODO: For now we just continue to the next neighbor
                //cout << "CHANNELS ARE MAXED OUT\n";
                continue;
            }

            //If our neighbor is the desired node, AND we're at the correct path length, save this path!
            if(neighbor == destNode && currentHop < hops) {

                (*p[currentPath]).edges[(*p[currentPath]).index] = &edgeList[edgeListIndex[currentNode]];


                for(int i = 0; i <= (*p[currentPath]).index; ++i) {
                    (*p[currentPath]).cost += 1;

                    //(*p[currentPath]).channelNum[i] = (*(*p[currentPath]).edges[i]).load;
                    int en = (*(*p[currentPath]).edges[i]).edgeNum;
                    int ind = 0;
                    (*p[currentPath]).channelNum[i] = -1;
                    for(int f = 0; f < MAX_CHANNELS; ++f) {
                        if(channels[en][f].numBackups == 0) {
                            (*p[currentPath]).channelNum[i] = f;
                            goto GOOD;
                        }
                    }
                    cout <<"MAJOR ERRORRRRRRRRR\n";
                    GOOD:
                    cout <<"";
                }

                (*p[currentPath]).primary = true;
                (*p[currentPath]).sourceNode = sourceNode;
                (*p[currentPath]).destNode = destNode;
                (*p[currentPath]).hops = hops;

                //Copy the whole path up until the dest node to the next path in the array.
                //Note that we don't copy the COST from the current primary path, as the cost is computed
                //independently for each primary path.
                (*p[currentPath+1]).sourceNode = sourceNode;
                (*p[currentPath+1]).destNode = destNode;
                (*p[currentPath+1]).hops = hops;
                (*p[currentPath+1]).index = (*p[currentPath]).index-1;
                for(int i = 0; i < (*p[currentPath]).index; ++i) {
                    (*p[currentPath+1]).edges[i] = (*p[currentPath]).edges[i];
                }

                currentPath += 1;

                (*p[currentPath]).index += 1;
                ++edgeListIndex[currentNode];
                goto LOOP;
            }

            if(!visited[neighbor]) {

                (*p[currentPath]).edges[(*p[currentPath]).index] = &edgeList[edgeListIndex[currentNode]];
                (*p[currentPath]).index += 1;

                st.push(neighbor);
                visited[neighbor] = 1;
                currentHop++;

                //continue the while loop, but increment the ELI first.
                ++edgeListIndex[currentNode];
                goto LOOP;
            }
        }

        NEXT_NODE:
        currentHop--;

        //Once we've visited all of this node's neighbors, we reset it so that a
        //different path involving this node can be explored.
        visited[currentNode] = 0;
        (*p[currentPath]).index -= 1;

        edgeListIndex[currentNode] = vertexList[currentNode];
        st.pop();
    }
    return currentPath;
}


int computeAllBackupPaths(int vertexList[], Edge edgeList[2*N_EDGES], Path *primaryPath, int hops, Connection2 *p[MAX_PATHS], Channel channels[2*N_EDGES][MAX_CHANNELS]) {
    //cout << "Preparing to compute all backup paths\n";
    //cout << "PP: S: " << (*primaryPath).sourceNode << " D: " << (*primaryPath).destNode << "\n";

    //initialize arrays
    int visited[N_NODES]; //visited[i] is 1 if node i has been visited on this path, 0 otherwise.
    int currentPath = 0;

    //edgeListIndex[i] contains the index into edgeList[] (aka the compact adj list) for node i.
    int edgeListIndex[N_NODES];

    (*(*p[currentPath]).backupPath).index = 0;

    //Initialize our search components
    for(int i = 0; i < N_NODES; ++i) {
        visited[i] = 0;
        edgeListIndex[i] = vertexList[i];
    }

    stack <int> st;
    int currentNode;
    int neighbor;
    int currentHop = 0;

    st.push((*primaryPath).sourceNode);
    visited[(*primaryPath).sourceNode] = 1;

    LOOP:
    while(st.size() > 0) {
        currentNode = st.top();
        //for each neighbor of currentNode
        for(; edgeListIndex[currentNode] < vertexList[currentNode+1]; ++edgeListIndex[currentNode]) {
            neighbor = edgeList[edgeListIndex[currentNode]].v2;
            int firstOpenChannel = MAX_CHANNELS+1;


            //If we're too far away from our source node, backtrack.
            if(currentHop >= hops) {
                goto NEXT_NODE;
            }

            //Ensure that that we do not use ANY edges that the primary path uses.
            for(int k = 0; k <= (*primaryPath).index; ++k) {
                if(&edgeList[edgeListIndex[currentNode]] == (*primaryPath).edges[k]) {
                    goto END_OF_LOOP;
                }
            }

            (*(*p[currentPath]).backupPath).freeEdges[(*(*p[currentPath]).backupPath).index] = false;
            (*(*p[currentPath]).backupPath).channelNum[(*(*p[currentPath]).backupPath).index] = -1;

            //For every channel on this edge, if every backup path allocated for this channel
            //is disjoint from current backup path candidate, then this edge is free!
            //for each channel on this edge
            for(int ch = 0; ch < MAX_CHANNELS; ++ch) {

                CHANNEL_START:
                if(channels[edgeListIndex[currentNode]][ch].primary == true) {
                    continue;
                }

                //At this point, we know that there are no primary paths on this channel
                //Thus we must check and see if it is "free".


                //we COULD use this channel, but there may be a "free" one further down.
                if(channels[edgeListIndex[currentNode]][ch].numBackups == 0) {
                    if(ch < firstOpenChannel) {
                        firstOpenChannel = ch;
                    }
                    continue;
                }

                //Check every connection currently on protected on the channel
                for(int bup = 0; bup < channels[edgeListIndex[currentNode]][ch].numBackups; ++bup) {

                    //At this point, we know that there is at least one path protected on this channel.
                    //Technically, we should also know that it's not a primary path.

                    //for each edge of the protected connection's primary path
                    for(int e = 0; e <= (*(*channels[edgeListIndex[currentNode]][ch].backupsOnChannel[bup]).primaryPath).index; ++e) {

                        //see if its the same edge as used by our primary path.
                        for(int te = 0; te <=(*(*p[currentPath]).primaryPath).index; ++te ) {

                            if((*(*channels[edgeListIndex[currentNode]][ch].backupsOnChannel[bup]).primaryPath).edges[e] == (*(*p[currentPath]).primaryPath).edges[te]) {
                                //There is a non-disjoint primary path on this channel, so it is unusable.
                                goto CHANNEL_LOOP_END;
                            }
                        }
                    }
                }

                //If every protected connection is primary disjoint
                //Then this channel is free af.
                (*(*p[currentPath]).backupPath).freeEdges[(*(*p[currentPath]).backupPath).index] = true;
                (*(*p[currentPath]).backupPath).channelNum[(*(*p[currentPath]).backupPath).index] = ch;
                goto CHANNEL_END;

                CHANNEL_LOOP_END:
                cout <<"";
            }

            CHANNEL_END:
            //Here we check to see if we were unable to find a "free" path
            if((*(*p[currentPath]).backupPath).freeEdges[(*(*p[currentPath]).backupPath).index] == false) {
                if(firstOpenChannel < MAX_CHANNELS) {
                    (*(*p[currentPath]).backupPath).channelNum[(*(*p[currentPath]).backupPath).index] = firstOpenChannel;
                }else {
                    //We could not find any open channel.
                    goto END_OF_LOOP;
                }
            }


            //If our neighbor is the desired node, AND we're at the correct path length, save this path!
            if(neighbor == (*primaryPath).destNode && currentHop < hops) {

                (*(*p[currentPath]).backupPath).edges[(*(*p[currentPath]).backupPath).index] = &edgeList[edgeListIndex[currentNode]];

                (*(*p[currentPath]).backupPath).sourceNode = (*primaryPath).sourceNode;
                (*(*p[currentPath]).backupPath).destNode = (*primaryPath).destNode;
                (*(*p[currentPath]).backupPath).hops = hops;
                (*(*p[currentPath]).backupPath).primary = false;

                //Make sure we don't use the same path for the primary and backup paths.
                if(comparePath((*(*p[currentPath]).backupPath),*primaryPath)) {
                    goto NEXT_NODE;
                }

                //Copy the whole path up until the dest node to the next path in the array.
                (*(*p[currentPath+1]).backupPath).sourceNode = (*primaryPath).sourceNode;
                (*(*p[currentPath+1]).backupPath).destNode = (*primaryPath).destNode;
                (*(*p[currentPath+1]).backupPath).hops = hops;
                (*(*p[currentPath+1]).backupPath).index = (*(*p[currentPath]).backupPath).index-1;

                for(int i = 0; i < (*(*p[currentPath]).backupPath).index; ++i) {
                    (*(*p[currentPath+1]).backupPath).edges[i] = (*(*p[currentPath]).backupPath).edges[i];
                    (*(*p[currentPath+1]).backupPath).freeEdges[i] = (*(*p[currentPath]).backupPath).freeEdges[i];
                    (*(*p[currentPath+1]).backupPath).channelNum[i] = (*(*p[currentPath]).backupPath).channelNum[i];
                }

                currentPath += 1;

                (*(*p[currentPath]).backupPath).index += 1;
                ++edgeListIndex[currentNode];
                goto LOOP;
            }

            if(!visited[neighbor]) {

                (*(*p[currentPath]).backupPath).edges[(*(*p[currentPath]).backupPath).index] = &edgeList[edgeListIndex[currentNode]];
                (*(*p[currentPath]).backupPath).index += 1;

                st.push(neighbor);
                visited[neighbor] = 1;
                currentHop++;

                //continue the while loop, but increment the ELI first.
                ++edgeListIndex[currentNode];
                goto LOOP;
            }
            END_OF_LOOP:
            cout <<"";
        }

        NEXT_NODE:
        currentHop--;

        //Once we've visited all of this node's neighbors, we reset it so that a
        //different path involving this node can be explored.
        visited[currentNode] = 0;
        (*(*p[currentPath]).backupPath).index -= 1;

        edgeListIndex[currentNode] = vertexList[currentNode];
        st.pop();
    }
    return currentPath;
}

bool computeBackupPath(int vertexList[], Edge edgeList[2*N_EDGES], Connection conns[NUM_CONNECTIONS], int connectionNum, int hops) {
    //cout << "-----------COMPUTING BACKUP PATH ---------------\n";
    int visited[N_NODES]; //visited[i] is 1 if node i has been visited on this path, 0 otherwise.
    int edgeListIndex[N_NODES];
    Path hopefulPath;
    hopefulPath.index = 0;

    //Initialize our search components
    for(int i = 0; i < N_NODES; ++i) {
        visited[i] = 0;
        edgeListIndex[i] = vertexList[i];
    }

    stack <int> st;
    int currentNode;
    int neighbor;
    int currentHop = 0;

    st.push(conns[connectionNum].primaryPath.sourceNode);
    visited[conns[connectionNum].primaryPath.sourceNode] = 1;

    LOOP:
    while(st.size() > 0) {

        currentNode = st.top();
        //for each neighbor of currentNode
        for(; edgeListIndex[currentNode] < vertexList[currentNode+1]; ++edgeListIndex[currentNode]) {
            neighbor = edgeList[edgeListIndex[currentNode]].v2;

            //Ensure that that we do not use ANY edges that the primary path uses.
            for(int k = 0; k <= conns[connectionNum].primaryPath.index; ++k) {
                if(&edgeList[edgeListIndex[currentNode]] == conns[connectionNum].primaryPath.edges[k]) {
                    goto END_OF_LOOP;
                }
            }

            //If we're too far away from our source node, backtrack.
            if(currentHop >= hops) {
                goto NEXT_NODE;
            }

            //if this edge is at max capacity, we would want to check and see if we can share a channel with one of the paths.
            if(edgeList[edgeListIndex[currentNode]].load == MAX_CHANNELS) {
                //For now we just continue to the next neighbor
                //cout << "CHANNELS ARE MAXED OUT\n";
                continue;
            }

            if(neighbor == conns[connectionNum].primaryPath.destNode && currentHop != hops-1) {
                continue;
            }

            //If our neighbor is the desired node, print out the current path.
            if(neighbor == conns[connectionNum].primaryPath.destNode && currentHop == hops-1) {
                visited[neighbor] = 1;

                hopefulPath.edges[hopefulPath.index] = &edgeList[edgeListIndex[currentNode]];
                hopefulPath.sourceNode = conns[connectionNum].primaryPath.sourceNode;
                hopefulPath.destNode = conns[connectionNum].primaryPath.destNode;
                hopefulPath.hops = currentHop+1;

                //Make sure we don't use the same path for the primary and backup paths.
                if(comparePath(hopefulPath,conns[connectionNum].primaryPath)) {
                    goto NEXT_NODE;
                }

                conns[connectionNum].backupPath = hopefulPath;

                //Debugging
                for(int j = 0; j < hopefulPath.index; ++j) {
                    //cout << (*hopefulPath.edges[j]).v1 << " -> ";
                }
                //cout << (*hopefulPath.edges[hopefulPath.index]).v1 << " -> " << (*hopefulPath.edges[hopefulPath.index]).v2 << "\n\n";
                return true;
            }

            if(!visited[neighbor]) {

                st.push(neighbor);
                visited[neighbor] = 1;
                currentHop++;

                hopefulPath.edges[hopefulPath.index] = &edgeList[edgeListIndex[currentNode]];
                hopefulPath.index += 1;

                //continue the while loop, but increment the ELI first.
                ++edgeListIndex[currentNode];
                goto LOOP;
            }

            END_OF_LOOP:
            cout << "";
        }

        NEXT_NODE:
        currentHop--;

        //Once we've visited all of this node's neighbors, we reset it so that a
        //different path involving this node can be explored.
        visited[currentNode] = 0;
        hopefulPath.index -= 1;

        edgeListIndex[currentNode] = vertexList[currentNode];
        st.pop();
    }
    //cout << "COULD NOT COMPUTE A BACKUP PATH\n";
    //TODO: We need a graceful way to handle not being able to compute a backup path.
    return false;
}


/**
HELPER FUNCTIONS
**/
bool comparePath(const Path& p1, const Path& p2)
{
    //cout << p1.sourceNode << " = " << p2.sourceNode << "\n";
    if((p1).sourceNode != (p2).sourceNode) {
        return false;
    }
    //cout << p1.destNode << " = " << p2.destNode << "\n";
    if((p1).destNode != (p2).destNode) {
        return false;
    }
    //cout << p1.hops << " = " << p2.hops << "\n";
    if((p1).hops != (p2).hops || (p1).index != (p2).index) {
        return false;
    }

    for(int i = 0; i <= p1.index; ++i) {

        if(p1.edges[i] != p2.edges[i]) {
            return false;
        }
    }

    return true;
}

/**
INPUT FUNCTIONS
**/
void readGraph(int vertexList[], Edge compactEdgeList[2*N_EDGES]) {
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

    int counter = 0;
    for(int i = 0; i < N_NODES; ++i) {
        vertexList[i] = counter;
        for(int j = 0; j < N_NODES; ++j) {
            if(edgeList[i][j] != 0) {
                compactEdgeList[counter].v1 = i;
                compactEdgeList[counter].v2 = j;
                compactEdgeList[counter].load = 0;
                compactEdgeList[counter].edgeNum = counter;
                counter++;
            }
        }
    }
    vertexList[N_NODES] = 2*N_EDGES;
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

void randomConnections(int vertexList[],Edge edgeList[2*N_EDGES],int sampleNum) {
    Connection conns[NUM_CONNECTIONS];

    int numIncompleteConnections = 0; //Number of connections that couldn't be allocated completely.

    //Generate NUM_CONNECTIONS random connections
    for(int i = 0; i < NUM_CONNECTIONS; ++i) {

        //Choose random source and destination nodes
        int v1 = rand() % N_NODES;
        int v2 = rand() % N_NODES;

        //Make sure we're not trying to go to from v_i to v_i... because that's boring.
        if(v1 == v2) {
            i--;
            continue;
        }

        conns[i].sourceNode = v1;
        conns[i].destNode = v2;

        if(single_connection_N_hops(vertexList,edgeList,10,i,conns) == false) {
            numIncompleteConnections++;
        }
    }

    //Prints the current load for every edge of the graph.
    for(int i = 0; i < 2*N_EDGES; ++i) {
        cout << edgeList[i].v1 << " -> " << edgeList[i].v2 << " | LOAD: " << edgeList[i].load << "\n";
    }

    cout << "Number of incomplete connections: " << numIncompleteConnections << "\n";
    cout << "Number of complete connections: " << (NUM_CONNECTIONS - numIncompleteConnections) << "\n";

    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
        printConnection(&conns[i]);
    }
    //exportNetworkLoad(conns,edgeList,sampleNum,numIncompleteConnections);
}


bool single_connection_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int hops, int connectionNum, Connection conns[NUM_CONNECTIONS]) {
    printf("-------DEFINING_CONNECTION_%d-------\n",connectionNum);

    //Allocate the primary path (Make it the shortest path for now).
    //We start with paths that are 1 hop away, continuing with longer paths until we successfully find a valid one.
    for(int hop = 1; hop < hops; ++hop) {
        bool pathFound = computePrimaryPath(vertexList,edgeList,conns[connectionNum].sourceNode,conns[connectionNum].destNode,hop,&conns[connectionNum].primaryPath);
        if(pathFound == true) {
            //This becomes the primary path for this connection.
            conns[connectionNum].primaryPath.primary = true;
            conns[connectionNum].primaryPath.active = true;
            conns[connectionNum].validPrimary = true;
            break;
        }else {
            conns[connectionNum].validPrimary = false;//TODO: See comment below.
        }
    }

    //Allocate the backup path
    //We know that we didn't find a path < primaryPath.hops.
    for(int hop = conns[connectionNum].primaryPath.hops; hop < hops; ++hop) {
        bool pathFound = computeBackupPath(vertexList,edgeList,conns,connectionNum,hop);
        //TODO: Find a graceful way to handle not being able to compute a backup path.
        if(pathFound == true) {
            conns[connectionNum].backupPath.primary = false;
            conns[connectionNum].backupPath.active = false;
            conns[connectionNum].validBackup = true;
            break;
        }else {
            conns[connectionNum].validBackup = false; //TODO: I think this is sloppy. Consider making a default ctor for the struct
        }
    }

    if(conns[connectionNum].validPrimary == true && conns[connectionNum].validBackup == true) {
        return true;
    }

    return false;
}


/**
    Allocate primary path for the connection. Finds a path from SN -> DN of EXACTLY "hops" hops.
**/
bool computePrimaryPath(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode, int hops, Path *p) {
    //printf("-------SINGLE_PATH_%d_HOPS-------\n",hops);
    //printf("-------|%d -> %d|----------------\n",sourceNode,destNode);

    //initialize arrays
    int visited[N_NODES]; //visited[i] is 1 if node i has been visited on this path, 0 otherwise.

    //edgeListIndex[i] contains the index into edgeList[] (aka the compact adj list) for node i.
    int edgeListIndex[N_NODES];

    (*p).index = 0;

    //Initialize our search components
    for(int i = 0; i < N_NODES; ++i) {
        visited[i] = 0;
        edgeListIndex[i] = vertexList[i];
    }

    stack <int> st;
    int currentNode;
    int neighbor;
    int currentHop = 0;

    st.push(sourceNode);
    visited[sourceNode] = 1;

    LOOP:
    while(st.size() > 0) {

        currentNode = st.top();
        //for each neighbor of currentNode
        for(; edgeListIndex[currentNode] < vertexList[currentNode+1]; ++edgeListIndex[currentNode]) {
            neighbor = edgeList[edgeListIndex[currentNode]].v2;

            //If we're too far away from our source node, backtrack.
            if(currentHop >= hops) {
                goto NEXT_NODE;
            }

            //if this edge is at max capacity (i.e. has no free channels),
            // we would want to check and see if we can share a channel with one of the paths.
            if(edgeList[edgeListIndex[currentNode]].load == MAX_CHANNELS) {
                //TODO: For now we just continue to the next neighbor
                //cout << "CHANNELS ARE MAXED OUT\n";
                continue;
            }

            //Check to see if we're at the correct path length.
            if(neighbor == destNode && currentHop != hops-1) {
                continue;
            }

            //If our neighbor is the desired node, AND we're at the correct path length, save this path!
            if(neighbor == destNode && currentHop == hops-1) {
                visited[neighbor] = 1;

                (*p).edges[(*p).index] = &edgeList[edgeListIndex[currentNode]];

                //Now that we have the path set, increase the load on each edge.
                for(int i = 0; i <= (*p).index; ++i) {
                    //(*(*p).edges[i]).load++;
                }

                (*p).sourceNode = sourceNode;
                (*p).destNode = destNode;
                (*p).hops = hops;
                (*p).primary = true;
                return true;
            }

            if(!visited[neighbor]) {

                (*p).edges[(*p).index] = &edgeList[edgeListIndex[currentNode]];
                (*p).index += 1;

                st.push(neighbor);
                visited[neighbor] = 1;
                currentHop++;

                //continue the while loop, but increment the ELI first.
                ++edgeListIndex[currentNode];
                goto LOOP;
            }
        }

        NEXT_NODE:
        currentHop--;

        //Once we've visited all of this node's neighbors, we reset it so that a
        //different path involving this node can be explored.
        visited[currentNode] = 0;

        (*p).index -= 1;

        edgeListIndex[currentNode] = vertexList[currentNode];
        st.pop();
    }
    return false;
}



/**
OUTPUT FUNCTIONS
**/
void printConnection(Connection *connection) {
    if((*connection).validPrimary == true) {
        //Debug output for development.
        cout << "PRINTING CONNECTION\n";
        printf("PRIMARY - Hops: %d, Index: %d\n",(*connection).primaryPath.hops,(*connection).primaryPath.index);
        printf("Source: %d, Dest: %d\n",(*connection).sourceNode, (*connection).destNode);
        for(int j = 0; j < (*connection).primaryPath.index; ++j) {
            cout << (*(*connection).primaryPath.edges[j]).v1 << " -> ";
        }
        cout << (*(*connection).primaryPath.edges[(*connection).primaryPath.index]).v1 << " -> " << (*(*connection).primaryPath.edges[(*connection).primaryPath.index]).v2 << "\n\n";
    }else {
        cout << "COULDN'T ALLOCATE A VALID PRIMARY PATH\n";
    }


    if((*connection).validBackup == true) {
        printf("BACKUP - Hops: %d, Index: %d\n",(*connection).backupPath.hops,(*connection).backupPath.index);
        printf("Source: %d, Dest: %d\n",(*connection).sourceNode, (*connection).destNode);
        for(int j = 0; j < (*connection).backupPath.index; ++j) {
            cout << (*(*connection).backupPath.edges[j]).v1 << " -> ";
        }
        cout << (*(*connection).backupPath.edges[(*connection).backupPath.index]).v1 << " -> " << (*(*connection).backupPath.edges[(*connection).backupPath.index]).v2 << "\n\n";
    }else {
        cout << "UNABLE TO ALLOCATE A VALID BACKUP PATH\n";
    }
}

void printPath(Path *path) {
    if((*path).primary == true) {
        cout << "Primary Path ";
    }else {
        cout << "Backup Path ";
    }
    cout << "| Cost = " << (*path).cost << "\n";

    for(int j = 0; j < (*path).index; ++j) {
        cout << (*(*path).edges[j]).v1 << " -> ";
    }
    cout << (*(*path).edges[(*path).index]).v1 << " -> " << (*(*path).edges[(*path).index]).v2 << "\n\n";
}

/*
void exportNetworkLoad(Connection conns[NUM_CONNECTIONS],Edge edgeList[2*N_EDGES],int sampleNum, int numIncompleteConnections) {
    ofstream outputFile;
    string filename ("./load_per_edge/load_per_edge_");
    filename += to_string(sampleNum);
    filename += ".csv";

    cout << filename << "\n";
    outputFile.open(filename,ios_base::app);
    outputFile << numIncompleteConnections << "\n";
    for(int i = 0; i < 2*N_EDGES; ++i) {
        outputFile << i << "," << edgeList[i].load << "\n";

    }
    outputFile.close();
    for(int j = 0; j < NUM_CONNECTIONS; ++j) {
        conns[j] = Connection();
    }
}*/

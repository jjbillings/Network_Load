
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

#define NUM_CONNECTIONS 1000
#define MAX_CHANNELS 30
#define SAMPLES 1

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

int computeAllPrimaryPaths(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, Path *p[MAX_PATHS]);
int computeAllBackupPaths(int vertexList[], Edge edgeList[2*N_EDGES], Path *primaryPath, int hops, Connection2 *p[MAX_PATHS], Channel channels[2*N_EDGES][MAX_CHANNELS]);

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

    //TODO: Currently, some backups seem to be "creating" new channels. Restrict max channels
    for(int x = 0; x < 20; ++x) {
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

    }


    //print out each edge, with its load and the number of backups protected on each channel.
    for(int m = 0; m < 2*N_EDGES; ++m) {
        cout << "LOAD: " << edgeList[m].v1 << " -> " << edgeList[m].v2 << ": " << edgeList[m].load << " | ";
        if(edgeList[m].load > 0) {
            for(int c = 0; c < edgeList[m].load; ++c) {
                cout << "C" << c << ": " << channels[m][c].numBackups << " ";
                if(channels[m][c].primary == true) {
                    cout << "PRIMARY | ";
                }
            }
        }
        cout << "\n";

    }

    return 0;
}

void increaseNetworkLoad(Connection2 *connection, Channel channels[2*N_EDGES][MAX_CHANNELS]) {

    if((*(*connection).primaryPath).index <= 0) {
        cout << "Primary Path DNE?\n";
        return;
    }
    //Increment the network load; put the backup on the channels

    //Here we are incrementing the network load for the PRIMARY PATH
    for(int i = 0; i <= (*(*connection).primaryPath).index; ++i) {
        //Every edge in the primary path gets its load increased
        channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*(*connection).primaryPath).edges[i]).load].primary = true;
        channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*(*connection).primaryPath).edges[i]).load].backupsOnChannel[0] = connection;
        channels[(*(*(*connection).primaryPath).edges[i]).edgeNum][(*(*(*connection).primaryPath).edges[i]).load].numBackups += 1;
        (*(*(*connection).primaryPath).edges[i]).load += 1;
    }

    //Here we are increasing the network load for the BACKUP PATH
    for(int i = 0; i <= (*(*connection).backupPath).index; ++i) {
        //Temp
        Edge *e = (*(*connection).backupPath).edges[i];
        int cNum = (*(*connection).backupPath).channelNum[i];

        //first path to use this channel, or this is not a free edge for the backup path.
        if(channels[(*e).edgeNum][cNum].numBackups == 0 || (*(*connection).backupPath).freeEdges[i] == false) {
            (*e).load += 1;
        }

        //Marks that the connection is protected on this channel.
        int en = (*e).edgeNum;
        int numbs = channels[en][cNum].numBackups;
        channels[en][cNum].primary = false;
        channels[en][cNum].backupsOnChannel[numbs] = connection;
        channels[en][cNum].numBackups += 1;
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
    int k = computeAllPrimaryPaths(vertexList,edgeList,sourceNode,destNode,N_NODES,paths);
    if(k == 0) {
        //cout << "NO PRIMARY PATHS POSSIBLE";
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
        Path *cheapest;
        for(int backup = 0; backup < bp; ++backup) {
            if((*(*bps[backup]).backupPath).index > 0 && (*(*bps[backup]).backupPath).cost < minCost) {
                cheapest = (*bps[backup]).backupPath;
                minCost = (*cheapest).cost;
            }
        }

        (*conns[j]).backupPath = cheapest;
        (*conns[j]).validBackup = true;

        (*conns[j]).combinedCost = (*(*conns[j]).backupPath).cost + (*(*conns[j]).primaryPath).cost;
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

    Connection2 ret = cheapestCon;
    for(int i = 0; i < k; ++i) {
        if(i != ind) {
            delete conns[i];
            delete paths[i];
        }
    }
    return ret;
}



int computeAllPrimaryPaths(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, Path *p[MAX_PATHS]) {
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
            bool channelProb = false;
            int cNum = 0;


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

            //For every channel on this edge, if every backup path allocated for this channel
            //is disjoint from current backup path candidate, then this edge is free!
            //for each channel on this edge
            for(int ch = 0; ch < edgeList[edgeListIndex[currentNode]].load; ++ch) {

                CHANNEL_START:
                if(channels[edgeListIndex[currentNode]][ch].primary == true) {
                    ch+=1;
                    if(ch < edgeList[edgeListIndex[currentNode]].load) {
                        goto CHANNEL_START;
                    }else {
                        goto CHANNEL_END;
                    }
                }

                //If there are no paths (primary or backup) allocated on this channel, then it costs 1, since now
                //a primary cannot use this one.
                if(channels[edgeListIndex[currentNode]][ch].numBackups == 0) {
                    channelProb = true;
                    cNum = ch;
                    goto CHANNEL_END;
                }

                //Check every connection currently on protected on the channel
                for(int bup = 0; bup < channels[edgeListIndex[currentNode]][ch].numBackups; ++bup) {

                    //for each edge of the protected connection's primary path
                    for(int e = 0; e <= (*(*channels[edgeListIndex[currentNode]][ch].backupsOnChannel[bup]).primaryPath).index; ++e) {
                        //see if its the same edge as used by our primary path.
                        for(int te = 0; te <=(*(*p[currentPath]).primaryPath).index; ++te ) {

                            if((*(*channels[edgeListIndex[currentNode]][ch].backupsOnChannel[bup]).primaryPath).edges[e] == (*(*p[currentPath]).primaryPath).edges[te]) {
                                //cout << "NON-DISJOINT PRIMARY PATHS\n";
                                channelProb = true;
                                //If we run into a problem, go to the next channel.
                                ch+=1;
                                if(ch < edgeList[edgeListIndex[currentNode]].load) {
                                    goto CHANNEL_START;
                                }else {
                                    cNum = ch;
                                    goto CHANNEL_END;
                                }

                            }
                        }
                    }

                }
                //If every protected connection is primary disjoint
                //Then this channel is free af.
                cNum = ch;
                goto CHANNEL_END;
            }

            CHANNEL_END:
            if(channelProb == false) {
                (*(*p[currentPath]).backupPath).freeEdges[(*(*p[currentPath]).backupPath).index] = true;
                (*(*p[currentPath]).backupPath).channelNum[(*(*p[currentPath]).backupPath).index] = cNum; // save which channel we are using for this link
            }else {
                (*(*p[currentPath]).backupPath).freeEdges[(*(*p[currentPath]).backupPath).index] = false;
                (*(*p[currentPath]).backupPath).channelNum[(*(*p[currentPath]).backupPath).index] = cNum;
            }

            //if this edge is at max capacity (i.e. has no free channels),
            // we would want to check and see if we can share a channel with one of the paths.
            if(edgeList[edgeListIndex[currentNode]].load == MAX_CHANNELS) {
                //cout << "CHANNELS ARE MAXED OUT, CAN WE SHARE?\n";
                if(channelProb == true) {
                    continue;
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
                    (*(*p).edges[i]).load++;
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
}

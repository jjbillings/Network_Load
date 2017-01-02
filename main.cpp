
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

#define NUM_CONNECTIONS 800
#define MAX_CHANNELS 50



struct Edge {
    int v1;
    int v2;
    int load; //load <= maxChannels
    int maxChannels;
    //Path *channels[MAX_CHANNELS];
};

struct Path {
    //TODO: Re-evaluate if we need source/dest nodes since we have them in the Connection Struct.
    int sourceNode;
    int destNode;

    int hops;
    int index;
    //Every path that uses a particular edge just has a reference to it (not a copy), so they can each manipulate it.
    Edge *edges[N_NODES];
    bool primary;
    bool active;

};

struct Connection {
	int sourceNode;
	int destNode;
    bool validBackup;
    bool validPrimary;
	Path backupPath;
	Path primaryPath;
	int channel; //IF we want to stick to one channel.
};



void randomConnections(int vertexList[],Edge edgeList[2*N_EDGES],int sampleNum);
//bool single_path_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, int paths[NUM_CONNECTIONS][N_NODES], int connectionNum, Path structPaths[NUM_CONNECTIONS]);
bool computeBackupPath(int vertexList[], Edge edgeList[2*N_EDGES], Connection conns[NUM_CONNECTIONS], int connectionNum, int hops);

bool single_path_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode,int hops, Path *p);
void single_connection_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int hops, int connectionNum, Connection conns[NUM_CONNECTIONS]);


void readGraph(int vertexList[],Edge compactEdgeList[2*N_EDGES]);
void printPath(int sourceNode, int destNode, int parent[]);
void exportNetworkLoad(Connection conns[NUM_CONNECTIONS],Edge edgeList[2*N_EDGES],int sampleNum);

bool comparePath(const Path& p1, const Path& p2);

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

/*
 *
 */
int main(int argc, char** argv) {
    //Stores the starting index of each vertex into the edgeList/edgeWeights
    int vertexList[N_NODES+1];
    Edge edgeList[2*N_EDGES];

    //init random number generator
    srand(time(NULL));
    for(int i = 0; i < 40; ++i) {
        readGraph(vertexList,edgeList);
        randomConnections(vertexList,edgeList,i);
    }

    return 0;
}

void randomConnections(int vertexList[],Edge edgeList[2*N_EDGES],int sampleNum) {
    int paths[NUM_CONNECTIONS][N_NODES];
    Path structPaths[NUM_CONNECTIONS];
    Connection conns[NUM_CONNECTIONS];

    //Init our paths storage
    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
        structPaths[i].index = 0; //TODO: Auto-init to 0???
        for(int j = 0; j < N_NODES; ++j) {
            paths[i][j] = -1;
        }
    }



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

        single_connection_N_hops(vertexList,edgeList,10,i,conns);

        /*
        //First, compute the shortest path from v1 to v2.
        for(int hop = 3; hop < N_NODES; ++hop) {
            bool pathFound = single_path_N_hops(vertexList,edgeList,v1,v2,hop,paths,i,structPaths);
            if(pathFound == true) {
                cout << "YAY WE FOUND ONE\n";
                break;
            }else {
                cout << "WE DIDN't FIND A PATH\n";
            }
        }
        */
    }



    //Prints the current load for every edge of the graph.
    for(int i = 0; i < 2*N_EDGES; ++i) {
        cout << edgeList[i].v1 << " -> " << edgeList[i].v2 << " | LOAD: " << edgeList[i].load << "\n";
    }


    /*
    //Prints the paths/
    //TODO: probably going to get rid of the paths array pretty soon and rely solely on the structPaths.
    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
        for(int j = 0; j < N_NODES; ++j) {
            cout << paths[i][j] << " ";
        }
        cout << "\n";
    }

    //Print out the paths using the Path struct.
    for(int i = 0; i < NUM_CONNECTIONS; ++i) {
        cout << "PATH: " << structPaths[i].sourceNode << " -> " << structPaths[i].destNode << ". HOPS: " << structPaths[i].hops << "\n";
        for(int j = 0; j < structPaths[i].index; ++j) {
            cout << (*structPaths[i].edges[j]).v1 << " -> ";
        }
        cout << (*structPaths[i].edges[structPaths[i].index]).v1 << " -> " << (*structPaths[i].edges[structPaths[i].index]).v2 << "\n\n";
    }

    for(int i = 2; i < 10; ++i) {
        computeBackupPath(vertexList,edgeList,structPaths[1],structPaths[3],i);
    }
    */
}

void exportNetworkLoad(Connection conns[NUM_CONNECTIONS],Edge edgeList[2*N_EDGES],int sampleNum) {
    ofstream outputFile;
    string filename ("./load_per_edge/load_per_edge_");
    filename += to_string(sampleNum);
    filename += ".csv";

    cout << filename << "\n";
    outputFile.open(filename,ios_base::app);
    for(int i = 0; i < 2*N_EDGES; ++i) {
        outputFile << i << "," << edgeList[i].load << "\n";

    }
    outputFile.close();
    for(int j = 0; j < NUM_CONNECTIONS; ++j) {
        conns[j] = Connection();
    }
}

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
                counter++;
            }
        }
    }
    vertexList[N_NODES] = 2*N_EDGES;
}


void single_connection_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int hops, int connectionNum, Connection conns[NUM_CONNECTIONS]) {
    printf("-------DEFINING_CONNECTION_%d-------\n",connectionNum);

    //Allocate the primary path (Make it the shortest path for now).
    //We start with paths that are 1 hop away, continuing with longer paths until we successfully find a valid one.
    for(int hop = 1; hop < N_NODES; ++hop) {
        bool pathFound = single_path_N_hops(vertexList,edgeList,conns[connectionNum].sourceNode,conns[connectionNum].destNode,hop,&conns[connectionNum].primaryPath);
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
    for(int hop = conns[connectionNum].primaryPath.hops; hop < N_NODES; ++hop) {
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

    /*
    if(conns[connectionNum].validPrimary == true) {
        //Debug output for development.
        cout << "PRINTING CONNECTION " << connectionNum << "\n";
        printf("PRIMARY - Hops: %d, Index: %d\n",conns[connectionNum].primaryPath.hops,conns[connectionNum].primaryPath.index);
        printf("Source: %d, Dest: %d\n",conns[connectionNum].sourceNode, conns[connectionNum].destNode);
        for(int j = 0; j < conns[connectionNum].primaryPath.index; ++j) {
            cout << (*conns[connectionNum].primaryPath.edges[j]).v1 << " -> ";
        }
        cout << (*conns[connectionNum].primaryPath.edges[conns[connectionNum].primaryPath.index]).v1 << " -> " << (*conns[connectionNum].primaryPath.edges[conns[connectionNum].primaryPath.index]).v2 << "\n\n";
    }else {
        cout << "COULDN'T ALLOCATE A VALID PRIMARY PATH\n";
    }


    if(conns[connectionNum].validBackup == true) {
        printf("BACKUP - Hops: %d, Index: %d\n",conns[connectionNum].backupPath.hops,conns[connectionNum].backupPath.index);
        printf("Source: %d, Dest: %d\n",conns[connectionNum].sourceNode, conns[connectionNum].destNode);
        for(int j = 0; j < conns[connectionNum].backupPath.index; ++j) {
            cout << (*conns[connectionNum].backupPath.edges[j]).v1 << " -> ";
        }
        cout << (*conns[connectionNum].backupPath.edges[conns[connectionNum].backupPath.index]).v1 << " -> " << (*conns[connectionNum].backupPath.edges[conns[connectionNum].backupPath.index]).v2 << "\n\n";
    }else {
        cout << "UNABLE TO ALLOCATE A VALID BACKUP PATH\n";
    }*/

}


/**
    Allocate primary path for the connection. Finds a path from SN -> DN of EXACTLY "hops" hops.
**/
bool single_path_N_hops(int vertexList[], Edge edgeList[2*N_EDGES],int sourceNode, int destNode, int hops, Path *p) {
    //printf("-------SINGLE_PATH_%d_HOPS-------\n",hops);
    //printf("-------|%d -> %d|----------------\n",sourceNode,destNode);

    //cout << sourceNode << "," << destNode << "," << maxHops << "\n";
    int numPaths = 0;
    //initialize arrays
    int parent[N_NODES]; //parent[i] represents the parent of node i in this path, or -1 if no parent.
    int visited[N_NODES]; //visited[i] is 1 if node i has been visited on this path, 0 otherwise.

    //edgeListIndex[i] contains the index into edgeList[] (aka the compact adj list) for node i.
    int edgeListIndex[N_NODES];

    (*p).index = 0;

    //Initialize our search components
    for(int i = 0; i < N_NODES; ++i) {
        parent[i] = -1;
        visited[i] = 0;
        edgeListIndex[i] = vertexList[i];
    }

    stack <int> st;
    int currentNode;
    int neighbor;
    int currentHop = 0;

    st.push(sourceNode);
    visited[sourceNode] = 1;
    parent[sourceNode] = sourceNode;

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
            //TODO: Verify that, since we are allocating the primary path, if channels are maxed then there's nothing we can do.
            if(edgeList[edgeListIndex[currentNode]].load == MAX_CHANNELS) {
                //TODO: For now we just continue to the next neighbor
                cout << "CHANNELS ARE MAXED OUT\n";
                continue;
            }

            //Check to see if we're at the correct path length.
            if(neighbor == destNode && currentHop != hops-1) {
                continue;
            }

            //If our neighbor is the desired node, AND we're at the correct path length, save this path!
            if(neighbor == destNode && currentHop == hops-1) {
                parent[neighbor] = currentNode;
                visited[neighbor] = 1;
                numPaths++;

                (*p).edges[(*p).index] = &edgeList[edgeListIndex[currentNode]];

                //Now that we have the path set, increase the load on each edge.
                for(int i = 0; i <= (*p).index; ++i) {
                    (*(*p).edges[i]).load++;
                }

                (*p).sourceNode = sourceNode;
                (*p).destNode = destNode;
                (*p).hops = hops;

                return true;
            }

            if(!visited[neighbor]) {

                (*p).edges[(*p).index] = &edgeList[edgeListIndex[currentNode]];
                (*p).index += 1;

                st.push(neighbor);
                parent[neighbor] = currentNode;
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
        parent[currentNode] = -1;

        (*p).index -= 1;

        edgeListIndex[currentNode] = vertexList[currentNode];
        st.pop();
    }
    return false;
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
                //TODO: For now we just continue to the next neighbor
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
    cout << "COULD NOT COMPUTE A BACKUP PATH\n";
    //TODO: We need a graceful way to handle not being able to compute a backup path.
    return false;
}

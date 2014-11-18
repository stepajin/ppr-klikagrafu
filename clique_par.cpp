#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <stack>
#include <vector>
#include <algorithm>
#include <fstream>
#include "mpi.h"

#define DEBUG_WORK_EXCHANGE false
#define DEBUG_ADUV false
#define DEBUG_RESULTS false

#define MASTER_PROCESS 0
#define CHAR_FALSE 48 
#define CHAR_TRUE 49 
#define NONE 0
#define WHITE 1
#define BLACK 2

#define QUANTUM 100

#define MSG_TAG_INPUT_SIZE 1
#define MSG_TAG_INPUT_graph 2
#define MSG_TAG_WORK_DONATION 3
#define MSG_TAG_TOKEN 4
#define MSG_TAG_TERMINATION 5
#define MSG_TAG_RESULT 6
#define MSG_TAG_WORK_REQUEST 7

int PID, P;

using namespace std;

/*************
*
*    STATE 
*
**************/

struct State {
	vector<int> * nodes;
	State * next;
	State * prev;

	State() {
		nodes = new vector<int>();
		next = NULL;
	}

	State(vector<int> * nodes) {
		this->nodes = nodes;
	}

	~State() {
		delete nodes;
	}

	void print() {
		if (nodes->size() == 0) {
			cout << "empty" << endl;
			return;
		}

		for (unsigned int i = 0; i < nodes->size(); i++) {
			cout << nodes->at(i) << " ";
		}
		cout << endl;
	}

	bool contain(int node) {
		return find(nodes->begin(), nodes->end(), node) != nodes->end();
	}

	int * getArray(int length) {
		int * array = new int[length];
		array[0] = nodes->size();
		for (int i = 0; i < nodes->size(); i++) {
			array[i+1] = nodes->at(i);
		}

		return array;
	}
};

/***************************************
*
*    STACK EXTENDED BY POP FROM BOTTOM
*
***************************************/

class Stack {
	State * top;
	State * bottom;
	int size;

public:

	Stack() {
		top = NULL;
		bottom = NULL;
		size = 0;
	}

	State * getTop() {
		return top;
	}

	State * getBottom() {
		return bottom;
	}

	int getSize() {
		return size;
	}

	void push(State * state) {
		state->prev = top;
		if (top != NULL) {
			top->next = state;
		}
		top = state;

		if (bottom == NULL) {
			bottom = top;
		}

		size++;
	}

	State * pop() {
		State * state = top;
	
		if (top != NULL) {
			top = top->prev;
		}

		if (top != NULL) {
			top->next = NULL;
		} else {
			bottom = NULL;
		}

		size--;
		return state;
	}

	State * popBottom() {
		State * state = bottom;

		if (bottom != NULL) {
			bottom = bottom->next;
		}

		if (bottom != NULL) {
			bottom->prev = NULL;
		} else {
			top = NULL;
		}

		size--;
		return state;
	}

	bool empty() {
		return size == 0;
	}

	void print() {
		State * node = top;
		while (node != NULL) {
			node->print();
			node = node->prev;
		}

		if (PID == 0) cout << "-----" << endl;
	}	
};


/***************
*
*	FUNCTIONS
*
****************/

void graphviz(char ** graph, int size, int * result) {
	cout << "graph G {" << endl;
	cout << "rankdir=LR;"<< endl;

	vector<int> * v = new vector<int>();
	for (int i = 1; i <= result[0]; i++) {
		v->push_back(result[i]);
	}
	State maxClique(v);

	for (int i = 0; i < size; i++) {
		bool isInClique = false;

		if (maxClique.contain(i)) {
			cout << i << " [color=\"red\"]" << endl;
			isInClique = true;
		}

		for (int k = 0; k < size; k++) {
			if (graph[i][k] == 49 && k >= i) {
				cout << i << " -- " << k ;
	
				if (isInClique && maxClique.contain(k)) {
					cout << " [color=\"red\"];";
				}

				cout << endl;
			}
		}
	}

	cout << "}" << endl;
}

bool isClique(State * state, int addedNode, char ** graph) {
	for (unsigned int i = 0; i < state->nodes->size(); i++) {
		if (graph[state->nodes->at(i)][addedNode] == CHAR_FALSE) {
			return false;
		}
	}

	return true;
}

/******************
*
*	Handle input
*
******************/

int getInputSize(ifstream & infile) {
	int size;
	MPI_Status status;

	if (PID == MASTER_PROCESS) {
		if (infile.is_open()) {
			/* Read input size */
			infile >> size;
			infile.get(); /* newline */
		} else {
			size = 0;
		}
	
		/* Broadcast input size */
		for (int dest = 1; dest < P; dest++) {
		    MPI_Send(&size, 1, MPI_INT, dest, MSG_TAG_INPUT_SIZE, MPI_COMM_WORLD);
		}

	} else {
		/* Receive input size from master process */
	 	MPI_Recv (&size, 1, MPI_INT, MASTER_PROCESS, MSG_TAG_INPUT_SIZE, MPI_COMM_WORLD, &status);
	}

	return size;
}

char ** getInputgraph(int size, ifstream & infile) {
	int position = 0;
	MPI_Status status;

	/* Allocate input graph */
	char ** graph = new char * [size];
	for (int i = 0; i < size; i++)
		graph[i] = new char[size+1];
	
	int bufferLength = 2*size*size;
   	char * buffer = new char[bufferLength];
   	
   	/* Get input graph */	
	if (PID == MASTER_PROCESS) {

		/* Read input graph */
		for (int i = 0; i < size; i++) {
			infile.getline(graph[i], size+1);
		}

		/* Broadcast input graph */
		for (int i = 0; i < size; i++) {
	   		MPI_Pack(graph[i], size, MPI_CHAR, buffer, bufferLength, &position, MPI_COMM_WORLD);
   		}

		for (int dest = 1; dest < P; dest++) {
			MPI_Send(buffer, position, MPI_PACKED, dest, MSG_TAG_INPUT_graph, MPI_COMM_WORLD);
		}

	} else {
   		/* Receive input graph */
	    MPI_Recv(buffer, bufferLength, MPI_PACKED, MASTER_PROCESS, MSG_TAG_INPUT_graph, MPI_COMM_WORLD, &status);
      	
      	for (int i = 0; i < size; i++) {
		   	MPI_Unpack(buffer, bufferLength, &position, graph[i], size, MPI_CHAR, MPI_COMM_WORLD);
      	}	
	}

   	delete [] buffer;

   	return graph;
}

/*******************
*
*	Work donation
*
*******************/

bool donateWork(int applicant, Stack * stack, int size) {
	int * buffer;

	if (stack->getSize() <= 1) {
		buffer = new int[size + 1];
		buffer[0] = 0;
		MPI_Send(buffer, size + 1, MPI_INT, applicant, MSG_TAG_WORK_DONATION, MPI_COMM_WORLD);

		if (DEBUG_WORK_EXCHANGE) {
			cout << PID << " refused " << applicant << endl;
		}

		delete [] buffer;
		return false;

	} else {
		State * state = stack->popBottom();
		buffer = state->getArray(size);
		MPI_Send(buffer, size + 1, MPI_INT, applicant, MSG_TAG_WORK_DONATION, MPI_COMM_WORLD);

		if (DEBUG_WORK_EXCHANGE) {
			cout << PID << " donated " << applicant << " state: ";
			state->print();
		}

		delete state;
		delete [] buffer;
		return true;
	}
}

void requestWork(int donor) {
	MPI_Request request;
	MPI_Send(&PID, 1, MPI_INT, donor, MSG_TAG_WORK_REQUEST, MPI_COMM_WORLD);

	if (DEBUG_WORK_EXCHANGE) {
		cout << PID << " requested " << donor << endl;
	}
}

void receiveWork(Stack * stack, int graphSize) {
	MPI_Status status;
	int * buffer = new int[graphSize + 1];
	MPI_Recv (buffer, graphSize + 1, MPI_INT, MPI_ANY_SOURCE, MSG_TAG_WORK_DONATION, MPI_COMM_WORLD, &status);

	int stateSize = buffer[0];

	if (stateSize > 0) {
		vector<int> * vect = new vector<int>();
	
		for (int i = 1; i < stateSize+1; i++) {
			vect->push_back(buffer[i]);
		}

		State * state = new State(vect);
		stack->push(state);

		if (DEBUG_WORK_EXCHANGE) {
			cout << PID << " accepted state: ";
			state->print();
		}
	
	} else if (DEBUG_WORK_EXCHANGE) {
		cout << PID << " accepted refusal" << endl;
	}

	delete [] buffer;
}

/****************
*
*	ADUV
*
*****************/

void broadcastTermination() {
	int something = 0;
	for (int dest = 1; dest < P; dest++) {
		MPI_Send(&something, 1, MPI_INT, dest, MSG_TAG_TERMINATION, MPI_COMM_WORLD);
	}
}

void sendToken(int token) {
	int dest = (PID + 1) % P;
   	if (PID == MASTER_PROCESS) {
    	token = WHITE;
   	}

	MPI_Send(&token, 1, MPI_INT, dest, MSG_TAG_TOKEN, MPI_COMM_WORLD);

	if (DEBUG_ADUV) {
		string color = (token == BLACK)? "black" : "white";
		cout << PID << " sent " << color << " token" << endl;
	}
}

/****************
*
*	Find result
*
****************/

void sendLocalResult(int size, State * maxClique) {
	MPI_Status status;
	int * buffer = maxClique->getArray(size + 1);
	MPI_Send(buffer, size + 1, MPI_INT, MASTER_PROCESS, MSG_TAG_RESULT, MPI_COMM_WORLD);
	delete [] buffer;
}

int * getGlobalResult(int size, State * maxClique) {
	MPI_Status status;

	int * result = maxClique->getArray(size + 1);
		
	for (int src = 1; src < P; src++) {
		int * buffer = new int[size+1];
		MPI_Recv (buffer, size + 1, MPI_INT, src, MSG_TAG_RESULT, MPI_COMM_WORLD, &status);

		if (buffer[0] > result[0]) {
			delete [] result;
			result = buffer;
		} else {
			delete [] buffer;
		}
	}

	return result;
}

/*************
*
*    MAIN 
*
**************/

int main(int argc, char *argv[]) {
	int size, position = 0, counter, flag = 0, token, donor, applicant;
	int processColor;
	MPI_Status status;
	char ** graph;
	ifstream infile;
	Stack * stack = new Stack();
	State * maxClique;
	bool termination = false, waitingForReply = false, donated, tokenReceived;
	double time1, time2;

	/* Start up MPI */
  	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &PID);
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = MPI_Wtime();

	/* Handle input */
	if (PID == MASTER_PROCESS) {
		infile.open (argv[1], ifstream::in);
		if (!infile.is_open()) {
			if (PID == 0) cout << "Wrong parameters" << endl;
		}
	}

	size = getInputSize(infile);
	if (size == 0) {
		MPI_Finalize();
		return 1;
	}

	graph = getInputgraph(size, infile);
	
	if (PID == MASTER_PROCESS) {
		infile.close();
	}

   	/* Initialize stack */
   	for (int i = 0; i <= size - 1; i++) {
		if (i % P == PID) {
			vector<int> * v = new vector<int>;
			v->push_back(i);
			State * state = new State(v);
			stack->push(state);
		}
	}

	maxClique = new State(); //empty
	processColor = WHITE;
	counter = 0;
	token = (PID == MASTER_PROCESS)? BLACK : NONE;
	donor = (PID + 1) % P;
	donated = false;

	/******************
	*
	*	The main loop
	*
	*******************/

	while (!termination) {

		if (counter++ == QUANTUM) {
			counter = 0;

			/*********************
			*
			*	Handle messages
			*
			*********************/

		    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

		    if (flag) {
		    	switch (status.MPI_TAG) {
		    		case MSG_TAG_WORK_REQUEST:
    				 	MPI_Recv (&applicant, 1, MPI_INT, MPI_ANY_SOURCE, MSG_TAG_WORK_REQUEST, MPI_COMM_WORLD, &status);
		    			donated = donateWork(applicant, stack, size);
		    			if (donated && PID > applicant) {
		    				processColor = BLACK;
		    				donated = false;
		    			}
		    			break;

		    		case MSG_TAG_WORK_DONATION:
		    			receiveWork(stack, size);
		    			waitingForReply = false;
						donor = (donor + 1) % P;
				    	if (donor == PID) {
			    			donor = (donor + 1) % P;
						}
		    			break;
		    			
		    		case MSG_TAG_TOKEN:
    				 	MPI_Recv (&token, 1, MPI_INT, MPI_ANY_SOURCE, MSG_TAG_TOKEN, MPI_COMM_WORLD, &status);
    				 	
    				 	if (DEBUG_ADUV) {
			 				string color = (token == BLACK)? "black" : "white";
    				 		cout << PID << " accepted " << color << " token" << endl;
    				 	}
		    			
		    			if (processColor == BLACK) {
		    				if (DEBUG_ADUV) {
			    				cout << PID << " coloured token to black" << endl;  
							}
		    				token = BLACK;
		    			}

		    			if (PID == MASTER_PROCESS && token == WHITE) {
		    				broadcastTermination();
				    		termination = true;
				    	}
		    			break;

		    		case MSG_TAG_TERMINATION:
		    			int foo;
    				 	MPI_Recv (&foo, 1, MPI_INT, MASTER_PROCESS, MSG_TAG_TERMINATION, MPI_COMM_WORLD, &status);
		    			termination = true;
		    			break;

		    		default: 
		    			cout << "???" << endl;
		    			break;
		    	}
		    }

		    /**********************
			*
			*	Handle empty stack
			*
		    ***********************/

		    if (stack->empty() && !termination) {
		    	if (!waitingForReply && token == NONE) {
		    		requestWork(donor);
		    		waitingForReply = true;
		    	}

		    	if (token != NONE) {
			    	sendToken(token);
			    	token = NONE;
			    	processColor = WHITE;
		    	}
		    }
		}

		/*************************
		*
		*	Expand state space
		*
		**************************/

		if (!stack->empty()) {
			State * state = stack->getTop();
			stack->pop();

			if (state->nodes->size() > maxClique->nodes->size()) {
				delete maxClique;
				maxClique = state;
			}

			for (int node = state->nodes->back() + 1; node <= size - 1; node++) {
				if (isClique(state, node, graph)) {
					vector<int> * v = new vector<int>;
					*v = *(state->nodes);
					v->push_back(node);
					stack->push(new State(v));
				}
			}

			if (state != maxClique)
				delete state;
		}
	}

	/***************************
	*
	*	Find and print result
	*
	****************************/

	if (DEBUG_RESULTS) {
		usleep(1000);
		if (PID == MASTER_PROCESS) {
			cout << "---" << endl;
		}

		usleep(1000);
		cout << PID << ": ";
		maxClique->print();
	}

	if (PID != MASTER_PROCESS) {
		sendLocalResult(size, maxClique);
	} else {
		int * result = getGlobalResult(size, maxClique);
	
		/* Print result */
		if (DEBUG_RESULTS) {
			usleep(1000);
			cout << "---" << endl;
		}

		if (argc == 3 && argv[2][0] == 49) {
			graphviz(graph, size, result);
		} else {
			cout << "Size: " << result[0] << endl;
			cout << "Nodes: ";
		
			for (int i = 1; i <= result[0]; i++) {
				cout << result[i] << " ";
			}
			cout << endl;
		}

		delete [] result;
	}

	delete maxClique;

	/* Shut down MPI */
	MPI_Barrier(MPI_COMM_WORLD);
	time2 = MPI_Wtime();
	if (PID == MASTER_PROCESS && argc < 3)
	  	cout << "Time of computation: " << time2 - time1 << "s" << endl; 
  	MPI_Finalize();

  	return 0;
}
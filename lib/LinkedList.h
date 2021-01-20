#include <stdio.h>
#include <stdlib.h>
#ifndef _LinkedList_h_gf
#define _LinkedList_h_gf

typedef struct Node Node;
typedef struct LinkedList LinkedList;
//Basic singly linked data node
struct Node
{
  Node* next;
  void* data;
};

struct LinkedList
{
  Node* head;
  Node* tail;
  Node* curr;
};

//frees the node data
extern void freeNode(Node* n);
//frees the linkedlist, if freeNodes == 1 then nodes are deleted, if freeNodes == 0 then not
extern void freeLinkedList(LinkedList* this, int freeNodes);
//returns the data stored by a node
extern void* fetchData(Node* n);
//Creates a new node
extern Node* new_Node(void* data);
extern void append(LinkedList* this, void* data);
extern void push(LinkedList* this, void* data);

extern Node* pop(LinkedList* this);
//General linkedList search and deletion. Assumes node is in the list, otherwise will go out of bounds
extern void LinkedList_delete(LinkedList* this, Node* tgt);
//Takes two nodes, links if a is a valid node (can serve as an insertion or a deletion, simply replace b with NULL to delete)
extern int linkNodes(Node* a, Node* b);
extern void writeToNode(Node* a, void* data);
extern LinkedList* new_LinkedList();
//Create a linkedList iterator to make it easier to march down the list without directly accessing pointers
typedef struct LinkedListIterator LinkedListIterator;
struct LinkedListIterator
{
  LinkedList* target;
  Node* curr;
};
//Instantiates a new linkedlist iterator
extern LinkedListIterator* new_LinkedListIterator(LinkedList* this);
//returns the current node and moves to the next node
extern Node* LinkedListIterator_Next(LinkedListIterator* this);
//frees the memory and deletes the iterator
extern void LinkedListIterator_free(LinkedListIterator* this);
#endif

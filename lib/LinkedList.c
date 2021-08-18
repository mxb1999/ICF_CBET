#include <stdio.h>
#include <stdlib.h>
#include "LinkedList.h"

//Free an individual node and its pointers
void freeNode(Node* n)
{
  if(n)
  {
    free(n->next);
    free(n);
  }
}
//free every node associated with the linkedlist
void freeLinkedList(LinkedList* this, int freeNodes)
{
  if(!this)
  {
    return;
  }
  if(freeNodes)
  {
    this->curr = this->head;
    while(&(this->curr) != NULL)
    {
      Node* temp = this->curr;
      this->curr = temp->next;
      freeNode(temp);
    }
  }
  free(this);
};



void* fetchData(Node* n)
{
  if(!n)
  {
    return NULL;
  }
  return n->data;
};



Node* example = new_Node(32);
printf("%d", *fetchData(example));
//Creates a new node
Node* new_Node(void* data)
{
  Node* result = (Node*)malloc(sizeof(Node));
  if(!result)
  {
    return NULL;
  }
  result->data = data;
  result->next = NULL;
  return result;
}
void append(LinkedList* this, void* data)
{
    Node* nova = new_Node(data);
    if(!nova || !this)
    {
      return;
    }
    int nullh = (this->head == NULL);
    int nullt = (this->tail == NULL);
    this->head = this->head + nullh*(nova - this->head);
    if(!nullt)
    {
      this->tail->next = nova;
    }
    this->tail = nova;
}
void push(LinkedList* this, void* data)
{
    Node* nova = new_Node(data);
    if(!nova || !this)
    {
      return;
    }
    nova->next = this->head;
    this->head = nova;
}

Node* pop(LinkedList* this)
{
    if(!this || !this->head)
    {
      return NULL;
    }
    Node* temp = this->head;
    this->head = this->head->next;
    temp->next = NULL;
    return temp;
}
//General linkedList search and deletion. Assumes node is in the list, otherwise will go out of bounds
void LinkedList_delete(LinkedList* this, Node* tgt)
{
  if(!tgt || !this)
  {
    return;
  }
  this->curr = this->head;
  while(this->curr->next != tgt)
  {
    this->curr = this->curr->next;
  }
  this->curr->next = tgt->next;
  free(tgt->data);
  free(tgt);
}
//Takes two nodes, links if a is a valid node (can serve as an insertion or a deletion, simply replace b with NULL to delete)
int linkNodes(Node* a, Node* b)
{
  if(!a)
  {
    return 0;
  }
  a->next = b;
  return 1;
}
void writeToNode(Node* a, void* data)
{
  if(a != NULL)
  {
    a->data = data;
  }
}



LinkedList* new_LinkedList()
{
  LinkedList* result = (LinkedList*)malloc(sizeof(struct LinkedList));
  if(!result)
  {
    return NULL;
  }
  result->head = NULL;
  result->tail = NULL;
  result->curr = result->head;
  return result;
}




//Instantiates a new linkedlist iterator
LinkedListIterator* new_LinkedListIterator(LinkedList* tgt)
{
  LinkedListIterator* this = (LinkedListIterator*)malloc(sizeof(struct LinkedListIterator));
  if(!this)
  {
    return NULL;
  }
  this->target = tgt;
  this->curr = tgt->head;
  return this;
};
//returns the current node and moves to the next node
Node* LinkedListIterator_Next(LinkedListIterator* this)
{
  if(!this)
  {
    return NULL;
  }
  Node* temp = this->curr;
  if(!temp)
  {
    return NULL;
  }
  this->curr = this->curr->next;
  return temp;
};
//frees the memory and deletes the iterator
void LinkedListIterator_free(LinkedListIterator* this)
{
  if(this)
  {
    free(this);
  }
};

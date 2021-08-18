#include <iostream>

#ifndef TSQUEUE
#define TSQUEUE
  //Thread safe queue: Prevent multiple threads from pushing simultaneously
  //How to increase throughput?
  //Keep an array of insertions, start at location 0 and write until full

  template <typename T>


  struct Node
  {
  private:
    T val;
    Node* next;
  public:
    Node();
    Node(T t);
    T getVal();
    Node* next();
    ~Node();
  };
  Class tsqueue
  {
  private:
    int threads;//
    int size;
    T* priorArr;
    Node* first;
    Node* last;
  public:
    tsqueue(int t);
    void push(T val, int thread);
    T pop();
    ~tsqueue();

  };

#include <iostream>

#include "graph-utilities/semaphore.hpp"

using namespace std;


// Semaphore class (passive wait variant)

void
semaphore::post() {
  unique_lock<mutex> lck(mtx);
  atomic_fetch_add_explicit(&count, 1, memory_order_acq_rel);
  cv.notify_one();
}


void
semaphore::wait() {
  unique_lock<mutex> lck(mtx);
  while (count.load(memory_order_acquire) == 0)
    cv.wait(lck);

  atomic_fetch_add_explicit(&count, -1, memory_order_acq_rel);
}


// Spin-semaphore class (passive wait variant)

void
spin_semaphore::initialize(int memberCount) {
  unique_lock<mutex> lck(mtx);

  if (cvs) {
    cerr << "ERROR: Trying to initialize an already initialized spin-semaphore." << endl;
    exit(EXIT_FAILURE);
  }

  cvs = new condition_variable[memberCount];
  members = memberCount;
  maxIdx = 0;
  currIdx = 0;
  overflow = 0;
}


void
spin_semaphore::reset() {
  unique_lock<mutex> lck(mtx, defer_lock);

  if (cvs) {
    while (currIdx != maxIdx) {
      lck.lock();
      cvs[currIdx].notify_all();
      lck.unlock();
      currIdx = (currIdx + 1) % members;
    }
  }

  lck.lock();
  maxIdx = 0;
  currIdx = 0;
  members = 0;
  overflow = 0;

  if (cvs) {
    delete[] cvs;
    cvs = NULL;
  }

  lck.unlock();
}


void
spin_semaphore::post() {
  unique_lock<mutex> lck(mtx);

  if (currIdx == maxIdx) {
    overflow++;
  } else {
    cvs[currIdx].notify_all();
    currIdx = (currIdx + 1) % members;
  }
}


void
spin_semaphore::wait() {
  unique_lock<mutex> lck(mtx);

  if (overflow > 0) {
    overflow--;
  } else {
    int myIdx = maxIdx;

    maxIdx = (maxIdx + 1) % members;
    cvs[myIdx].wait(lck);
  }
}



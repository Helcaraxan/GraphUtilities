#include "queries.hpp"

using namespace std;


// ReachabilityQuery declarations

Vertex *
ReachabilityQuery::getSource() {
  return source.load(memory_order_acquire);
}


Vertex *
ReachabilityQuery::getTarget() {
  return target.load(memory_order_acquire);
}


bool
ReachabilityQuery::getAnswer() {
  return answer.load(memory_order_acquire);
}


bool
ReachabilityQuery::getError() {
  return error.load(memory_order_acquire);
}


SearchMethod
ReachabilityQuery::getMethod() {
  return method.load(memory_order_acquire);
}


void
ReachabilityQuery::setAnswer(bool value) {
  answer.store(value, memory_order_release);
}


void
ReachabilityQuery::setError(bool value) {
  error.store(value, memory_order_release);
}


void
ReachabilityQuery::setCancel(bool value) {
  cancel.store(value, memory_order_release);
}


void
ReachabilityQuery::setInternal(uint64_t value) {
  internal.store(value, memory_order_release);
}


void
ReachabilityQuery::setMethod(SearchMethod value) {
  method.store(value, memory_order_release);
}


bool
ReachabilityQuery::getCancel() {
  return cancel.load(memory_order_acquire);
}


uint64_t
ReachabilityQuery::getInternal() {
  return internal.load(memory_order_acquire);
}

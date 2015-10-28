#include "graph-utilities/queries.hpp"

using namespace std;


// ReachabilityQuery declarations

Vertex *
ReachabilityQuery::getSource() const {
  return source.load(memory_order_acquire);
}


Vertex *
ReachabilityQuery::getTarget() const {
  return target.load(memory_order_acquire);
}


bool
ReachabilityQuery::getAnswer() const {
  return answer.load(memory_order_acquire);
}


bool
ReachabilityQuery::getError() const {
  return error.load(memory_order_acquire);
}


SearchMethod
ReachabilityQuery::getMethod() const {
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
ReachabilityQuery::setInternal(ReachabilityQuery::InternalQueue * value) {
  internal.store(value, memory_order_release);
}


void
ReachabilityQuery::setMethod(SearchMethod value) {
  method.store(value, memory_order_release);
}


bool
ReachabilityQuery::getCancel() const {
  return cancel.load(memory_order_acquire);
}


ReachabilityQuery::InternalQueue *
ReachabilityQuery::getInternal() const {
  return internal.load(memory_order_acquire);
}

#ifndef CLUSTRVIZ_BASE_H
#define CLUSTRVIZ_BASE_H 1

#include <RcppEigen.h>
#include <vector>
#include <set>

#define CLUSTRVIZ_STATUS_UPDATE_TIME_SECS 0.1 // Every 0.1s
#define CLUSTRVIZ_STATUS_WIDTH_CHECK 20 // Every 20 status updates * 0.1s => every 2s

// Helper to determine if STL set contains an element
//
// In general, this is not efficient because one wants to do something
// with the element and/or its location, but here we really only need containment
template <typename T>
bool contains(const std::set<T>& container, T element){
  typename std::set<T>::const_iterator it = container.find(element);
  return it != container.end();
}

// Prototypes - utils.cpp
Eigen::MatrixXd MatrixProx(const Eigen::MatrixXd&,
                           double,
                           const Eigen::VectorXd&,
                           bool);

#endif

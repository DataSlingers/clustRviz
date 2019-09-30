#include "clustRviz.h"

Rcpp::List get_cluster_assignments_impl(const Eigen::MatrixXi& E,
                                        const Eigen::VectorXi& edge_ind,
                                        int n){

  // We use a simple (depth-first?) search to determine the connected components
  // of the graphs. Since we frequently need to check if a vertex is in a component,
  // we represent each component as a std::set<int>, and we store the components
  // in a std::vector
  std::vector<std::set<int> > components;
  std::set<int> all_vertices_seen;

  // Iterate over all possible edges - this big loop is a no-op where edge_ind == 0
  for(unsigned int i = 0; i < edge_ind.size(); i++){

    if(edge_ind(i) != 0){ // If the edge is present
      int edge_begin = E(i, 0);
      int edge_end   = E(i, 1);

      // We begin by looking for the first vertex (here called "begin") in each component
      bool found_component_begin = false;

      for(unsigned int j = 0; j < components.size(); j++){
        std::set<int>& component_j = components[j];

        if(contains(component_j, edge_begin)){
          // Once we found a component containing the "begin" vertex, let's see if
          // it contains the "end" vertex.
          found_component_begin = true;
          bool found_component_end = false;

          // If begin and end are already in the same component, we don't need
          // to do anything.
          if(contains(component_j, edge_end)){
            found_component_end = true;
            break; // Continue to next edge
          }

          // Now check other components
          for(unsigned int k = 0; k < components.size(); k++){
            if(k != j){ // We handled k == j above

              std::set<int>& component_k = components[k];

              if(contains(component_k, edge_end)){
                found_component_end = true;
                // This implies components j and k are connected, but weren't already
                //
                // First we copy all the elements of component_k into component_j
                component_j.insert(component_k.begin(), component_k.end());
                // Now we drop component k
                components.erase(components.begin() + k);
                break; // No need to check other components
              }
            }
          }

          // If we never found `edge_end` in any component, we add it to component J
          // since it is connected to `edge_begin.`
          if(!found_component_end){
            component_j.insert(edge_end);
            all_vertices_seen.insert(edge_end);
          }

          break; // No need to check other components,
                 // since edge_begin can only be in one component
        }
      }

      // If we didnt' find edge_begin in any component, we first check for edge_begin
      if(!found_component_begin){

        // First check if edge_end is anywhere:
        // If it is, then we add edge_begin to the same component
        bool found_component_end_inner = false;

        for(unsigned int j = 0; j < components.size(); j++){
          std::set<int>& component_j = components[j];
          if(contains(component_j, edge_end)){
            // We didn't find edge_begin, but we do have edge_end, so let's add
            // edge_begin to the same component
            component_j.insert(edge_begin);
            all_vertices_seen.insert(edge_begin);
            found_component_end_inner = true;
            break;
          }
        }

        // If we can't find edge_begin or edge_end anywhere, they are both new
        // and get there own new component
        if(!found_component_end_inner){
          std::set<int> new_component{edge_begin, edge_end};
          all_vertices_seen.insert(edge_begin);
          all_vertices_seen.insert(edge_end);
          components.push_back(new_component);
        }
      }
    }
  }

  // Sort components in decreasing size order

  // Now build things in a way that works for R
  //
  // Add singleton components for isolated vertices
  for(int i = 1; i <= n; i++){ // We're using R's (1-based) vertex labels
    if(!contains(all_vertices_seen, i)){
      std::set<int> new_component{i};
      components.push_back(new_component);
    }
  }

  // Sort components by smallest vertex index
  // This is independent of the order of the edge set / algorithm used
  std::sort(components.begin(),
            components.end(),
            [](const std::set<int>& left, const std::set<int>& right){
              return *left.begin() < *right.begin();
            });

  int num_components = components.size();

  Rcpp::IntegerVector component_sizes(num_components, 1);
  Rcpp::IntegerVector component_indicators(n, -1);

  // Assign labels - loop over components and then elements within components
  // Requires irregular access to component_indicators, but it's O(1) (=O(n) total) instead
  // of searching through all the sets repeatedly
  for(unsigned int i = 0; i < components.size(); i++){
    std::set<int>& component_i = components[i];
    component_sizes[i] = component_i.size();

    for(int j : component_i){
      component_indicators[j - 1] = i + 1; // j - 1 because edges are 1-indexed (coming form R)
                                           // Similarly, we start counting components at 1
    }
  }

  return Rcpp::List::create(Rcpp::Named("membership") = Rcpp::wrap(component_indicators),
                            Rcpp::Named("csize") = Rcpp::wrap(component_sizes),
                            Rcpp::Named("no") = Rcpp::wrap(num_components));
}

// Get cluster assignments
//
// Given the output of CARP/CBASS (in vectorized form), perform the actual cluster
// assignments
//
// Input: E - a two column matrix with the edge set used for clustering.
//            (The i-th row is (j, k) if there is an edge between j and k)
//        edge_ind - A 0/1 matrix with NCOL(edge_ind) == NROW(E)
//                   Element (i, j) is 1 if the j-th edge (as defined by E)
//                   was "active" in the i-th stage of clustering. (I.e., were those
//                   two elements fused together). Each row is one iteration of CARP/CBASS
//        n - the number of observations
//
// Output - A list of NROW(edge_ind) `Rcpp::List`s, each of which has the
//          following elements:
//            - membership - an `Rcpp::IntegerVector` with cluster labels
//            - csize      - The number of elements in each cluster
//            - no         - the number of clutsers
// TODO: Do we need all three of these?
// [[Rcpp::export]]
Rcpp::List get_cluster_assignments(const Eigen::MatrixXi& E,
                                   const Eigen::MatrixXi& edge_ind,
                                   int n){
  Rcpp::List return_object(edge_ind.rows());

  for(Eigen::Index i = 0; i < edge_ind.rows(); i++){
    return_object[i] = get_cluster_assignments_impl(E, edge_ind.row(i), n);
  }

  return return_object;
}

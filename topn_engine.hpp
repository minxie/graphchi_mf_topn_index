#ifndef __MX_TOPN_ENGINE
#define __MX_TOPN_ENGINE

#include <map>
#include <set>
#include <algorithm>
#include <iostream>

/* Parameters which should be defined in the main algorithm */
std::vector<VertexDataType> *latent_factors;
float (*pprediction_func_test)(const vertex_data&, const vertex_data&, const float, double &, void *) = NULL;

/* TopN rec parameters */
int n_top = -1;

bool sort_items_c(std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {
  return a.second > b.second;
}

/* General purpose TopN program */
struct GeneralTopNProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {

  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
    if (vertex.num_outedges() > 0) {
      /* Exclude neighbors */
      std::map<unsigned int, bool> h_neighbor;
      for (int e = 0; e < vertex.num_edges(); e++) {
        h_neighbor[vertex.edge(e)->vertex_id()] = true;
      }

      /* Calculate TopN recommendations */
      std::vector<std::pair<unsigned int, double> > rec_vec;
      vertex_data & vdata = (*latent_factors)[vertex.id()];
      for (unsigned int i = M; i < M + N; i++) {
      //  if (h_neighbor.find(i) == h_neighbor.end()) { // Not observed
          vertex_data & nbr_latent = latent_factors_inmem[i];
          double prediction;
          (*pprediction_func_test)(vdata, nbr_latent, 0, prediction, NULL);
          rec_vec.push_back(std::make_pair(i, prediction));
       // }
      }
      std::partial_sort(rec_vec.begin(), rec_vec.begin()+n_top, rec_vec.end(), sort_items_c);
    }
  }

  void before_iteration(int iteration, graphchi_context & gcontext) {
  }

  void after_iteration(int iteration, graphchi_context &gcontext) {
   gcontext.set_last_iteration(0);  
}

};

void run_general_topn_program(graphchi_engine<VertexDataType, EdgeDataType> *engine,
                              std::vector<vertex_data> *latent_factors_inmem,
                              float (*prediction_func)(const vertex_data & user, 
                                                       const vertex_data & movie, 
                                                       float rating, double & prediction, 
                                                       void * extra)) {
  latent_factors = latent_factors_inmem;
  pprediction_func_test = prediction_func;

  GeneralTopNProgram test_prog;
  engine->run(test_prog, 1);
}




#endif //__MX_TOPN_ENGINE

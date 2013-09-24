#ifndef __MX_TOPN_ENGINE
#define __MX_TOPN_ENGINE


struct GeneralTopNProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {

  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
  }

  void before_iteration(int iteration, graphchi_context & gcontext) {
  }

  void after_iteration(int iteration, graphchi_context &gcontext) {
  }

};

void run_general_topn_program(graphchi_engine<VertexDataType, EdgeDataType> *engine,
                              std::vector<vertex_data> *latent_factors_inmem) {
  
}




#endif //__MX_TOPN_ENGINE

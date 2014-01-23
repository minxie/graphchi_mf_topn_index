/**
 * @file
 * @author  Danny Bickson
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.

 *
 * @section DESCRIPTION
 *
 * Matrix factorization with the Stochastic Gradient Descent (SGD) algorithm.
 * Algorithm is described in the papers:
 * 1) Matrix Factorization Techniques for Recommender Systems Yehuda Koren, Robert Bell, Chris Volinsky. In IEEE Computer, Vol. 42, No. 8. (07 August 2009), pp. 30-37. 
 * 2) Takács, G, Pilászy, I., Németh, B. and Tikk, D. (2009). Scalable Collaborative Filtering Approaches for Large Recommender Systems. Journal of Machine Learning Research, 10, 623-656.
 *
 * 
 */


#include "eigen_wrapper.hpp"
#include "common.hpp"

double sgd_lambda = 1e-3; //sgd regularization
double sgd_gamma = 1e-3;  //sgd step size
double sgd_step_dec = 0.9; //sgd step decrement

struct vertex_data {
  vec pvec; //storing the feature vector

  vertex_data() {
    pvec = zeros(D);
  }
  void set_val(int index, float val){
    pvec[index] = val;
  }
  float get_val(int index){
    return pvec[index];
  }

};

#include "util.hpp"

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
typedef vertex_data VertexDataType;
typedef float EdgeDataType;  // Edges store the "rating" of user->movie pair

graphchi_engine<VertexDataType, EdgeDataType> * pengine = NULL; 
graphchi_engine<VertexDataType, EdgeDataType> * pvalidation_engine = NULL; 
std::vector<vertex_data> latent_factors_inmem;

#include "rmse.hpp"
#include "rmse_engine.hpp"
#include "topn_engine.hpp"
#include "topn_engine_kd2.hpp"
#include "io.hpp"
#include "kdtree.hpp"
#include "rtree.hpp"
/** compute a missing value based on SGD algorithm */
float sgd_predict(const vertex_data& user, 
    const vertex_data& movie, 
    const float rating, 
    double & prediction, 
    void * extra = NULL){


  prediction = dot_prod(user.pvec,movie.pvec);

  //truncate prediction to allowed values
  prediction = std::min((double)prediction, maxval);
  prediction = std::max((double)prediction, minval);
  //return the squared error
  float err = rating - prediction;
  assert(!std::isnan(err));
  return err*err; 

}


/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct SGDVerticesInMemProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
  mutex lock;

  /**
   * Called before an iteration is started.
   */
  void before_iteration(int iteration, graphchi_context &gcontext) {
    reset_rmse(gcontext.execthreads);
  }



  /**
   * Called after an iteration has finished.
   */
  void after_iteration(int iteration, graphchi_context &gcontext) {
    sgd_gamma *= sgd_step_dec;
    training_rmse(iteration, gcontext);
    run_validation(pvalidation_engine, gcontext);
   /* std::vector<double> lbo (D,0);
    std::vector<double> rbo (D,1);
    double copy[N];   
    for(int j = 0; j < D; j++)
    {
        for(unsigned int i = M; i < M + N; i++)
            copy[i - M] = latent_factors_inmem.at(i).pvec[j];
        std::sort(copy, copy+N);
        lbo.at(j) = copy[0];
        rbo.at(j) = copy[N-1];
    }
    std::cout << "left bound: ("; 
    for(int i = 0; i < D; i++)
	std::cout << lbo.at(i) << ", ";
    std::cout << ").  right bound: (";
    for(int i = 0; i < D; i++)
	std::cout << rbo.at(i) << ", ";
    std::cout << ")." << std::endl;*/
  }

  /**
   *  Vertex update function.
   */
  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
    //go over all user nodes
    if ( vertex.num_outedges() > 0){
      vertex_data & user = latent_factors_inmem[vertex.id()]; 
      //go over all ratings
      for(int e=0; e < vertex.num_edges(); e++) {
        float observation = vertex.edge(e)->get_data();                
        vertex_data & movie = latent_factors_inmem[vertex.edge(e)->vertex_id()];
        double estScore;
        rmse_vec[omp_get_thread_num()] += sgd_predict(user, movie, observation, estScore);
        double err = observation - estScore;
        if (std::isnan(err) || std::isinf(err))
          logstream(LOG_FATAL)<<"SGD got into numerical error. Please tune step size using --sgd_gamma and sgd_lambda" << std::endl;

        // lock.lock();
        // std::cout << err << std::endl;
        // lock.unlock();

        //NOTE: the following code is not thread safe, since potentially several
        //user nodes may updates this item gradient vector concurrently. However in practice it
        //did not matter in terms of accuracy on a multicore machine.
        //if you like to defend the code, you can define a global variable
        //mutex mymutex;
        //
        //and then do: mymutex.lock()
        movie.pvec += sgd_gamma*(err*user.pvec - sgd_lambda*movie.pvec);
        //and here add: mymutex.unlock();
        user.pvec += sgd_gamma*(err*movie.pvec - sgd_lambda*user.pvec);
      }
    }

  }



};



//dump output to file
void output_sgd_result(std::string filename) {
  MMOutputter_mat<vertex_data> user_mat(filename + "_U.mm", 0, M, "This file contains SGD output matrix U. In each row D factors of a single user node.", latent_factors_inmem);
  MMOutputter_mat<vertex_data> item_mat(filename + "_V.mm", M ,M+N,  "This file contains SGD  output matrix V. In each row D factors of a single item node.", latent_factors_inmem);

  logstream(LOG_INFO) << "SGD output files (in matrix market format): " << filename << "_U.mm" <<
                                                                           ", " << filename + "_V.mm " << std::endl;
}


int main(int argc, const char ** argv) {

  D = 2;
  M = 0;
  N = 100;

  latent_factors_inmem.resize(100);

  for (int i=0; i < (int)100; i++){
    for (int j=0; j<D; j++)
      latent_factors_inmem[i].pvec[j] = 1.0 * drand48();
  }
  
  RTree T;
  T.build_rtree(latent_factors_inmem);
  T.print_rtree();
  
  return 0;
}

#ifndef __TOPN_ENGINE_KD2_H_
#define __TOPN_ENGINE_KD2_H_

#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <queue>

#include "kdtree.hpp"
/* Parameters which should be defined in the main algorithm */
//std::vector<VertexDataType> *latent_factors;

/* TopN rec parameters */
//int n_top = -1;
kd_Node* root;

struct kdnodepq {
	kd_Node* rectangle;
	double dist;
	kdnodepq(kd_Node* node, double distance)
	{
	    rectangle = node; 
	    dist = distance;
	}
};

struct kdnodepqCmp {
	bool operator() (const kdnodepq& p1, const kdnodepq& p2)
	{
	    	return p1.dist <= p2.dist;
	}
};

struct vertexpq {
	unsigned int index;
	double dist;
	vertexpq(int item_index, double distance)
	{
	    index = item_index;
	    dist = distance;
	}
};

struct vertexpqCmp {
	bool operator()	(const vertexpq& p1, const vertexpq& p2)
	{
	    return p1.dist >= p2.dist;
	}
};


double distance(kd_Node* rectangle, const vertex_data& user)
{   
    //calculate the distance between a rectangle and a query user
    double dist=0;
    for(int i = 0; i < D; i++)
    {//calculate the norm of user vector and its distance from the rectangle 
        if(user.pvec[i] >= 0)
            dist += rectangle->rbound.at(i) * user.pvec[i];
        else
            dist += rectangle->lbound.at(i) * user.pvec[i];
    }
    return dist;
} 

void searchkdtree( const vertex_data& user, const std::map<unsigned int, bool>& h_neighbor, unsigned int uid, kd_Node* mnode,
						std::priority_queue<vertexpq, std::vector<vertexpq>, vertexpqCmp> &candi_pq,
						double upbound)
{
  if(upbound < distance(mnode, user))
  {
     if(mnode->lchild == NULL && mnode->rchild == NULL)
     {
        for(int i = 0; i < mnode->item.size(); i++)
        { //  exclude neighbors
            unsigned int item_index = mnode->item.at(i).first;
            if(h_neighbor.find(item_index) == h_neighbor.end())
            {
                double prediction = dot_prod(user.pvec, mnode->item.at(i).second.pvec);
                if(prediction > upbound)
                {
                    candi_pq.push(vertexpq(mnode->item.at(i).first, prediction));
                    if(candi_pq.size() > n_top)
                    {
                        candi_pq.pop();
                        upbound = candi_pq.top().dist;    //update upper bound
		    }
                    //upbound = candi_pq.top().dist;    //update upper bound
               }
            }
        }
    }
    else
    {
        double ldist = distance(mnode->lchild, user);
        double rdist = distance(mnode->rchild, user);   
        if(ldist <= rdist)
        {
            searchkdtree(user, h_neighbor, uid, mnode->lchild, candi_pq, upbound);
            searchkdtree(user, h_neighbor, uid, mnode->rchild, candi_pq, upbound);

	}
        else
        {
            searchkdtree(user, h_neighbor, uid, mnode->rchild, candi_pq, upbound);
            searchkdtree(user, h_neighbor, uid, mnode->lchild, candi_pq, upbound);

	}
    }

  }
	
    return ;
}

/* General purpose TopN program */
struct KdTopNProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    

  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
    if (vertex.num_outedges() > 0) {
     
      std::map<unsigned int, bool> h_neighbor;
      for (int e = 0; e < vertex.num_edges(); e++) {
        h_neighbor[vertex.edge(e)->vertex_id()] = true;
      }
      // Calculate TopN recommendations 
      std::vector<std::pair<unsigned int, double> > rec_vec;
      std::vector<std::pair<unsigned int, double> > reverse_rec;
      vertex_data & user = (*latent_factors)[vertex.id()];  
      std::priority_queue<vertexpq, std::vector<vertexpq>, vertexpqCmp> candi_pq;
      double upbound = 0;      
      searchkdtree(user, h_neighbor, vertex.id(), root, candi_pq, upbound);
	
      while (candi_pq.size() > 0)
      {
        reverse_rec.push_back(std::make_pair(candi_pq.top().index, candi_pq.top().dist));
        candi_pq.pop();
      }
      for(int i = n_top-1; i >= 0; i--)
      {
        rec_vec.push_back(reverse_rec.at(i));
      }
	/*std::ofstream ofs("../../result/rec_result_kd", std::ofstream::out | std::ofstream::app);
	for (int i = 0; i < n_top; i++)
	{
	  ofs << "<" << vertex.id() << ", " << rec_vec.at(i).first << ", " << rec_vec.at(i).second << ">";
	}	
	ofs << std::endl;
        ofs.close();*/
	// std::cout << "found top n items for user: "<< vertex.id() << std::endl;
    }  
  } 
  
  void before_iteration(int iteration, graphchi_context & gcontext) {

  }

  void after_iteration(int iteration, graphchi_context &gcontext) {
   gcontext.set_last_iteration(0);  
}


};

void run_kd_topn_program(graphchi_engine<VertexDataType, EdgeDataType> *engine,
                              std::vector<vertex_data> *latent_factors_inmem, 
                                                kd_Node* index_root) {
  latent_factors = latent_factors_inmem;
  root = index_root;

  KdTopNProgram test_prog;
  engine->run(test_prog, 1);
}

#endif //__TOPN_ENGINE_M_H_

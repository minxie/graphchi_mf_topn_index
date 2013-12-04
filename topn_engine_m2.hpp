#ifndef __TOPN_ENGINE_M_H_
#define __TOPN_ENGINE_M_H_

#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <queue>

#include "mtree.hpp"
/* Parameters which should be defined in the main algorithm */
//std::vector<VertexDataType> *latent_factors;

/* TopN rec parameters */
//int n_top = -1;
m_Tree* root;

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


double bound(m_Tree* circle, const vertex_data& user)
{   
    //calculate upper bound between a query user and a M-tree node
    double dist=0;
    for(int i = 0; i < D; i++)
      dist += user.pvec[i] * circle->center.at(i);
    dist += circle->radius * (sqrt(dot_prod(user.pvec, user.pvec)));
    return dist;
}

void searchmtree( const vertex_data& user, const std::map<unsigned int, bool>& h_neighbor, unsigned int uid, m_Tree* mnode,
						std::priority_queue<vertexpq, std::vector<vertexpq>, vertexpqCmp> &candi_pq,
						double upbound)
{
  if(upbound < bound(mnode, user))
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
			upbound = candi_pq.top().dist;
                        candi_pq.pop();
                    }
                    //upbound = candi_pq.top().dist;    //update upper bound
               }
            }
        }
    }
    else
    {
        double ldist=0, rdist=0;
        for(int i = 0; i < D; i++)
        {
            ldist += user.pvec[i] * mnode->lchild->center.at(i);
            rdist += user.pvec[i] * mnode->rchild->center.at(i);
        }
        if(ldist <= rdist)
        {
            searchmtree(user, h_neighbor, uid, mnode->lchild, candi_pq, upbound);
            searchmtree(user, h_neighbor, uid, mnode->rchild, candi_pq, upbound);

	}
        else
        {
            searchmtree(user, h_neighbor, uid, mnode->rchild, candi_pq, upbound);
            searchmtree(user, h_neighbor, uid, mnode->lchild, candi_pq, upbound);

	}
    }

  }
	
    return ;
}

/* General purpose TopN program */
struct MTopNProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    

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
      searchmtree(user, h_neighbor, vertex.id(), root, candi_pq, upbound);
      
      int size = candi_pq.size();
      while (candi_pq.size() > 0)
      {
        reverse_rec.push_back(std::make_pair(candi_pq.top().index, candi_pq.top().dist));
        candi_pq.pop();
      }
      for(int i = n_top-1; i >= 0; i--)
      {
        rec_vec.push_back(reverse_rec.at(i));
      }
	std::ofstream ofs("../../result/rec_result_m", std::ofstream::out | std::ofstream::app);
	for (int i = 0; i < n_top; i++)
	{
	  ofs << "<" << vertex.id() << ", " << rec_vec.at(i).first << ", " << rec_vec.at(i).second << ">";
	}	
	ofs << std::endl;
        ofs.close();
    }  
  } 
  
  void before_iteration(int iteration, graphchi_context & gcontext) {

  }

  void after_iteration(int iteration, graphchi_context &gcontext) {
   gcontext.set_last_iteration(0);  
}


};

void run_m_topn_program(graphchi_engine<VertexDataType, EdgeDataType> *engine,
                              std::vector<vertex_data> *latent_factors_inmem, 
                                                m_Tree* index_root) {
  latent_factors = latent_factors_inmem;
  root = index_root;
  //std::cout << "test program begin" << std::endl;
  MTopNProgram test_prog;
  engine->run(test_prog, 1);
}

#endif //__TOPN_ENGINE_M_H_

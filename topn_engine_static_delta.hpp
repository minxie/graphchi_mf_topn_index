#ifndef __MX_TOPN_ENGINE
#define __MX_TOPN_ENGINE

#include <map>
#include <set>
#include <algorithm>
#include <iostream>

/* Parameters which should be defined in the main algorithm */
std::vector<VertexDataType> *latent_factors;
float (*pprediction_func_test)(const vertex_data&, const vertex_data&, const float, double &, void *) = NULL;

std::vector<double> lbound (D, 0);
std::vector<double> rbound (D, 1);


/* TopN rec parameters */
int n_top = -1;
/* Top maintain items */
int num = 100;
int delta_max = -1;
int delta_min = 100;

bool sort_items_c(std::pair<unsigned int, double> a, std::pair<unsigned int, double> b) {
  return a.second > b.second;
}

struct util
{
	/* data */
	double utility;
	int index;
	util(double utility1, int index1)
	{
		utility = utility1;
		index = index1;
	}
};


bool utilcmp(util u1, util u2)
{
	return u1.utility > u2.utility;
}

double shortest_dis(const vertex_data& item, const vertex_data& user, double theta)
{	//	calculate the shortest distance of the item within the range of a user
	double distance = 0;
	for (int i = 0; i < D; i++)
	{
		if (item.pvec[i] > 0)
			distance += (user.pvec[i] - theta) * item.pvec[i];
		else
			distance += (user.pvec[i] + theta) * item.pvec[i];
	}
	return distance;
}

void find_top(const vertex_data& user, unsigned int u_index, std::vector<std::pair<unsigned int, double> > rec_vec)
{

	/* the adaptive range of the current user */
	double theta = 0.005;
	/* the number of items to be maintained besides top n */
	int delta = 0;
	/* indicator of the sign of each user vector dimension */
	int flag[D];
	
	std::vector<std::pair<unsigned int, double> > shortest;
        std::vector<util> Util;
        double uppervalue = 0;   
	vertex_data cur_user = user;   

	for (int i = 0; i < n_top; i++)
	{
		double dist = shortest_dis( (*latent_factors)[rec_vec.at(i).first], cur_user, theta);
		shortest.push_back(std::make_pair(rec_vec.at(i).first, dist));
	}
	double utility_value = 0;
	for (int i = 0; i < D; i++)
	{
		if(cur_user.pvec[i] > 1e-7)
		{
			flag[i] = 0;
			utility_value = rbound[i] * (cur_user.pvec[i] + theta) / cur_user.pvec[i];
			Util.push_back(util(utility_value, i));
		}
			
		else if (cur_user.pvec[i] < -1e-7)
		{
			flag[i] = 1;
			cur_user.pvec[i] = -cur_user.pvec[i];
			utility_value = (-lbound[i]) * (cur_user.pvec[i] - theta) / cur_user.pvec[i];
			Util.push_back(util(utility_value, i));
		}
		else
		{
			std::cout << "User item has zero dimension value: " << u_index << ", " << i << std::endl;
			flag[i] = 2;
			utility_value = 0;
			Util.push_back(util(utility_value, i));
		}
	}
	/* sort the utility in descending order */
	std::sort(Util.begin(), Util.begin()+D, utilcmp);
	int tmp = 1;

	while(tmp)
	{
		int selected = n_top + delta;
		unsigned int index = rec_vec.at(selected).first;
                vertex_data & cur_item = (*latent_factors)[index];	
		//	budget is the production of the last item selected and the current user
		double budget = dot_prod( (*latent_factors)[rec_vec.at(selected-1).first].pvec, cur_user.pvec);

		double value = 0;
		for(int i = 0; i < D; i++)
		{
			int dim = Util.at(i).index;						//	current dimension based on the order of utility

			if(flag[dim] == 0)									//	u[i] is initially positive
			{
				budget -= cur_user.pvec[dim] * lbound[dim];
				if(cur_user.pvec[dim] * (rbound[dim] - lbound[dim]) < budget)		// current dimension can be completely selected
				{
					value += (cur_user.pvec[dim] + theta) * rbound[dim];
					budget -= cur_user.pvec[dim] * (rbound[dim] - lbound[dim]);
				}

				else
				{
					value += (cur_user.pvec[dim] + theta) * budget / cur_user.pvec[dim];
					budget = 0;
					break;
				}
			}
			else												//	u[i] is initially negative
			{
				budget -= cur_user.pvec[dim] * (-rbound[dim]);	//	budget will increase here
				if(cur_user.pvec[dim] * (rbound[dim] - lbound[dim]) < budget)
				{
					value += (cur_user.pvec[dim] + theta) * rbound[dim];
					budget -= cur_user.pvec[dim] * (rbound[dim] - lbound[dim]);
				}
				else
				{
					budget = 0;
					break;
				}
			}
		}
		int cur_num = 0;
		for(int i = 0; i < shortest.size(); i++)
		{
			if(shortest.at(i).second >= value)
				cur_num++;
		}
		if(cur_num >= n_top)		//	found enough items to maintain
		{
			if(delta > delta_max)
			    delta_max = delta;
			if(delta < delta_min)
			    delta_min = delta;
			uppervalue = value;
			tmp = 0;
		}
		else
		{
			double dist = shortest_dis(cur_item, cur_user, theta);
			shortest.push_back(std::make_pair(index, dist));
			delta++;	
		}
	}

	std::ofstream ofs("../../result/rec_result_inc", std::ofstream::out | std::ofstream::app);
	ofs << "User " << u_index << ", Value " << uppervalue << ", Maintain items: " << n_top+delta << std::endl;
	
	for(int i = 0; i < n_top+delta; i++)
	{
	    ofs << "<" << u_index << ", " << rec_vec.at(i).first << ", " << rec_vec.at(i).second << ", " << shortest.at(i).second << ">";
	}
	ofs << std::endl;
	ofs.close();

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

      vertex_data & cur_user = (*latent_factors)[vertex.id()];
      for (unsigned int i = M; i < M + N; i++) {
        if (h_neighbor.find(i) == h_neighbor.end()) { // Not observed
          vertex_data & nbr_latent = latent_factors_inmem[i];
          double prediction;
          (*pprediction_func_test)(cur_user, nbr_latent, 0, prediction, NULL);
          rec_vec.push_back(std::make_pair(i, prediction));
        }
      }
      //std::partial_sort(rec_vec.begin(), rec_vec.begin()+n_top, rec_vec.end(), sort_items_c);
	std::partial_sort(rec_vec.begin(), rec_vec.begin()+num+n_top, rec_vec.end(), sort_items_c);
	find_top(cur_user, vertex.id(), rec_vec);

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
  double copy[N];
  for(int j = 0; j < D; j++)
  {
      for(unsigned int i = M; i < M + N; i++)
          copy[i - M] = latent_factors_inmem->at(i).pvec[j];
       std::sort(copy, copy+N);
       lbound.at(j) = copy[0];
       rbound.at(j) = copy[N-1];
  }


  GeneralTopNProgram test_prog;
  engine->run(test_prog, 1);

  std::cout << "Max delta: " << delta_max << " Min delta: " << delta_min << std::endl;
}



#endif //__MX_TOPN_ENGINE

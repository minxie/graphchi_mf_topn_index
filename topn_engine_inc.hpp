#ifndef __MX_TOPN_ENGINE
#define __MX_TOPN_ENGINE

#include <map>
#include <set>
#include <algorithm>
#include <iostream>

/* Parameters which should be defined in the main algorithm */
std::vector<VertexDataType> *latent_factors;
std::vector<VertexDataType> *lf_previous;
float (*pprediction_func_test)(const vertex_data&, const vertex_data&, const float, double &, void *) = NULL;
std::vector<double> *lbound ;
std::vector<double> *rbound ;
std::vector<double> *ldelta ;
std::vector<double> *rdelta ;
/* Parameters  */
std::multimap<unsigned int, unsigned int> item_top;
std::vector<unsigned int> top_one;

int cur_iter;
unsigned int num_buffer;
/* TopN rec parameters */
int n_top = -1;
/* Top maintain items */
int num = 100;
int delta_max = -1;
int delta_min = 100;
int delta_sum;
double theta = 0.05;

struct rec_buf
{
    bool valid;
    unsigned int index;         //      first item points to it
    int num;
    std::vector<unsigned int> item;
    double uppervalue;

    rec_buf() {}
};

/*struct user_node {
//  unsigned int index;
  rec_buf * pbuffer;
  user_node(rec_buf * y)
  {
        pbuffer = y;
  }
};*/


std::vector<rec_buf *> user_buffer;

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

double update_uppervalue(double theta, const vertex_data& user)
{
	double update = 0;
	for (int i = 0; i < D; i++)
	{
	    if ( (*ldelta)[i] * (user.pvec[i]-theta) > update)
		update = (*ldelta)[i] * (user.pvec[i]-theta);
	    if ( (*ldelta)[i] * (user.pvec[i]+theta) > update)
		update = (*ldelta)[i] * (user.pvec[i]+theta);
	    if ( (*rdelta)[i] * (user.pvec[i]-theta) > update)
		update = (*rdelta)[i] * (user.pvec[i]-theta);
	    if ( (*rdelta)[i] * (user.pvec[i]+theta) > update)
		update = (*rdelta)[i] * (user.pvec[i]+theta);  
	}
	return update;
}

unsigned int update_topone(const vertex_data& user, std::vector<unsigned int> itemset)
{
	unsigned int max_index;
	double max_dist=0;
	for(int i = 0; i < itemset.size(); i++)
	{
	    if(dot_prod(user.pvec, (*latent_factors)[itemset.at(i)].pvec) > max_dist)
	    {
		max_dist = dot_prod(user.pvec, (*latent_factors)[itemset.at(i)].pvec);
		max_index = itemset.at(i);
	    }
	}
	return max_index;
}
	
void find_top(const vertex_data& user, unsigned int u_index, std::vector<std::pair<unsigned int, double> > rec_vec)
{
	/* the adaptive range of the current user */
	theta = 0.005;
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
		/* Insert the top1 item-user pair into map */
		//item_top.insert(std::pair<unsigned int, unsigned int> (rec_vec.at(i).first, u_index));	
	}

	double utility_value = 0;
	for (int i = 0; i < D; i++)
	{
		if(cur_user.pvec[i] > 1e-7)
		{
			flag[i] = 0;
			utility_value = (*rbound)[i] * (cur_user.pvec[i] + theta) / cur_user.pvec[i];
			Util.push_back(util(utility_value, i));
		}
			
		else if (cur_user.pvec[i] < -1e-7)
		{
			flag[i] = 1;
			cur_user.pvec[i] = -cur_user.pvec[i];
			utility_value = (-(*lbound)[i]) * (cur_user.pvec[i] - theta) / cur_user.pvec[i];
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
	{	//std::cout << "rec_Vec size: " << rec_vec.size() << std::endl;
		int selected = n_top + delta;
		if(selected >= rec_vec.size())
			std::cout << "*************Overflow at line 176*************" << std::endl;
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
				budget -= cur_user.pvec[dim] * (*lbound)[dim];
				if(cur_user.pvec[dim] * ( (*rbound)[dim] - (*lbound)[dim]) < budget)		// current dimension can be completely selected
				{
					value += (cur_user.pvec[dim] + theta) * (*rbound)[dim];
					budget -= cur_user.pvec[dim] * ((*rbound)[dim] - (*lbound)[dim]);
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
				budget -= cur_user.pvec[dim] * (-( (*rbound)[dim]) );	//	budget will increase here
				if(cur_user.pvec[dim] * ( (*rbound)[dim] - (*lbound)[dim] ) < budget)
				{
					value += (cur_user.pvec[dim] + theta) * (*rbound)[dim];
					budget -= cur_user.pvec[dim] * ( (*rbound)[dim] - (*lbound)[dim]);
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
	/* Create new buffer */
	rec_buf * tmp_buffer = new rec_buf();
	tmp_buffer->index = u_index;
	for(int i = 0; i < n_top+delta; i++)
	{
	    	if(i >= rec_vec.size())
			std::cout << "*************Overflow at line 243*************" << std::endl;
	    	tmp_buffer->item.push_back(rec_vec.at(i).first);
	}
	tmp_buffer->uppervalue = uppervalue;
	tmp_buffer->num = 1;
	num_buffer++;
	delta_sum += delta;
	if(u_index >= user_buffer.size())
		std::cout << "*************Overflow at line 253*************" << std::endl;
        if(cur_iter == 0) 
	{
		user_buffer.at(u_index) = tmp_buffer;
	}
	else
	{/* Check whether the previous buffer has to be released */
		std::cout << "User " << u_index << "current iteration-------" << cur_iter << std::endl;
		
		user_buffer.at(u_index)->num--;
		if( user_buffer.at(u_index)->num == 0)
		{
		    rec_buf * tmp_pointer = new rec_buf();
		    tmp_pointer = user_buffer.at(u_index);
		    user_buffer.at(u_index) = tmp_buffer;
		    free(tmp_pointer);
		    tmp_pointer = NULL;
		}
		else
		{
		    user_buffer.at(u_index) = tmp_buffer;
		}
	}
	free(tmp_buffer);
	tmp_buffer = NULL;

	std::cout << "Create buffer for usr " << u_index << std::endl;

	std::ofstream ouf("../../result/rec_result_inc", std::ofstream::out | std::ofstream::app);
	ouf << "User " << u_index << ", Value " << uppervalue << ", Maintain items: " << n_top+delta << std::endl;
	
	for(int i = 0; i < n_top+delta; i++)
	{
	    ouf << "<" << u_index << ", " << rec_vec.at(i).first << ", " << rec_vec.at(i).second << ", " << shortest.at(i).second << ">";
	}
	ouf << std::endl;
	ouf.close();

}


/* General purpose TopN program */
struct GeneralTopNProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {

  void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
    if (vertex.id() < M) {
	std::cout << "begin finding items for user " << vertex.id() << std::endl;
      /* Exclude neighbors */
      std::map<unsigned int, bool> h_neighbor;
	std::cout << "User " << vertex.id()  << " neighbors " << vertex.num_edges() <<std::endl;
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
	std::cout << "excluded neighbors for user " << vertex.id() << std::endl;
      //std::partial_sort(rec_vec.begin(), rec_vec.begin()+n_top, rec_vec.end(), sort_items_c);
	/* iteration >= 1, needs to search the maintainted items first */
	if(cur_iter > 0)
	{
		std::cout << "User " << vertex.id() << "********** iteration " << cur_iter << std::endl;
	  	std::cout << " User " <<vertex.id()   << "lf_previous.size " << lf_previous->size() << std::endl;

	  	vertex_data & pre_user = (*lf_previous)[vertex.id()];
		/* First check if the updated user is still within the theta range */
		std::cout << "User " << vertex.id() << "user_buffer size " << user_buffer.size() <<std::endl;
		rec_buf * pbuff = user_buffer.at(vertex.id());
		std::vector<unsigned int> & pre_vec = pbuff->item;
	  	std::cout << "User " << vertex.id() << "previous item set size " << pre_vec.size() << std::endl;
	  	int compare = 0;
	  	double new_dis = 0;
	  	double new_up = pbuff->uppervalue;
	  	new_up += update_uppervalue(theta, pre_user);	
	  	for (int i = 0; i < pre_vec.size(); i++)
	  	{
	  		new_dis = dot_prod(cur_user.pvec, latent_factors_inmem[pre_vec.at(i)].pvec);
			if(new_dis >= new_up)
		    		compare++; 
	  	}
	  	/* Chekc if updated user is still within the theta range */
	 	std::cout << "Check if user " <<vertex.id()  << " is still in range " << compare <<  std::endl;
	  	if(compare < n_top)
	  	{
			std::cout << "User " << vertex.id() << " Not in range " << top_one.size() << std::endl;
			unsigned int top1 = top_one.at(vertex.id());
			/* Check map to search for users with same top1 item */	
			int count = item_top.count(top1);
			std::cout << "Check if there is other user with same top 1 item as user " << vertex.id() << std::endl;
			if (count == 0)
			{
	  			std::partial_sort(rec_vec.begin(), rec_vec.begin()+num+n_top, rec_vec.end(), sort_items_c);
	  			/* Update top one item */
	  			top_one.at(vertex.id()) = rec_vec.at(0).first;
				/* Insert map */
				item_top.insert(std::pair<unsigned int, unsigned int>(rec_vec.at(0).first, vertex.id()));
				/* Update buffer pointer */
	  			find_top(cur_user, vertex.id(), rec_vec);
			}
 			else
			{	/* Check with every user who has the same top1 item */	
				std::multimap<unsigned int, unsigned int>::iterator it;
    				for (it=item_top.equal_range(top1).first; it!=item_top.equal_range(top1).second; ++it)
				{
					/* Check if the user still points to this buffer */
					std::cout << "User " << vertex.id() << " Check if the user in map is still valid" << std::endl;
					std::cout << "User " << vertex.id() << " " << (*it).second << std::endl;
					if(top_one.at( (*it).second ) != top1)
						continue;
					pre_user = (*lf_previous)[ (*it).second ];
					/* Check if the buffer satisfied the condition */
					pbuff = user_buffer.at( (*it).second);
					pre_vec = pbuff->item;
					compare = 0;
					new_dis = 0;
					new_up = pbuff->uppervalue;
					for (int i = 0; i < pre_vec.size(); i++)
	  				{
	  					new_dis = dot_prod(cur_user.pvec, latent_factors_inmem[pre_vec.at(i)].pvec);
	  					new_up += update_uppervalue(theta, cur_user);	
						if(new_dis >= new_up)
		    				compare++; 
	  				}
					std::cout << "User " << vertex.id() << " Check if the buffer can be accepted" << std::endl;
					if(compare >= n_top)
					{	/* Accept current buffer */
						user_buffer.at(vertex.id())->num--;
						std::cout << "User " << vertex.id() << " Check if old buffer needs to be released" << std::endl;
						if(user_buffer.at(vertex.id())->num == 0)
						{
							rec_buf * tmp_p = new rec_buf();
							tmp_p = user_buffer.at(vertex.id());
							user_buffer.at(vertex.id())= pbuff;
							free(tmp_p);
							tmp_p = NULL;
						}
						else
						{
							user_buffer.at(vertex.id()) = pbuff;
						}
						pbuff->num++;
						/* Update top1 item */
						top_one.at(vertex.id()) = update_topone(cur_user, pbuff->item);
						/* Insert map */
						item_top.insert(std::pair<unsigned int, unsigned int>(top_one.at(vertex.id()), vertex.id() )  );
						continue;
					}
					else
					{
	  					std::partial_sort(rec_vec.begin(), rec_vec.begin()+num+n_top, rec_vec.end(), sort_items_c);
	  					/* Update top one item */
	  					top_one.at(vertex.id()) = rec_vec.at(0).first;
						/* Insert map */
						item_top.insert(std::pair<unsigned int, unsigned int>(rec_vec.at(0).first, vertex.id()) );
						/* Update buffer pointer */
	  					find_top(cur_user, vertex.id(), rec_vec);
					}				
				}			
			}		

	  	}
	  	else
			std::cout << "User " << vertex.id() << "is still in range" << std::endl;

	}
	/* first iteration */
	else
	{	std::cout << "iteration " << cur_iter << " user " << vertex.id() << std::endl;
	  	std::partial_sort(rec_vec.begin(), rec_vec.begin()+num+n_top, rec_vec.end(), sort_items_c);
	  	/* Record top one item */
		if (vertex.id() >= top_one.size())
			std::cout << "*************Overflow at line 418*************" << std::endl;
	  	top_one.at(vertex.id()) = rec_vec.at(0).first;
	  	std::cout << rec_vec.at(0).first << ", " << vertex.id() << std::endl;
	  	item_top.insert(std::pair<unsigned int, unsigned int>(rec_vec.at(0).first, vertex.id()) );
	  	std::cout << "entering find_top function for user " << vertex.id() << std::endl;
	  	find_top(cur_user, vertex.id(), rec_vec);
	  	std::cout << "found topn for user " << vertex.id() << std::endl;
	}

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
			      std::vector<vertex_data> *llf_previous,
                              float (*prediction_func)(const vertex_data & user, 
                                                       const vertex_data & movie, 
                                                       float rating, double & prediction, 
                                                       void * extra),
			      int iter, 
                              std::vector<double> *lldelta, std::vector<double> *rrdelta,
                              std::vector<double> *llbound, std::vector<double> *rrbound) {
  std::cout << "running testing program" << std::endl;
  latent_factors = latent_factors_inmem;
  lf_previous = llf_previous;
  pprediction_func_test = prediction_func;
  ldelta = lldelta;
  rdelta = rrdelta;
  lbound = llbound;
  rbound = rrbound;
  cur_iter = iter;  
  num_buffer = 0;
  delta_sum = 0;
  if(cur_iter == 0)
  {
	rec_buf * tmp_buf = NULL;
	for(unsigned int i = 0; i < M; i++)
	{
	  top_one.push_back(0);
	  user_buffer.push_back(tmp_buf);
	}
  }  


  GeneralTopNProgram test_prog;
  engine->run(test_prog, 1);
  std::ofstream ofs("../../result/result_inc", std::ofstream::out | std::ofstream::app);

  ofs << "Iteration: " << cur_iter << " Buffer created: " << num_buffer << std::endl;
  ofs << "Total delta: " << delta_sum <<"Max delta: " << delta_max << " Min delta: " << delta_min << std::endl;
  ofs.close();

}



#endif //__MX_TOPN_ENGINE

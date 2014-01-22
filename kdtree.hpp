#ifndef _KDTREE_H__
#define _KDTREE_H__


#include <vector>
#include <list>
#include <queue>

struct kd_Node
{
	kd_Node* lchild;
	kd_Node* rchild;
	int llevel;
	std::vector<double> lbound;
	std::vector<double> rbound;
	std::list<kd_Node*> child;
	std::vector<std::pair<unsigned int, VertexDataType> > item;
	kd_Node() {}
};


int level, range;


kd_Node* build_kdtree(int level, std::vector<double> lbound, std::vector<double> rbound, std::vector<std::pair<unsigned int, vertex_data> > tids)
{
	kd_Node *root = new kd_Node();
	root->lbound = lbound;
	root->rbound = rbound;
	root->llevel = level;
	if(tids.size() <= 50 || level == range)
        {
                root->lchild = NULL;
                root->rchild = NULL;
                for(unsigned int i = 0; i < tids.size(); i++)
                        root->item.push_back(tids.at(i));	//the type of tids.at(i) is also a pair
                //set the range of the kd tree
                range = level;
                return root;
        }
	
	//building kd tree
	int att = level % D;
	double * copy = new double[tids.size()];
	double median;
	for(unsigned int i = 0; i < tids.size(); i++)
		copy[i] = tids.at(i).second.pvec[att];
	std::sort(copy, copy+tids.size());
	if(tids.size() % 2 == 1)
		median = copy[tids.size()/2 + 1];
	else
		median = (copy[tids.size()/2 + 1] + copy[tids.size()/2]) / 2;

	std::vector<std::pair<unsigned int, vertex_data> > ltids;
	std::vector<std::pair<unsigned int, vertex_data> > rtids;
	for(unsigned int i = 0; i < tids.size(); i++)
	{
		if(tids.at(i).second.pvec[att] <= median)
		{
			ltids.push_back(tids.at(i));
		}
		else
		{
			rtids.push_back(tids.at(i));
		}
	}
	delete [] copy;
	double tmp_val = rbound.at(att);
	rbound.at(att) = median;
	root->lchild = build_kdtree(level+1, lbound, rbound, ltids);
	lbound.at(att) = median;
	rbound.at(att) = tmp_val;
	root->rchild = build_kdtree(level+1, lbound, rbound, rtids);
	//child nodes list
	root->child.push_back(root->lchild);
	root->child.push_back(root->rchild);
	return root;
}

void print_kd (kd_Node* root)
{
	std::cout << "kd node level:  " << root->llevel << std::endl;
	std::cout << "lbound: < ";
	for(int i = 0; i < D; i++)
	    std::cout << root->lbound.at(i) << ", ";
	std::cout << " >" << std::endl;
	std::cout << "rbound: < ";
        for(int j = 0; j < D; j++)
            std::cout << root->rbound.at(j) << ", ";
        std::cout << " >" << std::endl;

	if(root -> lchild)
	    print_kd(root->lchild);
	if(root -> rchild)	
	    print_kd(root->rchild);
	if(root->lchild == NULL && root->rchild == NULL)
	{
	    std::cout << "number of vertex in leaf node: "<< root->item.size() << ". They are: ";
	    for(unsigned int i=0; i < root->item.size(); i++)
		std::cout << root->item.at(i).first << ", ";
	    std::cout << std::endl;
	}
	return ;
}


kd_Node* init_kdtree(std::vector<vertex_data> *latent_factors_inmem)
{
	std::vector<std::pair<unsigned int, vertex_data> > item;
	double copy[N];
	for(unsigned int i = M; i < M + N; i++)
	{
		item.push_back(std::make_pair(i, latent_factors_inmem->at(i)));
		copy[i - M] = latent_factors_inmem->at(i).pvec[0];
	}
	std::sort(copy, copy+N);
	level = 0;
	range = -1;

	std::vector<double> lbound (D, 0);
	std::vector<double> rbound (D, 1);
	lbound.at(0) = copy[0];
	rbound.at(0) = copy[N-1];
	for(int j = 1; j < D; j++)
        {
            for(unsigned int i = M; i < M + N; i++)
                copy[i - M] = latent_factors_inmem->at(i).pvec[j];
            std::sort(copy, copy+N);
            lbound.at(j) = copy[0];
            rbound.at(j) = copy[N-1];
        }
	
	kd_Node* init = build_kdtree(level, lbound, rbound, item);
//	print_kd(init);

	return init;
}


#endif

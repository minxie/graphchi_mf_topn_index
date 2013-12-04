#ifndef _MTREE_H__
#define _MTREE_H__
  
#include <vector>
#include <list>
#include <queue>
#include <cmath>

struct m_Tree
{
    m_Tree* lchild;
    m_Tree* rchild;
    std::vector<double> center;
    double radius;
    int llevel;

    std::list<m_Tree*> child;
    std::vector<std::pair<unsigned int, VertexDataType> > item; 
    m_Tree() {}
};

int level, range;

std::vector<double> find_mean(std::vector<std::pair<unsigned int, vertex_data> > tids)
{
    std::vector<double> mean (D, 0);
    for(int j = 0; j < D; j++)
    {
        for (unsigned int i = 0; i < tids.size(); i++)
        {
            mean.at(j) += tids.at(i).second.pvec[j];
        }
        mean.at(j) = mean.at(j)/tids.size();
    }
    return mean;
}

double find_radius(std::vector<double> center, std::vector<std::pair<unsigned int, vertex_data> > tids)
{
    double radius = 0;
    double distance;
    std::vector<double> diff;
    for (unsigned int i = 0; i < tids.size(); i++)
    {
        distance = 0;
        for(int j = 0; j < D; j++)
        {
            distance += (center.at(j) - tids.at(i).second.pvec[j]) * (center.at(j) - tids.at(i).second.pvec[j]); 
        }
	distance = sqrt(distance);
        if(radius < distance)
            radius = distance;
    }
    return radius;
}


std::pair<vec, double> SplitMetricTree(std::vector<std::pair<unsigned int, vertex_data> > tids)
{
    unsigned int num = tids.size();
    srand((unsigned int )time(NULL));
    unsigned int rand_node = rand() % num;
    unsigned int anode, bnode;
    double distance;
    double max = -0xffff;
    for (unsigned int i = 0; i < num; i++)
    {     
        distance = 0;
        for(int j = 0; j < D; j++)
        {
            distance += (tids.at(rand_node).second.pvec[j] - tids.at(i).second.pvec[j]) * (tids.at(rand_node).second.pvec[j] - tids.at(i).second.pvec[j]);
        }
        if(max < distance)
        {
            max = distance;
            anode = i;
        }
    }
    max = -0xffff;
    for (unsigned int i = 0; i < num; i++)
    {
        distance = 0;
        for(int j = 0; j < D; j++)
        {
            distance += (tids.at(anode).second.pvec[j] - tids.at(i).second.pvec[j]) * (tids.at(anode).second.pvec[j] - tids.at(i).second.pvec[j]);
        }
        if(max < distance)
        {
            max = distance;
            bnode = i;
        }
    }
    double len_a = dot_prod(tids.at(anode).second.pvec, tids.at(anode).second.pvec);
    double len_b = dot_prod(tids.at(bnode).second.pvec, tids.at(bnode).second.pvec);
    vec w(D);
    for(int i = 0; i < D; i++)
	w[i] = (tids.at(bnode).second.pvec[i] - tids.at(anode).second.pvec[i]);
    double b = (len_a - len_b) / 2;
    return std::make_pair(w, b);

}

m_Tree* build_mtree (int level, std::vector<std::pair<unsigned int, vertex_data> > tids)
{

    m_Tree* root = new m_Tree();
    root->center = find_mean(tids);
    root->radius = find_radius(root->center, tids);
    root->llevel = level;
    if (tids.size() <= 50)
    {
        root->lchild = NULL;
        root->rchild = NULL;
        range = level;
        for(unsigned int i = 0; i < tids.size(); i++)
            root->item.push_back(tids.at(i));
        return root;
    }
    else
    {
	std::pair<vec, double> split = SplitMetricTree(tids);
        std::vector<std::pair<unsigned int, vertex_data> > ltids;
        std::vector<std::pair<unsigned int, vertex_data> > rtids;
	vec w (split.first);
	double b = split.second;
        for (unsigned int i = 0; i < tids.size(); i++)
        {
	    double tmp = dot_prod(w, tids.at(i).second.pvec);
            if (tmp + b <= 0 )
                ltids.push_back(tids.at(i));
            else
                rtids.push_back(tids.at(i));
        }
        root->lchild = build_mtree(level+1, ltids);
        root->rchild = build_mtree(level+1, rtids);
        root->child.push_back(root->lchild);
        root->child.push_back(root->rchild);
	return root;
    }
}

void print_m (m_Tree* root)
{
    std::cout << "m node level:  " << root->llevel << std::endl;
    std::cout << "center: <";
    for(int i = 0; i < D; i++)
        std::cout << root->center.at(i) << ", ";
    std::cout << ">;  ";
    std::cout << "radius: " << root->radius << std::endl;
    if(root->lchild)
	print_m(root->lchild);
    if(root->rchild)
	print_m(root->rchild);

    if(root->lchild == NULL && root->rchild == NULL)
    {
        std::cout << "number of vertex in leaf node: "<< root->item.size() << ". They are: ";
        for(unsigned int i=0; i < root->item.size(); i++)
            std::cout << root->item.at(i).first << ", ";
        std::cout << std::endl;
    } 
    return ;
}

m_Tree* init_mtree(std::vector<vertex_data> *latent_factors_inmem)
{
    std::vector<std::pair<unsigned int, vertex_data> > item;
    for(unsigned int i = M; i < M + N; i++)
    {
        item.push_back(std::make_pair(i, latent_factors_inmem->at(i)));
    }
    level = 0;
    range = -1;
    m_Tree* init = build_mtree(level, item);
//    print_m(init);
    return init;
}

#endif                                                                                                                                    

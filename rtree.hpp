#ifndef __RTREE_H
#define __RTREE_H

#include <algorithm>
#include <cmath>
#include <vector>
using std::vector;

bool sort_items_rtree(std::pair<double, unsigned int> a, std::pair<double, unsigned int> b) {
  return a.first > b.second;
}

class RTreeNode {
public:
  RTreeNode() {
  }
  ~RTreeNode() {
  }
public:
  vector<RTreeNode*> _children;
  vec _lbound;
  vec _rbound;
  int _count;
  vector<int> _tids;
};

class RTree {
public:
  RTree() {
  }

  ~RTree() {
    // Release memory for rtree nodes
    recursive_delete(_root);
  }

private:
  void recursive_delete(RTreeNode* node) {
    if (node->_children.size() != 0) {
      for (unsigned int i = 0; i < node->_children.size(); i++) {
        recursive_delete(node->_children[i]);
      }
    }
    delete node;
  }

public:
  void print_rtree() {
    print_rtree_node(0, _root);
  }

  void print_rtree_node(int level, RTreeNode *node) {
    for (int i = 0; i < level; i++)
      std::cout << " ";
    std::cout << node->_count << std::endl;
    for (int i = 0; i < level; i++)
      std::cout << " ";
    for (int i = 0; i < D; i++)
      std::cout << node->_lbound[i] << " ";
    std::cout << std::endl;
    for (int i = 0; i < level; i++)
      std::cout << " ";
    for (int i = 0; i < D; i++)
      std::cout << node->_rbound[i] << " ";
    std::cout << std::endl;
    if (node->_tids.size() != 0) {
      for (unsigned int i = 0; i < node->_tids.size(); i++)
        std::cout << node->_tids[i] << " ";
      std::cout << std::endl;
    }
  }

  void build_rtree(vector<vertex_data> *latent_factors_inmem) {
    _data = latent_factors_inmem;
    
    vec lbound(D), rbound(D);
    for (int i = 0; i < D; i++) {
      lbound[i] = 10000;
      rbound[i] = -10000;
    }

    double max_val = -10000;
    int max_idx = 0;
    for (int i = 0; i < D; i++) {
      if (rbound[i] - lbound[i] > max_val) {
        max_idx = i;
      }
    }

    int NODE_SIZE = 50;
    vector<std::pair<double, unsigned int> > items;
    for(unsigned int i = M; i < M + N; i++) {
      items.push_back(std::make_pair(latent_factors_inmem->at(i).pvec[max_idx], i));
    }
    std::sort(items.begin(), items.end(), sort_items_rtree);
    
    vector<RTreeNode*> *cnode_list = new vector<RTreeNode*>();
    vector<RTreeNode*> *nnode_list = new vector<RTreeNode*>();

    // Initiate leaf nodes
    int num_nodes = (int)(N / (float)(NODE_SIZE) + 0.5);
    for (int ni = 0; ni < num_nodes; ni++) {
      RTreeNode *node = new RTreeNode();

      node->_lbound = zeros(D);
      node->_rbound = zeros(D);
      for (int i = 0; i < D; i++) {
        node->_lbound[i] = 10000;
        node->_rbound[i] = -10000;
      }

      node->_count = 0;
      int tmpN = (int)N;
      for (int i = ni * NODE_SIZE; i < std::min(ni * NODE_SIZE + NODE_SIZE + 1, tmpN); i++) {
        node->_tids.push_back(items[i].second);
        ++node->_count;
        for (int j = 0; j < D; j++) {
          if (latent_factors_inmem->at(items[i].second).pvec[j] < node->_lbound[j])
            node->_lbound[j] = latent_factors_inmem->at(items[i].second).pvec[j];
          if (latent_factors_inmem->at(items[i].second).pvec[j] > node->_rbound[j])
            node->_rbound[j] = latent_factors_inmem->at(items[i].second).pvec[j];
        }
      }
      
      cnode_list->push_back(node);
    }

    // Initiate intermediate nodes
    while (cnode_list->size() > NODE_SIZE) {
      num_nodes = (int)(cnode_list->size() / (float)(NODE_SIZE) + 0.5);
      for (int ni = 0; ni < num_nodes; ni++) {
        RTreeNode *node = new RTreeNode();

        node->_lbound = zeros(D);
        node->_rbound = zeros(D);
        for (int i = 0; i < D; i++) {
          node->_lbound[i] = 10000;
          node->_rbound[i] = -10000;
        }

        node->_count = 0;
        tmpN = (int)(cnode_list->size());
        for (int i = ni * NODE_SIZE; i < min(ni * NODE_SIZE + NODE_SIZE + 1, tmpN); i++) {
          node->_children.push_back((*cnode_list)[i]);
          node->_count += (*cnode_list)[i]->_count;
          for (int j = 0; j < D; j++) {
            if ((*cnode_list)[i]->_lbound[j] < node->_lbound[j])
              node->_lbound[j] = (*cnode_list)[i]->_lbound[j];
            if ((*cnode_list)[i]->_rbound[j] > node->_rbound[j])
              node->_rbound[j] = (*cnode_list)[i]->_rbound[j];
          }
        }

        nnode_list->push_back(node);
      }
      cnode_list->clear();
      vector<RTreeNode*> *tmp_list = cnode_list;
      cnode_list = nnode_list;
      nnode_list = tmp_list;
    }

    // Final root setup
    _root = new RTreeNode();
    _root->_lbound = zeros(D);
    _root->_rbound = zeros(D);
    for (int i = 0; i < D; i++) {
      _root->_lbound[i] = 10000;
      _root->_rbound[i] = -10000;
    }
    
    _root->_count = 0;        
    for (int i = 0; i < cnode_list->size(); i++) {
      _root->_children.push_back((*cnode_list)[i]);
      _root->_count += (*cnode_list)[i]->_count;
      for (int j = 0; j < D; j++) {
        if ((*cnode_list)[i]->_lbound[j] < _root->_lbound[j])
          _root->_lbound[j] = (*cnode_list)[i]->_lbound[j];
        if ((*cnode_list)[i]->_rbound[j] > _root->_rbound[j])
          _root->_rbound[j] = (*cnode_list)[i]->_rbound[j];
      }
    }

    delete cnode_list;
    delete nnode_list;
  }

public:
  RTreeNode* _root;
  vector<vertex_data> * _data;
};




#endif

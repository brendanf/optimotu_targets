#include <Rcpp.h>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <forward_list>

struct cluster {
  int min_d, index;
  std::vector<int>* members;
  cluster *parent;
  std::unordered_set<cluster*> *children;
};

struct cluster_compare{
  bool operator() (const cluster* c1, const cluster* c2) const {
    if (c1->parent == nullptr) return false;
    if (c2->parent == nullptr) return true;
    return c1->parent->min_d < c2->parent->min_d;
  }
};

// this only reassigns within the clust_array,
// it does not touch the clusters themselves
cluster* reassign(cluster *c_old, cluster *c_new, cluster** clust_array,
                  const int m, const int n, const int min_d, int max_d = -1) {
  cluster * c2p = c_old->parent;
  if (max_d == -1) max_d = c2p == nullptr ? m : c_old->parent->min_d;
  if (c_old->members == nullptr) {
    // above the tip level, we only know our children, not our tips.
    // so recurse to children
    for (auto c : *c_old->children) {
      reassign(c, c_new, clust_array, m, n, min_d, max_d);
    }
  } else {
    for (int tip : *c_old->members) {
      for (int i = min_d; i < max_d; i++) {
        clust_array[i * n + tip] = c_new;
      }
    }
  }
  return c2p;
}

// this only reassigns within the clust_array,
// it does not touch the clusters themselves, except that it deletes c_old
cluster* reassign_all(cluster* c_old, cluster* c_new, cluster** clust_array,
                             const int m, const int n) {
  cluster * c2p = reassign(c_old, c_new, clust_array, m, n, c_old->min_d);
  delete c_old->members;
  delete c_old;
  return c2p;
}

// [[Rcpp::export]]
Rcpp::List single_linkage(const std::string file, const int n, const float dmin, const float dmax, const float dstep) {
  const int m = (int) floorf((dmax - dmin)/dstep) + 1;
  // pointer to which cluster a sequence belongs to at each threshold
  cluster **clust_array = new cluster*[m*n];
  cluster *c, *c1, *c2;
  int i, j, j1, j2;
  for (j = 0; j < n; j++) {
    c = new cluster;
    c->min_d = 0;
    c->index = j;
    c->members = new std::vector<int>{j};
    c->parent = nullptr;
    c->children = nullptr;
    for (i = 0; i < m; i++) {
      clust_array[i*n + j] = c;
    }
  }
  std::ifstream infile(file);
  int seq1, seq2;
  float dist;
  while(infile >> seq1 >> seq2 >> dist) {
    if (seq1 == seq2) continue;
    i = std::max((int) floor((dist - dmin) / dstep), 0);
    c1 = clust_array[i * n + seq1];
    c2 = clust_array[i * n + seq2];
    while (c1 != c2) {
      /* we need to merge c1 and c2, and all their parents, starting at i and
       going up to m. */
      // clusters are not strictly sorted, but be sure the smallest member
      // comes first.
      j1 = c1->index;
      j2 = c2->index;
      if (j1 > j2) {
        j = j2;
        j2 = j1;
        j1 = j;
        c = c2;
        c2 = c1;
        c1 = c;
      }
      cluster *cnew;
      cluster *c1p = c1->parent;
      cluster *c2p = c2->parent;
      if (i == 0) {
        // we are merging at the tip level
        // here we need to actually modify members
        cnew = c1;
        std::copy(
          c2->members->end(),
          c2->members->begin(),
          c1->members->end()
        );
        reassign_all(c2, c1, clust_array, m, n);
      } else {
        if (i == c1->min_d) {
          // we don't need to create a new cluster for c1, we can just modify this
          cnew = c1;
          cnew->children->insert(c2);
          c2p->children->erase(c2);
          c2p->children->insert(cnew);
          if (i == c2->min_d) {
            // oldchild is not needed anymore, so we need to update the parent of its children
            cnew->children->reserve(cnew->children->size() + c2->children->size());
            for (cluster * child : *c2->children) {
              child->parent = cnew;
              cnew->children->insert(child);
            }
            reassign_all(c2, cnew, clust_array, m, n);
          } else {
            cnew->children->insert(c2);
            reassign(c2, cnew, clust_array, m, n, i);
          }
        } else if (i == c2->min_d) {
          cnew = c2;
          cnew->index = j1;
          cnew->children->insert(c1);
          c1p->children->erase(c1);
          c1p->children->insert(cnew);
          reassign(c1, cnew, clust_array, m, n, i);
        } else {
          //create a new cluster as parent of current c1 and c2
          cnew = new cluster;
          cnew->min_d = i;
          cnew->index = j1;
          cnew->children = new std::unordered_set<cluster*> {{c1, c2}};
          cnew->members = nullptr;
          c1p->children->erase(c1);
          c1p->children->insert(cnew);
          reassign(c1, cnew, clust_array, m, n, i);
          c2p->children->erase(c2);
          c2p->children->insert(cnew);
          reassign(c2, cnew, clust_array, m, n, i);
        }
        if (c1p == nullptr) {
          if (c2p == nullptr) {
            cnew->parent = nullptr;
            i = m;
          } else {
            cnew->parent = c1p;
            c1 = cnew;
            c2 = c2p;
            i = c2p->min_d;
          }
        } else if (c2p == nullptr || c1p->min_d <= c2p->min_d) {
          cnew->parent = c2p;
          c1 = c1p;
          c2 = cnew;
          i = c1p->min_d;
        } else {
          cnew->parent = c1p;
          c1 = cnew;
          c2 = c2p;
          i = c2p->min_d;
        }
      }
      // members are not tracked during the cluster construction phase;
    // this gets populated when we create the output.
    }
  }
  auto tip_clusters = std::set<cluster*, cluster_compare>(clust_array, clust_array + n);
  int n_clusters = tip_clusters.size();
  auto clusters = std::forward_list<cluster*>(tip_clusters.begin(), tip_clusters.end());
  Rcpp::List out = Rcpp::List::create(m);
  for (i = 0; i < m; i++) {
    Rcpp::List rclusters = Rcpp::List::create();
    for (const auto& c : clusters) {
      Rcpp::IntegerVector members = Rcpp::IntegerVector::import(c->members->begin(), c->members->end());
      rclusters.push_back(members);
    }
    out[m] = rclusters;
  }
  return out;
}

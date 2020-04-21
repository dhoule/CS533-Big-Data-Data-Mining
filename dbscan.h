/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   Files: mpi_main.cpp clusters.cpp  clusters.h utils.h utils.cpp          */
/*        dbscan.cpp dbscan.h kdtree2.cpp kdtree2.hpp          */
/*      geometric_partitioning.h geometric_partitioning.cpp  */
/*                                 */
/*   Description: an mpi implementation of dbscan clustering algorithm       */
/*        using the disjoint set data structure        */
/*                                                                           */
/*   Author:  Md. Mostofa Ali Patwary                                        */
/*            EECS Department, Northwestern University                       */
/*            email: mpatwary@eecs.northwestern.edu                          */
/*                                                                           */
/*   Copyright, 2012, Northwestern University                                */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/*   Please cite the following publication if you use this package       */
/*                       */
/*   Md. Mostofa Ali Patwary, Diana Palsetia, Ankit Agrawal, Wei-keng Liao,  */
/*   Fredrik Manne, and Alok Choudhary, "A New Scalable Parallel DBSCAN      */
/*   Algorithm Using the Disjoint Set Data Structure", Proceedings of the    */
/*   International Conference on High Performance Computing, Networking,     */
/*   Storage and Analysis (Supercomputing, SC'12), pp.62:1-62:11, 2012.      */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _DBSCAN_
#define _DBSCAN_

#include "utils.h"
#include "clusters.h"

namespace NWUClustering {

  class ClusteringAlgo : public Clusters {
  public:
    ClusteringAlgo(){ }
    virtual ~ClusteringAlgo();

    void set_dbscan_params(double eps, int minPts, double seed_percentage);
    void  get_clusters_distributed();
    void  writeCluster_distributed(string outfilename);

    void  trivial_compression(vector <int>* data, vector < vector <int> >* parser);
    void  trivial_decompression(vector <int>* data);

    void modify_status_vectors(int pid, kdtree2_result_vector &ne, kdtree2_result_vector &ne_outer);
    int binarySearch(vector<int>& dirty, int l, int r, int needle, int& flag);
    void getSeeds();

  public:
    
    double  m_epsSquare; // AKA radius. It is the square of the "radius" given by the user
    int   m_minPts; // The minimum number of points, given by the user, to start a cluster
    int   m_compression;
    double m_perc_of_dataset; // The percentage of the points each node is to use, given by the user. Default value is 1.0(all points)

    vector <int> m_parents; // Elements hold the pointers of the clustering tree
    vector <int> m_parents_pr; // Elements hold the pointers for which node the point is in
    vector <int> m_child_count; // number of "children" a possible centroid has

    vector <int> neededIndices; // Values are KdTree indices of points to be used.

    vector <int> m_member; // Values are either 0 or 1. It's size = size_of(m_pts.m_i_num_points). Used to determine if a border point or not.
    vector <int> m_corepoint; // Values are either 0 or 1. It's size = size_of(m_pts.m_i_num_points). Used to determine center points.

    vector <int> triage; // local points that have been found. Deals with `ne` vector. The "growing set."
    vector <int> assessed; // local points that have been checked as center points, or have been found in the intersection of 2/+ neighborhoods.
    vector <int> assessed_outer; // remote points that have been seen already. Deals with the `ne_outer` vector.
  };  

  void run_dbscan_algo_uf_mpi_interleaved(ClusteringAlgo& dbs); // union find dbscan algorithm using mpi with interleaved communication
  void get_neighborhood_points(ClusteringAlgo& dbs, kdtree2_result_vector &ne, kdtree2_result_vector &ne_outer, int pid);
  void unionize_neighborhood(ClusteringAlgo& dbs, kdtree2_result_vector &ne, kdtree2_result_vector &ne_outer, int pid, vector < vector <int > >* p_cur_insert);
};

#endif

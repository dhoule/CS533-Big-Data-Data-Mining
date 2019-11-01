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

#include "dbscan.h"

namespace NWUClustering {

  /*
    eps = epsilon/radious
    minPts = minimum number of points need to make a cluster
  */
  void ClusteringAlgo::set_dbscan_params(double eps, int minPts) {
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    if(rank == proc_of_interest) cout << "in dbscan line: 36" << " in ClusteringAlgo::set_dbscan_params" << endl;
    m_epsSquare =  eps * eps;
    m_minPts =  minPts;
    m_messages_per_round = -1; // always -1
    m_compression = 0; // can set to 1 if want to compress specailly in the first round of communication
    if(rank == proc_of_interest) cout << "in dbscan line: 41" << " eps: " << eps << endl;
    if(rank == proc_of_interest) cout << "in dbscan line: 42" << " minPts: " << minPts << endl;
    if(rank == proc_of_interest) cout << "in dbscan line: 43" << " m_epsSquare: " << m_epsSquare << endl;
    if(rank == proc_of_interest) cout << "in dbscan line: 44" << " m_minPts: " << m_minPts << endl;
    if(rank == proc_of_interest) cout << "in dbscan line: 45" << " m_messages_per_round: " << m_messages_per_round << endl;
    if(rank == proc_of_interest) cout << "in dbscan line: 46" << " m_compression: " << m_compression << endl;
  }

  // Destructor
  ClusteringAlgo::~ClusteringAlgo() {
    m_noise.clear();
    m_visited.clear();
    m_parents.clear();
    m_parents_pr.clear();
    m_child_count.clear();

    m_corepoint.clear();
    m_member.clear();
  }

  /*
    Function keeps track of the time it takes to complete for some reason.
    called in run_dbscan_algo_uf_mpi_interleaved()
  */
  void ClusteringAlgo::trivial_decompression(vector <int>* data, int nproc, int rank, int round, double& dcomtime) {
    double start = MPI_Wtime();
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    if(rank == proc_of_interest) cout << "in dbscan line: 69" << " in ClusteringAlgo::trivial_decompression" << endl;
    vector <int> parser;
    // allocates a MINIMUM amount of memory, the size of 'data'
    parser.reserve((*data).size());
    // assign "data" to another variable
    parser = (*data); if(rank == proc_of_interest) cout << "in dbscan line: 74" << " copying parser = (*data);" << endl;
    // removes all elements, destroying them. The size becomes 0.
    (*data).clear(); if(rank == proc_of_interest) cout << "in dbscan line: 76" << " clearing (*data).clear();" << endl;
    // declaring an INT automatically initializes it to 0...
    int pid_count = parser[0], pos, i, j, pid, npid, npid_count;
    if(rank == proc_of_interest) cout << "in dbscan line: 79" << " pid_count: " << pid_count << endl;
    pos++;// the value becomes 1
    // TODO I have no idea what this does. Need to have outputs to determine what is going on;
      // besides rebuilding the 'data' variable
    while(pid_count > 0) {
      pid = parser[pos++];
      npid_count = parser[pos++];
      if(rank == proc_of_interest) cout << "in dbscan line: 86" << " pid: " << pid << endl;
      if(rank == proc_of_interest) cout << "in dbscan line: 87" << " npid_count: " << npid_count << endl;
      for(j = 0; j < npid_count; j++) {
        // push_back() adds elements to the "back" of the vector
        (*data).push_back(pid); 
        (*data).push_back(parser[pos++]);
        if(rank == proc_of_interest) cout << "in dbscan line: 92" << " pid: " << pid << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 93" << " parser[pos - 1]: " << parser[pos - 1] << endl;
      }

      pid_count--;
    }
    // removes all elements, destroying them. The size becomes 0.
    parser.clear();

    double stop = MPI_Wtime();
    // increment the time counter
    dcomtime += (stop - start);
  }
  
  // called in run_dbscan_algo_uf_mpi_interleaved()
  void ClusteringAlgo::trivial_compression(vector <int>* data, vector < vector <int> >* parser, int nproc, int rank, int round, double& comtime, double& sum_comp_rate) {
    // get the starting time before doing anything in this function
    double start = MPI_Wtime();
    double org = 0, comp = 0;
    int pairs, pid, npid, i, j, pid_count, npid_count;
    int irank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &irank); if(irank == proc_of_interest) cout << "in dbscan line: 113" << " in ClusteringAlgo::trivial_compression" << endl;
    if(irank == proc_of_interest) cout << "in dbscan line: 114" << " nproc: " << nproc << endl;
    if(irank == proc_of_interest) cout << "in dbscan line: 115" << " rank: " << rank << endl;
    if(irank == proc_of_interest) cout << "in dbscan line: 116" << " round: " << round << endl;
    if(irank == proc_of_interest) cout << "in dbscan line: 117" << " comtime: " << comtime << endl;
    if(irank == proc_of_interest) cout << "in dbscan line: 118" << " sum_comp_rate: " << sum_comp_rate << endl;
    // The number of "pairs" in the "data" vector
    pairs = (*data).size()/2; // TODO this must be why the data must be a factor of 2???
    if(irank == proc_of_interest) cout << "in dbscan line: 121" << " pairs: " << pairs << endl;
    // TODO don't know what "org" is supposed to stand for, 'original' maybe
    org = (*data).size();   
    if(irank == proc_of_interest) cout << "in dbscan line: 124" << " org: " << org << endl;
    // loop over 'data' and add elements to the back of a vector 'parser'[pid]
    for(i = 0; i < pairs; i++) {
      pid = (*data)[2 * i];
      npid = (*data)[2 * i + 1];
      if(irank == proc_of_interest) cout << "in dbscan line: 129" << " pid: " << pid << endl;
      if(irank == proc_of_interest) cout << "in dbscan line: 130" << " npid: " << npid << endl;
      (*parser)[pid].push_back(npid);
    }
    // empty the 'data' vector, and set the size to 0
    (*data).clear();
    // inititalize 'pid_count' to 0, and add it to the back of the 'data' vector
    pid_count = 0;
    (*data).push_back(pid_count); // uniques pids, should update later
    // 'm_pts' is the current cluster's struct object???
    // Loop rebuilds the 'data' vector and clears out dimensions of the 'parser' vector
    if(irank == proc_of_interest) cout << "in dbscan line: 140" << " m_pts->m_i_num_points: " << m_pts->m_i_num_points << endl;
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      npid_count = (*parser)[i].size();
      if(irank == proc_of_interest) cout << "in dbscan line: 143" << " npid_count: " << npid_count << endl;
      if(npid_count > 0) {
        (*data).push_back(i);
        (*data).push_back(npid_count);
        if(irank == proc_of_interest) cout << "in dbscan line: 147" << " i: " << i << endl;
        if(irank == proc_of_interest) cout << "in dbscan line: 148" << " npid_count: " << npid_count << endl;
        for(j = 0; j < npid_count; j++) {
          (*data).push_back((*parser)[i][j]);  
          if(irank == proc_of_interest) cout << "in dbscan line: 151" << " (*parser)[i][j]: " << (*parser)[i][j] << endl;        
        }
        pid_count++;
        if(irank == proc_of_interest) cout << "in dbscan line: 154" << " pid_count: " << pid_count << endl; 
        (*parser)[i].clear();
      }
    }
    // assign the first element the value of unique cluster IDs
    (*data)[0] = pid_count; if(irank == proc_of_interest) cout << "in dbscan line: 159" << " pid_count: " << pid_count << endl; 
    // "computed" is the new size of the 'data' vector
    comp = (*data).size(); if(irank == proc_of_interest) cout << "in dbscan line: 159" << " comp: " << comp << endl; 
    // Get the stopping time of the function, increase the incrimenter
    double stop = MPI_Wtime();
    comtime += (stop - start);
    // increase the "speed up" incrementer
    sum_comp_rate += (comp / org);
  }

  // called in run_dbscan_algo_uf_mpi_interleaved(), but is commented out.
    // the comments above it talk about "compression"
    // seems to be older code of trivial_compression(), but wasn't taken out.
  /*
  void ClusteringAlgo::convert(vector < vector <int> >* data, int nproc, int rank, int round) {
    int j, tid, size, pid, v1, v2, pairs, count;
    vector < vector <int> > parser;
    vector <int> init, verify;
    int min, max;

    for(tid = 0; tid < nproc; tid++) {
      pairs = (*data)[tid].size()/2;
      
      verify.resize(2 * pairs, -1);

      if(pairs == 0)
        continue;

      min = m_pts->m_i_num_points;
      max = -1;
      for(pid = 0; pid < pairs; pid++) {
        if((*data)[tid][2 * pid] < min)
          min = (*data)[tid][2 * pid];
        
        if((*data)[tid][2 * pid] > max)
          max  = (*data)[tid][2 * pid];

        verify[2 * pid] = (*data)[tid][2 * pid];
        verify[2 * pid + 1] = (*data)[tid][2 * pid + 1];
      }

      init.clear();
      parser.resize(max - min + 1, init);
      
      for(pid = 0; pid < pairs; pid++) {
        v2 = (*data)[tid].back();
        (*data)[tid].pop_back();

        v1 = (*data)[tid].back();
                    (*data)[tid].pop_back();

        parser[v1 - min].push_back(v2);
      }

      count = 0;
      (*data)[tid].push_back(-1); // insert local root count later

      for(pid = min; pid <= max; pid++) {
        size = parser[pid - min].size();
        if(size > 0) {
          count++;
          (*data)[tid].push_back(pid);
          (*data)[tid].push_back(size);

          for(j = 0; j < size; j++) {
            (*data)[tid].push_back(parser[pid - min].back());
            parser[pid - min].pop_back();
          }
        }
      }

      (*data)[tid][0] = count;
      parser.clear();

      count = (*data)[tid][0];
      int k = 1, size, u = 0;
      for(pid = 0; pid < count; pid++) {   
        v2 = (*data)[tid][k++];
        size = (*data)[tid][k++];

        for(j = k; j < size; j++, k++) {
          v1 = (*data)[tid][k++];

          if(v2 != verify[u++])
            cout << "SOMETHING IS WRONG" << endl;

          if(v1 != verify[u++])
            cout << "SOMETHING IS WRONG" << endl;
        } 
      }
    }
  } */

  // called in mpi_main.cpp.
    // The function merges Points from other nodes
    // TODO this needs to be looked at in more detail
  void ClusteringAlgo::get_clusters_distributed() {
    // Determine the current node's rank within the cluster, and the size of the cluster itself
    int rank, nproc, i;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if(rank == proc_of_interest) cout << "in dbscan line: 260" << " in ClusteringAlgo::get_clusters_distributed" << endl;
    // get the root all the local points first
    // store in a message buffer in case need to ask other processors

    vector < vector <int > > merge_received;
    vector < vector <int > > merge_send1;
    vector < vector <int > > merge_send2;
    vector <int> init;
    // resize the vectors to the size of the number of nodes,
      // initialized with the value of the 'init' vector.
    merge_received.resize(nproc, init);
    merge_send1.resize(nproc, init);
    merge_send2.resize(nproc, init);

    int pid;
    // loop over the other dimensions of the vectors, 
      // resizing them to the minimum length of the number of points that 'm_pts' has.
    for(pid = 0; pid < nproc; pid++) {
      if(rank == proc_of_interest) cout << "in dbscan line: 278" << " pid: " << pid << endl;
      merge_received[pid].reserve(m_pts->m_i_num_points);
      merge_send1[pid].reserve(m_pts->m_i_num_points);
      merge_send2[pid].reserve(m_pts->m_i_num_points);
    }
    // These vectors, together, seems to be used to perform a bubble sort type of operation
      // 'p_cur_send' and 'p_cur_insert' have actual uses.
    vector < vector <int > >* pswap;
    vector < vector <int > >* p_cur_send;
    vector < vector <int > >* p_cur_insert;
    // assign pointers???
    p_cur_send = &merge_send1;
    p_cur_insert = &merge_send2;
    if(rank == proc_of_interest) cout << "in dbscan line: 291" << " p_cur_send: " << p_cur_send << endl;
    if(rank == proc_of_interest) cout << "in dbscan line: 292" << " p_cur_insert: " << p_cur_insert << endl;
    // I don't know why reserve() and resize() are called together???
    m_child_count.reserve(m_pts->m_i_num_points);
    m_child_count.resize(m_pts->m_i_num_points, 0);   
    if(rank == proc_of_interest) cout << "in dbscan line: 296" << " m_pts->m_i_num_points: " << m_pts->m_i_num_points << endl;
    int root, local_continue_to_run = 0, global_continue_to_run;
    // loop over the points
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      // find the point containing i
      root = i;
      if(rank == proc_of_interest) cout << "in dbscan line: 302" << " root: " << root << endl;
      // continue the loop as long as the node's ID matches the parent's ID???
      while(m_parents_pr[root] == rank) {
        if(m_parents[root] == root) {
          if(rank == proc_of_interest) cout << "in dbscan line: 306" << " root: " << root << endl;
          break;
        }
        root = m_parents[root];
        if(rank == proc_of_interest) cout << "in dbscan line: 310" << " root: " << root << endl;
      }
      
      if(m_parents[root] == root && m_parents_pr[root] == rank) { // root is a local root
        if(rank == proc_of_interest) cout << "in dbscan line: 314" << " m_parents[root] == root && m_parents_pr[root] == rank: TRUE" << endl;
        // set the root of i directly to root
        m_parents[i] = root;
        m_child_count[root] = m_child_count[root] + 1; // increase the child count by one
        if(rank == proc_of_interest) cout << "in dbscan line: 318" << " m_parents[i]: " << m_parents[i] << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 319" << " m_child_count[root]: " << m_child_count[root] << endl;
        //m_parents_pr[i] = rank; // NO NEED TO SET THIS AS IT        
      } else {
        // set up info to request data from other nodes
        (*p_cur_insert)[m_parents_pr[root]].push_back(0); // flag: 0 means query and 1 means a reply
        (*p_cur_insert)[m_parents_pr[root]].push_back(m_parents[root]);
        (*p_cur_insert)[m_parents_pr[root]].push_back(i);
        (*p_cur_insert)[m_parents_pr[root]].push_back(rank);
        local_continue_to_run++;
        if(rank == proc_of_interest) cout << "in dbscan line: 328" << " m_parents[root]: " << m_parents[root] << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 329" << " i: " << i << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 330" << " rank: " << rank << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 331" << " local_continue_to_run: " << local_continue_to_run << endl;
      }     
    }
    
    // MAY BE REMOVED
    //MPI_Barrier(MPI_COMM_WORLD);
    // TODO 'round' is the name of a function in C++. Might need to change it, just to be on the safe side.
    int pos, round = 0, quadraples, scount, tid, tag = 0, rtag, rsource, rcount, isend[nproc], irecv[nproc], flag;
    // TODO change to malloc() for each of the arrays...
    MPI_Request s_req_recv[nproc], s_req_send[nproc], d_req_send[nproc], d_req_recv[nproc]; // better to malloc the memory
    MPI_Status  s_stat, d_stat_send[nproc], d_stat;
    int target_point, source_point, source_pr;
    if(rank == proc_of_interest) cout << "in dbscan line: 343" << " Start communications..." << endl;
    while(1) {
      global_continue_to_run = 0;
      // Combines values from all processes and distributes the result back to all processes
      MPI_Allreduce(&local_continue_to_run, &global_continue_to_run, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(rank == proc_of_interest) cout << "in dbscan line: 348" << " calling MPI_Allreduce()" << endl;
      if(rank == proc_of_interest) cout << "in dbscan line: 349" << " local_continue_to_run: " << local_continue_to_run << endl;
      if(rank == proc_of_interest) cout << "in dbscan line: 350" << " global_continue_to_run: " << global_continue_to_run << endl;
      if(global_continue_to_run == 0)
        break;
      // bubble sort operation
      pswap = p_cur_insert;
      p_cur_insert = p_cur_send;
      p_cur_send = pswap;
      // wipe out the current vectors of 'p_cur_insert'...
      for(tid = 0; tid < nproc; tid++)
        (*p_cur_insert)[tid].clear();
  
      scount = 0;
      // loop over the messages that need to be sent, and use a non-blocking operation
      for(tid = 0; tid < nproc; tid++) {
        isend[tid] = (*p_cur_send)[tid].size();
        if(rank == proc_of_interest) cout << "in dbscan line: 365" << " isend[tid]: " << isend[tid] << endl;
        if(isend[tid] > 0) {
          MPI_Isend(&(*p_cur_send)[tid][0], isend[tid], MPI_INT, tid, tag + 1, MPI_COMM_WORLD, &d_req_send[scount]);
          scount++;
          if(rank == proc_of_interest) cout << "in dbscan line: 365" << " scount: " << scount << endl;
        }
      }
      // TODO need a MPI_Wait()
      if(rank == proc_of_interest) cout << "in dbscan line: 373" << " isend[0]: " << isend[0] << endl;
      MPI_Alltoall(&isend[0], 1, MPI_INT, &irecv[0], 1, MPI_INT, MPI_COMM_WORLD);

      rcount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(irecv[tid] > 0) {
          merge_received[tid].clear();
          merge_received[tid].assign(irecv[tid], -1);
          MPI_Irecv(&merge_received[tid][0], irecv[tid], MPI_INT, tid, tag + 1, MPI_COMM_WORLD, &d_req_recv[rcount]);
          rcount++;
          if(rank == proc_of_interest) cout << "in dbscan line: 383" << " irecv[tid]: " << irecv[tid] << endl;
          if(rank == proc_of_interest) cout << "in dbscan line: 384" << " rcount: " << rcount << endl;
        }
      }

      local_continue_to_run = 0;
      // loop over the received messages
      for(tid = 0; tid < rcount; tid++) {
        // wait for any replay
        MPI_Waitany(rcount, &d_req_recv[0], &pos, &d_stat);
        if(rank == proc_of_interest) cout << "in dbscan line: 393" << " rcount: " << rcount << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 394" << " d_req_recv[0]: " << d_req_recv[0] << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 395" << " pos: " << pos << endl;
        rtag = d_stat.MPI_TAG;
        rsource = d_stat.MPI_SOURCE;
        if(rank == proc_of_interest) cout << "in dbscan line: 398" << " rtag: " << rtag << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 399" << " rsource: " << rsource << endl;
        // 'tag' is incremented after this loop
        if(rtag == tag + 1) {
          // TODO find out the purpose of dividing by 4
          quadraples = merge_received[rsource].size()/4;
          if(rank == proc_of_interest) cout << "in dbscan line: 404" << " quadraples: " << quadraples << endl;
          for(pid = 0; pid < quadraples; pid++) {
            // TODO - this is dependent on the elements being in the correct order!!!
            // get the quadraple
            // back() returns the last element of the vector
            source_pr = merge_received[rsource].back();
            // pop_back() removes the last element of the vector
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 412" << " source_pr: " << source_pr << endl;
            source_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 415" << " source_point: " << source_point << endl;
            target_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 418" << " target_point: " << target_point << endl;
            flag = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 418" << " flag: " << flag << endl;
            if(flag == 0) {    
              root = target_point;
              if(rank == proc_of_interest) cout << "in dbscan line: 424" << " root: " << root << endl;
              if(rank == proc_of_interest) cout << "in dbscan line: 425" << " m_parents_pr[root]: " << m_parents_pr[root] << endl;
              while(m_parents_pr[root] == rank) {
                if(m_parents[root] == root)
                  break;
                root = m_parents[root];
                if(rank == proc_of_interest) cout << "in dbscan line: 430" << " root: " << root << endl;
              }
              if(rank == proc_of_interest) cout << "in dbscan line: 432" << " m_parents[root] == root && m_parents_pr[root] == rank: " << m_parents[root] == root && m_parents_pr[root] == rank << endl;
              if(m_parents[root] == root && m_parents_pr[root] == rank) { // root is a local root
              
                m_child_count[root] = m_child_count[root] + 1; // increase the child count by one
                // have to return the child about root
                if(rank == proc_of_interest) cout << "in dbscan line: 437" << " m_child_count[root]: " << m_child_count[root] << endl;
                (*p_cur_insert)[source_pr].push_back(1);
                (*p_cur_insert)[source_pr].push_back(source_point);
                (*p_cur_insert)[source_pr].push_back(m_parents[root]);
                (*p_cur_insert)[source_pr].push_back(m_parents_pr[root]);
                local_continue_to_run++;
              } else {
                (*p_cur_insert)[m_parents_pr[root]].push_back(0);
                (*p_cur_insert)[m_parents_pr[root]].push_back(m_parents[root]);
                (*p_cur_insert)[m_parents_pr[root]].push_back(source_point);
                (*p_cur_insert)[m_parents_pr[root]].push_back(source_pr);
                local_continue_to_run++;
              }
              if(rank == proc_of_interest) cout << "in dbscan line: 450" << " local_continue_to_run: " << local_continue_to_run << endl;
            } else {
              // got a reply, so just set the parent
              m_parents[target_point] = source_point;
              m_parents_pr[target_point] = source_pr;
              if(rank == proc_of_interest) cout << "in dbscan line: 455" << " m_parents[target_point]: " << m_parents[target_point] << endl;
              if(rank == proc_of_interest) cout << "in dbscan line: 456" << " m_parents_pr[target_point]: " << m_parents_pr[target_point] << endl;
            }
          }
        }
      }

      tag++;
      if(rank == proc_of_interest) cout << "in dbscan line: 463" << " tag: " << tag << endl;
      round++;
      if(rank == proc_of_interest) cout << "in dbscan line: 465" << " round: " << round << endl;
      if(rank == proc_of_interest) cout << "in dbscan line: 466" << " scount: " << scount << endl;
      if(scount > 0)
        MPI_Waitall(scount, &d_req_send[0], &d_stat_send[0]); // wait for all the sending operation
    }

    // MAY BE REMOVED
    //MPI_Barrier(MPI_COMM_WORLD);

    int final_cluster_root = 0, total_final_cluster_root = 0;

    int points_in_cluster_final = 0, total_points_in_cluster_final = 0;
    // determine the number of points in the current cluster, as well as the total system???
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      if(m_parents[i] == i && m_parents_pr[i] == rank && m_child_count[i] > 1) {
        points_in_cluster_final += m_child_count[i];
        final_cluster_root++;
        if(rank == proc_of_interest) cout << "in dbscan line: 482" << " points_in_cluster_final: " << points_in_cluster_final << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 483" << " final_cluster_root: " << final_cluster_root << endl;
      }
    }
    
    MPI_Allreduce(&points_in_cluster_final, &total_points_in_cluster_final, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&final_cluster_root, &total_final_cluster_root, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int total_points = 0;
    MPI_Allreduce(&m_pts->m_i_num_points, &total_points, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(rank == proc_of_interest) cout << "in dbscan line: 493 Points in clusters " << total_points_in_cluster_final << " Noise " << total_points - total_points_in_cluster_final << " Total points " << total_points << endl;

    if(rank == proc_of_interest) cout << "in dbscan line: 495 Total number of clusters " << total_final_cluster_root << endl;

    vector<int> global_roots;
    global_roots.resize(nproc, 0);
    // Gathers data from all processes
    MPI_Allgather(&final_cluster_root, sizeof(int), MPI_BYTE, &global_roots[0], sizeof(int), MPI_BYTE, MPI_COMM_WORLD); 
  

    int cluster_offset = 0;

    for(i = 0; i <= rank; i++) {
      cluster_offset += global_roots[i];
      if(rank == proc_of_interest) cout << "in dbscan line: 507 cluster_offset: " << cluster_offset << endl;
    }

    m_pid_to_cid.clear();
    m_pid_to_cid.resize(m_pts->m_i_num_points, -1);
    if(rank == proc_of_interest) cout << "in dbscan line: 512 m_pts->m_i_num_points: " << m_pts->m_i_num_points << endl;
    // assign for the global roots only
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      if(m_parents[i] == i && m_parents_pr[i] == rank) {
        if(m_child_count[i] > 1) {
          m_pid_to_cid[i] = cluster_offset;
          cluster_offset++;
          if(rank == proc_of_interest) cout << "in dbscan line: 519 m_pid_to_cid[i]: " << m_pid_to_cid[i] << endl;
          if(rank == proc_of_interest) cout << "in dbscan line: 520 cluster_offset: " << cluster_offset << endl;
        } else {
          m_pid_to_cid[i] = 0; // noise point
          if(rank == proc_of_interest) cout << "in dbscan line: 523 noise point: " << i << endl;
        }
      }
    }
      
    for(i = 0; i < m_pts->m_i_num_points; i++) {
      if(m_parents_pr[i] == rank) {
        if(m_parents[i] != i) { //skip the noise points
          m_pid_to_cid[i] = m_pid_to_cid[m_parents[i]];
          if(rank == proc_of_interest) cout << "in dbscan line: 532 m_pid_to_cid[i]: " << m_pid_to_cid[i] << endl;
        }
      } else {
        // ask the outer to to send back the clusterID
        (*p_cur_insert)[m_parents_pr[i]].push_back(0);
        (*p_cur_insert)[m_parents_pr[i]].push_back(m_parents[i]);
        (*p_cur_insert)[m_parents_pr[i]].push_back(i);
        (*p_cur_insert)[m_parents_pr[i]].push_back(rank);
        if(rank == proc_of_interest) cout << "in dbscan line: 540 m_parents[i]: " << m_parents[i] << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 541 i: " << i << endl;
      }
    }
    
    // MAY BE REMOVED
    //MPI_Barrier(MPI_COMM_WORLD);

    /*
      TODO The following section begins to merge the info from other nodes
    */

    tag++;
    if(rank == proc_of_interest) cout << "in dbscan line: 553 tag: " << tag << endl;
    int later_count;
    for(later_count = 0; later_count < 2; later_count++) {
      pswap = p_cur_insert;
      p_cur_insert = p_cur_send;
      p_cur_send = pswap;

      for(tid = 0; tid < nproc; tid++)
        (*p_cur_insert)[tid].clear();
  
      scount = 0;
      if(rank == proc_of_interest) cout << "in dbscan line: 564 scount: " << scount << endl;
      for(tid = 0; tid < nproc; tid++) {
        isend[tid] = (*p_cur_send)[tid].size();
        if(rank == proc_of_interest) cout << "in dbscan line: 567 isend[tid]: " << isend[tid] << endl;
        if(isend[tid] > 0) {
          MPI_Isend(&(*p_cur_send)[tid][0], isend[tid], MPI_INT, tid, tag + 1, MPI_COMM_WORLD, &d_req_send[scount]);
          scount++;
          if(rank == proc_of_interest) cout << "in dbscan line: 571 scount: " << scount << endl;
        }
      }
      // TODO need a MPI_Wait()
      MPI_Alltoall(&isend[0], 1, MPI_INT, &irecv[0], 1, MPI_INT, MPI_COMM_WORLD);
      if(rank == proc_of_interest) cout << "in dbscan line: 576 isend[0]: " << isend[0] << endl;
      rcount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(irecv[tid] > 0) {
          if(rank == proc_of_interest) cout << "in dbscan line: 580 irecv[tid]: " << irecv[tid] << endl;
          merge_received[tid].clear();
          merge_received[tid].assign(irecv[tid], -1);
          MPI_Irecv(&merge_received[tid][0], irecv[tid], MPI_INT, tid, tag + 1, MPI_COMM_WORLD, &d_req_recv[rcount]);
          rcount++;
          if(rank == proc_of_interest) cout << "in dbscan line: 585 rcount: " << rcount << endl;
        }
      }

      for(tid = 0; tid < rcount; tid++) {
        MPI_Waitany(rcount, &d_req_recv[0], &pos, &d_stat);
        if(rank == proc_of_interest) cout << "in dbscan line: 591 rcount: " << rcount << endl;
        rtag = d_stat.MPI_TAG;
        rsource = d_stat.MPI_SOURCE;
        if(rank == proc_of_interest) cout << "in dbscan line: 594 rtag: " << rtag << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 595 rsource: " << rsource << endl;
        if(rank == proc_of_interest) cout << "in dbscan line: 596 rtag == tag + 1: " << rtag == tag + 1 << endl;
        if(rtag == tag + 1) {
          quadraples = merge_received[rsource].size()/4;
          if(rank == proc_of_interest) cout << "in dbscan line: 599 quadraples: " << quadraples << endl;
          for(pid = 0; pid < quadraples; pid++) {
            // get the quadraple
            source_pr = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 604 source_pr: " << source_pr << endl;
            source_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 607 source_point: " << source_point << endl;
            target_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 610 target_point: " << target_point << endl;
            flag = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(rank == proc_of_interest) cout << "in dbscan line: 613 flag: " << flag << endl;
            if(flag == 0) {         
              (*p_cur_insert)[source_pr].push_back(1);
              (*p_cur_insert)[source_pr].push_back(source_point);
              (*p_cur_insert)[source_pr].push_back(m_pid_to_cid[m_parents[target_point]]);
              (*p_cur_insert)[source_pr].push_back(-1); // One extra INT, may be needed in future
            } else {
              // got a reply, so just set the parent
              m_pid_to_cid[target_point] = source_point; // this assigns the clusterID
            }
          }
        }
      }

      if(scount > 0)
        MPI_Waitall(scount, &d_req_send[0], &d_stat_send[0]); // wait for all the sending operation
    
      //MPI_Barrier(MPI_COMM_WORLD); // MAY NEED TO ACTIVATE THIS
      tag++;
    }

    merge_received.clear();
    merge_send1.clear();
    merge_send2.clear();
    init.clear();
    global_roots.clear();
  }

  // called in mpi_main.cpp.
  // writes data to file
  void ClusteringAlgo::writeCluster_distributed(string outfilename) {
    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if(rank == proc_of_interest) cout << "in dbscan line: 647" << " in ClusteringAlgo::writeCluster_distributed" << endl;
    int i;

    // get the total number of points
    int total_points = 0;
    MPI_Allreduce(&m_pts->m_i_num_points, &total_points, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
    vector<int> point_count;
    point_count.resize(nproc, 0);
    MPI_Allgather(&m_pts->m_i_num_points, sizeof(int), MPI_BYTE, &point_count[0], sizeof(int), MPI_BYTE, MPI_COMM_WORLD);

    int ret, ncfile;
    
    string outfilename_dis = outfilename;
    outfilename_dis = outfilename_dis; //"_clusters.nc";
  
    // create the file, if exists, open the file.
    ret = ncmpi_create(MPI_COMM_WORLD, outfilename_dis.c_str(), NC_CLOBBER, MPI_INFO_NULL, &ncfile);  
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    } 

    MPI_Offset num_particles = total_points;
    int num_particles_id;

    ret = ncmpi_def_dim(ncfile, "num_particles", num_particles, &num_particles_id);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    string column_name_initial = "position_col_X";
    stringstream column_name_id;
    string column_name;

    int j, ncolumn = 1, col_id = 0, varid[m_pts->m_i_dims + 1]; //number of column is 1, col_id is 0 as we use the first one

    // write the columns
    for(j = 0; j < m_pts->m_i_dims; j++) {
      column_name_id.str("");
      column_name_id << j;
  
      column_name = column_name_initial + column_name_id.str(); 
          
      ret = ncmpi_def_var(ncfile, column_name.c_str(), NC_FLOAT, ncolumn, &col_id, &varid[j]);          
      if (ret != NC_NOERR) {
        handle_error(ret, __LINE__);
        return;
      }        
    }

    column_name = "cluster_id";

    ret = ncmpi_def_var(ncfile, column_name.c_str(), NC_INT, ncolumn, &col_id, &varid[j]);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    ret = ncmpi_enddef(ncfile);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    MPI_Offset start[2], count[2];
      
    start[0] = 0;
    for(i = 0; i < rank; i++) {
      start[0] += point_count[i]; 
    }

    count[0] = point_count[rank];
    start[1] = 0; // this to satisfy PnetCDF requirement
    count[1] = 1;//dim_sizes[dimids[1]];
    
    // allocate memory  
    float *data = new float[count[0] * count[1]];

    // write the data columns
    for(j = 0; j < m_pts->m_i_dims; j++) {
      // get the partial column data
      for(i = 0; i < m_pts->m_i_num_points; i++)
        data[i] = m_pts->m_points[i][j];

      // write the data
      ret = ncmpi_put_vara_float_all(ncfile, varid[j], start, count, data);
      if (ret != NC_NOERR) {
        handle_error(ret, __LINE__);
        return;
      }
    }
    delete [] data;

    int *data_id = new int[count[0] * count[1]];    

    //write the cluster_ids
    for(i = 0; i < m_pts->m_i_num_points; i++)
      data_id[i] = m_pid_to_cid[i];

    ret = ncmpi_put_vara_int_all(ncfile, varid[m_pts->m_i_dims], start, count, data_id);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    delete [] data_id;  

    // close the file
    ret = ncmpi_close(ncfile);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }
  
    //cout << "rank " << rank << " AT the end of wrting PnetCDF file" << endl;
  }

  // "uf" == "Union Find"
  // called in mpi_main.cpp
    // Function gets the union of 2 tress, to create a larger cluster
  void run_dbscan_algo_uf_mpi_interleaved(ClusteringAlgo& dbs) {
    double start = MPI_Wtime();     
    int i, pid, j, k, npid, prID;
    int rank, nproc, mpi_namelen;
    kdtree2_result_vector ne;
    kdtree2_result_vector ne_outer;
    if(rank == proc_of_interest) cout << "in dbscan line: 775" << " in run_dbscan_algo_uf_mpi_interleaved" << endl;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // initialize some parameters
    dbs.m_clusters.clear();
      
    // assign parent to itestf
    dbs.m_parents.resize(dbs.m_pts->m_i_num_points, -1);
    dbs.m_parents_pr.resize(dbs.m_pts->m_i_num_points, -1);

    int total_points = 0, points_per_pr[nproc], start_pos[nproc];

    // getting the total number of local points and assigning postions  
    MPI_Allgather(&dbs.m_pts->m_i_num_points, 1, MPI_INT, &points_per_pr[0], 1, MPI_INT, MPI_COMM_WORLD);
    
    for(i = 0; i < nproc; i++) {
      start_pos[i] = total_points;
      total_points += points_per_pr[i];
      if(rank == proc_of_interest) cout << "in dbscan line: 794" << " start_pos[i]: " << start_pos[i] << endl;
      if(rank == proc_of_interest) cout << "in dbscan line: 795" << " total_points: " << total_points << endl;
    }

    // assign proc IDs
    vector <int> vec_prID;
    vec_prID.resize(total_points, -1);

    k = 0;
    for(i = 0; i < nproc; i++) {
      for(j = 0; j < points_per_pr[i]; j++) {
        if(rank == proc_of_interest) cout << "in dbscan line: 805" << " vec_prID[k]: " << vec_prID[k] << endl;
        vec_prID[k++] = i;
        if(rank == proc_of_interest) cout << "in dbscan line: 807" << " vec_prID[k]: " << vec_prID[k] << endl;
      }
    }

    // restting the membership and corepoints values
    dbs.m_member.resize(dbs.m_pts->m_i_num_points, 0);
    dbs.m_corepoint.resize(dbs.m_pts->m_i_num_points, 0);

    vector<int>* ind = dbs.m_kdtree->getIndex();
    vector<int>* ind_outer = dbs.m_kdtree_outer->getIndex();

    // setting paretns to itself and corresponding proc IDs
    for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
      pid = (*ind)[i];
      if(rank == proc_of_interest) cout << "in dbscan line: 821" << " pid: " << pid << endl;
      dbs.m_parents[pid] = pid;
      dbs.m_parents_pr[pid] = rank;
      if(rank == proc_of_interest) cout << "in dbscan line: 821" << " dbs.m_parents[pid]: " << dbs.m_parents[pid] << endl;
      if(rank == proc_of_interest) cout << "in dbscan line: 821" << " dbs.m_parents_pr[pid]: " << dbs.m_parents_pr[pid] << endl;
    }

    vector < vector <int > > merge_received;
    vector < vector <int > > merge_send1;
    vector < vector <int > > merge_send2;
    vector <int> init;
    int rtag, rsource, tag = 0, pos = 0, scount, rcount, isend[nproc], irecv[nproc];
    
    merge_received.resize(nproc, init);
    merge_send1.resize(nproc, init);
    merge_send2.resize(nproc, init);
    
    // reserving communication buffer memory
    for(pid = 0; pid < nproc; pid++) {
      merge_received[pid].reserve(dbs.m_pts->m_i_num_points * nproc);
      merge_send1[pid].reserve(dbs.m_pts->m_i_num_points * nproc);
      merge_send2[pid].reserve(dbs.m_pts->m_i_num_points * nproc);
    }

    int root, root1, root2, tid;

    vector < vector <int > >* pswap;
    vector < vector <int > >* p_cur_send;
    vector < vector <int > >* p_cur_insert;

    p_cur_send = &merge_send1;
    p_cur_insert = &merge_send2;

    if(rank == proc_of_interest) cout << "Init time " << MPI_Wtime() - start << endl; 

    MPI_Barrier(MPI_COMM_WORLD);
    // TODO continue here
    // the main part of the DBSCAN algorithm (called local computation)
    start = MPI_Wtime();
    for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
      pid = (*ind)[i];
      
      // getting the local neighborhoods of local point
      ne.clear();
      dbs.m_kdtree->r_nearest_around_point(pid, 0, dbs.m_epsSquare, ne);
      
      ne_outer.clear();
      vector<float> qv(dbs.m_pts->m_i_dims);

      for (int u = 0; u < dbs.m_pts->m_i_dims; u++)
        qv[u] = dbs.m_kdtree->the_data[pid][u];

      // getting the remote neighborhood of the local point
      if(dbs.m_pts_outer->m_i_num_points > 0)
        dbs.m_kdtree_outer->r_nearest(qv, dbs.m_epsSquare, ne_outer);
    
      qv.clear();
      
      if(ne.size() + ne_outer.size() >= dbs.m_minPts) {
        // pid is a core point
        root = pid;

        dbs.m_corepoint[pid] = 1;
        dbs.m_member[pid] = 1;
        
        // traverse the rmote neighbors and add in the communication buffers  
        for(j = 0; j < ne_outer.size(); j++) {
          npid = ne_outer[j].idx;

          (*p_cur_insert)[dbs.m_pts_outer->m_prIDs[npid]].push_back(pid);
          (*p_cur_insert)[dbs.m_pts_outer->m_prIDs[npid]].push_back(dbs.m_pts_outer->m_ind[npid]);
        }
        
        //traverse the local neighbors and perform union operation
        for (j = 0; j < ne.size(); j++) {
          npid = ne[j].idx;

          // get the root containing npid
          root1 = npid;
          root2 = root;

          if(dbs.m_corepoint[npid] == 1 || dbs.m_member[npid] == 0) {
            dbs.m_member[npid] = 1;

            // REMS algorithm to (union) merge the trees
            while(dbs.m_parents[root1] != dbs.m_parents[root2]) {
              if(dbs.m_parents[root1] < dbs.m_parents[root2]) {
                if(dbs.m_parents[root1] == root1) {
                  dbs.m_parents[root1] = dbs.m_parents[root2];
                  root = dbs.m_parents[root2];
                  break;
                }

                // splicing comression technique
                int z = dbs.m_parents[root1];
                dbs.m_parents[root1] = dbs.m_parents[root2];
                root1 = z;
              } else {
                if(dbs.m_parents[root2] == root2) {
                  dbs.m_parents[root2] = dbs.m_parents[root1];
                  root = dbs.m_parents[root1];
                  break;
                }

                // splicing compressio technique
                int z = dbs.m_parents[root2];
                dbs.m_parents[root2] = dbs.m_parents[root1];                  
                root2 = z;
              }
            }
          }
        }
      }
    }
      
    MPI_Barrier(MPI_COMM_WORLD);

    int v1, v2, par_proc, triples, local_count, global_count;
    double temp_inter_med, inter_med, stop = MPI_Wtime();

    if(rank == proc_of_interest) cout << "Local computation took " << stop - start << endl;

    inter_med = MPI_Wtime();


    start = stop;
    i = 0;
      
    MPI_Request s_req_recv[nproc], s_req_send[nproc], d_req_send[nproc], d_req_recv[nproc]; // better to malloc the memory
    MPI_Status  s_stat, d_stat_send[nproc], d_stat;

    start = MPI_Wtime();

    local_count = 0;
  
    // performing additional compression for the local points that are being sent 
    // this steps identifies the points that actually going to connect the trees in other processors
    // this step will eventually helps further compression before the actual communication happens
    for(tid = 0; tid < nproc; tid++) {
      triples = (*p_cur_insert)[tid].size()/2;
      local_count += triples;

      for(pid = 0; pid < triples; pid++) {
        v1 = (*p_cur_insert)[tid][2 * pid];

        root1 = v1;
        while(dbs.m_parents[root1] != root1)
          root1 = dbs.m_parents[root1];

        while(dbs.m_parents[v1] != root1) {
          int tmp = dbs.m_parents[v1];
          dbs.m_parents[v1] = root1;
          v1 = tmp;
        }

        (*p_cur_insert)[tid][2 * pid] = root1;
      }
    }

    local_count = local_count/nproc;
          
    global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //message_per_round
    int uv, uf, um, ul, ucount;
    int local_continue_to_run, global_continue_to_run;
    double dcomtime = 0, comtime = 0, sum_comp_rate = 0;
    vector <vector <int> > parser;
    vector <int> init_ex;
    parser.resize(dbs.m_pts->m_i_num_points, init_ex);
    
    while(1) {
      pswap = p_cur_insert;
      p_cur_insert = p_cur_send;  
      p_cur_send = pswap;

      for(tid = 0; tid < nproc; tid++)
        (*p_cur_insert)[tid].clear();

      // Uncommend the following if you want to compress in the first round, where compression ratio could be high
      //if(dbs.m_compression == 1 && i == 0)
      //{
      //  dbs.trivial_compression(p_cur_send, &parser, nproc, rank, i, comtime, sum_comp_rate);
      //}
      
      // send all the data
      // compress the data before send
      //dbs.convert(p_cur_send, nproc, rank, i);
    
      scount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(dbs.m_compression == 1 && i == 0 && (*p_cur_send)[tid].size() > 0) {
          dbs.trivial_compression(&(*p_cur_send)[tid], &parser, nproc, rank, i, comtime, sum_comp_rate);
        }

        isend[tid] = (*p_cur_send)[tid].size();
        if(isend[tid] > 0) {
          MPI_Isend(&(*p_cur_send)[tid][0], isend[tid], MPI_INT, tid, tag + 1, MPI_COMM_WORLD, &d_req_send[scount]);
          scount++;
        }
      }
      // TODO need a MPI_Wait()
      MPI_Alltoall(&isend[0], 1, MPI_INT, &irecv[0], 1, MPI_INT, MPI_COMM_WORLD);

      rcount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(irecv[tid] > 0) {
          merge_received[tid].clear();
          merge_received[tid].assign(irecv[tid], -1);
          MPI_Irecv(&merge_received[tid][0], irecv[tid], MPI_INT, tid, tag + 1, MPI_COMM_WORLD, &d_req_recv[rcount]);
          rcount++;
        }
      }

      local_count = 0;
      
      //get the data and process them
      for(tid = 0; tid < rcount; tid++) {
        MPI_Waitany(rcount, &d_req_recv[0], &pos, &d_stat);
      
        rtag = d_stat.MPI_TAG;
        rsource = d_stat.MPI_SOURCE;
  
        if(rtag == tag + 1) {
          // process received the data now
          if(dbs.m_messages_per_round == -1 && i == 0) {
            if(dbs.m_compression == 1) {
              // call the decompression function
              dbs.trivial_decompression(&merge_received[rsource], nproc, rank, i, dcomtime);

              triples = merge_received[rsource].size()/2;
              par_proc = rsource;
            } else {
              triples = merge_received[rsource].size()/2;
              par_proc = rsource;
            }
          } else {
            triples = merge_received[rsource].size()/3;
          }

          for(pid = 0; pid < triples; pid++) {
            // get the pair
            v1 = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            
            if((dbs.m_messages_per_round == -1 && i > 0) || (dbs.m_messages_per_round != -1)) {
              par_proc = merge_received[rsource].back();
              merge_received[rsource].pop_back();
            }
        
            v2 = merge_received[rsource].back();
            merge_received[rsource].pop_back();
        
            int con = 0;
            
            if(i > 0)
              con = 1;
            else if (i == 0 && (dbs.m_corepoint[v1] == 1 || dbs.m_member[v1] == 0)) { 
              dbs.m_member[v1] = 1;
              con = 1;
            }

            if(con == 1) {
        
              root1 = v1;
              // this will find the boundary vertex or the root if the root is in this processor
              while(dbs.m_parents_pr[root1] == rank) {
                if(dbs.m_parents[root1] == root1)
                  break;

                root1 = dbs.m_parents[root1];
              }
        
              // compress the local path
              while(v1 != root1 && vec_prID[v1] == rank) {
                int tmp = dbs.m_parents[v1];
                dbs.m_parents[v1] = root1;
                v1 = tmp;
              }
  
              if(dbs.m_parents[root1] == v2 && dbs.m_parents_pr[root1] == par_proc) {
                //same_set++;
                continue;
              }           
                
              if(par_proc == rank) {
                if(dbs.m_parents[root1] == dbs.m_parents[v2])
                  continue;
              }
                
              if(dbs.m_parents[root1] == root1 && dbs.m_parents_pr[root1] == rank) { // root1 is a local root
                if(start_pos[rank] + root1 < start_pos[par_proc] + v2) {
                  // do union
                  dbs.m_parents[root1] = v2;
                  dbs.m_parents_pr[root1] = par_proc;
                  continue;
                } else {
                  // ask the parent of v2
                  (*p_cur_insert)[par_proc].push_back(root1);
                  (*p_cur_insert)[par_proc].push_back(dbs.m_parents_pr[root1]);
                  (*p_cur_insert)[par_proc].push_back(v2);

                  local_count++;
                }
              } else {
                // root1 is not local
                if(start_pos[dbs.m_parents_pr[root1]] + root1 < start_pos[par_proc] + v2) {
                  // ask the parent of root1
                  (*p_cur_insert)[dbs.m_parents_pr[root1]].push_back(v2);
                  (*p_cur_insert)[dbs.m_parents_pr[root1]].push_back(par_proc);
                  (*p_cur_insert)[dbs.m_parents_pr[root1]].push_back(dbs.m_parents[root1]);
                    
                  local_count++;
                } else {
                  // ask the parent of v2
                  (*p_cur_insert)[par_proc].push_back(dbs.m_parents[root1]);
                  (*p_cur_insert)[par_proc].push_back(dbs.m_parents_pr[root1]);
                  (*p_cur_insert)[par_proc].push_back(v2);

                  local_count++;
                }
              }
            }
          }
          
          merge_received[rsource].clear();
        } else {
          cout << "rank " << rank << " SOMETHING IS WRONG" << endl;
        }
      }
        
      if(scount > 0)
        MPI_Waitall(scount, &d_req_send[0], &d_stat_send[0]); 

      tag += 2; // change the tag value although not important
      
      local_continue_to_run = 0;
      local_count = 0;
      for(tid = 0; tid < nproc; tid++) {
        local_count += (*p_cur_insert)[tid].size()/3; // TODO why divide by 3
        if((*p_cur_insert)[tid].size() > 0)
          local_continue_to_run = 1;
      }
      
      local_count = local_count / nproc;

      global_count = 0;
      global_continue_to_run = 0;
            
      MPI_Allreduce(&local_continue_to_run, &global_continue_to_run, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      if(global_continue_to_run == 0)
        break;
      i++;
    }

    stop = MPI_Wtime(); 
    if(rank == proc_of_interest) cout << "Merging took " << stop - start << endl;

    pswap = NULL;
    p_cur_insert = NULL;
    p_cur_send = NULL;

    for(tid = 0; tid < nproc; tid++) {
      merge_received[tid].clear();
      merge_send1[tid].clear();
      merge_send2[tid].clear();
    }

    merge_received.clear();
    merge_send1.clear();
    merge_send2.clear();
    ind = NULL;
    ind_outer = NULL;

    vec_prID.clear();
    ne.clear();
    ne_outer.clear();
    parser.clear();
    init_ex.clear();
    init.clear();
  }
};


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

#include "geometric_partitioning.h"

namespace NWUClustering {
  void get_extra_points(ClusteringAlgo& dbs) {
    #ifdef _DEBUG_GP
    MPI_Barrier(MPI_COMM_WORLD);
    double end, start = MPI_Wtime();
    #endif

    int k, rank, nproc, i, j;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    #ifdef _DEBUG
    if(rank == proc_of_interest) cout << "extra point time part 0 strating " << endl;
    #endif
    
    interval* local_box = new interval[dbs.m_pts->m_i_dims];
    compute_local_bounding_box(dbs, local_box);

    // extend the box
    float eps = sqrt(dbs.m_epsSquare);

    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      local_box[i].upper += eps;
      local_box[i].lower -= eps; 
    }
  
    // all together all the extending bounding box
    interval* gather_local_box = new interval[dbs.m_pts->m_i_dims * nproc];

    // gather the local bounding box first
    MPI_Allgather(local_box, sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, gather_local_box,
    sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, MPI_COMM_WORLD);

    bool if_inside, overlap;
    int count = 0, gcount;

    vector <float> empty;
    vector <vector <float> > send_buf;
    vector <vector <float> > recv_buf;
    send_buf.resize(nproc, empty);
    recv_buf.resize(nproc, empty);

    vector <int> empty_i;
    vector <vector <int> > send_buf_ind;
    vector <vector <int> > recv_buf_ind;
    send_buf_ind.resize(nproc, empty_i);
    recv_buf_ind.resize(nproc, empty_i);

    for(k = 0; k < nproc; k++) {
      if(k == rank) // self
        continue;

      // check the two extended bounding box of proc rank and k. If they don't overlap, there must be no points 
      // SHOULD SUBTRACT EPS      

      overlap = true;
      for(j = 0; j < dbs.m_pts->m_i_dims; j++) {
        if(gather_local_box[rank * dbs.m_pts->m_i_dims + j].lower < gather_local_box[k * dbs.m_pts->m_i_dims +j].lower) {
          //if(gather_local_box[rank * dbs.m_pts->m_i_dims + j].upper < gather_local_box[k * dbs.m_pts->m_i_dims + j].lower)
          if(gather_local_box[rank * dbs.m_pts->m_i_dims + j].upper - gather_local_box[k * dbs.m_pts->m_i_dims + j].lower < eps) {
            overlap = false;
            break;
          }
        } else {
          //if(gather_local_box[k * dbs.m_pts->m_i_dims + j].upper < gather_local_box[rank * dbs.m_pts->m_i_dims + j].lower)
          if(gather_local_box[k * dbs.m_pts->m_i_dims + j].upper - gather_local_box[rank * dbs.m_pts->m_i_dims + j].lower < eps) {
            overlap = false;
            break;
          }
        }
      }

      // the two bouding boxes are different, so continue to the next processors
      if(overlap == false)
        continue;

      // get the overlapping regions
      for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
        if_inside = true;
        for(j = 0; j < dbs.m_pts->m_i_dims; j++) {
          if(dbs.m_pts->m_points[i][j] < gather_local_box[k * dbs.m_pts->m_i_dims + j].lower || 
            dbs.m_pts->m_points[i][j] > gather_local_box[k * dbs.m_pts->m_i_dims + j].upper) {
            if_inside = false;
            break;
          }
        }
        
        if(if_inside == true) {
          for(j = 0; j < dbs.m_pts->m_i_dims; j++)
            send_buf[k].push_back(dbs.m_pts->m_points[i][j]);
          
          send_buf_ind[k].push_back(i);
          count++;
        }
      }
    }

    #ifdef _DEBUG_GP
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
      #ifdef _DEBUG
      if(rank == proc_of_interest) cout << "extra point time part 2: " << end - start << endl;        
      #endif
    start = end;
    #endif

    // now send buf have all the points. Send the size first to everyone
    vector <int> send_buf_size, recv_buf_size;
    send_buf_size.resize(nproc, 0);
    recv_buf_size.resize(nproc, 0);

    for(i = 0; i < nproc; i++)
      send_buf_size[i] = send_buf[i].size();
  
    MPI_Alltoall(&send_buf_size[0], 1, MPI_INT, &recv_buf_size[0], 1, MPI_INT, MPI_COMM_WORLD);

    //return;
    int tag = 200, send_count, recv_count;
    MPI_Request req_send[2 * nproc], req_recv[2 * nproc];
    MPI_Status stat_send[2 * nproc], stat_recv;

    recv_count = 0;
    for(i = 0; i < nproc; i++) {
      if(recv_buf_size[i] > 0) {
        recv_buf[i].resize(recv_buf_size[i], 0);
        recv_buf_ind[i].resize(recv_buf_size[i] / dbs.m_pts->m_i_dims, -1);

        MPI_Irecv(&recv_buf[i][0], recv_buf_size[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &req_recv[recv_count++]);
        MPI_Irecv(&recv_buf_ind[i][0], recv_buf_size[i] / dbs.m_pts->m_i_dims, MPI_INT, i, tag + 1, MPI_COMM_WORLD, &req_recv[recv_count++]);
      }
    }

    send_count = 0;
    for(i = 0; i < nproc; i++) {
      if(send_buf_size[i] > 0) {
        MPI_Isend(&send_buf[i][0], send_buf_size[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &req_send[send_count++]);
        MPI_Isend(&send_buf_ind[i][0], send_buf_size[i] / dbs.m_pts->m_i_dims, MPI_INT, i, tag + 1, MPI_COMM_WORLD, &req_send[send_count++]);
      }
    }

    int rtag, rsource, rpos;

    #if _GET_PARTION_STAT == 0
    dbs.allocate_outer(dbs.m_pts->m_i_dims);
    #endif

    for(i = 0; i < recv_count; i++) {
      MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
    
      rtag = stat_recv.MPI_TAG;
      rsource = stat_recv.MPI_SOURCE;
    
      if(rtag == tag) {
        // process the request
        //cout << "proc " << rank << " add points called " << endl;
        #if _GET_EXTRA_POINT_STAT == 0  // WHY THIS IS HERE??????????????????????????????????????????????????????????????
        dbs.addPoints(rsource, recv_buf_size[rsource], dbs.m_pts->m_i_dims, recv_buf[rsource]);
        #endif
        recv_buf[rsource].clear();
      } else if(rtag == tag + 1) {
        // postpond this computation and call update points later
        // processing immediately might lead to invalid computation
      }
    } 
    
    // MAY NOT NEED THIS
    if(send_count > 0)
      MPI_Waitall(send_count, &req_send[0], &stat_send[0]);

    #ifdef _DEBUG_GP
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
      #ifdef _DEBUG
      if(rank == proc_of_interest) cout << "extra point time part 3: " << end - start << endl;        
      #endif
    start = end;
    #endif

    // got all the points
    // now update the indices of the outer points

    #if _GET_EXTRA_POINT_STAT == 0 // WHY THIS IS HERE??????????????????????????????????????????????????????????????
    dbs.updatePoints(recv_buf_ind);
    #endif

    MPI_Reduce(&count, &gcount, 1, MPI_INT, MPI_SUM, proc_of_interest, MPI_COMM_WORLD);

    empty.clear();
    send_buf.clear();
    recv_buf.clear();
    send_buf_size.clear();
    recv_buf_size.clear();
    send_buf_ind.clear();
    recv_buf_ind.clear();

    delete [] gather_local_box;
    delete [] local_box;
  }

  void start_partitioning(ClusteringAlgo& dbs) {
    int r_count, s_count, rank, nproc, i, j, k;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // compute the local bouding box for each dimention
    interval* box = new interval[dbs.m_pts->m_i_dims];
    
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) { 
      box[i].upper = dbs.m_pts->m_box[i].upper;
      box[i].lower = dbs.m_pts->m_box[i].lower;
    }

    // compute the global bouding box for each dimention
    interval* gbox = new interval[dbs.m_pts->m_i_dims];
    compute_global_bounding_box(dbs, box, gbox, nproc);

    #ifdef _DEBUG
    if(rank == proc_of_interest) cout << "Partitioning: Pos 1" << endl;
    #endif

    // find the loop count for nproc processors
    int internal_nodes, partner_rank, loops, b, color, sub_rank, d, max, sub_nprocs;

    MPI_Status status;

    loops = 0;
    i = nproc;
    internal_nodes = 1;
    while((i = i >> 1) > 0) {
      loops++;
      internal_nodes = internal_nodes << 1;
    }
    
    internal_nodes = internal_nodes << 1;
    
    //gbox for each node in the tree [ONLY upto to reaching each processor]
    interval** nodes_gbox = new interval*[internal_nodes];
    for(i = 0; i < internal_nodes; i++)
      nodes_gbox[i] = new interval[dbs.m_pts->m_i_dims];
    
    copy_global_box_to_each_node(dbs, nodes_gbox, gbox, internal_nodes);
  
    vector <float> send_buf;
    vector <int>   invalid_pos_as;
    vector <float> recv_buf;
        
    int pow2_i;
    float median;

    for(i = 0; i < loops; i++) {
      pow2_i = POW2(i);
      b  = nproc - (int) (nproc / pow2_i);
      color = (int)((rank & b) / POW2(loops - i ));
      partner_rank = rank ^ (int)(nproc/POW2(i + 1));

      MPI_Comm new_comm;
      MPI_Comm_split(MPI_COMM_WORLD, color, rank, &new_comm);
      MPI_Comm_rank(new_comm, &sub_rank);

      if(sub_rank == 0) {
        d = 0;
        for(j = 1; j < dbs.m_pts->m_i_dims; j++) {
          if(nodes_gbox[pow2_i + color][j].upper - nodes_gbox[pow2_i + color][j].lower > 
              nodes_gbox[pow2_i + color][d].upper - nodes_gbox[pow2_i + color][d].lower)
            d = j;
        }
      } 
  
      MPI_Bcast(&d, 1, MPI_INT, 0, new_comm);

      // compute the median in this dimension
      float median  = get_median(dbs, d, new_comm);   

      s_count = get_points_to_send(dbs, send_buf, invalid_pos_as, median, d, rank, partner_rank);

      if (rank < partner_rank) {
        MPI_Sendrecv(&s_count, 1, MPI_INT, partner_rank, 4, &r_count, 1, MPI_INT, partner_rank, 5, MPI_COMM_WORLD, &status);
        recv_buf.resize(r_count * dbs.m_pts->m_i_dims, 0.0);
        MPI_Sendrecv(&send_buf[0], s_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 2,
              &recv_buf[0], r_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD, &status);
        send_buf.clear();
      } else {
        MPI_Sendrecv(&s_count, 1, MPI_INT, partner_rank, 5, &r_count, 1, MPI_INT, partner_rank, 4, MPI_COMM_WORLD, &status);
        recv_buf.resize(r_count * dbs.m_pts->m_i_dims, 0.0);
        MPI_Sendrecv(&send_buf[0], s_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 3, 
              &recv_buf[0], r_count * dbs.m_pts->m_i_dims, MPI_FLOAT, partner_rank, 2, MPI_COMM_WORLD, &status);
        send_buf.clear();
      }

      update_points(dbs, s_count, invalid_pos_as, recv_buf);
      recv_buf.clear();
  
      copy_box(dbs, nodes_gbox[LOWER(pow2_i+color)], nodes_gbox[pow2_i+color]);
      nodes_gbox[LOWER(pow2_i+color)][d].upper =  median;
      copy_box(dbs, nodes_gbox[UPPER(pow2_i+color)], nodes_gbox[pow2_i+color]);
      nodes_gbox[UPPER(pow2_i+color)][d].lower =  median; 

      MPI_Comm_free(&new_comm);
    }

    // free the allocated memory
    for(i = 0; i < nproc; i++)
      delete [] nodes_gbox[i];

    delete [] nodes_gbox;
    delete [] gbox;
    delete [] box;

  }

  void update_points(ClusteringAlgo& dbs, int s_count, vector <int>& invalid_pos_as, vector <float>& recv_buf) {
    int i, j, k, l, r_count = recv_buf.size() / dbs.m_pts->m_i_dims;

    if(r_count >= s_count) {
      //invalid_pos_as.reserve(dbs.m_pts->m_i_num_points + r_count - s_count);
      invalid_pos_as.resize(dbs.m_pts->m_i_num_points + r_count - s_count, 1);

      //allocate memory for the points
      dbs.m_pts->m_points.resize(dbs.m_pts->m_i_num_points + r_count - s_count);
      for(int ll = 0; ll < dbs.m_pts->m_i_num_points + r_count - s_count; ll++)
        dbs.m_pts->m_points[ll].resize(dbs.m_pts->m_i_dims);

      j = 0;
      for(i = 0; i < invalid_pos_as.size(); i++) {
        if(invalid_pos_as[i] == 1) {
          for(k = 0; k < dbs.m_pts->m_i_dims; k++)
            dbs.m_pts->m_points[i][k] = recv_buf[j++];
        }
      }     

      dbs.m_pts->m_i_num_points = dbs.m_pts->m_i_num_points + r_count - s_count;
    } else {
      j = 0;
      i = 0;  
      if(recv_buf.size() > 0) {
        for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
          if(invalid_pos_as[i] == 1) {
            for(k = 0; k < dbs.m_pts->m_i_dims; k++)
              dbs.m_pts->m_points[i][k] = recv_buf[j++];
          
            if(j == recv_buf.size()) {
              i++;
              break;
            }
          }
        }
      }
      
      l = dbs.m_pts->m_i_num_points;
      for( ; i < invalid_pos_as.size(); i++) {
        if(invalid_pos_as[i] == 1) {
          while(l > i) {
            l--;
            if(invalid_pos_as[l] == 0)
              break;
          }

          if(invalid_pos_as[l] == 0)  
            for(k = 0; k < dbs.m_pts->m_i_dims; k++)
              dbs.m_pts->m_points[i][k] = dbs.m_pts->m_points[l][k];
        }
      }

      //allocate memory for the points
      dbs.m_pts->m_points.resize(dbs.m_pts->m_i_num_points + r_count - s_count);
      for(int ll = 0; ll < dbs.m_pts->m_i_num_points + r_count - s_count; ll++)
        dbs.m_pts->m_points[ll].resize(dbs.m_pts->m_i_dims);

      dbs.m_pts->m_i_num_points = dbs.m_pts->m_i_num_points + r_count - s_count;
    }   
  }

  int get_points_to_send(ClusteringAlgo& dbs, vector <float>& send_buf, vector <int>& invalid_pos_as, float median, int d, int rank, int partner_rank) {
    int i, count = 0, j;
    send_buf.reserve(dbs.m_pts->m_i_num_points * dbs.m_pts->m_i_dims);
    invalid_pos_as.clear();
    invalid_pos_as.resize(dbs.m_pts->m_i_num_points, 0);

    for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
      if (rank < partner_rank) {
        if(dbs.m_pts->m_points[i][d] > median) {
          invalid_pos_as[i] = 1;
          count++;
          for(j = 0; j < dbs.m_pts->m_i_dims; j++)
            send_buf.push_back(dbs.m_pts->m_points[i][j]);
        }
      } else {
        if(dbs.m_pts->m_points[i][d] <= median) {
          invalid_pos_as[i] = 1;
          count++;
          for(j = 0; j < dbs.m_pts->m_i_dims; j++)
            send_buf.push_back(dbs.m_pts->m_points[i][j]);
        }
      }
    }

    return count;
  } 

  float get_median(ClusteringAlgo& dbs, int d, MPI_Comm& new_comm) { 

    // ADDITIONAL CODE
    float median;
    
    vector <float> data;
    data.reserve(dbs.m_pts->m_i_num_points);
    data.resize(dbs.m_pts->m_i_num_points, 0);

    for (int k=0; k < dbs.m_pts->m_i_num_points; k++)
      data[k] = dbs.m_pts->m_points[k][d];

    median = findKMedian(data, data.size()/2);
    data.clear();

    int proc_count;
    MPI_Comm_size(new_comm, &proc_count);

    vector <float> all_medians;
    all_medians.resize(proc_count, 0);

    MPI_Allgather(&median, sizeof(int), MPI_BYTE, &all_medians[0], sizeof(int), MPI_BYTE, new_comm);  

    median = findKMedian(all_medians, all_medians.size()/2); 
    all_medians.clear();
    
    return median;  
  }

  void compute_local_bounding_box(ClusteringAlgo& dbs, interval* box) {
    int i, j;

    //we assume each processor has at least one point
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      box[i].upper = dbs.m_pts->m_points[0][i];
      box[i].lower = dbs.m_pts->m_points[0][i];
    }
  
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      for(j = 1; j < dbs.m_pts->m_i_num_points; j++) {
        if(box[i].lower > dbs.m_pts->m_points[j][i])
          box[i].lower = dbs.m_pts->m_points[j][i];
        else if(box[i].upper < dbs.m_pts->m_points[j][i])
          box[i].upper = dbs.m_pts->m_points[j][i];
      }
    }
  }

  void compute_global_bounding_box(ClusteringAlgo& dbs, interval* box, interval* gbox, int nproc) {
    int i, j, k;
  
    interval* gather_local_box = new interval[dbs.m_pts->m_i_dims * nproc];
  
    // gather the local bounding box first
    MPI_Allgather(box, sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, gather_local_box, 
          sizeof(interval) * dbs.m_pts->m_i_dims, MPI_BYTE, MPI_COMM_WORLD);

    // compute the global bounding box
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      gbox[i].lower = gather_local_box[i].lower;
      gbox[i].upper = gather_local_box[i].upper;
      
      k = i;
      for(j = 0; j < nproc; j++, k += dbs.m_pts->m_i_dims) {
        if(gbox[i].lower > gather_local_box[k].lower)
          gbox[i].lower = gather_local_box[k].lower;
        
        if(gbox[i].upper < gather_local_box[k].upper)
          gbox[i].upper = gather_local_box[k].upper;
      }
    }
    
    delete [] gather_local_box;
  }

  void copy_global_box_to_each_node(ClusteringAlgo& dbs, interval** nodes_gbox, interval* gbox, int internal_nodes) {
    int i, j;
    for(i = 0; i < internal_nodes; i++) {
      for(j = 0; j < dbs.m_pts->m_i_dims; j++) {
        nodes_gbox[i][j].upper = gbox[j].upper;
        nodes_gbox[i][j].lower = gbox[j].lower;
      }
    }
  }
  
  void copy_box(ClusteringAlgo& dbs, interval* target_box, interval* source_box) {
    for(int j = 0; j < dbs.m_pts->m_i_dims; j++) {
      target_box[j].upper = source_box[j].upper;
      target_box[j].lower = source_box[j].lower;
    }
  }
};


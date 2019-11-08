
#include "clusters.h"

namespace NWUClustering {

  // destructor 
    // can have 4 parameters: m_pts,m_kdtree,m_pts_outer,m_kdtree_outer
    // Clears, and deletes Points, Points_Outer, & kdtree2 objects
  Clusters::~Clusters() { // Pure virtual destructor 
    if(m_pts) {
      m_pts->m_points.clear();
      delete m_pts;
      m_pts = NULL;
    }

    if(m_kdtree) {
      delete m_kdtree;
      m_kdtree = NULL;
    }

    if(m_pts_outer) {
      m_pts_outer->m_prIDs.clear();
      m_pts_outer->m_ind.clear();
      m_pts_outer->m_points.clear();
      delete m_pts_outer;
      m_pts_outer = NULL;
    }

    if(m_kdtree_outer) {
      delete m_kdtree_outer;
      m_kdtree_outer = NULL;
    }
  }

  // Sets up the Points_Outer objects
    // Called in geometric_partitioning.cpp
  bool Clusters::allocate_outer(int dims) {
    // int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if(rank == proc_of_interest) cout << "in Clusters line: 38" << " in Clusters::allocate_outer" << endl;
    if(m_pts_outer == NULL) {
      m_pts_outer = new Points_Outer;
      m_pts_outer->m_prIDs.clear();
      m_pts_outer->m_ind.clear();
      m_pts_outer->m_i_dims = dims;
      m_pts_outer->m_i_num_points = 0;
    }
  }

  // Adds points to a cluster object
    // Called in geometric_partitioning.cpp
  bool Clusters::addPoints(int source, int buf_size, int dims, vector<float>& raw_data) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if(rank == proc_of_interest) cout << "in Clusters line: 53" << " in Clusters::updatePoints" << endl;
    // if(rank == proc_of_interest) cout << "in Clusters line: 54" << " source: " << source << endl;
    // if(rank == proc_of_interest) cout << "in Clusters line: 55" << " buf_size: " << buf_size << endl;
    // if(rank == proc_of_interest) cout << "in Clusters line: 56" << " dims: " << dims << endl;
    // used as itterables
    int i, j, k; 
    // pos = original num of points, then used as a key in 'm_pts_outer->m_points'
    int pos;
    // num_points = new number of points for the cluster
    int num_points = buf_size / dims;
    // TODO This is for tracing
    //cout << "add points called" << endl; 
    // if(rank == proc_of_interest) cout << "in Clusters line: 65" << " num_points: " << num_points << endl;
    // incorrect dimension
      // the number of dimensions/attributes can never change
    if(m_pts_outer->m_i_dims != dims)
      return false;

    // Save the original number of points, then increase by the new "num_points" value
    pos = m_pts_outer->m_i_num_points;
    m_pts_outer->m_i_num_points += num_points;
    //m_pts_outer->m_points.resize(extents[m_pts_outer->m_i_num_points][dims]);
    // if(rank == proc_of_interest) cout << "in Clusters line: 75" << " m_pts_outer->m_i_num_points: " << m_pts_outer->m_i_num_points << endl;
    //allocate memory for the points
      // resize() - Resizes the container so that it contains n elements
    m_pts_outer->m_points.resize(m_pts_outer->m_i_num_points);
    // resize each new point object
    for(int ll = 0; ll < m_pts_outer->m_i_num_points; ll++)
      m_pts_outer->m_points[ll].resize(dims);
    // resize the point IDs
    m_pts_outer->m_prIDs.resize(m_pts_outer->m_i_num_points, -1);

    k = 0;
    // loop over the new number of points
    for(i = 0; i < num_points; i++) {
      // loop over the number of dimensions/attributes
      for(j = 0; j < dims; j++) {
        // assign the new 'raw_data' elements to the newest points
        m_pts_outer->m_points[pos][j] = raw_data[k++]; // TODO possible buffer overflow????????
        // if(rank == proc_of_interest) cout << "in Clusters line: 92" << " m_pts_outer->m_points[pos][j]: " << m_pts_outer->m_points[pos][j] << endl;
      }
      // assign the cluster ID for the point
      m_pts_outer->m_prIDs[pos] = source; // TODO check that source is a cluster ID
      // increment the counter for the key...
      pos++;
    }

    //cout << "outer " << m_pts_outer->m_i_num_points << endl;

    return true;
  }

  // Updates OUTER points' cluster IDs
    // Called in geometric_partitioning.cpp
  bool Clusters::updatePoints(vector< vector<int> >& raw_ind) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if(rank == proc_of_interest) cout << "in Clusters line: 110" << " in Clusters::updatePoints" << endl;
    // used as itterables
    int i, j = -1;
    // I believe these are cluster IDs
    int source = -1,  prev_source = -1;
    // resize 'm_ind', by 'm_i_num_points' spaces, new elements are initialized as copies of -1
    m_pts_outer->m_ind.resize(m_pts_outer->m_i_num_points, -1);
    // loop over the Outer points, and update cluster IDs
    for(i = 0; i < m_pts_outer->m_i_num_points; i++) {
      source = m_pts_outer->m_prIDs[i];
      // if(rank == proc_of_interest) cout << "in Clusters line: 120" << " source: " << source << endl;
      if(source != prev_source)
        j = 0;

      m_pts_outer->m_ind[i] = raw_ind[source][j++];
      // if(rank == proc_of_interest) cout << "in Clusters line: 125" << " m_pts_outer->m_ind[i]: " << m_pts_outer->m_ind[i] << endl;
      prev_source = source;
    }

    return true;
  }

  /*
    Called from mpi_main.cpp
    Reads the binary file. 
    Input data points are equally partioned, each core reading their corresponding part of the file. 
    Later the points will be partioned geometrically.
    Points are created here, based off of the data points.
    According to the README, 'num_points' & 'dims' HAVE to be the 1st 2 things in the file (each 4 bytes).
  */
  int Clusters::read_file(char* infilename, int isBinaryFile) {
    ssize_t numBytesRead; // TODO this isn't even used
    // used as itterables
    int i, j;
    // rank = current node's ID, nproc = total number of nodes in the system
    int rank, nproc;
    int num_points, dims;
    // Get the current node's 'rank'
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the total number of nodes in the system
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    // if(rank == proc_of_interest) cout << "in Clusters line: 151" << " in Clusters::read_file" << endl;
    // only supports binary data file
    if(isBinaryFile == 1) {
      // NOTE: the input data points are equally partioned and each core read their corresponding part
      // later the points will be partioned geometrically

      ifstream file (infilename, ios::in|ios::binary);
      if(file.is_open()) {
        // size of the dataset
        file.read((char*)&num_points, sizeof(int));
        // number of dimensions for each data point/record
        file.read((char*)&dims, sizeof(int));

        // #ifdef _DEBUG 
        // if(rank == proc_of_interest) cout << "in Clusters line: 165" << "Points: " << num_points << " dims: " << dims << endl;
        // #endif

        // compute the respective segments of the file to read.
          // sch = section or search?
          // lower = the starting position in the file to start reading from
          // upper = the last position in the file to read from
        long long sch, lower, upper;
        // determine the size of the section of file to read
        if(num_points % nproc == 0)
          sch = num_points/nproc; // the 2 are evenly divisible
        else
          sch = num_points/nproc + 1; // taking care of the remainder, by adding a whole 1

        // #ifdef _DEBUG
        if(rank == proc_of_interest) {
          // cout << "in Clusters line: 181 Segment size: " << sch << endl;
          // cout << "in Clusters line: 182 Points per process on average: " << num_points/nproc << endl;
        }
        // #endif

        lower = sch * rank; // the pointer for the starting location
        upper = sch * (rank + 1); // the pointer for the end location
        // if the "end location" exceeds the number of points, reduce it to the number of points
        if(upper > num_points)
          upper = num_points;
        
        if(rank == proc_of_interest) {
          // cout << "in Clusters line: 193" << " upper: " << upper << endl;
          // cout << "in Clusters line: 194" << " lower: " << lower << endl;
        }
        // allocate memory for points
        m_pts = new Points;
        m_pts->m_i_dims = dims; // number of dimensions/attributes
        m_pts->m_i_num_points = upper - lower; // the number of points

        // interval struct defined in kdtree2.hpp
          // m_box, is of type interval*
        m_pts->m_box = new interval[m_pts->m_i_dims]; 

        //allocate memory for the points
        //m_pts->m_points.resize(extents[num_points][dims]);
        //m_pts->m_points.resize(extents[m_pts->m_i_num_points][dims]);
        // if(rank == proc_of_interest) cout << "in Clusters line: 208 m_i_dims: " << m_pts->m_i_dims << endl;
        // if(rank == proc_of_interest) cout << "in Clusters line: 209 m_i_num_points: " << m_pts->m_i_num_points << endl;

        //allocate memory for the points
        m_pts->m_points.resize(m_pts->m_i_num_points);
        for(int ll = 0; ll < m_pts->m_i_num_points; ll++)
          m_pts->m_points[ll].resize(dims);

        // point_coord_type = typedef float point_coord_type; in 'utils.h'
        point_coord_type* pt;         
        // initializes 'pt' variable
        pt = (point_coord_type*) malloc(dims * sizeof(point_coord_type));

        // fseek to the respective position of the file
        file.seekg(lower * dims * sizeof(point_coord_type), ios::cur);
        // loop over the area of the file for the node
          // TODO possible speed up: put 'upper - lower' into a local variable...
        for (i = 0; i < upper - lower; i++) {
          // signature: istream& read (char* s, streamsize n);
          // Extracts n characters from the stream and stores them in the array pointed to by s.
          file.read((char*)pt, dims * sizeof(point_coord_type));
          
          for (j = 0; j < dims; j++) {
            m_pts->m_points[i][j] = pt[j];
            // if(rank == proc_of_interest) cout << "in Clusters line: 232 m_pts->m_points[i][j]: " << m_pts->m_points[i][j] << endl;
            if(i == 0) { 
              m_pts->m_box[j].upper = m_pts->m_points[i][j];
              m_pts->m_box[j].lower = m_pts->m_points[i][j];
            } else {
              if(m_pts->m_box[j].lower > m_pts->m_points[i][j])
                m_pts->m_box[j].lower = m_pts->m_points[i][j];
              else if(m_pts->m_box[j].upper < m_pts->m_points[i][j])
                m_pts->m_box[j].upper = m_pts->m_points[i][j];
            }
            // if(rank == proc_of_interest) cout << "in Clusters line: 242 m_pts->m_box[j].lower: " << m_pts->m_box[j].lower << endl;
            // if(rank == proc_of_interest) cout << "in Clusters line: 243 m_pts->m_box[j].upper: " << m_pts->m_box[j].upper << endl;
          }
        }

        free(pt);
        pt = NULL;

        file.close();
      } else {
        cout << "rank " << rank << " Error: no such file: " << infilename << endl;
        return -1;
      }
    } else {
      cout << "Only supports binary data: Failed to read data" << endl;
      return -1;
    }
    return 0;
  }
  
  // Called from mpi_main...
  int Clusters::build_kdtree() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if(rank == proc_of_interest) cout << "in Clusters line: 266 " << " in Clusters::build_kdtree" << endl;
    // if the Point objects weren't created, don't bother going further...
    if(m_pts == NULL) {
      cout << "Point set is empty" << endl;
      return -1;
    }

    //m_kdtree = new kdtree2(m_pts->m_points, true);
    m_kdtree = new kdtree2(m_pts->m_points, false);

    if(m_kdtree == NULL) {
      cout << "Falied to allocate new kd tree for orginial points" << endl;
      return -1;
    }

    return 0;   
  } 
  // Called from mpi_main...
  int Clusters::build_kdtree_outer() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if(rank == proc_of_interest) cout << "in Clusters line: 287 " << " in Clusters::build_kdtree_outer" << endl;
    if(m_pts_outer == NULL) {
      cout << "Outer point set is empty" << endl;
      return -1;
    }

    if(m_pts_outer->m_i_num_points > 0) {

      //m_kdtree_outer = new kdtree2(m_pts->m_points_outer, true);
      m_kdtree_outer = new kdtree2(m_pts_outer->m_points, false);

      if(m_kdtree_outer == NULL) {
        cout << "Falied to allocate new kd tree for outer points" << endl;
        return -1;
      }
    }

    return 0;   
  } 
}

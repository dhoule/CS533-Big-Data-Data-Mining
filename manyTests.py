
import subprocess
import os





datasets = ['clus50k.bin', 'part64.bin', 'texture17695.bin', 't8.8k.bin', 't7.10k.bin', 't5.8k.bin', 't4.8k.bin']#, 'edge17695.bin']#, 'part128.bin']
# datasets = ['part64.bin']

# This block is used to create a new list of incremental values.
  # The Rock Cluster doesn't have the `pip` command, and I don't
  # have authorization to add to it. The values must there for
  # be hardcoded before being sent to the cluster. 
# increments = reversed(np.arange(0.90,1.01, 0.01))
# temp = ""
# for increment in increments:
#   temp += str(round(increment,2)) + ", "
# print temp[:-2]

gprof = 0 # 0 = not looking for the profiler. 1 = looking for the profiler files.

# loop over datasets
for dataset in datasets:

  # print "\n\nStarting Dataset: " + dataset + "\n\n"
  # Need to determine the epsilon depending on the dataset used
  if(dataset == 'clus50k.bin'):
    eps = '25'
    mp = '5'
  elif(dataset == 'part64.bin'):
    eps = '0.01'
    mp = '5'
  elif(dataset == 'part128.bin'):
    eps = '0.008'
    mp = '5'
  elif(dataset == 'texture17695.bin'):
    eps = '3'
    mp = '2'
  elif(dataset == 't8.8k.bin'):
    eps = '10'
    mp = '10'
  elif(dataset == 't7.10k.bin'):
    eps = '10'
    mp = '12'
  elif(dataset == 't5.8k.bin'):
    eps = '8'
    mp = '21'
  elif(dataset == 't4.8k.bin'):
    eps = '10'
    mp = '20'
  elif(dataset == 'edge17695.bin'):
    eps = '3'
    mp = '2'

  # The only number of nodes allowed are: 2^x
  nodes = [2,4,8,16]
  # nodes = [2,4]
  datasetFileName = dataset + '_total.txt'
  gDatasetFileName = 'gprof_' + dataset + '_total.txt'

  # loop over allowed nodes
  for node in nodes:
    # print "\tStarting Node Size: " + str(node) + "\n"
    # Uncomment whichever list is needed at the time, before sending to the cluster
    # 1.0-0.01
    # increments = [1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82, 0.81, 0.8, 0.79, 0.78, 0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.7, 0.69, 0.68, 0.67, 0.66, 0.65, 0.64, 0.63, 0.62, 0.61, 0.6, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52, 0.51, 0.5, 0.49, 0.48, 0.47, 0.46, 0.45, 0.44, 0.43, 0.42, 0.41, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.32, 0.31, 0.3, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
    # 1.0-0.9
    # increments = [1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9]
    # 1.0-0.1
    # increments = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    increments = [1.0]
    # The Node's file name is made based on the loop options: dataset_node
    nodeFileName = dataset + '_' + str(node) + '.txt'
    gNodeFileName = 'gprof_' + dataset + '_' + str(node) + '.txt'
    # loop over percentage of dataset to use
    for increment in increments:
      # print "\t\tStarting Increment Size: " + str(increment) + "\n"
      # Create a file to hold the output to CMD line. The name is made based on the loop options: dataset_node_increment
      fileName = dataset + '_' + str(node) + '_' + str(increment) + '.txt'
      gFileName = 'gprof_' + dataset + '_' + str(node) + '_' + str(increment) + '.txt'
      with open(fileName,'w+') as fout:
        # do each test 10 times
        for i in range(10):
          # Each command line argument, must be its own string element. No spaces within the strings.
          # args = ['mpiexec','-n','16','./mpi_dbscan','-i','part64.bin', '-b', '-m','5','-e','0.01', '-k','0.75']
          args = ['mpiexec','-n',str(node),'./mpi_dbscan','-i',dataset, '-b', '-m',mp,'-e',eps, '-p',str(increment)]
          # Syntax: subprocess.call(args, *, stdin=None, stdout=None, stderr=None, shell=False)
          # Output is diverted to file
          subprocess.call(args,stdout=fout)
          if(gprof == 1):
            os.system('gprof -s ./mpi_dbscan gmon.out-*')
            if os.path.exists(gFileName):
              append_write = 'a' # append if already exists
            else:
              append_write = 'w' # make a new file if not
            with open(gFileName,append_write) as gfout:
              os.system('gprof ./mpi_dbscan gmon.sum > temp.txt')
              with open('temp.txt','r') as gin:
                gfout.write("*************  " + str(i) + "  *************\n\n")
                fl = gin.readlines()
                for l in fl:
                  gfout.write("\t" + l)
                gfout.write("\n\n")
            os.system('rm gmon*')
            os.system('rm temp.txt')
      
      # Setting the starting values
      minpts = None
      dimensions = None
      minPtsInClts = None 
      maxPtsInClts = None 
      meanPtsInClts = None 
      minNoise = None
      maxNoise = None
      meanNoise = None
      minNumOfClusters = None
      maxNumOfClusters = None
      meanNumOfClusters = None
      totalPts = None
      minTime = None
      maxTime = None
      meanTime = None
      # Compile the above file info.
      with open(fileName,'r') as fread:
        fl = fread.readlines()
        for l in fl:
          # From the file, need to grab dataset, nodes, eps, minpts, perctOfDB, dimensions, ptsInClts, noise, totalPts, numOfClusters
            # Get min, max, and mean. Min and max need all info preserved. Delete file afterwards.
            # The program output I use has been modified. You will need to update accordingly.
          if (minpts == None) and ("MinPts" in l):
            minpts = l.split(' ')[3].strip()
          
          if (dimensions == None) and ("Dimensions" in l):
            dimensions = l.split(':')[-1].strip()
          
          if "Points" in l:
            t = l.split(' ')

            cluster = int(t[3].strip().strip(','))
            if None == minPtsInClts:
              minPtsInClts = cluster
              maxPtsInClts = cluster
              meanPtsInClts = 0
            if cluster < minPtsInClts:
              minPtsInClts = cluster
            if cluster > maxPtsInClts:
              maxPtsInClts = cluster
            meanPtsInClts += cluster

            noise = int(t[5].strip().strip(','))
            if None == minNoise:
              minNoise = noise
              maxNoise = noise
              meanNoise = 0
            if noise < minNoise:
              minNoise = noise 
            if noise > maxNoise:
              maxNoise = noise 
            meanNoise += noise

            if totalPts == None:
              totalPts = int(t[8].strip())

          if "Total number" in l:
            cluster = int(l.split(' ')[4])
            if None == minNumOfClusters:
              minNumOfClusters = cluster
              maxNumOfClusters = cluster
              meanNumOfClusters = 0
            if cluster < minNumOfClusters:
              minNumOfClusters = cluster
            if cluster > maxNumOfClusters:
              maxNumOfClusters = cluster
            meanNumOfClusters += cluster

          if "DBSCAN" in l:
            time = float(l.split(' ')[8])
            if None == minTime:
              minTime = time
              maxTime = time
              meanTime = 0
            if time < minTime:
              minTime = time
            if time > maxTime:
              maxTime = time
            meanTime += time
    
      # Get the averages
      meanPtsInClts = meanPtsInClts/10.0
      meanNoise = meanNoise /10.0
      meanNumOfClusters = meanNumOfClusters/10.0
      meanTime = meanTime/10.0

      # Delete the unneeded file
      os.remove(fileName)
      flag = False
      if os.path.exists(nodeFileName):
        append_write = 'a' # append if already exists
      else:
        append_write = 'w' # make a new file if not

      # Create/append the Node's file to hold the computed values. 
      with open(nodeFileName,append_write) as fout:
        fout.write(str(node) + "," + str(increment) + "," + str(eps) + "," + str(minpts) + "," + str(dimensions) + "," + str(minPtsInClts) + "," + str(maxPtsInClts) + "," + str(meanPtsInClts) + "," + str(minNoise) + "," + str(maxNoise) + "," + str(meanNoise) + "," + str(minNumOfClusters) + "," + str(maxNumOfClusters) + "," + str(meanNumOfClusters) + "," + str(totalPts) + "," + str(minTime) + "," + str(maxTime) + "," + str(meanTime) + "\n")
    
      if(gprof == 1):
        if os.path.exists(gNodeFileName):
          append_write = 'a' # append if already exists
        else:
          append_write = 'w' # make a new file if not
        with open(gNodeFileName,append_write) as gfout:
          gfout.write("*************  " + str(increment) + "  *************\n\n")
          with open(gFileName,'r') as fin:
            fl = fin.readlines()
            for l in fl:
              gfout.write("\t" + l)
        os.remove(gFileName)

    if os.path.exists(datasetFileName):
      append_write = 'a' # append if already exists
    else:
      flag = True
      append_write = 'w' # make a new file if not

    # Create/append the Node's file to hold the computed values. 
    with open(datasetFileName,append_write) as fout:
      if flag:
        fout.write("Nodes,Percent_used,Epsilon,MinPts,Dims,MinPts_in_clts,MaxPts_in_clts,MeanPts_in_clts,MinNoise,MaxNoise,MeanNoise,MinNumOfClusters,MaxNumOfClusters,MeanNumOfClusters,TotalPts,MinTime,MaxTime,MeanTime\n")
      with open(nodeFileName,'r') as fread:
        fl = fread.readlines()
        for l in fl:
          fout.write(l)

    if(gprof == 1):
      if os.path.exists(gDatasetFileName):
        append_write = 'a' # append if already exists
      else:
        flag = True
        append_write = 'w' # make a new file if not
      with open(gDatasetFileName,append_write) as fout:
        fout.write("*************  " + str(node) + "  *************\n\n")
        with open(gNodeFileName,'r') as fread:
          fl = fread.readlines()
          for l in fl:
            fout.write("\t" + l)
          fout.write("\n\n")
      os.remove(gNodeFileName)

    # The individual node file is no longer needed, as it has been joined into the dataset file
    os.remove(nodeFileName)
    # print "\tEnded Node Size: " + str(node) + "\n"

  # print "\n\nEnded Dataset: " + dataset + "\n\n"


# Output of ./mpi_dbscan. NOTE: This has been changed from the original output.
# Dataset used: clus50k.bin
# Number of process cores: 2
# Epsilon: 25 MinPts: 5 Percent_of_dataset_used: 1
# Dimensions of each point: 10
# Parallel DBSCAN (init, local computation, and merging) took 19.6082 seconds

# Points in clusters: 46914, Noise: 3086, Total points: 50000
# Total number of clusters: 51

# Contents of dataset_node.txt files
# Nodes: 2
# Percent_used: 0.96
# Epsilon: 0.01
# MinPts: 5
# Dims: 3
# MinPts_in_clts: 1259
# MaxPts_in_clts: 1286
# MeanPts_in_clts: 1275.4
# MinNoise: 60154
# MaxNoise: 60181
# MeanNoise: 60164.6
# MinNumOfClusters: 113
# MaxNumOfClusters: 116
# MeanNumOfClusters: 114.6
# TotalPts: 61440
# MinTime: 0.421928
# MaxTime: 1.42898
# MeanTime: 0.6606667







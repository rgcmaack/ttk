/// \ingroup vtk
/// \class ttkMergeTreeVisualization
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Visualization module for merge trees.

#pragma once

#include <FTMTree.h>
#include <MergeTreeVisualization.h>

#include <ttkAlgorithm.h>

// VTK Includes
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

class ttkMergeTreeVisualization : public ttk::MergeTreeVisualization {
private:
  // redefine types
  using idNode = ttk::ftm::idNode;
  using SimplexId = ttk::SimplexId;
  using FTMTree_MT = ttk::ftm::FTMTree_MT;

  // Visualization parameters
  bool PlanarLayout = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  bool OutputSegmentation = false;
  int ShiftMode = 0; // 0: Star ; 1: Star Barycenter ; 2: Line ; 3: Double Line
  int MaximumImportantPairs = 0;
  int MinimumImportantPairs = 0;

  // Offset
  int iSampleOffset = 0;
  int noSampleOffset = 0;
  double prevXMaxOffset = 0;

  // Print only one tree
  int printTreeId = -1; // -1 for all
  int printClusterId = -1; // -1 for all

  // Barycenter position according alpha
  bool BarycenterPositionAlpha = false;
  double Alpha = 0.5;

  // Used for critical type and for point coordinates
  std::vector<vtkUnstructuredGrid *> treesNodes;
  std::vector<std::vector<int>>
    treesNodeCorrMesh; // used to acces treesNodes given input trees

  // Segmentation
  std::vector<vtkDataSet *> treesSegmentation;

  // Clustering output
  std::vector<int> clusteringAssignment;
  std::vector<std::vector<std::vector<std::tuple<idNode, idNode, double>>>>
    outputMatchingBarycenter;

  // Barycenter output
  std::vector<std::vector<float>> allBaryPercentMatch;

  // Temporal Subsampling Output
  std::vector<bool> interpolatedTrees;

  // Output
  vtkUnstructuredGrid *vtkOutputNode{};
  vtkUnstructuredGrid *vtkOutputArc{};
  vtkUnstructuredGrid *vtkOutputSegmentation{};

  // Matching output
  vtkUnstructuredGrid *vtkOutputNode1{}, *vtkOutputNode2{}; // input data
  std::vector<std::vector<SimplexId>> nodeCorr1, nodeCorr2;
  vtkUnstructuredGrid *vtkOutputMatching{}; // output

  // Filled by the algorithm
  std::vector<std::vector<SimplexId>> nodeCorr;
  std::vector<double> clusterShift;
  double prevXMax = 0;

public:
  ttkMergeTreeVisualization(){};
  ~ttkMergeTreeVisualization(){};

  // ==========================================================================
  // Getter / Setter
  // ==========================================================================
  // Visualization parameters
  void setPlanarLayout(bool b) {
    PlanarLayout = b;
  }
  void setDimensionSpacing(double d) {
    DimensionSpacing = d;
  }
  void setDimensionToShift(int i) {
    DimensionToShift = i;
  }
  void setOutputSegmentation(bool b) {
    OutputSegmentation = b;
  }
  void setShiftMode(int mode) {
    ShiftMode = mode;
  }
  void setMaximumImportantPairs(int maxPairs) {
    MaximumImportantPairs = maxPairs;
  }
  void setMinimumImportantPairs(int minPairs) {
    MinimumImportantPairs = minPairs;
  }

  // Offset
  void setISampleOffset(int offset) {
    iSampleOffset = offset;
  }
  void setNoSampleOffset(int offset) {
    noSampleOffset = offset;
  }
  void setPrevXMaxOffset(double offset) {
    prevXMaxOffset = offset;
  }

  // Print only one tree
  void setPrintTreeId(int id) {
    printTreeId = id;
  }
  void setPrintClusterId(int id) {
    printClusterId = id;
  }

  // Barycenter position according alpha
  void setBarycenterPositionAlpha(bool pos) {
    BarycenterPositionAlpha = pos;
  }
  void setAlpha(double alpha) {
    Alpha = alpha;
  }

  // Used for critical type and if not planar layout for point coordinates
  void setTreesNodes(std::vector<vtkUnstructuredGrid *> &nodes) {
    treesNodes = nodes;
  }
  void setTreesNodeCorrMesh(std::vector<std::vector<int>> &nodeCorrMesh) {
    treesNodeCorrMesh = nodeCorrMesh;
  }
  void setTreesNodes(vtkUnstructuredGrid *nodes) {
    treesNodes.clear();
    treesNodes.emplace_back(nodes);
  }
  void setTreesNodeCorrMesh(std::vector<int> &nodeCorrMesh) {
    treesNodeCorrMesh.clear();
    treesNodeCorrMesh.emplace_back(nodeCorrMesh);
  }

  // Segmentation
  void setTreesSegmentation(std::vector<vtkDataSet *> &segmentation) {
    treesSegmentation = segmentation;
  }
  void setTreesSegmentation(vtkDataSet *segmentation) {
    treesSegmentation.clear();
    treesSegmentation.emplace_back(segmentation);
  }

  // Clustering output
  void setClusteringAssignment(std::vector<int> &asgn) {
    clusteringAssignment = asgn;
  }
  void setOutputMatchingBarycenter(
    std::vector<std::vector<std::vector<std::tuple<idNode, idNode, double>>>>
      &matching) {
    outputMatchingBarycenter = matching;
  }

  // Barycenter output
  std::vector<std::vector<float>> getAllBaryPercentMatch() {
    return allBaryPercentMatch;
  }
  void
    setAllBaryPercentMatch(std::vector<std::vector<float>> &baryPercentMatch) {
    allBaryPercentMatch = baryPercentMatch;
  }

  // Temporal Subsampling Output
  void setInterpolatedTrees(std::vector<bool> &isInterpolatedTrees) {
    interpolatedTrees = isInterpolatedTrees;
  }

  // Output
  void setVtkOutputNode(vtkUnstructuredGrid *vtkNode) {
    vtkOutputNode = vtkNode;
  }
  void setVtkOutputArc(vtkUnstructuredGrid *vtkArc) {
    vtkOutputArc = vtkArc;
  }
  void setVtkOutputSegmentation(vtkUnstructuredGrid *vtkSegmentation) {
    vtkOutputSegmentation = vtkSegmentation;
  }

  // Matching output
  void setVtkOutputNode1(vtkUnstructuredGrid *vtkNode1) {
    vtkOutputNode1 = vtkNode1;
  }
  void setVtkOutputNode2(vtkUnstructuredGrid *vtkNode2) {
    vtkOutputNode2 = vtkNode2;
  }
  void setNodeCorr1(std::vector<std::vector<SimplexId>> &nodeCorrT) {
    nodeCorr1 = nodeCorrT;
  }
  void setNodeCorr2(std::vector<std::vector<SimplexId>> &nodeCorrT) {
    nodeCorr2 = nodeCorrT;
  }
  void setVtkOutputMatching(vtkUnstructuredGrid *vtkMatching) {
    vtkOutputMatching = vtkMatching;
  }
  void setOutputMatching(
    std::vector<std::tuple<idNode, idNode, double>> &matching) {
    outputMatchingBarycenter = std::vector<
      std::vector<std::vector<std::tuple<idNode, idNode, double>>>>(1);
    outputMatchingBarycenter[0]
      = std::vector<std::vector<std::tuple<idNode, idNode, double>>>(1);
    outputMatchingBarycenter[0][0] = matching;
  }

  // Filled by the algorithm
  std::vector<std::vector<SimplexId>> getNodeCorr() {
    return nodeCorr;
  }
  std::vector<double> getClusterShift() {
    return clusterShift;
  }
  double getPrevXMax() {
    return prevXMax;
  }

  // ==========================================================================
  // Matching Visualization
  // ==========================================================================
  template <class dataType>
  void makeMatchingOutput(FTMTree_MT *tree1, FTMTree_MT *tree2) {
    std::vector<FTMTree_MT *> trees{tree1, tree2};
    std::vector<FTMTree_MT *> barycenters;

    makeMatchingOutput<dataType>(trees, barycenters);
  }

  template <class dataType>
  void makeMatchingOutput(std::vector<FTMTree_MT *> &trees,
                          std::vector<FTMTree_MT *> &barycenters) {
    int numInputs = trees.size();
    int NumberOfBarycenters = barycenters.size();
    bool clusteringOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters
      = std::max(NumberOfBarycenters, 1); // to always enter the outer loop
    if(not clusteringOutput)
      numInputs = 1;

    vtkNew<vtkUnstructuredGrid> vtkMatching{};
    vtkNew<vtkPoints> pointsM{};

    // Fields
    vtkNew<vtkIntArray> matchingID{};
    matchingID->SetName("MatchingID");
    vtkNew<vtkIntArray> matchingType{};
    matchingType->SetName("MatchingType");
    vtkNew<vtkDoubleArray> matchPers{};
    matchPers->SetName("MeanMatchedPersistence");
    vtkNew<vtkDoubleArray> costArray{};
    costArray->SetName("Cost");
    vtkNew<vtkIntArray> tree1NodeIdField{};
    tree1NodeIdField->SetName("tree1NodeId");
    vtkNew<vtkIntArray> tree2NodeIdField{};
    tree2NodeIdField->SetName("tree2NodeId");

    vtkNew<vtkFloatArray> matchingPercentMatch{};
    matchingPercentMatch->SetName("MatchingPercentMatch");

    // Iterate through clusters and trees
    printMsg(
      "// Iterate through clusters and trees", ttk::debug::Priority::VERBOSE);
    int count = 0;
    for(int c = 0; c < NumberOfBarycenters; ++c) {
      for(int i = 0; i < numInputs; ++i) {
        if((printTreeId == -1 and printClusterId != -1 and c != printClusterId)
           or (printTreeId != -1 and printClusterId == -1 and i != printTreeId)
           or (printTreeId != -1 and printClusterId != -1
               and (c != printClusterId or i != printTreeId)))
          continue;
        for(std::tuple<idNode, idNode, double> match :
            outputMatchingBarycenter[c][i]) {
          vtkIdType pointIds[2];
          idNode tree1NodeId = std::get<0>(match);
          idNode tree2NodeId = std::get<1>(match);
          double cost = std::get<2>(match);
          FTMTree_MT *tree1;
          FTMTree_MT *tree2 = trees[i];
          if(not clusteringOutput) {
            tree1 = trees[0];
            tree2 = trees[1];
          } else
            tree1 = barycenters[c];

          // Get first point
          printMsg("// Get first point", ttk::debug::Priority::VERBOSE);
          SimplexId pointToGet1 = clusteringOutput ? nodeCorr2[c][tree1NodeId]
                                                   : nodeCorr1[0][tree1NodeId];
          double *point1 = vtkOutputNode2->GetPoints()->GetPoint(pointToGet1);
          const SimplexId nextPointId1 = pointsM->InsertNextPoint(point1);
          pointIds[0] = nextPointId1;

          // Get second point
          printMsg("// Get second point", ttk::debug::Priority::VERBOSE);
          SimplexId pointToGet2 = clusteringOutput ? nodeCorr1[i][tree2NodeId]
                                                   : nodeCorr1[1][tree2NodeId];
          double *point2 = vtkOutputNode1->GetPoints()->GetPoint(pointToGet2);
          const SimplexId nextPointId2 = pointsM->InsertNextPoint(point2);
          pointIds[1] = nextPointId2;

          // Add cell
          printMsg("// Add cell", ttk::debug::Priority::VERBOSE);
          vtkMatching->InsertNextCell(VTK_LINE, 2, pointIds);

          // Add arc matching percentage
          printMsg(
            "// Add arc matching percentage", ttk::debug::Priority::VERBOSE);
          if(allBaryPercentMatch.size() != 0)
            matchingPercentMatch->InsertNextTuple1(
              allBaryPercentMatch[c][tree1NodeId]);

          // Add tree1 and tree2 node ids
          printMsg(
            "// Add tree1 and tree2 node ids", ttk::debug::Priority::VERBOSE);
          tree1NodeIdField->InsertNextTuple1(pointToGet1);
          tree2NodeIdField->InsertNextTuple1(pointToGet2);

          // Add matching ID
          matchingID->InsertNextTuple1(count);

          // Add matching type
          printMsg("// Add matching type", ttk::debug::Priority::VERBOSE);
          int thisType = 0;
          int tree1NodeDown
            = tree1->getNode(tree1NodeId)->getNumberOfDownSuperArcs();
          int tree1NodeUp
            = tree1->getNode(tree1NodeId)->getNumberOfUpSuperArcs();
          int tree2NodeDown
            = tree2->getNode(tree2NodeId)->getNumberOfDownSuperArcs();
          int tree2NodeUp
            = tree2->getNode(tree2NodeId)->getNumberOfUpSuperArcs();
          if(tree1NodeDown != 0 and tree1NodeUp != 0 and tree2NodeDown != 0
             and tree2NodeUp != 0)
            thisType = 1; // Saddle to Saddle
          if(tree1NodeDown == 0 and tree1NodeUp != 0 and tree2NodeDown == 0
             and tree2NodeUp != 0)
            thisType = 2; // Leaf to leaf
          if(tree1NodeDown != 0 and tree1NodeUp == 0 and tree2NodeDown != 0
             and tree2NodeUp == 0)
            thisType = 3; // Root to root
          matchingType->InsertNextTuple1(thisType);

          // Add mean matched persistence
          printMsg(
            "// Add mean matched persistence", ttk::debug::Priority::VERBOSE);
          double tree1Pers = tree1->getNodePersistence<dataType>(tree1NodeId);
          double tree2Pers = tree2->getNodePersistence<dataType>(tree2NodeId);
          double meanPersistence = (tree1Pers + tree2Pers) / 2;
          matchPers->InsertNextTuple1(meanPersistence);

          // Add cost
          costArray->InsertNextTuple1(cost);

          count++;
        }
      }
    }
    vtkMatching->SetPoints(pointsM);
    vtkMatching->GetCellData()->AddArray(matchingType);
    vtkMatching->GetCellData()->AddArray(matchPers);
    vtkMatching->GetCellData()->AddArray(matchingID);
    vtkMatching->GetCellData()->AddArray(costArray);
    vtkMatching->GetCellData()->AddArray(tree1NodeIdField);
    vtkMatching->GetCellData()->AddArray(tree2NodeIdField);
    if(allBaryPercentMatch.size() != 0)
      vtkMatching->GetCellData()->AddArray(matchingPercentMatch);
    vtkOutputMatching->ShallowCopy(vtkMatching);
  }

  // ==========================================================================
  // Trees Visualization
  // ==========================================================================
  template <class dataType>
  void makeTreesOutput(FTMTree_MT *tree1) {
    std::vector<FTMTree_MT *> trees{tree1};

    makeTreesOutput<dataType>(trees);
  }

  template <class dataType>
  void makeTreesOutput(FTMTree_MT *tree1, FTMTree_MT *tree2) {
    std::vector<FTMTree_MT *> trees{tree1, tree2};

    makeTreesOutput<dataType>(trees);
  }

  template <class dataType>
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees) {
    std::vector<FTMTree_MT *> barycenters;
    clusteringAssignment = std::vector<int>(trees.size(), 0);

    makeTreesOutput<dataType>(trees, barycenters);
  }

  template <class dataType>
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees,
                       std::vector<FTMTree_MT *> &barycenters) {
    int numInputs = trees.size();
    int numInputsOri = numInputs;
    int NumberOfBarycenters = barycenters.size();
    bool clusteringOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters
      = std::max(NumberOfBarycenters, 1); // to always enter the outer loop

    // Bounds
    printMsg("Bounds and branching", ttk::debug::Priority::VERBOSE);
    std::vector<std::tuple<double, double, double, double, double, double>>
      allBounds(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      if(OutputSegmentation) {
        double *tBounds = treesSegmentation[i]->GetBounds();
        allBounds[i] = std::make_tuple(tBounds[0], tBounds[1], tBounds[2],
                                       tBounds[3], tBounds[4], tBounds[5]);
      } else if(treesNodes[i] != nullptr)
        allBounds[i]
          = getRealBounds(treesNodes[i], trees[i], treesNodeCorrMesh[i]);
      else
        allBounds[i] = allBounds[0];
      /*if(PlanarLayout){
        // TODO correctly manage bounds when planar layout
        std::vector<double> tBound = tupleToVector(allBounds[i]);
        int level = getTreeDepth(trees[i]);
        for(int j = 0; j < tBound.size(); ++j)
          tBound[j] = tBound[j]*level/5;
        allBounds[i] = vectorToTuple(tBound);
      }*/
    }
    std::vector<std::tuple<double, double, double, double, double, double>>
      allBaryBounds(barycenters.size());
    std::vector<std::vector<idNode>> allBaryBranching(barycenters.size());
    std::vector<std::vector<int>> allBaryBranchingID(barycenters.size());
    for(size_t c = 0; c < barycenters.size(); ++c) {
      allBaryBounds[c] = getMaximalBounds(allBounds, clusteringAssignment, c);
      barycenters[c]->getTreeBranching(
        allBaryBranching[c], allBaryBranchingID[c]);
    }
    if(not clusteringOutput)
      allBaryBounds.emplace_back(
        getMaximalBounds(allBounds, clusteringAssignment, 0));

    // ----------------------------------------------------------------------
    // Make Trees Output
    // ----------------------------------------------------------------------
    printMsg("--- Make Trees Output", ttk::debug::Priority::VERBOSE);
    std::vector<FTMTree_MT *> treesOri(trees);
    if(ShiftMode == 1) { // Star Barycenter
      trees.clear();
      clusteringAssignment.clear();
      for(unsigned int j = 0; j < barycenters.size(); ++j) {
        trees.emplace_back(barycenters[j]);
        clusteringAssignment.emplace_back(j);
      }
      numInputs = trees.size();
    }
    // - Declare VTK arrays
    vtkNew<vtkUnstructuredGrid> vtkArcs{};
    vtkNew<vtkPoints> points{};

    // Node fields
    vtkNew<vtkIntArray> criticalType{};
    criticalType->SetName("CriticalType");
    vtkNew<vtkFloatArray> persistenceNode{};
    persistenceNode->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDNode{};
    clusterIDNode->SetName("ClusterID");
    vtkNew<vtkIntArray> isDummyNode{};
    isDummyNode->SetName("isDummyNode");
    vtkNew<vtkIntArray> branchNodeID{};
    branchNodeID->SetName("BranchNodeID");
    vtkNew<vtkFloatArray> scalar{};
    scalar->SetName("Scalar");
    vtkNew<vtkIntArray> isImportantPairsNode{};
    isImportantPairsNode->SetName("isImportantPair");
    vtkNew<vtkIntArray> nodeID{};
    nodeID->SetName("NodeId");

    vtkNew<vtkIntArray> treeIDNode{};
    treeIDNode->SetName("TreeID");
    vtkNew<vtkIntArray> branchBaryNodeID{};
    branchBaryNodeID->SetName("BranchBaryNodeID");
    vtkNew<vtkIntArray> isInterpolatedTreeNode{};
    isInterpolatedTreeNode->SetName("isInterpolatedTree");

    vtkNew<vtkFloatArray> percentMatch{};
    percentMatch->SetName("PercentMatchNode");
    vtkNew<vtkFloatArray> persistenceBaryNode{};
    persistenceBaryNode->SetName("PersistenceBarycenter");

    // Arc fields
    vtkNew<vtkFloatArray> persistenceArc{};
    persistenceArc->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDArc{};
    clusterIDArc->SetName("ClusterID");
    vtkNew<vtkIntArray> isImportantPairsArc{};
    isImportantPairsArc->SetName("isImportantPair");
    vtkNew<vtkIntArray> isDummyArc{};
    isDummyArc->SetName("isDummyArc");
    vtkNew<vtkIntArray> branchID{};
    branchID->SetName("BranchID");
    vtkNew<vtkIntArray> upNodeId{};
    upNodeId->SetName("upNodeId");
    vtkNew<vtkIntArray> downNodeId{};
    downNodeId->SetName("downNodeId");

    vtkNew<vtkIntArray> treeIDArc{};
    treeIDArc->SetName("TreeID");
    vtkNew<vtkIntArray> branchBaryID{};
    branchBaryID->SetName("BranchBaryID");
    vtkNew<vtkIntArray> isInterpolatedTreeArc{};
    isInterpolatedTreeArc->SetName("isInterpolatedTree");

    vtkNew<vtkFloatArray> percentMatchArc{};
    percentMatchArc->SetName("PercentMatchArc");
    vtkNew<vtkFloatArray> persistenceBaryArc{};
    persistenceBaryArc->SetName("PersistenceBarycenter");

    // Segmentation
    vtkNew<vtkAppendFilter> appendFilter{};

    // Internal data
    int cellCount = 0;
    int pointCount = 0;
    bool foundOneInterpolatedTree = false;
    nodeCorr = std::vector<std::vector<SimplexId>>(numInputs);
    clusterShift = std::vector<double>(NumberOfBarycenters, 0);
    allBaryPercentMatch = std::vector<std::vector<float>>(NumberOfBarycenters);

    // --------------------------------------------------------
    // Iterate through all clusters
    // --------------------------------------------------------
    printMsg("Iterate through all clusters", ttk::debug::Priority::VERBOSE);
    double importantPairsOriginal = importantPairs_;
    for(int c = 0; c < NumberOfBarycenters; ++c) {

      // Get radius
      printMsg("// Get radius", ttk::debug::Priority::VERBOSE);
      double delta_max = std::numeric_limits<double>::lowest();
      int noSample = 0 + noSampleOffset;
      for(int i = 0; i < numInputsOri; ++i) {
        delta_max = std::max(
          (std::get<3>(allBounds[i]) - std::get<2>(allBounds[i])), delta_max);
        delta_max = std::max(
          (std::get<1>(allBounds[i]) - std::get<0>(allBounds[i])), delta_max);
        if(clusteringAssignment[i] != c)
          continue;
        noSample += 1;
      }
      double radius = delta_max * 2 * DimensionSpacing;
      int iSample = 0 + iSampleOffset - 1;

      if(c < NumberOfBarycenters - 1)
        clusterShift[c + 1] = radius * 4 + clusterShift[c];

      // Line/Double line attributes
      prevXMax = 0 + prevXMaxOffset;
      std::vector<double> allPrevXMax;
      double prevYMax = std::numeric_limits<double>::lowest();

      // ------------------------------------------
      // Iterate through all trees of this cluster
      // ------------------------------------------
      printMsg("Iterate through all trees of this cluster",
               ttk::debug::Priority::VERBOSE);
      for(int i = 0; i < numInputs; ++i) {
        if(clusteringAssignment[i] != c)
          continue;

        iSample += 1;

        if((printTreeId == -1 and printClusterId != -1 and c != printClusterId)
           or (printTreeId != -1 and printClusterId == -1 and i != printTreeId)
           or (printTreeId != -1 and printClusterId != -1
               and (c != printClusterId or i != printTreeId)))
          continue;

        // Manage important pairs threshold
        importantPairs_ = importantPairsOriginal;
        if(MaximumImportantPairs > 0 or MinimumImportantPairs > 0) {
          std::vector<std::tuple<idNode, idNode, dataType>> pairs;
          trees[i]->getPersistencePairsFromTree(pairs, false);
          if(MaximumImportantPairs > 0) {
            double tempThreshold
              = 0.999 * std::get<2>(pairs[pairs.size() - MaximumImportantPairs])
                / std::get<2>(pairs[pairs.size() - 1]);
            tempThreshold *= 100;
            importantPairs_ = std::max(importantPairs_, tempThreshold);
          }
          if(MinimumImportantPairs > 0) {
            double tempThreshold
              = 0.999 * std::get<2>(pairs[pairs.size() - MinimumImportantPairs])
                / std::get<2>(pairs[pairs.size() - 1]);
            tempThreshold *= 100;
            importantPairs_ = std::min(importantPairs_, tempThreshold);
          }
        }

        // Get is interpolated tree (temporal subsampling)
        bool isInterpolatedTree = false;
        if(interpolatedTrees.size() != 0)
          isInterpolatedTree = interpolatedTrees[i];
        foundOneInterpolatedTree |= isInterpolatedTree;

        // Get branching
        printMsg("// Get branching", ttk::debug::Priority::VERBOSE);
        std::vector<idNode> treeBranching;
        std::vector<int> treeBranchingID;
        trees[i]->getTreeBranching(treeBranching, treeBranchingID);

        // Get shift
        printMsg("// Get shift", ttk::debug::Priority::VERBOSE);
        double angle = 360.0 / noSample * iSample;
        double pi = 3.14159265359;
        double diff_x = 0, diff_y = 0;
        double alphaShift
          = BarycenterPositionAlpha ? (-radius + 2 * radius * Alpha) * -1 : 0;
        switch(ShiftMode) {
          case 0: // Star
            diff_x
              = -1 * radius * std::cos(-1 * angle * pi / 180) + clusterShift[c];
            diff_y = -1 * radius * std::sin(-1 * angle * pi / 180);
            break;
          case 1: // Star Barycenter
            diff_x = clusterShift[c] + alphaShift;
            diff_y = 0;
            break;
          case 2: // Line
            diff_x = prevXMax + radius;
            break;
          case 3: // Double Line
            diff_x = prevXMax + radius;
            if(i >= numInputs / 2) {
              diff_y = -(prevYMax + radius / 2);
              diff_x = allPrevXMax[i - int(numInputs / 2)] + radius;
            } else
              allPrevXMax.emplace_back(prevXMax);
            break;
          default:
            break;
        }

        // Get dimension shift
        printMsg("// Get dimension shift", ttk::debug::Priority::VERBOSE);
        double diff_z = PlanarLayout ? 0 : -std::get<4>(allBounds[i]);
        // TODO DimensionToShift for Planar Layout
        if(not PlanarLayout)
          if(DimensionToShift != 0) { // is not X
            if(DimensionToShift == 2) // is Z
              diff_z = diff_x;
            else if(DimensionToShift == 1) // is Y
              diff_y = diff_x;
            diff_x = -std::get<0>(allBounds[i]);
          }

        // Planar layout
        printMsg("// Planar Layout", ttk::debug::Priority::VERBOSE);
        std::vector<float> layout;
        if(PlanarLayout) {
          double refPersistence;
          if(clusteringOutput)
            refPersistence = barycenters[0]->getNodePersistence<dataType>(
              barycenters[0]->getRoot());
          else
            refPersistence
              = trees[0]->getNodePersistence<dataType>(trees[0]->getRoot());
          treePlanarLayout<dataType>(
            trees[i], allBaryBounds[c], refPersistence, layout);
        }

        // Internal arrays
        printMsg("// Internal arrays", ttk::debug::Priority::VERBOSE);
        int cptNode = 0;
        nodeCorr[i] = std::vector<SimplexId>(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeSimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeDummySimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> layoutCorr(trees[i]->getNumberOfNodes());
        std::vector<idNode> treeMatching(trees[i]->getNumberOfNodes(), -1);
        if(clusteringOutput)
          for(auto match : outputMatchingBarycenter[c][i])
            treeMatching[std::get<1>(match)] = std::get<0>(match);
        // _ m[i][j] contains the node in treesOri[j] matched to the node i in
        // the barycenter
        std::vector<std::vector<idNode>> baryMatching(
          trees[i]->getNumberOfNodes(), std::vector<idNode>(numInputsOri, -1));
        if(ShiftMode == 1) {
          for(size_t j = 0; j < outputMatchingBarycenter[c].size(); ++j)
            for(auto match : outputMatchingBarycenter[c][j])
              baryMatching[std::get<0>(match)][j] = std::get<1>(match);
          allBaryPercentMatch[c]
            = std::vector<float>(trees[i]->getNumberOfNodes(), 100);
        }

        // ----------------------------
        // Tree traversal
        // ----------------------------
        printMsg("// Tree traversal", ttk::debug::Priority::VERBOSE);
        std::queue<idNode> queue;
        queue.emplace(trees[i]->getRoot());
        while(!queue.empty()) {
          idNode node = queue.front();
          queue.pop();
          idNode nodeOrigin = trees[i]->getNode(node)->getOrigin();

          // Push children to the queue
          printMsg(
            "// Push children to the queue", ttk::debug::Priority::VERBOSE);
          std::vector<idNode> children;
          trees[i]->getChildren(node, children);
          for(auto child : children)
            queue.emplace(child);

          // --------------
          // Insert point
          // --------------
          printMsg("// Get and insert point", ttk::debug::Priority::VERBOSE);
          int nodeMesh = -1;
          int nodeMeshTreeIndex = -1;
          double noMatched = 0.0;
          double point[3] = {0, 0, 0};
          if(ShiftMode == 1) { // Star barycenter
            for(int j = 0; j < numInputsOri; ++j) {
              if((int)baryMatching[node][j] != -1) {
                nodeMesh = treesNodeCorrMesh[j][baryMatching[node][j]];
                double *pointTemp
                  = treesNodes[j]->GetPoints()->GetPoint(nodeMesh);
                for(int k = 0; k < 3; ++k)
                  point[k] += pointTemp[k];
                noMatched += 1;
                nodeMeshTreeIndex = j;
              }
            }
            for(int k = 0; k < 3; ++k)
              point[k] /= noMatched;
          } else if(not isInterpolatedTree) {
            nodeMesh = treesNodeCorrMesh[i][node];
            double *pointP = treesNodes[i]->GetPoints()->GetPoint(nodeMesh);
            for(int p = 0; p < 3; ++p)
              point[p] = pointP[p];
          }
          if(PlanarLayout) {
            layoutCorr[node] = cptNode;
            point[0] = layout[cptNode];
            point[1] = layout[cptNode + 1];
            point[2] = 0;
            cptNode += 2;
          }
          point[0] += diff_x;
          point[1] += diff_y;
          point[2] += diff_z;

          // Get x Max and y Min for next iteration if needed (double line mode)
          prevXMax = std::max(prevXMax, point[0]);
          if(ShiftMode == 3) { // Double line
            if(i < numInputs / 2)
              prevYMax = std::max(prevYMax, point[1]);
            if(i == int(numInputs / 2) - 1)
              prevXMax = 0;
          }

          // TODO too many dummy nodes are created
          bool dummyNode = PlanarLayout and not branchDecompositionPlanarLayout_
                           and !trees[i]->isRoot(node)
            /*and !isLeaf(trees[i], node)
            and isBranchOrigin(trees[i], node)*/
            ;
          if(dummyNode)
            treeDummySimplexId[node] = points->InsertNextPoint(
              point); // will be modified when processing son
          SimplexId nextPointId = points->InsertNextPoint(point);
          treeSimplexId[node] = nextPointId;
          nodeCorr[i][node] = nextPointId;
          if(dummyNode)
            nodeCorr[i][node] = treeDummySimplexId[node];

          // --------------
          // Insert cell connecting parent
          // --------------
          printMsg(
            "// Add cell connecting parent", ttk::debug::Priority::VERBOSE);
          if(!trees[i]->isRoot(node)) {
            vtkIdType pointIds[2];
            pointIds[0] = treeSimplexId[node];

            idNode nodeParent = trees[i]->getParentSafe(node);
            // TODO too many dummy cells are created
            bool dummyCell = PlanarLayout
                             and not branchDecompositionPlanarLayout_
                             and treeBranching[node] == nodeParent
                             and !trees[i]->isRoot(nodeParent);
            if(PlanarLayout and branchDecompositionPlanarLayout_) {
              pointIds[1] = treeSimplexId[treeBranching[node]];
            } else if(dummyCell) {
              double dummyPoint[3]
                = {point[0], layout[layoutCorr[nodeParent] + 1] + diff_y,
                   0. + diff_z};
              SimplexId dummyPointId = treeDummySimplexId[nodeParent];
              points->SetPoint(dummyPointId, dummyPoint);
              vtkIdType dummyPointIds[2];
              dummyPointIds[0] = dummyPointId;
              dummyPointIds[1] = treeSimplexId[nodeParent];
              vtkArcs->InsertNextCell(VTK_LINE, 2, dummyPointIds);
              pointIds[1] = dummyPointId;
            } else
              pointIds[1] = treeSimplexId[nodeParent];

            vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);

            // --------------
            // Arc field
            // --------------
            int toAdd = (dummyCell ? 2 : 1);
            for(int toAddT = 0; toAddT < toAdd; ++toAddT) {
              // Add arc matching percentage
              if(ShiftMode == 1) // Star Barycenter
                percentMatchArc->InsertNextTuple1(
                  allBaryPercentMatch[c][allBaryBranching[c][node]]);

              // Add branch bary ID
              printMsg(
                "// Push arc bary branch id", ttk::debug::Priority::VERBOSE);
              if(clusteringOutput and ShiftMode != 1) {
                int tBranchID = -1;
                if(treeMatching[node] < allBaryBranchingID[c].size()) {
                  tBranchID = allBaryBranchingID[c][treeMatching[node]];
                  if(!trees[i]->isLeaf(node)
                     && treeMatching[nodeOrigin] < allBaryBranchingID[c].size())
                    tBranchID = allBaryBranchingID[c][treeMatching[nodeOrigin]];
                }
                branchBaryID->InsertNextTuple1(tBranchID);
              }

              // Add branch ID
              int tBranchID = treeBranchingID[node];
              branchID->InsertNextTuple1(tBranchID);

              // Add up and down nodeId
              upNodeId->InsertNextTuple1(treeSimplexId[nodeParent]);
              downNodeId->InsertNextTuple1(treeSimplexId[node]);

              // Add arc persistence
              printMsg(
                "// Push arc persistence", ttk::debug::Priority::VERBOSE);
              idNode nodeToGetPers = treeBranching[node];
              if(PlanarLayout and branchDecompositionPlanarLayout_)
                nodeToGetPers = node;
              double persToAdd
                = trees[i]->getNodePersistence<dataType>(nodeToGetPers);
              persistenceArc->InsertNextTuple1(persToAdd);

              // Add arc persistence barycenter
              if(clusteringOutput and ShiftMode != 1) {
                idNode nodeToGet = treeBranching[node];
                if(PlanarLayout and branchDecompositionPlanarLayout_)
                  nodeToGet = node;
                if(treeMatching[nodeToGet] < allBaryBranchingID[c].size())
                  persistenceBaryArc->InsertTuple1(
                    cellCount, barycenters[c]->getNodePersistence<dataType>(
                                 treeMatching[nodeToGet]));
              }

              // Add arc cluster ID
              clusterIDArc->InsertNextTuple1(clusteringAssignment[i]);

              // Add arc tree ID
              treeIDArc->InsertNextTuple1(i);

              // Add isImportantPair
              bool isImportant = false;
              idNode nodeToGetImportance = treeBranching[node];
              if(PlanarLayout and branchDecompositionPlanarLayout_)
                nodeToGetImportance = node;
              isImportant = trees[i]->isImportantPair<dataType>(
                nodeToGetImportance, importantPairs_);
              isImportantPairsArc->InsertNextTuple1(isImportant);

              // Add isDummyArc
              bool isDummy = toAdd == 2 and toAddT == 0;
              isDummyArc->InsertNextTuple1(isDummy);

              // Add isInterpolatedTree
              isInterpolatedTreeArc->InsertNextTuple1(isInterpolatedTree);

              cellCount++;
            }
          }

          // --------------
          // Node field
          // --------------
          int toAdd = (dummyNode ? 2 : 1);
          for(int toAddT = 0; toAddT < toAdd; ++toAddT) {
            // Add node id
            nodeID->InsertNextTuple1(treeSimplexId[node]);

            // Add node scalar
            scalar->InsertNextTuple1(trees[i]->getValue<dataType>(node));

            // Add criticalType
            printMsg("// Add criticalType", ttk::debug::Priority::VERBOSE);
            int criticalTypeT = -1;
            if(not isInterpolatedTree) {
              if(ShiftMode == 1) {
                if(nodeMeshTreeIndex != -1) {
                  auto array
                    = treesNodes[nodeMeshTreeIndex]->GetPointData()->GetArray(
                      "CriticalType");
                  criticalTypeT = array->GetTuple1(nodeMesh);
                }
              } else {
                auto array
                  = treesNodes[i]->GetPointData()->GetArray("CriticalType");
                criticalTypeT = array->GetTuple1(nodeMesh);
              }
            } else {
              // TODO
            }
            criticalType->InsertNextTuple1(criticalTypeT);

            // Add node matching percentage
            if(ShiftMode == 1) { // Star Barycenter
              float percentMatchT = noMatched * 100 / numInputs;
              percentMatch->InsertNextTuple1(percentMatchT);
              allBaryPercentMatch[c][node] = percentMatchT;
            }

            // Add node branch bary id
            printMsg(
              "// Add node bary branch id", ttk::debug::Priority::VERBOSE);
            if(clusteringOutput and ShiftMode != 1) {
              int tBranchID = -1;
              if(treeMatching[node] < allBaryBranchingID[c].size()) {
                tBranchID = allBaryBranchingID[c][treeMatching[node]];
                if(!trees[i]->isLeaf(node)
                   && treeMatching[nodeOrigin] < allBaryBranchingID[c].size())
                  tBranchID = allBaryBranchingID[c][treeMatching[nodeOrigin]];
              }
              branchBaryNodeID->InsertNextTuple1(tBranchID);
            }

            // Add node branch bary id
            int tBranchID = treeBranchingID[node];
            if(not trees[i]->isLeaf(node))
              tBranchID = treeBranchingID[nodeOrigin];
            branchNodeID->InsertNextTuple1(tBranchID);

            // Add node persistence
            printMsg("// Push node persistence", ttk::debug::Priority::VERBOSE);
            persistenceNode->InsertNextTuple1(
              trees[i]->getNodePersistence<dataType>(node));

            // Add node persistence barycenter
            if(clusteringOutput and ShiftMode != 1)
              if(treeMatching[node] < allBaryBranchingID[c].size())
                persistenceBaryNode->InsertTuple1(
                  pointCount, barycenters[c]->getNodePersistence<dataType>(
                                treeMatching[node]));

            // Add node clusterID
            clusterIDNode->InsertNextTuple1(clusteringAssignment[i]);

            // Add arc tree ID
            treeIDNode->InsertNextTuple1(i);

            // Add isDummyNode
            bool isDummy
              = toAdd == 2 and toAddT == 1 and !trees[i]->isRoot(node);
            isDummyNode->InsertNextTuple1(isDummy);

            // Add isInterpolatedTree
            isInterpolatedTreeNode->InsertNextTuple1(isInterpolatedTree);

            // Add isImportantPair
            bool isImportant = false;
            isImportant
              = trees[i]->isImportantPair<dataType>(node, importantPairs_);
            isImportantPairsNode->InsertNextTuple1(isImportant);

            pointCount++;
          }

          printMsg("end loop", ttk::debug::Priority::VERBOSE);
        } // end tree traversal

        // --------------
        // Manage segmentation
        // --------------
        printMsg("// Shift segmentation", ttk::debug::Priority::VERBOSE);
        if(OutputSegmentation and not PlanarLayout) {
          vtkNew<vtkUnstructuredGrid> iTreesSegmentationCopy{};
          iTreesSegmentationCopy->DeepCopy(treesSegmentation[i]);
          auto iVkOutputSegmentationTemp
            = vtkUnstructuredGrid::SafeDownCast(iTreesSegmentationCopy);
          for(int p = 0;
              p < iVkOutputSegmentationTemp->GetPoints()->GetNumberOfPoints();
              ++p) {
            double *point = iVkOutputSegmentationTemp->GetPoints()->GetPoint(p);
            point[0] += diff_x;
            point[1] += diff_y;
            point[2] += diff_z;
            iVkOutputSegmentationTemp->GetPoints()->SetPoint(p, point);
          }
          appendFilter->AddInputData(iVkOutputSegmentationTemp);
        }
      }
    }
    for(int i = persistenceBaryNode->GetNumberOfTuples(); i < pointCount; ++i)
      persistenceBaryNode->InsertNextTuple1(0);
    for(int i = persistenceBaryArc->GetNumberOfTuples(); i < cellCount; ++i)
      persistenceBaryArc->InsertNextTuple1(0);

    // --- Add VTK arrays to output
    printMsg("// Add VTK arrays to output", ttk::debug::Priority::VERBOSE);
    // Manage node output
    vtkOutputNode->SetPoints(points);
    vtkOutputNode->GetPointData()->AddArray(criticalType);
    vtkOutputNode->GetPointData()->AddArray(persistenceNode);
    vtkOutputNode->GetPointData()->AddArray(clusterIDNode);
    vtkOutputNode->GetPointData()->AddArray(treeIDNode);
    vtkOutputNode->GetPointData()->AddArray(isDummyNode);
    vtkOutputNode->GetPointData()->AddArray(branchNodeID);
    vtkOutputNode->GetPointData()->AddArray(nodeID);
    vtkOutputNode->GetPointData()->AddArray(isImportantPairsNode);
    if(not branchDecompositionPlanarLayout_)
      vtkOutputNode->GetPointData()->AddArray(scalar);
    if(clusteringOutput and ShiftMode != 1) {
      vtkOutputNode->GetPointData()->AddArray(branchBaryNodeID);
      vtkOutputNode->GetPointData()->AddArray(persistenceBaryNode);
    }
    if(foundOneInterpolatedTree)
      vtkOutputNode->GetPointData()->AddArray(isInterpolatedTreeNode);
    if(ShiftMode == 1) // Star Barycenter
      vtkOutputNode->GetPointData()->AddArray(percentMatch);

    // Manage arc output
    vtkArcs->SetPoints(points);
    vtkArcs->GetCellData()->AddArray(persistenceArc);
    vtkArcs->GetCellData()->AddArray(clusterIDArc);
    vtkArcs->GetCellData()->AddArray(treeIDArc);
    vtkArcs->GetCellData()->AddArray(isImportantPairsArc);
    vtkArcs->GetCellData()->AddArray(isDummyArc);
    vtkArcs->GetCellData()->AddArray(branchID);
    vtkArcs->GetCellData()->AddArray(upNodeId);
    vtkArcs->GetCellData()->AddArray(downNodeId);
    if(clusteringOutput and ShiftMode != 1) {
      vtkArcs->GetCellData()->AddArray(branchBaryID);
      vtkArcs->GetCellData()->AddArray(persistenceBaryArc);
    }
    if(foundOneInterpolatedTree)
      vtkArcs->GetCellData()->AddArray(isInterpolatedTreeArc);
    if(ShiftMode == 1) // Star Barycenter
      vtkArcs->GetCellData()->AddArray(percentMatchArc);
    if(vtkOutputArc == vtkOutputNode)
      vtkArcs->GetPointData()->ShallowCopy(vtkOutputNode->GetPointData());
    vtkOutputArc->ShallowCopy(vtkArcs);

    // Manage segmentation output
    if(OutputSegmentation and not PlanarLayout) {
      appendFilter->Update();
      vtkOutputSegmentation->ShallowCopy(appendFilter->GetOutput());
    }

    //
    if(ShiftMode == 1) // Star Barycenter
      trees = treesOri;
  }

  // ==========================================================================
  // Bounds Utils
  // ==========================================================================
  std::tuple<double, double, double, double, double, double>
    getRealBounds(vtkUnstructuredGrid *treeNodes,
                  FTMTree_MT *tree,
                  std::vector<int> &nodeCorrT) {
    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();
    std::queue<idNode> queue;
    queue.emplace(tree->getRoot());
    while(!queue.empty()) {
      idNode node = queue.front();
      queue.pop();
      double *point = treeNodes->GetPoints()->GetPoint(nodeCorrT[node]);
      x_min = std::min(x_min, point[0]);
      x_max = std::max(x_max, point[0]);
      y_min = std::min(y_min, point[1]);
      y_max = std::max(y_max, point[1]);
      z_min = std::min(z_min, point[2]);
      z_max = std::max(z_max, point[2]);
      std::vector<idNode> children;
      tree->getChildren(node, children);
      for(auto child : children)
        queue.emplace(child);
    }
    return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
  }

  std::tuple<double, double, double, double, double, double>
    getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree) {
    std::vector<int> nodeCorrT(tree->getNumberOfNodes());
    for(size_t i = 0; i < nodeCorrT.size(); ++i)
      nodeCorrT[i] = i;

    return getRealBounds(treeNodes, tree, nodeCorrT);
  }
};

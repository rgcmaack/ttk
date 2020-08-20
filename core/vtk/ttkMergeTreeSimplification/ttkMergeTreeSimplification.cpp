#include <ttkMergeTreeSimplification.h>

#include <vtkInformation.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMergeTreeSimplification);

ttkMergeTreeSimplification::ttkMergeTreeSimplification(){
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(2);
}

ttkMergeTreeSimplification::~ttkMergeTreeSimplification(){}

int ttkMergeTreeSimplification::FillInputPortInformation(int port, vtkInformation* info) {
    switch(port){
        case 0:
            info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
            return 1;
        case 1:
            info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
            return 1;
        default:
            return 0;
    }
}

int ttkMergeTreeSimplification::FillOutputPortInformation(int port, vtkInformation* info) {
    switch(port){
        case 0:
            info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
            return 1;
        case 1:
            info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 1);
            return 1;
        default:
            return 0;
    }
}

int ttkMergeTreeSimplification::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    ttk::Timer timer;
    this->printMsg("Retrieving Input Data",0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

    // Get the input
    auto i_MergeTree = vtkUnstructuredGrid::GetData( inputVector[0] );
    auto i_Segmentation = vtkDataSet::GetData( inputVector[1] );

    const size_t nPoints = i_MergeTree->GetNumberOfPoints();
    const size_t nEdges = i_MergeTree->GetNumberOfCells();

    // Get scalar arrays
    auto i_mtScalars = this->GetInputArrayToProcess(0, inputVector);
    auto i_mtScalarsData = (float*) ttkUtils::GetVoidPointer(i_mtScalars);

    // Get arrays with fixed name
    auto i_mtOrderFieldData = (int*) ttkUtils::GetVoidPointer( i_MergeTree->GetPointData()->GetArray("OrderField") );
    auto i_mtBranchIds = vtkIntArray::SafeDownCast(i_MergeTree->GetPointData()->GetArray("BranchId"));
    auto i_mtBranchIdsData = (int*) ttkUtils::GetVoidPointer( i_mtBranchIds );
    auto i_mtVertexIdsData = (int*) ttkUtils::GetVoidPointer( i_MergeTree->GetPointData()->GetArray("VertexId") );
    auto i_mtPairIdsData = (int*) ttkUtils::GetVoidPointer( i_MergeTree->GetPointData()->GetArray("PairId") );

    auto i_mtConnectivityData = (vtkIdType*) ttkUtils::GetVoidPointer( i_MergeTree->GetCells()->GetConnectivityArray() );
    auto i_mtOffsetsData = (vtkIdType*) ttkUtils::GetVoidPointer( i_MergeTree->GetCells()->GetOffsetsArray() );

    std::unordered_map<int,int> vertexIdToIndexMap;
    for(size_t i=0; i<nPoints; i++)
        vertexIdToIndexMap.emplace(i_mtVertexIdsData[i], i);

    auto o_mtBranchIds = vtkSmartPointer<vtkIntArray>::New();
    o_mtBranchIds->DeepCopy(i_mtBranchIds);
    auto o_mtBranchIdsData = (int*) ttkUtils::GetVoidPointer( o_mtBranchIds );

    std::vector<std::pair<int,int>> maxima;
    for(size_t i=0; i<nPoints; i++){
        o_mtBranchIdsData[i] = -1;

        // if not maximum continue
        if(i_mtBranchIdsData[i]!=i_mtVertexIdsData[i])
            continue;

        maxima.push_back({i_mtOrderFieldData[i], i});
    }
    std::sort(maxima.begin(),maxima.end());

    std::vector<int> indexToParentIndex(nPoints,-1);

    for(size_t i=0; i<nEdges; i++)
        indexToParentIndex[i_mtConnectivityData[i*2]] = i_mtConnectivityData[i*2+1];

    int nMaxima = maxima.size();

    for(int i=nMaxima-1; i>=0; i--){

        const auto& v = maxima[i].second;

        int parent = v;

        while(parent!=-1 && o_mtBranchIdsData[parent]==-1){
            o_mtBranchIdsData[parent] = v;
            parent = indexToParentIndex[parent];
        }
    }


    // std::map<int,int> parentMap


    // auto i_sSegmentationIds = vtkIntArray::SafeDownCast(i_Segmentation->GetPointData()->GetArray("SegmentationId"));

    // auto o_sBranchIds = vtkSmartPointer<vtkIntArray>::New();
    // o_sBranchIds->DeepCopy(i_sSegmentationIds);
    // o_sBranchIds->SetName("BranchId");
    // auto o_sBranchIdsData = (int*) ttkUtils::GetVoidPointer(o_sBranchIds);

    // auto o_sLeafIds = vtkSmartPointer<vtkIntArray>::New();
    // o_sLeafIds->DeepCopy(i_sSegmentationIds);
    // o_sLeafIds->SetName("LeafId");
    // auto o_sLeafIdsData = (int*) ttkUtils::GetVoidPointer(o_sLeafIds);


    // this->printMsg("Retrieving Input Data",1,timer.getElapsedTime());
    // timer.reStart();

    // this->printMsg("Simplifying Merge Tree",0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);
    // const size_t nPoints = i_MergeTree->GetNumberOfPoints();


    // // sort pairs
    // std::unordered_map<int,int> simplifiedBranchToPersistentBranchMap;

    // std::vector<std::pair<float,int>> pairs;
    // for(size_t i=0; i<nPoints; i++){

    //     // if not maximum continue
    //     if(i_mtBranchIdsData[i]!=i_mtVertexIdsData[i])
    //         continue;

    //     // at first each branch maps to itself
    //     simplifiedBranchToPersistentBranchMap.emplace(i_mtBranchIdsData[i],i_mtBranchIdsData[i]);

    //     // add persistence-maxima pair
    //     pairs.push_back({
    //         i_mtScalarsData[i] - i_mtScalarsData[vertexIdToIndexMap[ i_mtPairIdsData[i] ]],
    //         i
    //     });
    // }
    // std::sort(pairs.begin(),pairs.end());
    // const int nPairs = pairs.size();

    // float persistenceThreshold = (float) this->PersistenceThreshold;
    // const int globalMaxIndex = pairs[nPairs-1].second;
    // const int globalBranchId = i_mtBranchIdsData[globalMaxIndex];

    // for(int i = nPairs-2; i>=0; i--){ // -2: largest branch can never be simplified
    //     const int maxIndex = pairs[i].second;
    //     const int saddleId = i_mtPairIdsData[maxIndex];
    //     const int saddleIndex = vertexIdToIndexMap[ saddleId ];

    //     const float persistence = i_mtScalarsData[ maxIndex ] - i_mtScalarsData[ saddleIndex ];

    //     if(persistence<=persistenceThreshold){
    //         auto it = simplifiedBranchToPersistentBranchMap.find(i_mtBranchIdsData[maxIndex]);
    //         if(persistence<=0){
    //             it->second = globalBranchId;
    //         } else {
    //             it->second = simplifiedBranchToPersistentBranchMap[i_mtBranchIdsData[ saddleIndex ]];
    //         }
    //     }
    // }

    // for(size_t i=0; i<nPoints; i++){
    //     o_mtBranchIdsData[i] = simplifiedBranchToPersistentBranchMap[ i_mtBranchIdsData[i] ];
    // }

    // const size_t nSegmentationPoints = i_Segmentation->GetNumberOfPoints();
    // for(size_t i=0; i<nSegmentationPoints; i++){
    //     o_sBranchIdsData[i] = o_mtBranchIdsData[
    //         vertexIdToIndexMap[o_sBranchIdsData[i]]
    //     ];
    // }

    // Get output arrays
    auto o_MergeTree = vtkUnstructuredGrid::GetData( outputVector, 0 );
    // auto o_Segmentation = vtkDataSet::GetData( outputVector, 1 );

    o_MergeTree->ShallowCopy( i_MergeTree );
    o_MergeTree->GetPointData()->AddArray( o_mtBranchIds );

    // o_Segmentation->ShallowCopy( i_Segmentation );
    // o_Segmentation->GetPointData()->AddArray( o_sBranchIds );

    // this->printMsg("Simplifying Merge Tree",1,timer.getElapsedTime());

    return 1;
}
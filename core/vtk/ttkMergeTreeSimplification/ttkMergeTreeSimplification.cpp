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
#include <vtkUnsignedCharArray.h>

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

    const size_t nMTPoints = i_MergeTree->GetNumberOfPoints();
    const size_t nSPoints = i_Segmentation->GetNumberOfPoints();
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

    auto i_sSegmentationIds = vtkIntArray::SafeDownCast(i_Segmentation->GetPointData()->GetArray("SegmentationId"));

    this->printMsg("Retrieving Input Data",1,timer.getElapsedTime(),this->threadNumber_);
    timer.reStart();

    this->printMsg("Simplifying Merge Tree",0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

    std::unordered_map<int,int> vertexIdToIndexMap;
    for(size_t i=0; i<nMTPoints; i++)
        vertexIdToIndexMap.emplace(i_mtVertexIdsData[i], i);

    auto o_mtBranchIds = vtkSmartPointer<vtkIntArray>::New();
    o_mtBranchIds->DeepCopy(i_mtBranchIds);
    auto o_mtBranchIdsData = (int*) ttkUtils::GetVoidPointer( o_mtBranchIds );

    auto o_mtRegionType = vtkSmartPointer<vtkUnsignedCharArray>::New();
    o_mtRegionType->SetName("RegionType");
    o_mtRegionType->SetNumberOfTuples(nMTPoints);
    auto o_mtRegionTypeData = (unsigned char*) ttkUtils::GetVoidPointer( o_mtRegionType );

    auto o_mtPointType = vtkSmartPointer<vtkUnsignedCharArray>::New();
    o_mtPointType->SetName("PointType");
    o_mtPointType->SetNumberOfTuples(nMTPoints);
    auto o_mtPointTypeData = (unsigned char*) ttkUtils::GetVoidPointer( o_mtPointType );
    /*
    0 regular
    1 maximum
    2 saddle
    3 simplified
    */

    std::vector<std::pair<int,int>> maxima;
    for(size_t i=0; i<nMTPoints; i++){
        o_mtBranchIdsData[i] = -1;
        o_mtRegionTypeData[i] = 0;

        if(i_mtBranchIdsData[i]==i_mtVertexIdsData[i]){
            maxima.push_back({i_mtOrderFieldData[i], i});
            o_mtPointTypeData[i] = 1;
        } else {
            o_mtPointTypeData[i] = 2;
        }
    }
    std::sort(maxima.begin(),maxima.end());

    std::vector<int> indexToParentIndex(nMTPoints,-1);
    std::vector<int> childrenCounter(nMTPoints,0);

    for(size_t i=0; i<nEdges; i++){
        indexToParentIndex[i_mtConnectivityData[i*2]] = i_mtConnectivityData[i*2+1];

        childrenCounter[i_mtConnectivityData[i*2+1]]++;
    }

    // for(size_t i=0; i<nMTPoints; i++){
    //     if(childrenCounter[i]>2){
    //         this->printWrn(std::to_string(childrenCounter[i]));
    //     }
    // }

    int nMaxima = maxima.size();

    // compute branch decomposition for all arcs
    for(int i=nMaxima-1; i>=0; i--){
        const auto& vIndex = maxima[i].second;
        const auto& vBranchId = i_mtVertexIdsData[vIndex];

        int cIndex = vIndex;
        while(cIndex!=-1 && o_mtBranchIdsData[cIndex]==-1){
            o_mtBranchIdsData[cIndex] = vBranchId;
            cIndex = indexToParentIndex[cIndex];
        }
    }

    // override branch id for all non-persistent arcs
    float persistenceThreshold = (float) this->PersistenceThreshold;
    for(int i=nMaxima-2; i>=0; i--){ // -2: main branch is never simplified
        const auto& vIndex = maxima[i].second;
        const auto& uId = i_mtPairIdsData[vIndex];
        const auto& uIndex = vertexIdToIndexMap[uId];
        const auto& uBranchId = o_mtBranchIdsData[uIndex];

        const float persistence = i_mtScalarsData[ vIndex ] - i_mtScalarsData[ uIndex ];

        if(persistence<=persistenceThreshold){
            int cIndex = vIndex;

            while(o_mtBranchIdsData[cIndex]!=uBranchId){
                o_mtBranchIdsData[cIndex] = uBranchId;
                o_mtPointTypeData[cIndex] = 3;
                cIndex = indexToParentIndex[cIndex];
            }

            childrenCounter[cIndex]--;
            if(childrenCounter[cIndex]<2)
                o_mtPointTypeData[cIndex] = 3;
        }
    }

    // determine type
    for(int i=nMaxima-1; i>=0; i--){ // -2: main branch is never simplified
        const auto& vIndex = maxima[i].second;
        const auto& vBranchId = o_mtBranchIdsData[vIndex];

        int cIndex = vIndex;
        if(vBranchId==i_mtVertexIdsData[vIndex]){
            // if maxima was not simplified then mark region until first true saddle
            while(cIndex!=-1 && o_mtPointTypeData[cIndex]!=2){
                o_mtRegionTypeData[cIndex] = 1;
                cIndex = indexToParentIndex[cIndex];
            }
        } else {
            // otherwise check if simplified maxima reaches leaf region
            bool reachedLeafRegion = false;
            while(cIndex!=-1 && o_mtPointTypeData[cIndex]!=2){
                if(o_mtRegionTypeData[cIndex]==1){
                    reachedLeafRegion = true;
                    break;
                }
                cIndex = indexToParentIndex[cIndex];
            }

            cIndex = vIndex;
            const signed char label = reachedLeafRegion ? 1 : 0;
            while(cIndex!=-1 && o_mtPointTypeData[cIndex]!=2){
                o_mtRegionTypeData[cIndex] = label;
                cIndex = indexToParentIndex[cIndex];
            }
        }
    }

    this->printMsg("Simplifying Merge Tree",1,timer.getElapsedTime(),this->threadNumber_);
    timer.reStart();

    this->printMsg("Mapping Segmentation",0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

    // initilize output segmentation arrays
    auto o_sBranchIds = vtkSmartPointer<vtkIntArray>::New();
    o_sBranchIds->DeepCopy(i_sSegmentationIds);
    o_sBranchIds->SetName("BranchId");
    auto o_sBranchIdsData = (int*) ttkUtils::GetVoidPointer(o_sBranchIds);

    auto o_sRegionId = vtkSmartPointer<vtkIntArray>::New();
    o_sRegionId->DeepCopy(i_sSegmentationIds);
    o_sRegionId->SetName("RegionId");
    auto o_sRegionIdData = (int*) ttkUtils::GetVoidPointer(o_sRegionId);

    // map segmentation to simplification
    for(size_t i=0; i<nSPoints; i++){
        const auto& vIndex = vertexIdToIndexMap[ o_sBranchIdsData[i] ];
        o_sBranchIdsData[i] = o_mtBranchIdsData[ vIndex ];
        o_sRegionIdData[i] = o_mtRegionTypeData[ vIndex ]==1 ? o_mtBranchIdsData[ vIndex ] : -1;
    }

    // Get output arrays
    auto o_MergeTree = vtkUnstructuredGrid::GetData( outputVector, 0 );
    auto o_Segmentation = vtkDataSet::GetData( outputVector, 1 );

    o_MergeTree->ShallowCopy( i_MergeTree );
    o_MergeTree->GetPointData()->AddArray( o_mtBranchIds );
    o_MergeTree->GetPointData()->AddArray( o_mtRegionType );
    o_MergeTree->GetPointData()->AddArray( o_mtPointType );

    o_Segmentation->ShallowCopy( i_Segmentation );
    o_Segmentation->GetPointData()->AddArray( o_sBranchIds );
    o_Segmentation->GetPointData()->AddArray( o_sRegionId );

    this->printMsg("Mapping Segmentation",1,timer.getElapsedTime(),this->threadNumber_);

    return 1;
}
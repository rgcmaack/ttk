#include <ttkMergeTree.h>

#include <vtkInformation.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMergeTree);

ttkMergeTree::ttkMergeTree(){
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(3);
}

ttkMergeTree::~ttkMergeTree(){}

int ttkMergeTree::FillInputPortInformation(int port, vtkInformation* info) {
    if (port==0){
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
        return 1;
    }
    return 0;
}

int ttkMergeTree::FillOutputPortInformation(int port, vtkInformation* info) {
    switch(port){
        case 0:
            info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
            return 1;
        case 1:
            info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
            return 1;
        case 2:
            info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
            return 1;
        default:
            return 0;
    }
}

template<typename T>
int fillPointData(
    std::vector<ttk::Propagation<int>>& propagations,
    float* pointCoordData,
    int* vertexIdData,
    int* branchIdData,
    int* mtOrderFieldData,
    int* orderFieldData,

    vtkDataArray* outputScalars,
    vtkDataArray* inputScalars,
    std::unordered_map<int,std::pair<int,int>>& indexToIdMap,
    const ttk::Triangulation* triangulation
){
    auto inputScalarData = (T*) ttkUtils::GetVoidPointer( inputScalars );
    auto outputScalarData = (T*) ttkUtils::GetVoidPointer( outputScalars );

    for( const auto& it : indexToIdMap ) {
        const auto& v = it.first;

        mtOrderFieldData[it.second.first] = orderFieldData[v];
        vertexIdData[it.second.first] = v;
        outputScalarData[it.second.first] = inputScalarData[v];
        branchIdData[it.second.first] = propagations[it.second.second].criticalPoints[0];

        int offset = it.second.first*3;
        triangulation->getVertexPoint(
            v,
            pointCoordData[offset],
            pointCoordData[offset+1],
            pointCoordData[offset+2]
        );
    }

    return 1;
}

int ttkMergeTree::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    // Get the input
    auto input = vtkDataSet::GetData( inputVector[0] );
    size_t nVertices = input->GetNumberOfPoints();

    // Get triangulation of the input object
    auto triangulation = ttkAlgorithm::GetTriangulation( input );

    // Precondition triangulation
    this->PreconditionTriangulation( triangulation );

    // Get input array
    auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

    // Init order field
    auto orderField = vtkSmartPointer<vtkIntArray>::New();
    orderField->SetName( "OrderField" );
    orderField->SetNumberOfComponents(1);
    orderField->SetNumberOfTuples( nVertices );
    auto orderFieldData = (int*) ttkUtils::GetVoidPointer(orderField);

    // Compute order field
    {
        int status = -1;
        switch(inputScalars->GetDataType()) {
            vtkTemplateMacro(
                (status = this->computeOrderField<VTK_TT, int>(
                    orderFieldData,

                    nVertices,
                    (VTK_TT*) ttkUtils::GetVoidPointer(inputScalars)
                ))
            );
        }
        if(!status)
            return 0;
    }

    // Init segmentation
    auto segmentationIds = vtkSmartPointer<vtkIntArray>::New();
    segmentationIds->SetName( "SegmentationId" );
    segmentationIds->SetNumberOfComponents(1);
    segmentationIds->SetNumberOfTuples( nVertices );

    // Compute merge tree segmentation
    std::vector<ttk::Propagation<int>> propagations;
    {
        int status = 0;
        ttkVtkTemplateMacro(
            triangulation->getType(),
            orderField->GetDataType(),
            (
                status = this->computeMergeTreeSegmentation<int, TTK_TT>(
                    (int*) ttkUtils::GetVoidPointer(segmentationIds),
                    propagations,

                    (TTK_TT *) triangulation->getData(),
                    orderFieldData
                )
            )
        );
        if(!status)
            return 0;
    }

    // Finalize Output
    {
        ttk::Timer timer;
        this->printMsg(
            "Finalizing Output",0,0,
            ttk::debug::LineMode::REPLACE
        );

        // Create segmentation output
        {
            auto segmentation = vtkDataSet::GetData( outputVector, 2 );
            segmentation->ShallowCopy( input );

            auto segmentationPD = segmentation->GetPointData();
            segmentationPD->AddArray( orderField );
            segmentationPD->AddArray( segmentationIds );
        }

        // Compute merge tree output
        {
            const int nPropagations = propagations.size();

            std::unordered_map<int,std::pair<int,int>> indexToIdMap;
            for(int p=0, i=0; p<nPropagations; p++){
                const auto& criticalPoints = propagations[p].criticalPoints;
                const int nCriticalPoints = criticalPoints.size();

                for(int c=0; c<nCriticalPoints; c++){
                    auto it = indexToIdMap.find( criticalPoints[c] );
                    if(it==indexToIdMap.end())
                        indexToIdMap.insert({criticalPoints[c], {i++, p}});
                    else if(
                        orderFieldData[propagations[it->second.second].criticalPoints[0]]
                        <
                        orderFieldData[criticalPoints[0]]
                    )
                        it->second.second = p;
                }
            }

            const int nPoints = indexToIdMap.size();

            auto pointCoords = vtkSmartPointer<vtkFloatArray>::New();
            pointCoords->SetNumberOfComponents(3);
            pointCoords->SetNumberOfTuples(nPoints);
            auto pointCoordData = (float*) ttkUtils::GetVoidPointer( pointCoords );

            auto vertexIds = vtkSmartPointer<vtkIntArray>::New();
            vertexIds->SetName("VertexId");
            vertexIds->SetNumberOfTuples(nPoints);
            auto vertexIdData = (int*) ttkUtils::GetVoidPointer( vertexIds );

            auto branchIds = vtkSmartPointer<vtkIntArray>::New();
            branchIds->SetName( "BranchId" );
            branchIds->SetNumberOfTuples(nPoints);
            auto branchIdData = (int*) ttkUtils::GetVoidPointer( branchIds );

            auto outputScalars = vtkSmartPointer<vtkDataArray>::Take( inputScalars->NewInstance() );
            outputScalars->SetName( inputScalars->GetName() );
            outputScalars->SetNumberOfTuples(nPoints);

            auto mtOrderField = vtkSmartPointer<vtkIntArray>::New();
            mtOrderField->SetName( "OrderField" );
            mtOrderField->SetNumberOfTuples(nPoints);
            auto mtOrderFieldData = (int*) ttkUtils::GetVoidPointer( mtOrderField );

            auto pairId = vtkSmartPointer<vtkIntArray>::New();
            pairId->SetName( "PairId" );
            pairId->SetNumberOfTuples(nPoints);
            auto pairIdData = (int*) ttkUtils::GetVoidPointer( pairId );

            switch( inputScalars->GetDataType() ){
                vtkTemplateMacro(
                    (
                        fillPointData<VTK_TT>(
                            propagations,
                            pointCoordData,
                            vertexIdData,
                            branchIdData,
                            mtOrderFieldData,
                            orderFieldData,

                            outputScalars,
                            inputScalars,
                            indexToIdMap,
                            triangulation
                        )
                    )
                );
            }

            // critical points
            {
                auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
                offsets->SetNumberOfTuples( nPoints +1 );
                auto offsetsData = (vtkIdType*) ttkUtils::GetVoidPointer( offsets );
                for(int i=0; i<nPoints; i++)
                    offsetsData[i] = i;
                offsetsData[nPoints] = nPoints;

                auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
                connectivity->SetNumberOfTuples( nPoints );
                auto connectivityData = (vtkIdType*) ttkUtils::GetVoidPointer( connectivity );
                for(int i=0; i<nPoints; i++)
                    connectivityData[i] = i;

                auto cells = vtkSmartPointer<vtkCellArray>::New();
                cells->SetData(offsets, connectivity);

                auto points = vtkSmartPointer<vtkPoints>::New();
                points->SetData( pointCoords );

                auto criticalPoints = vtkUnstructuredGrid::GetData(outputVector, 0);
                criticalPoints->SetPoints(points);
                criticalPoints->SetCells(VTK_VERTEX, cells);

                auto criticalPointsPD = criticalPoints->GetPointData();
                criticalPointsPD->AddArray( vertexIds );
                criticalPointsPD->AddArray( branchIds );
                criticalPointsPD->AddArray( mtOrderField );
                criticalPointsPD->AddArray( outputScalars );
                criticalPointsPD->AddArray( pairId );
            }

            // edges
            {
                const int nEdges = nPoints - 1;

                auto cellBranchIds = vtkSmartPointer<vtkIntArray>::New();
                cellBranchIds->SetName( "BranchId" );
                cellBranchIds->SetNumberOfTuples(nEdges);
                auto cellBranchIdsData = (int*) ttkUtils::GetVoidPointer( cellBranchIds );

                auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
                offsets->SetNumberOfTuples( nEdges +1 );
                auto offsetsData = (vtkIdType*) ttkUtils::GetVoidPointer( offsets );
                for(int i=0; i<nEdges; i++)
                    offsetsData[i] = i*2;
                offsetsData[nEdges] = nEdges*2;

                auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
                connectivity->SetNumberOfTuples( nEdges*2 );
                auto connectivityData = (vtkIdType*) ttkUtils::GetVoidPointer( connectivity );

                for(int i=0, j=0, k=0; i<nPropagations; i++){
                    const auto& criticalPoints = propagations[i].criticalPoints;
                    const int nCriticalPoints = criticalPoints.size();

                    const auto& s0 = criticalPoints[0];
                    const auto& s1 = criticalPoints[nCriticalPoints-1];

                    pairIdData[ indexToIdMap[ s0 ].first ] = s1;
                    pairIdData[ indexToIdMap[ s1 ].first ] = s0;

                    for(int c=0; c<nCriticalPoints-1; c++){
                        cellBranchIdsData[k++] = propagations[indexToIdMap[criticalPoints[c]].second].criticalPoints[0];
                        connectivityData[j++] = indexToIdMap[ criticalPoints[c  ] ].first;
                        connectivityData[j++] = indexToIdMap[ criticalPoints[c+1] ].first;
                    }
                }

                auto cells = vtkSmartPointer<vtkCellArray>::New();
                cells->SetData(offsets, connectivity);

                auto points = vtkSmartPointer<vtkPoints>::New();
                points->SetData( pointCoords );

                auto edges = vtkUnstructuredGrid::GetData(outputVector, 1);
                edges->SetPoints(points);
                edges->SetCells(VTK_LINE, cells);

                auto edgesPD = edges->GetPointData();
                edgesPD->AddArray( vertexIds );
                edgesPD->AddArray( branchIds );
                edgesPD->AddArray( mtOrderField );
                edgesPD->AddArray( outputScalars );
                edgesPD->AddArray( pairId );

                edges->GetCellData()->AddArray( cellBranchIds );
            }
        }

        this->printMsg(
            "Finalizing Output",1,timer.getElapsedTime()
        );
        this->printMsg(ttk::debug::Separator::L1);
    }


    return 1;
}
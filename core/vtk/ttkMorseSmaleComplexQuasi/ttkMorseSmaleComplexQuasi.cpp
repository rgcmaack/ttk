#include <ttkMorseSmaleComplexQuasi.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>

//unused?
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkMorseSmaleComplexQuasi);

ttkMorseSmaleComplexQuasi::ttkMorseSmaleComplexQuasi() {
  this->setDebugMsgPrefix("MorseSmaleComplexQuasi");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(4);
}

ttkMorseSmaleComplexQuasi::~ttkMorseSmaleComplexQuasi() {
}

int ttkMorseSmaleComplexQuasi::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMorseSmaleComplexQuasi::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  } else if(port == 3) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkMorseSmaleComplexQuasi::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkPolyData::GetData(outputVector, 0);
  auto outputSeparatrices1 = vtkPolyData::GetData(outputVector, 1);
  auto outputSeparatrices2 = vtkPolyData::GetData(outputVector, 2);
  auto outputMorseComplexes = vtkDataSet::GetData(outputVector, 3);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(input->GetNumberOfPoints() == 0) {
    this->printErr("Input has no point.");
    return -1;
  }
  if(!outputCriticalPoints or !outputSeparatrices1 or !outputSeparatrices2
     or !outputMorseComplexes) {
    this->printErr("Output pointers are NULL.");
    return -1;
  }
#endif

  vtkDataArray *inputArray = this->GetOrderArray(input, 0);

#ifndef TTK_ENABLE_KAMIKAZE  
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }
#endif

  this->setInputOrderField(ttkUtils::GetVoidPointer(inputArray));

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE    
  if(triangulation == nullptr) {
    this->printErr("Triangulation is null");
    return 0;
  }
#endif

  this->preconditionTriangulation(triangulation);
  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();

#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfVertices) {
    this->printErr("Input has no vertices.");
    return -1;
  }
#endif

  vtkNew<ttkSimplexIdTypeArray> ascendingManifold{};
  vtkNew<ttkSimplexIdTypeArray> descendingManifold{};
  vtkNew<ttkSimplexIdTypeArray> morseSmaleManifold{};

#ifndef TTK_ENABLE_KAMIKAZE
  if(!ascendingManifold || !descendingManifold || !morseSmaleManifold) {
    this->printErr("Manifold vtkDataArray allocation problem.");
    return -1;
  }
#endif

  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(numberOfVertices);
  ascendingManifold->SetName("AscendingManifold");

  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(numberOfVertices);
  descendingManifold->SetName("DescendingManifold");

  morseSmaleManifold->SetNumberOfComponents(1);
  morseSmaleManifold->SetNumberOfTuples(numberOfVertices);
  morseSmaleManifold->SetName("MorseSmaleManifold");

  void *ascendingManifoldPtr = nullptr;
  void *descendingManifoldPtr = nullptr;
  void *morseSmaleManifoldPtr = nullptr;
  if(ComputeAscendingSegmentation)
    ascendingManifoldPtr = ttkUtils::GetVoidPointer(ascendingManifold);
  if(ComputeDescendingSegmentation)
    descendingManifoldPtr = ttkUtils::GetVoidPointer(descendingManifold);
  if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
     and ComputeFinalSegmentation)
    morseSmaleManifoldPtr = ttkUtils::GetVoidPointer(morseSmaleManifold);

  this->setOutputMorseComplexes(
    ascendingManifoldPtr, descendingManifoldPtr, morseSmaleManifoldPtr);

  // If all checks pass then log which array is going to be processed.
  this->printMsg("MSC Order computation starting...");
  this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

  int status = 0;

  ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
                     (status = dispatch<VTK_TT, TTK_TT>(
                      inputArray, outputCriticalPoints,
                      outputSeparatrices1,
                      outputSeparatrices2,
                      *static_cast<TTK_TT *>(triangulation->getData())
                      )));

  if(status != 0) {
    return 0;
  }
  
  outputMorseComplexes->ShallowCopy(input);
  // morse complexes
  if(ComputeAscendingSegmentation or ComputeDescendingSegmentation) {
    vtkPointData *pointData = outputMorseComplexes->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("Segmentation output has no point data.");
      return -1;
    }
#endif

    if(ComputeDescendingSegmentation)
      pointData->AddArray(descendingManifold);
    if(ComputeAscendingSegmentation)
      pointData->AddArray(ascendingManifold);
    if(ComputeAscendingSegmentation &&
      ComputeDescendingSegmentation &&
      ComputeFinalSegmentation)
      pointData->AddArray(morseSmaleManifold);
  }

  // return success
  return 1;
}

template <typename vtkArrayType, typename vectorType>
void setArray(vtkArrayType &vtkArray, vectorType &vector) {
  ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
}

template <typename dataType, typename triangulationType>
int ttkMorseSmaleComplexQuasi::dispatch(vtkDataArray *const inputArray,
                                        vtkPolyData *const outputCriticalPoints,
                                        vtkPolyData *const outputSeparatrices1,
                                        vtkPolyData *const outputSeparatrices2,
                                        const triangulationType &triangulation)
                                        {
  const auto scalars
    = static_cast<dataType *>(ttkUtils::GetVoidPointer(inputArray));
  
  // critical points
  criticalPoints_points.clear();
  criticalPoints_points_cellDimensions.clear();
  criticalPoints_points_cellIds.clear();

  // 1-separatrices
  SimplexId s1_numberOfPoints{};
  SimplexId s1_numberOfCells{};
  separatrices1_points.clear();
  separatrices1_points_smoothingMask.clear();
  separatrices1_points_cellDimensions.clear();
  separatrices1_points_cellIds.clear();
  separatrices1_cells_connectivity.clear();
  separatrices1_cells_sourceIds.clear();
  separatrices1_cells_destinationIds.clear();
  separatrices1_cells_separatrixIds.clear();
  separatrices1_cells_separatrixTypes.clear();
  separatrices1_cells_isOnBoundary.clear();
  std::vector<ttk::SimplexId> s1_separatrixFunctionMaxima{};
  std::vector<ttk::SimplexId> s1_separatrixFunctionMinima{};

  // 2-separatrices
  SimplexId s2_numberOfPoints{};
  SimplexId s2_numberOfCells{};
  separatrices2_points.clear();
  separatrices2_cells_connectivity.clear();

  this->setOutputCriticalPoints(
    &criticalPoints_points, &criticalPoints_points_cellDimensions,
    &criticalPoints_points_cellIds);

  this->setOutputSeparatrices1(
    &s1_numberOfPoints, &separatrices1_points,
    &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
    &separatrices1_points_cellIds, &s1_numberOfCells,
    &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
    &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
    &separatrices1_cells_separatrixTypes, &s1_separatrixFunctionMaxima,
    &s1_separatrixFunctionMinima, &separatrices1_cells_isOnBoundary);

  this->setOutputSeparatrices2(
    &s2_numberOfPoints, &separatrices2_points,
    &s2_numberOfCells, &separatrices2_cells_connectivity,
    &separatrices2_cells_caseTypes);

  const int ret = this->execute<dataType, triangulationType>(
    triangulation, SeparatricesManifold, ComputeSeparatrices, SeperateExtrema,
    FastSeparatrices
  );

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    this->printErr("MorseSmaleComplexOrder.execute() error");
    return -1;
  }
#endif

  // critical points
  {
    vtkNew<vtkPoints> points{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkSmartPointer<vtkDataArray> cellScalars{inputArray->NewInstance()};
    const auto nPoints = criticalPoints_points.size();

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !cellDimensions || !cellIds || !cellScalars) {
      this->printErr("Critical points vtkDataArray allocation problem.");
      return -1;
    }
#endif

    points->SetNumberOfPoints(nPoints);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");
    setArray(cellDimensions, criticalPoints_points_cellDimensions);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");
    setArray(cellIds, criticalPoints_points_cellIds);

    cellScalars->SetNumberOfComponents(1);
    cellScalars->SetName(inputArray->GetName());
    cellScalars->SetNumberOfTuples(nPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nPoints; ++i) {
      points->SetPoint(i, criticalPoints_points[i].data());
      cellScalars->SetTuple1(i, scalars[criticalPoints_points_cellIds[i]]);
    }

    outputCriticalPoints->SetPoints(points);

    vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(nPoints + 1);
    connectivity->SetNumberOfComponents(1);
    connectivity->SetNumberOfTuples(nPoints);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nPoints; ++i) {
      offsets->SetTuple1(i, i);
      connectivity->SetTuple1(i, i);
    }
    offsets->SetTuple1(nPoints, nPoints);
    vtkNew<vtkCellArray> cells{};
    cells->SetData(offsets, connectivity);
    outputCriticalPoints->SetVerts(cells);

    auto pointData = outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("outputCriticalPoints has no point data.");
      return -1;
    }
#endif

    pointData->AddArray(cellDimensions);
    pointData->AddArray(cellIds);
    pointData->AddArray(cellScalars);
  }

  // 1-separatrices
  {
    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<vtkSignedCharArray> smoothingMask{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> destinationIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkDoubleArray> separatrixFunctionMaxima{};
    vtkNew<vtkDoubleArray> separatrixFunctionMinima{};
    vtkNew<vtkDoubleArray> separatrixFunctionDiffs{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !smoothingMask || !cellDimensions || !cellIds
       || !sourceIds || !destinationIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("1-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices1_points);

    smoothingMask->SetNumberOfComponents(1);
    smoothingMask->SetName(ttk::MaskScalarFieldName);
    setArray(smoothingMask, separatrices1_points_smoothingMask);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");
    setArray(cellDimensions, separatrices1_points_cellDimensions);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");
    setArray(cellIds, separatrices1_points_cellIds);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");
    setArray(sourceIds, separatrices1_cells_sourceIds);

    destinationIds->SetNumberOfComponents(1);
    destinationIds->SetName("DestinationId");
    setArray(destinationIds, separatrices1_cells_destinationIds);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");
    setArray(separatrixIds, separatrices1_cells_separatrixIds);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");
    setArray(separatrixTypes, separatrices1_cells_separatrixTypes);

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");
    setArray(isOnBoundary, separatrices1_cells_isOnBoundary);

    vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(s1_numberOfCells + 1);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices1_cells_connectivity);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < s1_numberOfCells + 1; ++i) {
      offsets->SetTuple1(i, 2 * i);
    }

    vtkNew<vtkPoints> points{};
    points->SetData(pointsCoords);
    outputSeparatrices1->SetPoints(points);
    vtkNew<vtkCellArray> cells{};
#ifndef TTK_ENABLE_64BIT_IDS
    cells->Use32BitStorage();
#endif // TTK_ENABLE_64BIT_IDS
    cells->SetData(offsets, connectivity);
    outputSeparatrices1->SetLines(cells);

    auto pointData = outputSeparatrices1->GetPointData();
    auto cellData = outputSeparatrices1->GetCellData();

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData || !cellData) {
      this->printErr("outputSeparatrices1 has no point or no cell data.");
      return -1;
    }
#endif
  }

  // 2-separatrices
  {
    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<ttkSimplexIdTypeArray> caseTypes{};
    
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !caseTypes) {
      this->printErr("2-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices2_points);

    vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(s2_numberOfCells + 1);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices2_cells_connectivity);

    caseTypes->SetNumberOfComponents(1);
    caseTypes->SetName("CaseType");
    setArray(caseTypes, separatrices2_cells_caseTypes);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < s2_numberOfCells + 1; ++i) {
      offsets->SetTuple1(i, 3 * i);
    }

    vtkNew<vtkPoints> points{};
    points->SetData(pointsCoords);
    outputSeparatrices2->SetPoints(points);
    vtkNew<vtkCellArray> cells{};
#ifndef TTK_ENABLE_64BIT_IDS
    cells->Use32BitStorage();
#endif // TTK_ENABLE_64BIT_IDS
    cells->SetData(offsets, connectivity);
    outputSeparatrices2->SetPolys(cells);

    auto cellData = outputSeparatrices2->GetCellData();
    cellData->AddArray(caseTypes);
  }

  return 0;
}
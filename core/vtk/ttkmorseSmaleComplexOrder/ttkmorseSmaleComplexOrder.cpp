#include <ttkmorseSmaleComplexOrder.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>

vtkStandardNewMacro(ttkmorseSmaleComplexOrder);

ttkmorseSmaleComplexOrder::ttkmorseSmaleComplexOrder() {
  this->setDebugMsgPrefix("OrderMorseSmaleComplex");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(4);
}

ttkmorseSmaleComplexOrder::~ttkmorseSmaleComplexOrder() {
}

int ttkmorseSmaleComplexOrder::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkmorseSmaleComplexOrder::FillOutputPortInformation(int port,
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

int ttkmorseSmaleComplexOrder::RequestData(vtkInformation *request,
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

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);

#ifndef TTK_ENABLE_KAMIKAZE    
  if(triangulation == nullptr) {
    this->printErr("Triangulation is null");
    return 0;
  }
#endif

  this->preconditionTriangulation(triangulation);

  // If all checks pass then log which array is going to be processed.
  this->printMsg("MSC Order computation starting...");
  this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

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

  // Create an output array that has the same data type as the input array
  vtkSmartPointer<vtkDataArray> outputDescArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputDescArray->SetName("Descending Manifold"); // set array name
  outputDescArray->SetNumberOfComponents(1); // only one component per tuple
  outputDescArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  vtkSmartPointer<vtkDataArray> outputAscArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputAscArray->SetName("Ascending Manifold"); // set array name
  outputAscArray->SetNumberOfComponents(1); // only one component per tuple
  outputAscArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  this->setInputOrderField(ttkUtils::GetVoidPointer(inputArray));

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
      this->printErr("outputMorseComplexes has no point data.");
      return -1;
    }
#endif

    if(ComputeDescendingSegmentation)
      pointData->AddArray(descendingManifold);
    if(ComputeAscendingSegmentation)
      pointData->AddArray(ascendingManifold);
    if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
       and ComputeFinalSegmentation)
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
int ttkmorseSmaleComplexOrder::dispatch(vtkDataArray *const inputArray,
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

  if(ComputeCriticalPoints) {
    this->setOutputCriticalPoints(
      &criticalPoints_points, &criticalPoints_points_cellDimensions,
      &criticalPoints_points_cellIds);
  } else {
    this->setOutputCriticalPoints(
      nullptr, nullptr, nullptr);
  }

  const int ret = this->execute<dataType>(triangulation);

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
      cellScalars->SetTuple1(i, scalars[i]);
    }

    outputCriticalPoints->SetPoints(points);

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
  return 0;
}

// Templatize over the different input array data types and call the base code
  /*int status = 0; // this integer checks if the base code returns an error
  ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
                      (status = this->computeDescendingManifold<VTK_TT, TTK_TT>(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputDescArray),
                         (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         (TTK_TT *)triangulation->getData())));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
                      (status = this->computeAscendingManifold<VTK_TT, TTK_TT>(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputAscArray),
                         (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         (TTK_TT *)triangulation->getData())));
                        
  if(status != 1)
    return 0;
    
  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(input);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(outputDescArray);
  outputDataSet->GetPointData()->AddArray(outputAscArray);
  */
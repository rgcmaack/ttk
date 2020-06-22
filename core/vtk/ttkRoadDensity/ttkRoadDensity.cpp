#include <ttkRoadDensity.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <sstream>
#include <string>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <Triangulation.h>
#include <vtkAlgorithm.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

using namespace std;

vtkStandardNewMacro(ttkRoadDensity);

ttkRoadDensity::ttkRoadDensity() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkRoadDensity::~ttkRoadDensity() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkRoadDensity::FillInputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(
        vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
      break;
    // case 1:
    //   info->Set(
    //     vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    //   break;
    default:
      return 0;
  }
  return 1;
}

// see ttkAlgorithm::FillOutputPortInformation for details about this method
int ttkRoadDensity::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    default:
      return 0;
  }
  return 1;
}

int ttkRoadDensity::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  // get parameter
  std::string KernelFunction = this->GetKernelFunction();
  float KernelBandwidth = this->GetKernelBandwidth();

  // Get the input
  auto inputRefinedRoad = vtkUnstructuredGrid::GetData(inputVector[0]);
  // auto inputEventwithClosetPoint =
  // vtkUnstructuredGrid::GetData(inputVector[1]);

  int refinedRoadPointsNum = inputRefinedRoad->GetNumberOfPoints();
  // int eventPointsNum = inputEventwithClosetPoint->GetNumberOfPoints();

  // std::vector<int> rRoadPointtoEventSampleDictionary(refinedRoadPointsNum);

  // this->CalculateRRoadPointsEventSamples<vtkIdTypeArray>(
  //   rRoadPointtoEventSampleDictionary, eventPointsNum,
  //   (vtkIdTypeArray *)inputEventwithClosetPoint->GetPointData()->GetArray(
  //     "closestPoints"));

  auto refinedRoadPointWeight = vtkSmartPointer<vtkFloatArray>::New();
  refinedRoadPointWeight->SetName("pointWeight");
  refinedRoadPointWeight->SetNumberOfComponents(1);
  refinedRoadPointWeight->SetNumberOfTuples(refinedRoadPointsNum);

  // get neighbors for all the vertics in the rRoad network as a vetor such
  // that
  // could be send to core for density calculataion.
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputRefinedRoad);
  triangulation->preconditionVertexNeighbors();

  // auto nVertices = triangulation -> getNumberOfVertices();
  // for(int i = 0; i < nVertices; i++) {
  //   auto nNeighbours = triangulation->getVertexNeighborNumber(i);
  //   // string line = "Neighbours of " + to_string(i) + ": ";
  //   for(int j = 0; j < nNeighbours; j++) {
  //     ttk::SimplexId neighbourId = 0;
  //     triangulation->getVertexNeighbor(i, j, neighbourId);
  //     line += to_string(neighbourId) + " ";
  //   }
  //   // this->printMsg(line);
  //   // if(i > 10)
  //   //   break;
  // }

  // auto testingRefinedRoads = (float *)inputRefinedRoad->GetPointData()
  //                              ->GetArray("eventSample2rRoadPoint")
  //                              ->GetVoidPointer(0);

  // for(int i = 0; i < refinedRoadPointsNum; i++) {
  //   std::cout << "testingRefinedRoads: " << testingRefinedRoads[i] << std::endl;
  // }

  auto status_to_get_roadDensity = this->CalculateRoadWeight<float>(
    (float *)inputRefinedRoad->GetPointData()
      ->GetArray("eventSample2rRoadPoint")
      ->GetVoidPointer(0),
    refinedRoadPointsNum, triangulation,
    (float *)(inputRefinedRoad->GetPoints()->GetVoidPointer(0)),
    (float *)refinedRoadPointWeight->GetVoidPointer(0), KernelBandwidth);

  // finalize output
  {
    // for(int i = 0; i < nEdges * 3; i += 3)
    //   cout << cellIds[i] << " " << cellIds[i + 1] << " " << cellIds[i + 2]
    //        << endl;

    auto output = vtkUnstructuredGrid::GetData(outputVector);

    output->ShallowCopy(inputRefinedRoad);

    output->GetPointData()->AddArray(refinedRoadPointWeight);

    // output->Print(std::cout);
  } // Get the output

  return 1;
}

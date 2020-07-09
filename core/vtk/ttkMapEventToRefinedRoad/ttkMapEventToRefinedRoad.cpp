#include <ttkMapEventToRefinedRoad.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <sstream>
#include <string>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

#include <ttkUtils.h>

using namespace std;

vtkStandardNewMacro(ttkMapEventToRefinedRoad);

ttkMapEventToRefinedRoad::ttkMapEventToRefinedRoad() {
  this->SetThresholdtoDefineCloseness(500);

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkMapEventToRefinedRoad::~ttkMapEventToRefinedRoad() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkMapEventToRefinedRoad::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(
        vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
      break;
    case 1:
      info->Set(
        vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
      break;
    default:
      return 0;
  }
  return 1;
}

// see ttkAlgorithm::FillOutputPortInformation for details about this method
int ttkMapEventToRefinedRoad::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    default:
      return 0;
  }
  return 1;
}

int ttkMapEventToRefinedRoad::RequestData(vtkInformation *request,
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  float threshold4closeness = this->GetThresholdtoDefineCloseness();
  this->printMsg("threshold4closeness: " + to_string(threshold4closeness));
  // Get the input
  auto inputRefinedRoad = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto inputEvent = vtkUnstructuredGrid::GetData(inputVector[1]);

  int eventPointsNum = inputEvent->GetNumberOfPoints();
  int refinedRoadPointsNum = inputRefinedRoad->GetNumberOfPoints();
  std::cout << "eventPointsNum: " << eventPointsNum << std::endl;
  std::cout << "refinedRoadPointsNum: " << refinedRoadPointsNum << std::endl;

  auto eventSample2rRoadPoint = vtkSmartPointer<vtkFloatArray>::New();
  eventSample2rRoadPoint->SetName("eventSample2rRoadPoint");
  eventSample2rRoadPoint->SetNumberOfComponents(1);
  eventSample2rRoadPoint->SetNumberOfTuples(refinedRoadPointsNum);

  this->calculateClosetPointofRoadforEventData<float>(
    (float *)ttkUtils::GetVoidPointer(inputEvent->GetPoints()),
    (float *)ttkUtils::GetVoidPointer(inputRefinedRoad->GetPoints()),
    (float *)ttkUtils::GetVoidPointer(eventSample2rRoadPoint),
    // (float *)(inputEvent->GetPoints()->GetVoidPointer(0)),
    // (float *)(inputRefinedRoad->GetPoints()->GetVoidPointer(0)),
    // (float *)eventSample2rRoadPoint->GetVoidPointer(0),
    eventPointsNum, refinedRoadPointsNum, threshold4closeness);

  // finalize output
  {
    auto output = vtkUnstructuredGrid::GetData(outputVector);

    output->ShallowCopy(inputRefinedRoad);

    output->GetPointData()->AddArray(eventSample2rRoadPoint);

    output->Print(std::cout);
  } // Get the output

  return 1;
}

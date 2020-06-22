#include <ttkRoadRefinement.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <sstream>
#include <string>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkIntArray.h>
#include <vtkPointData.h>

#include <ttkUtils.h>

using namespace std;

vtkStandardNewMacro(ttkRoadRefinement);

ttkRoadRefinement::ttkRoadRefinement() {
  this->SetRoadSegmentUnit(10);

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkRoadRefinement::~ttkRoadRefinement() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkRoadRefinement::FillInputPortInformation(int port,
                                                vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(
        vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
      break;
    default:
      return 0;
  }
  return 1;
}

// see ttkAlgorithm::FillOutputPortInformation for details about this method
int ttkRoadRefinement::FillOutputPortInformation(int port,
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

int ttkRoadRefinement::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  float segmentUnit = this->GetRoadSegmentUnit();
  this->printMsg("segmentUnit: " + to_string(segmentUnit));

  // Get the input
  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);

  int nRefinedPoints = 0;
  int nRefinedEdges = 0;
//   this->getNumberOfRefinedPointsAndEdges<vtkIdType>(
//     (float *)ttkUtils::GetVoidPointer(input->GetPoints()),
//     input->GetCells()->GetPointer(),
//     input->GetNumberOfCells(),
//     nRefinedPoints,
//     nRefinedEdges,
//     segmentUnit
//   );

  std::cout << "nPoints " << nRefinedPoints << std::endl;
  std::cout << "nEdges " << nRefinedEdges << std::endl;

  int nPoints = nRefinedPoints + input->GetNumberOfPoints();
  int nEdges = nRefinedEdges + input->GetNumberOfCells();

  std::cout << "nPoints " << nPoints << std::endl;
  std::cout << "nEdges " << nEdges << std::endl;
  // TODO: initialize buffers
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nPoints);

  auto cells = vtkSmartPointer<vtkCellArray>::New();
//   auto cellIds
//     = (long long int *)cells->WritePointer(nEdges, 3 * nEdges);

  // TODO: fill buffers
  // vtkIdType = int ; int32; int64
//   int status = this->refineRoads<vtkIdType>(
//     (float *)  ttkUtils::GetVoidPointer(input->GetPoints()),
//     input->GetCells()->GetPointer(), segmentUnit,
//     (float *)ttkUtils::GetVoidPointer(points->GetVoidPointer),
    // cellIds,
    // input->GetNumberOfCells(),
//     input->GetNumberOfPoints()
// );

  // finalize output
  {
    // for(int i = 0; i < nEdges * 3; i += 3)
    //   cout << cellIds[i] << " " << cellIds[i + 1] << " " << cellIds[i + 2]
    //        << endl;

    auto output = vtkUnstructuredGrid::GetData(outputVector);
    output->SetPoints(points);
    output->SetCells(VTK_LINE, cells);

    // output->Print(std::cout);
  } // Get the output

  // auto output = vtkUnstructuredGrid::GetData(outputVector);

  // output->ShallowCopy(input);

  // output->Print(cout);

  return 1;
}

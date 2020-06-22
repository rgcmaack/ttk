#include <ttkRoadDataConverter.h>

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

vtkStandardNewMacro(ttkRoadDataConverter);

ttkRoadDataConverter::ttkRoadDataConverter() {
  this->SetFilepath("");

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkRoadDataConverter::~ttkRoadDataConverter() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkRoadDataConverter::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  return 0;
}

// see ttkAlgorithm::FillOutputPortInformation for details about this method
int ttkRoadDataConverter::FillOutputPortInformation(int port,
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

int ttkRoadDataConverter::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  std::string path = this->GetFilepath();
  this->printMsg("path: " + path);

  // Get the output
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  // TODO: compute buffer sizes
  int *pointsandnedges = this->getNumberOfPointsandEdges(path);

  // this->printMsg("pointandedges " + *pointsandnedges);

  int npoints = pointsandnedges[0];
  int nedges = pointsandnedges[1];


  // TODO: initialize buffers
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(npoints);

  auto cells = vtkSmartPointer<vtkCellArray>::New();
//   auto cellIds = (long long int*) cells->WritePointer(nedges, 3 * nedges);

  std::cout << "npoints buffer size:" << npoints << std::endl;
  std::cout << "nedges buffer size:" << 3 * nedges << std::endl;

  // // TODO: fill buffers
  // vtkIdType = int ; int32; int64
//   int status = this->parsePointCoords(
//     path,
//     (float *)points->GetVoidPointer(0),
//     cellIds
//   );

  // finalize output
  {

    // for(int i = 0; i < nedges * 3; i += 3)
    //   cout << cellIds[i] << " " << cellIds[i + 1] << " " << cellIds[i + 2]
    //        << endl;

    auto output = vtkUnstructuredGrid::GetData(outputVector);
    output->SetPoints(points);
    output->SetCells(VTK_LINE, cells);

    output->Print(std::cout);
  }

  return 1;
}

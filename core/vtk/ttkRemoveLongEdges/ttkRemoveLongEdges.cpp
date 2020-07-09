#include <ttkRemoveLongEdges.h>

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

vtkStandardNewMacro(ttkRemoveLongEdges);

ttkRemoveLongEdges::ttkRemoveLongEdges() {
  this->SetEdgeDistanceThreshold(1000);

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkRemoveLongEdges::~ttkRemoveLongEdges() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkRemoveLongEdges::FillInputPortInformation(int port,
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
int ttkRemoveLongEdges::FillOutputPortInformation(int port,
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

int ttkRemoveLongEdges::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  float edgeThreshold = this->GetEdgeDistanceThreshold();
  this->printMsg("edgeThreshold: " + to_string(edgeThreshold));

  // Get the input
  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);

  int nNewEdges = 0;

  // get buffer size
  this->getNumberOfnewEdges<vtkIdType>(
    (float *)ttkUtils::GetVoidPointer(input->GetPoints()),
    (vtkIdType *)ttkUtils::GetVoidPointer(
      input->GetCells()->GetConnectivityArray()),
    input->GetNumberOfCells(), nNewEdges, edgeThreshold);

  std::cout << "nNewEdges " << nNewEdges << std::endl;

  // TODO: initialize buffers
  auto points = input->GetPoints();

  auto cellOffsets = vtkSmartPointer<vtkIdTypeArray>::New();
  cellOffsets->SetNumberOfTuples(nNewEdges + 1);
  auto cellOffsetsData = (vtkIdType *)ttkUtils::GetVoidPointer(cellOffsets);
  for(int i = 0; i < nNewEdges + 1; i++)
    cellOffsetsData[i] = i * 2;

  auto cellConnectivity = vtkSmartPointer<vtkIdTypeArray>::New();
  cellConnectivity->SetNumberOfTuples(2 * nNewEdges);
  auto cellConnectivityData
    = (vtkIdType *)ttkUtils::GetVoidPointer(cellConnectivity);

  // TODO: fill buffers
  //// vtkIdType = int;int32;int64
  int status = this->removeEdges<vtkIdType>(
    (float *)ttkUtils::GetVoidPointer(input->GetPoints()),
    (vtkIdType *)ttkUtils::GetVoidPointer(
      input->GetCells()->GetConnectivityArray()),
    edgeThreshold, cellConnectivityData, input->GetNumberOfCells());

  // finalize output
  {
    auto output = vtkUnstructuredGrid::GetData(outputVector);
    output->SetPoints(points);

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetData(cellOffsets, cellConnectivity);
    output->SetCells(VTK_LINE, cellArray);

    // output->Print(std::cout);
  } // Get the output

  return 1;
}

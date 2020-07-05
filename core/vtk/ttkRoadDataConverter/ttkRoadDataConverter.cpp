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
#include <vtkUnsignedCharArray.h>

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
  if(path.empty()) {
    return 0;
  }

  int npoints = 0, nedges = 0;

  // TODO: compute buffer sizes
  int getBuffer = this->getNumberOfPointsandEdges(path, npoints, nedges);

  // TODO: initialize buffers
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(npoints);

  auto cellOffsets = vtkSmartPointer<vtkIdTypeArray>::New();
  cellOffsets->SetNumberOfTuples(nedges + 1);
  auto cellOffsetsData = (vtkIdType *)ttkUtils::GetVoidPointer(cellOffsets);
  for(int i = 0; i < nedges + 1; i++)
    cellOffsetsData[i] = i * 2;

  auto cellConnectivity = vtkSmartPointer<vtkIdTypeArray>::New();
  cellConnectivity->SetNumberOfTuples(nedges * 2);
  auto cellConnectivityData
    = (vtkIdType *)ttkUtils::GetVoidPointer(cellConnectivity);

  auto category4edge = vtkSmartPointer<vtkUnsignedCharArray>::New();
  category4edge->SetName("Category4Edge");
  category4edge->SetNumberOfComponents(1);
  category4edge->SetNumberOfTuples(nedges);

  std::vector<std::string> categoryDictionary;

  // TODO: fill buffers
  // vtkIdType = int ; int32; int64
  int status = this->parsePointCoords(
    path, (float *)ttkUtils::GetVoidPointer(points), cellConnectivityData,
    (unsigned char *)ttkUtils::GetVoidPointer(category4edge),
    categoryDictionary);

  // finalize output
  {
    // add the file and property to the output
    auto categoryDictionaryArray = vtkSmartPointer<vtkStringArray>::New();
    categoryDictionaryArray->SetName("CategoryDictionary");
    categoryDictionaryArray->SetNumberOfComponents(1);
    categoryDictionaryArray->SetNumberOfTuples(categoryDictionary.size());

    for(int i = 0; i < categoryDictionary.size(); i++)
      categoryDictionaryArray->SetValue(i, categoryDictionary[i]);

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetData(cellOffsets, cellConnectivity);

    auto output = vtkUnstructuredGrid::GetData(outputVector);
    output->SetPoints(points);
    output->SetCells(VTK_LINE, cellArray);

    // output->GetCellData()->AddArray(category4edge);
    // output->GetFieldData()->AddArray(categoryDictionaryArray);

    output->Print(std::cout);
  }

  return 1;
}

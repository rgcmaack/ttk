#include <ttkEventDataConverter.h>

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
#include <vtkStringArray.h>

using namespace std;

vtkStandardNewMacro(ttkEventDataConverter);

ttkEventDataConverter::ttkEventDataConverter() {
  this->SetFilepath("");

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkEventDataConverter::~ttkEventDataConverter() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkEventDataConverter::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  return 0;
}

// see ttkAlgorithm::FillOutputPortInformation for details about this method
int ttkEventDataConverter::FillOutputPortInformation(int port,
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

int ttkEventDataConverter::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  std::string path = this->GetFilepath();
  this->printMsg("path: " + path);

  // Get the output
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  // TODO: compute buffer sizes
  int nPoints = this->getNumberOfPoints(path);
  this->printMsg("nPoints: " + std::to_string(nPoints));

  // TODO: initialize buffers
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nPoints);

  auto categoryIndex = vtkSmartPointer<vtkUnsignedCharArray>::New();
  categoryIndex->SetName("CategoryIndex");
  categoryIndex->SetNumberOfComponents(1);
  categoryIndex->SetNumberOfTuples(nPoints);

  std::vector<string> categoryDictionary;

  // TODO: fill buffers
  // int status = this->parsePointCoords(
  //   path, (float *)points->GetVoidPointer(0),
  //   (unsigned char *)categoryIndex->GetVoidPointer(0), categoryDictionary);

  // finalize output
  // {
  //   // TODO: Cells
  //   auto cells = vtkSmartPointer<vtkCellArray>::New();
  //   auto cellIds = cells->WritePointer(nPoints, 2 * nPoints);
  //   for(int i = 0; i < nPoints; i++) {
  //     cellIds[2 * i] = 1;
  //     cellIds[2 * i + 1] = i;
  //   }

  //   // for(int i = 0; i < nPoints * 2; i += 2)
  //   //   cout << cellIds[i] << " " << cellIds[i + 1] << endl;

  //   auto output = vtkUnstructuredGrid::GetData(outputVector);
  //   output->SetPoints(points);
  //   output->SetCells(VTK_VERTEX, cells);

  //   auto categoryDictionaryArray = vtkSmartPointer<vtkStringArray>::New();
  //   categoryDictionaryArray->SetName("CategoryDictionary");
  //   categoryDictionaryArray->SetNumberOfComponents(1);
  //   categoryDictionaryArray->SetNumberOfTuples(categoryDictionary.size());

  //   for(int i = 0; i < categoryDictionary.size(); i++)
  //     categoryDictionaryArray->SetValue(i, categoryDictionary[i]);

  //   output->GetPointData()->AddArray(categoryIndex);

  //   output->GetFieldData()->AddArray(categoryDictionaryArray);

  //   // output->Print(std::cout);
  // }

  return 1;
}

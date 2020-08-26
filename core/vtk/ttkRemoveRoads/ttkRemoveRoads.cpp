#include <ttkRemoveRoads.h>

#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkUnsignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkRemoveRoads);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkRemoveRoads::ttkRemoveRoads() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkRemoveRoads::~ttkRemoveRoads() {
}

int ttkRemoveRoads::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRemoveRoads::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

int ttkRemoveRoads::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  vtkUnstructuredGrid *inputNetwork = vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *inputSelection = vtkUnstructuredGrid::GetData(inputVector[1]);

  vtkDataArray *inputArraySelection = this->GetInputArrayToProcess(0, inputSelection);

  if(!inputArraySelection){ return 0; }

  size_t nEdges = inputNetwork->GetNumberOfCells();
  size_t nSelectedEdges = inputArraySelection->GetNumberOfTuples();

  ttk::Timer timer;
  this->printMsg("Removing "+std::to_string(nSelectedEdges)+" of "+std::to_string(nEdges)+" edges.",0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

  auto selectionData = (int*) ttkUtils::GetVoidPointer(inputArraySelection);

  auto mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
  mask->SetName("MASK");
  mask->SetNumberOfTuples( nEdges );
  auto maskData = (unsigned char*) ttkUtils::GetVoidPointer(mask);

  for(size_t i=0; i<nEdges; i++){
      maskData[i] = 1;
  }

  for(size_t i=0; i<nSelectedEdges; i++){
      maskData[ selectionData[i] ] = 0;
  }

  auto temp = vtkSmartPointer<vtkUnstructuredGrid>::New();
  temp->ShallowCopy(inputNetwork);
  temp->GetCellData()->AddArray(mask);

  auto threshold = vtkSmartPointer<vtkThreshold>::New();
  threshold->SetInputData( temp );
  threshold->SetInputArrayToProcess(0,0,0,1,"MASK");
  threshold->ThresholdByUpper(1);
  threshold->Update();

  vtkUnstructuredGrid *outputNetwork = vtkUnstructuredGrid::GetData(outputVector, 0);

  outputNetwork->ShallowCopy(threshold->GetOutput());

  this->printMsg("Removing "+std::to_string(nSelectedEdges)+" of "+std::to_string(nEdges)+" edges.",1,timer.getElapsedTime(),this->threadNumber_);

  return 1;
}

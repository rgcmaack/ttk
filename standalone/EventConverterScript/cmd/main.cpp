/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include <vtkSmartPointer.h>
// #include <ttkHelloWorld.h>
// #include <ttkIcoSphere.h>
#include <ttkEventDataConverter.h>
// #include <ttkProgramBase.h>

#include <ttkCinemaWriter.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkXMLUnstructuredGridReader.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  ttk::globalDebugLevel_ = 3;

  // // debug eventDataConverter
  auto eventDataConverter = vtkSmartPointer<ttkEventDataConverter>::New();
  eventDataConverter->SetFilepath(
    "/home/ray/Documents/projects/hotspot/data/SeattleCrimes.geojson");

  auto writer = vtkSmartPointer<ttkCinemaWriter>::New();
  writer->SetDatabasePath(
    "/home/ray/Documents/projects/hotspot/data/SeattleCrimes.cdb");
  writer->SetInputConnection(0, eventDataConverter->GetOutputPort(0));
  writer->Update();

  return 1;
}

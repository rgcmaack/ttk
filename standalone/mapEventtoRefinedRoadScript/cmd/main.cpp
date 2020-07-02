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
#include <ttkMapEventToRefinedRoad.h>
#include <ttkRoadDataConverter.h>
#include <ttkRoadRefinement.h>
// #include <ttkRoadDensity.h>

#include <ttkCinemaWriter.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkXMLUnstructuredGridReader.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  ttk::globalDebugLevel_ = 3;

  // get event data
  auto eventDataConverter = vtkSmartPointer<ttkEventDataConverter>::New();
  eventDataConverter->SetFilepath(
    "/home/ray/Documents/projects/hotspot/data/SeattleCrimes.geojson");

  // get refined road data
  auto roadDataConverter = vtkSmartPointer<ttkRoadDataConverter>::New();
  roadDataConverter->SetFilepath(
    "/home/ray/Documents/projects/hotspot/data/SeattleRoads.geojson");
  auto filterPoints = vtkSmartPointer<vtkUnstructuredGridGeometryFilter>::New();
  filterPoints->SetInputConnection(0, roadDataConverter->GetOutputPort(0));

  auto roadSegments = vtkSmartPointer<ttkRoadRefinement>::New();
  roadSegments->SetRoadSegmentUnit(20);
  roadSegments->SetInputConnection(0, filterPoints->GetOutputPort(0));

  // map event to refined data
  auto mapEvent2RRoad = vtkSmartPointer<ttkMapEventToRefinedRoad>::New();
  mapEvent2RRoad->SetInputConnection(0, roadSegments->GetOutputPort(0));
  mapEvent2RRoad->SetInputConnection(1, eventDataConverter->GetOutputPort(0));

  auto writer = vtkSmartPointer<ttkCinemaWriter>::New();
  writer->SetDatabasePath(
    "/home/ray/Documents/projects/hotspot/data/road20units_500closeness.cdb");
  writer->SetInputConnection(0, mapEvent2RRoad->GetOutputPort(0));
  writer->Update();

  return 1;
}

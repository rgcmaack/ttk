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
#include <ttkRoadDensity.h>

#include <ttkArrayEditor.h>
#include <ttkCinemaWriter.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkXMLUnstructuredGridReader.h>

using namespace std;
using namespace ttk;

int main(int argc, char **argv) {

  ttk::globalDebugLevel_ = 3;

  // // debug eventDataConverter
  // auto eventDataConverter = vtkSmartPointer<ttkEventDataConverter>::New();
  // eventDataConverter->SetFilepath(
  //   "/home/local/ASUAD/rzhan100/ttk-data/data/PurdueEvent.json");

  // auto roadDataConverter = vtkSmartPointer<ttkRoadDataConverter>::New();
  // roadDataConverter->SetFilepath(
  //   "/home/local/ASUAD/rzhan100/ttk-data/data/purdue_road.geojson");
  // auto filterPoints =
  // vtkSmartPointer<vtkUnstructuredGridGeometryFilter>::New();
  // filterPoints->SetInputConnection(0, roadDataConverter->GetOutputPort(0));

  // ==============================================================
  // Load Data From Database
  // ==============================================================
  //   auto cinemaReader = vtkSmartPointer<ttkCinemaReader>::New();
  //   cinemaReader->SetDatabasePath(
  //     "/home/local/ASUAD/rzhan100/ttk-data/purdue_road_unique.cdb/");
  //   cinemaReader->SetFilePathColumnNames("FILE");

  //   auto cinemaQuery = vtkSmartPointer<ttkCinemaQuery>::New();
  //   cinemaQuery->SetSQLStatement( "SELECT * FROM vtkTable0 LIMIT 1" );
  //   cinemaQuery->SetInputConnection(0, cinemaReader->GetOutputPort(0));

  //   auto cinemaProductReader =
  //   vtkSmartPointer<ttkCinemaProductReader>::New();
  //   cinemaProductReader->SetFilepathColumnName("FILE");
  //   cinemaProductReader->SetInputConnection( 0,
  //   cinemaQuery->GetOutputPort(0)
  //   );

  // get the refinedRoad data
  // auto readerRoadUniqe = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  // readerRoadUniqe->SetFileName(
  //   "/home/local/ASUAD/rzhan100/ttk-data/"
  //   "purdue_road_unique.cdb/data/289383/289383_0.vtu");

  // auto roadSegments = vtkSmartPointer<ttkRoadRefinement>::New();
  // roadSegments->SetRoadSegmentUnit(20);
  // roadSegments->SetInputConnection(0, readerRoadUniqe->GetOutputPort(0));

  // get the event data with categroty and refinedroad data
  // auto readerEvent = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  // readerEvent->SetFileName(
  //   "/home/local/ASUAD/rzhan100/ttk-data/"
  //   "PurdueEventDatawithCategory.cdb/data/289383/289383_0.vtu");

  // auto readerRefinedRoad = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  // readerRefinedRoad->SetFileName(
  //   "/home/local/ASUAD/rzhan100/ttk-data/"
  //   "purdue_road_segment_20.cdb/data/289383/289383_0.vtu");

  // // send the two data to e2rr modeule
  // auto mapEvent2RRoad = vtkSmartPointer<ttkMapEventToRefinedRoad>::New();
  // mapEvent2RRoad->SetInputConnection(0, readerRefinedRoad->GetOutputPort(0));
  // mapEvent2RRoad->SetInputConnection(1, readerEvent->GetOutputPort(0));

  // auto writer = vtkSmartPointer<ttkCinemaWriter>::New();
  // writer->SetDatabasePath(
  //   "/home/local/ASUAD/rzhan100/ttk-data/rRoadPoints_20seg_100threashold_with_eventSamples.cdb");
  // writer->SetInputConnection(0, mapEvent2RRoad->GetOutputPort(0));
  // writer->Update();

  //   // write the event data with the near points of road to file
  //   auto addFD = vtkSmartPointer<ttkArrayEditor>::New();
  //   addFD->SetFieldDataString("Name,eventDatawithNearPoints");
  //   addFD->SetInputConnection(0, mapEvent2RRoad->GetOutputPort(0));

  // ////Code for Road Density Module
  // // get the refinedRoad data
  // auto readerRefinedRoad = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  // readerRefinedRoad->SetFileName(
  //   "/home/local/ASUAD/rzhan100/ttk-data/"
  //   "purdue_road_segment_20.cdb/data/289383/289383_0.vtu");

  // // get the event data with closet points
  auto readerRroadPointswithEventSamples
    = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  readerRroadPointswithEventSamples->SetFileName(
    "/home/local/ASUAD/rzhan100/ttk-data/"
    "rRoadPoints_20seg_100threashold_with_eventSamples.cdb/"
    "data/purdueRoadSegments_20.vtu");

  // send the two data to e2rr modeule
  auto RoadDensity = vtkSmartPointer<ttkRoadDensity>::New();
  RoadDensity->SetKernelFunction("Gaussian");
  RoadDensity->SetKernelBandwidth(500);
  RoadDensity->SetInputConnection(
    0, readerRroadPointswithEventSamples->GetOutputPort(0));

  //   // write the event data with the near points of road to file
  //   auto addFD = vtkSmartPointer<ttkArrayEditor>::New();
  //   addFD->SetFieldDataString("Name,eventDatawithNearPoints");
  //   addFD->SetInputConnection(0, mapEvent2RRoad->GetOutputPort(0));

  auto writer = vtkSmartPointer<ttkCinemaWriter>::New();
  writer->SetDatabasePath(
    "/home/local/ASUAD/rzhan100/ttk-data/roadDensity_20seg_100threshold4closeness_500band_linearKernel.cdb");
  writer->SetInputConnection(0, RoadDensity->GetOutputPort(0));
  writer->Update();

  return 1;
}

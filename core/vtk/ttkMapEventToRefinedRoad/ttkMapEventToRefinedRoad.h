/// \ingroup vtk
/// \class ttkMapEventToRefinedRoad
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkMapEventToRefinedRoadModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK Base Includes
#include <MapEventToRefinedRoad.h>

class TTKMAPEVENTTOREFINEDROAD_EXPORT ttkMapEventToRefinedRoad
  : public ttkAlgorithm,
    public ttk::MapEventToRefinedRoad {
public:
  static ttkMapEventToRefinedRoad *New();
  vtkTypeMacro(ttkMapEventToRefinedRoad, ttkAlgorithm);

  vtkGetMacro(ThresholdtoDefineCloseness, double);
  vtkSetMacro(ThresholdtoDefineCloseness, double);
  // std::string array2string(std::string *arr);

protected:
  ttkMapEventToRefinedRoad();
  ~ttkMapEventToRefinedRoad();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  double ThresholdtoDefineCloseness;
};
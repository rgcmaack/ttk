/// \ingroup vtk
/// \class ttkRoadRefinement
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkRoadRefinementModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK Base Includes
#include <RoadRefinement.h>

class TTKROADREFINEMENT_EXPORT ttkRoadRefinement
  : public ttkAlgorithm,
    public ttk::RoadRefinement {
public:
  static ttkRoadRefinement *New();
  vtkTypeMacro(ttkRoadRefinement, ttkAlgorithm);

  vtkGetMacro(RoadSegmentUnit, double);
  vtkSetMacro(RoadSegmentUnit, double);

  // std::string array2string(std::string *arr);

protected:
  ttkRoadRefinement();
  ~ttkRoadRefinement();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  double RoadSegmentUnit;
};
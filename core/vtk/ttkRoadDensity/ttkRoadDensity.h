/// \ingroup vtk
/// \class ttkRoadDensity
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkRoadDensityModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK Base Includes
#include <RoadDensity.h>

class TTKROADDENSITY_EXPORT ttkRoadDensity
  : public ttkAlgorithm,
    public ttk::RoadDensity {
public:
  static ttkRoadDensity *New();
  vtkTypeMacro(ttkRoadDensity, ttkAlgorithm);

  vtkGetMacro(KernelFunction, std::string);
  vtkSetMacro(KernelFunction, std::string);
  vtkGetMacro(KernelBandwidth, float);
  vtkSetMacro(KernelBandwidth, float);
  // std::string array2string(std::string *arr);

protected:
  ttkRoadDensity();
  ~ttkRoadDensity();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  std::string KernelFunction;
  float KernelBandwidth;
};
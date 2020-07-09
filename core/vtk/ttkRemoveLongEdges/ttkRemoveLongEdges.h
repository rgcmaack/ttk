/// \ingroup vtk
/// \class ttkRemoveLongEdges
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkRemoveLongEdgesModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK Base Includes
#include <RemoveLongEdges.h>

class TTKREMOVELONGEDGES_EXPORT ttkRemoveLongEdges
  : public ttkAlgorithm,
    public ttk::RemoveLongEdges {
public:
  static ttkRemoveLongEdges *New();
  vtkTypeMacro(ttkRemoveLongEdges, ttkAlgorithm);

  vtkGetMacro(EdgeDistanceThreshold, double);
  vtkSetMacro(EdgeDistanceThreshold, double);

  // std::string array2string(std::string *arr);

protected:
  ttkRemoveLongEdges();
  ~ttkRemoveLongEdges();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  double EdgeDistanceThreshold;
};
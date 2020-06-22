/// \ingroup vtk
/// \class ttkRoadDataConverter
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkRoadDataConverterModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK Base Includes
#include <RoadDataConverter.h>

class TTKROADDATACONVERTER_EXPORT ttkRoadDataConverter
  : public ttkAlgorithm,
    public ttk::RoadDataConverter {
public:
  static ttkRoadDataConverter *New();
  vtkTypeMacro(ttkRoadDataConverter, ttkAlgorithm);

  vtkGetMacro(Filepath, std::string);
  vtkSetMacro(Filepath, std::string);

  // std::string array2string(std::string *arr);

protected:
  ttkRoadDataConverter();
  ~ttkRoadDataConverter();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  std::string Filepath;
};
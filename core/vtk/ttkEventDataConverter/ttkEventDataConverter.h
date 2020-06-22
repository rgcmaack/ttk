/// \ingroup vtk
/// \class ttkEventDataConverter
/// \author Jonas Lukasczyk <jl@jluk.der>
/// \date 01.09.2019.
///
/// TODO

#pragma once

// VTK Module
#include <ttkEventDataConverterModule.h>

// VTK Includes
#include <vtkInformation.h>
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <EventDataConverter.h>

class TTKEVENTDATACONVERTER_EXPORT ttkEventDataConverter
    : public ttkAlgorithm
    , public ttk::EventDataConverter
{
    public:
        static ttkEventDataConverter *New();
        vtkTypeMacro(ttkEventDataConverter, ttkAlgorithm);

        vtkGetMacro(Filepath, std::string);
        vtkSetMacro(Filepath, std::string);

        // std::string array2string(std::string *arr);

      protected:
        ttkEventDataConverter();
        ~ttkEventDataConverter();

        int FillInputPortInformation(int port, vtkInformation* info) override;
        int FillOutputPortInformation(int port, vtkInformation* info) override;

        int RequestData(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector
        ) override;

        std::string Filepath;
};
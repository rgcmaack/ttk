#pragma once

// VTK Module
#include <ttkMergeTreeSimplificationModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <MergeTreeSimplification.h>

class TTKMERGETREESIMPLIFICATION_EXPORT ttkMergeTreeSimplification
    : public ttkAlgorithm    // we inherit from the generic ttkAlgorithm class
    , public ttk::MergeTreeSimplification // and we inherit from the base class
{
    private:
        double PersistenceThreshold{0};

    public:
        vtkSetMacro(PersistenceThreshold, double);
        vtkGetMacro(PersistenceThreshold, double);

        static ttkMergeTreeSimplification *New();
        vtkTypeMacro(ttkMergeTreeSimplification, ttkAlgorithm);

    protected:
        ttkMergeTreeSimplification();
        ~ttkMergeTreeSimplification() override;

        int FillInputPortInformation(int port, vtkInformation* info) override;
        int FillOutputPortInformation(int port, vtkInformation* info) override;
        int RequestData(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector
        ) override;
};
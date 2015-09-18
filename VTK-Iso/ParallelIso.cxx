/*
 * PVTKDemo.cpp
 *
 *  Created on: Aug 12, 2013
 *      Author: ollie
 */

#define SUPERNOVA

#include <vector>
#include <cmath>

#include "vtkNew.h"
#include "vtkImageData.h"
#include "vtkContourFilter.h"
#include "vtkMarchingCubes.h"
#include "vtkProcessIdScalars.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCompositeRenderManager.h"
#include "vtkTimerLog.h"

#include "vtkMPIController.h"

#include "vtkNonMergingPointLocator.h"

#ifdef SUPERNOVA

#include "vtkPNrrdReader.h"
#include "vtkImageResample.h"

#else

#include "vtkSampleFunction.h"
#include "vtkTangle.h"

static int EXTENT = 512;

#endif


#include "Stats.h"

static const float ISO_START=0.045;
static const float ISO_STEP=0.005;
static const int ISO_NUM=10;

Timer timer;

void process(vtkMultiProcessController* controller, void* vtkNotUsed(arg))
{
    int myId = controller->GetLocalProcessId();

    // Setup data source and filter pipeline

#ifdef SUPERNOVA
    vtkNew<vtkPNrrdReader> reader;
    std::string file = "/home/csewell/CGABenchmarks/data/normal_1349.nhdr";
    reader->SetFileName(file.c_str());

    vtkNew<vtkImageResample> resample;
    double resampleSize = 1.0f;
    resample->SetInputConnection(reader->GetOutputPort());
    resample->SetAxisMagnificationFactor(0,resampleSize);
    resample->SetAxisMagnificationFactor(1,resampleSize);
    resample->SetAxisMagnificationFactor(2,resampleSize);
    resample->UpdateInformation();
    resample->SetUpdateExtent(controller->GetLocalProcessId(),
       controller->GetNumberOfProcesses(), 0 );

    resample->Update();
#else
    vtkNew<vtkTangle> tangle;

    int rank = controller->GetLocalProcessId();
    vtkNew<vtkSampleFunction> sampler;
    sampler->SetImplicitFunction(tangle.GetPointer());
    sampler->SetModelBounds(-1, 1, -1, 1, -1, 1);
    sampler->SetSampleDimensions(EXTENT, EXTENT, EXTENT);
    sampler->ComputeNormalsOff();
    sampler->UpdateInformation();
    sampler->SetUpdateExtent(controller->GetLocalProcessId(),
       controller->GetNumberOfProcesses(), 0 );
    sampler->Update();

    if (rank == 1) {
        std::cout << controller->GetNumberOfProcesses() << std::endl;
        std::cout << sampler->GetOutput()->GetNumberOfCells()
                  << std::endl;
        }
#endif

    vtkNew<vtkMarchingCubes> isosurface;
   
#ifdef SUPERNOVA
    isosurface->SetInputData(resample->GetOutput());
    isosurface->SetValue(0, ISO_START);
#else
    isosurface->SetInputData(sampler->GetOutput());
    isosurface->SetValue(0, 0.3);
#endif
    
    isosurface->ComputeScalarsOff();
    isosurface->ComputeGradientsOff();
    isosurface->ComputeNormalsOn();

    vtkNonMergingPointLocator* simpleLocator = vtkNonMergingPointLocator::New();
    isosurface->SetLocator(simpleLocator);

    controller->Barrier();
    const size_t MAX_ITERATIONS = 1;
    std::vector<double> samples;
    samples.reserve(MAX_ITERATIONS);
    timer.Reset();
 
#ifdef SUPERNOVA
    for (int j = 0; j < ISO_NUM; ++j)
    {
        isosurface->SetValue(0, isosurface->GetValue(0) + ISO_STEP);
        isosurface->Update();
        controller->Barrier();
        vtkPolyData* output = isosurface->GetOutput();
        std::cout << (isosurface->GetValue(0)) << " " << output->GetNumberOfPoints() << std::endl;
    }
#else
    isosurface->Update();
    controller->Barrier();
#endif

    if (myId == 0)
    {
      samples.push_back(timer.GetElapsedTime());

      std::sort(samples.begin(), samples.end());
      stats::Winsorize(samples, 5.0);
      std::cout << "Benchmark \'VTK MPI Isosurface\' results:\n"
          << "\tmedian = " << stats::PercentileValue(samples, 50.0) << "s\n"
          << "\tmedian abs dev = " << stats::MedianAbsDeviation(samples) << "s\n"
          << "\tmean = " << stats::Mean(samples) << "s\n"
          << "\tstd dev = " << stats::StandardDeviation(samples) << "s\n"
          << "\tmin = " << samples.front() << "s\n"
          << "\tmax = " << samples.back() << "s\n"
          << "\t# of runs = " << samples.size() << "\n";
    }

#ifdef RENDER_OUTPUT

    vtkNew<vtkProcessIdScalars> procId;
    procId->SetInputData(isosurface->GetOutput());
    procId->SetController(controller);

    // Rendering objects.
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(procId->GetOutputPort());
    mapper->ImmediateModeRenderingOff();
    mapper->CreateDefaultLookupTable();
    mapper->SetScalarRange(0, controller->GetNumberOfProcesses());
    mapper->SetNumberOfPieces(controller->GetNumberOfProcesses());
    mapper->SetPiece(controller->GetLocalProcessId());

    // IMPORTANT: We need to call Update() after setting up NumberOfPieces and
    //            Piece correctly
    mapper->Update();

#ifndef SUPERNOVA
    std::cout << "Actual vtkImageData memory size on process "
	      << myId << ": "
	      << sampler->GetOutput()->GetActualMemorySize()
	      << " KiB" << std::endl;
#endif
    std::cout << "Number of polygon cells on process " << myId << ": "
	      << procId->GetOutput()->GetNumberOfCells() << endl;

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper.GetPointer());

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor.GetPointer());

    vtkNew<vtkRenderWindow> renWin;
    renWin->AddRenderer(renderer.GetPointer());
    renWin->SetSize(640, 480);

    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetRenderWindow(renWin.GetPointer());

    vtkNew<vtkCompositeRenderManager> renderManager;
    renderManager->SetController(controller);
    renderManager->SetRenderWindow(renWin.GetPointer());
    renderManager->InitializePieces();
    renderManager->InitializeRMIs();

    if (myId == 0)
        renderManager->ResetAllCameras();

    // Only the root process will have an active interactor. All
    // the other render windows will be slaved to the root.
    renderManager->StartInteractor();
#endif
}

int main(int argc, char* argv[])
{
    vtkNew<vtkMPIController> controller;
    controller->Initialize(&argc, &argv);

    // Execute the function named "process" on all processes
    controller->SetSingleMethod(process, 0);
    controller->SingleMethodExecute();

    // Clean-up and exit
    controller->Finalize();

    return 0;
}


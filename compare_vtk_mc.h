
#include <vtkm/cont/DeviceAdapter.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/DynamicArrayHandle.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/Pair.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkContourFilter.h>
#include <vtkFlyingEdges3D.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <vector>

namespace vtk
{
static void RunContourFilter(vtkImageData* image, int MAX_NUM_TRIALS)
{
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {

    vtkNew<vtkContourFilter> marching;
    marching->SetInputConnection(producer->GetOutputPort());

    vtkm::cont::Timer<> timer;

    marching->ComputeGradientsOff();
    marching->ComputeNormalsOff();
    marching->ComputeScalarsOn();
    marching->SetNumberOfContours(1);
    marching->SetValue(0, ISO_VALUE);

    marching->Update();

    double time = timer.GetElapsedTime();
    std::cout << "num cells: " << marching->GetOutput()->GetNumberOfCells() << std::endl;
    std::cout << "VTK Contour,Serial," << time << "," << i << std::endl;
    }
}

#ifdef VTK_HAS_FLYING_EDGES

static void RunFlyingEdges(vtkImageData* image, int MAX_NUM_TRIALS)
{
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {

    vtkNew<vtkFlyingEdges3D> marching;
    marching->SetInputConnection(producer->GetOutputPort());

    vtkm::cont::Timer<> timer;

    marching->ComputeGradientsOff();
    marching->ComputeNormalsOff();
    marching->ComputeScalarsOn();
    marching->SetNumberOfContours(1);
    marching->SetValue(0, ISO_VALUE);

    marching->Update();

    double time = timer.GetElapsedTime();
    std::cout << "num cells: " << marching->GetOutput()->GetNumberOfCells() << std::endl;
    std::cout << "VTK Flying Edges,Serial," << time << "," << i << std::endl;
    }
}

#endif

}
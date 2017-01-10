#include "gui.h"
#include <vtkPNGReader.h>
#include <vtkPlaneSource.h>
#include <vtkTextureMapToPlane.h>
#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkImageData.h>
#include <vtkImageExport.h>
#include <vtkCamera.h>
#include <vtkImageCanvasSource2D.h>
using namespace imr::gui;

Gui::Gui(vtkSmartPointer<vtkPNGReader> reader)
{
  RawPoints initial;
  RawPoints scan;
  Point pos;
//     vtkSmartPointer<vtkPNGReader> reader =  vtkSmartPointer<vtkPNGReader>::New();
//     reader->SetFileName("../data/belgioso-map5.png");

  vtkSmartPointer<vtkImageData> imageData;
     imageData = reader->GetOutput();
     reader->Update();
     int dims[3];
     imageData->GetDimensions (dims);
     std::cout << "DIMS: " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;



//      for(size_t r=0;r<dims[1];r++) {
//        for(size_t c=0;c<1000/*dims[1]*/;c++) {
//          std::cout << (*static_cast<unsigned char*>(imageData->GetScalarPointer(c, r, 0))==0 ? 'X' : ' ');
//        }
//        std::cout << std::endl;
//      }




    // Create a plane
//     double K = 43.36;
//     double Y = -42.86;
//     double X = -12.85;
//    plane->SetPoint1(X+1444.0/K/*width*/,Y,0.0);
//    plane->SetPoint2(X,Y+4265.0/K/*height*/,0.0);

    double K = 51.03;
    double X = -16.96;
    double Y = -43.25;


    vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
    plane->SetOrigin(X,Y,0);
    plane->SetCenter(0.0, 0.0, 0.0);
    plane->SetNormal(0.0, 0.0, 1.0);
    plane->SetPoint1(X+1872.0/K/*width*/,Y,0.0);
    plane->SetPoint2(X,Y+5015.0/K/*height*/,0.0);

    plane->SetResolution(1,1);



    // Apply the texture
    vtkSmartPointer<vtkTexture> texture = vtkSmartPointer<vtkTexture>::New();
    texture->SetInputConnection(reader->GetOutputPort());
    vtkSmartPointer<vtkTextureMapToPlane> texturePlane = vtkSmartPointer<vtkTextureMapToPlane>::New();
    texturePlane->SetInputConnection(plane->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    planeMapper->SetInputConnection(texturePlane->GetOutputPort());

    vtkSmartPointer<vtkActor> textureActor = vtkSmartPointer<vtkActor>::New();
    textureActor->SetMapper(planeMapper);
    textureActor->SetTexture(texture);




    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    RawPoints2vtkPoints(initial,points);
    map = vtkSmartPointer<vtkPolyData>::New();
    map->SetPoints(points);

    
    particles = vtkSmartPointer<vtkPolyData>::New();
    particles->SetPoints(points);
    
    
    points = vtkSmartPointer<vtkPoints>::New();
    RawPoints2vtkPoints(scan, points);
    scanMeasurement = vtkSmartPointer<vtkPolyData>::New();
    scanMeasurement->SetPoints(points);

    points = vtkSmartPointer<vtkPoints>::New();
    RawPoints2vtkPoints(pos,points);
    position = vtkSmartPointer<vtkPolyData>::New();
    position->SetPoints(points);

    points = vtkSmartPointer<vtkPoints>::New();
    RawPoints2vtkPoints(pos,points);
    icpPosition = vtkSmartPointer<vtkPolyData>::New();
    icpPosition->SetPoints(points);

    mapFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        mapFilter->AddInput(map);
    #else
        mapFilter->AddInputData(map);
    #endif
    mapFilter->Update();

    mapMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapMapper->SetInputConnection(mapFilter->GetOutputPort());

    mapActor = vtkSmartPointer<vtkActor>::New();
    mapActor->SetMapper(mapMapper);
    mapActor->GetProperty()->SetColor(0,1,1);
    mapActor->GetProperty()->SetPointSize(3);

    
    particlesFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        particlesFilter->AddInput(particles);
    #else
        particlesFilter->AddInputData(particles);
    #endif
    particlesFilter->Update();

    particlesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    particlesMapper->SetInputConnection(particlesFilter->GetOutputPort());

    particlesActor = vtkSmartPointer<vtkActor>::New();
    particlesActor->SetMapper(particlesMapper);
    particlesActor->GetProperty()->SetColor(1,0,0);
    particlesActor->GetProperty()->SetPointSize(3);
    
    
    
    
    scanFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        scanFilter->AddInput(scanMeasurement);
    #else
        scanFilter->AddInputData(scanMeasurement);
    #endif
    scanFilter->Update();

    scanMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    scanMapper->SetInputConnection(scanFilter->GetOutputPort());

    scanActor = vtkSmartPointer<vtkActor>::New();
    scanActor->SetMapper(scanMapper);
    scanActor->GetProperty()->SetColor(1,0,0);
    scanActor->GetProperty()->SetPointSize(3);

    positionFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        positionFilter->AddInput(position);
    #else
        positionFilter->AddInputData(position);
    #endif
    positionFilter->Update();

    positionMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    positionMapper->SetInputConnection(positionFilter->GetOutputPort());

    positionActor = vtkSmartPointer<vtkActor>::New();
    positionActor->SetMapper(positionMapper);
    positionActor->GetProperty()->SetColor(0,1,0);
    positionActor->GetProperty()->SetPointSize(5);

    icpPositionFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        icpPositionFilter->AddInput(icpPosition);
    #else
        icpPositionFilter->AddInputData(icpPosition);
    #endif
    icpPositionFilter->Update();

    icpPositionMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    icpPositionMapper->SetInputConnection(icpPositionFilter->GetOutputPort());

    icpPositionActor = vtkSmartPointer<vtkActor>::New();
    icpPositionActor->SetMapper(icpPositionMapper);
    icpPositionActor->GetProperty()->SetColor(0,0,1);
    icpPositionActor->GetProperty()->SetPointSize(5);


    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    interactorStyle = vtkSmartPointer<vtkInteractorStyleImage>::New();
    renderWindowInteractor->SetInteractorStyle(interactorStyle);

    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    interactorStyle = vtkSmartPointer<vtkInteractorStyleImage>::New();
    renderWindowInteractor->SetInteractorStyle(interactorStyle);


    renderer->AddActor(textureActor);
    renderer->AddActor(particlesActor);
    renderer->AddActor(mapActor);
    renderer->AddActor(scanActor);
    renderer->AddActor(positionActor);
//    renderer->AddActor(icpPositionActor);
    renderer->SetBackground(.5,0.5,0.5);
    renderWindow->Render();
}

void Gui::RawPoints2vtkPoints(const RawPoints &pts, vtkSmartPointer<vtkPoints> vtkPts)
{
    for(int i=0; i<pts.size(); i++) {
        vtkPts->InsertNextPoint(pts[i].x, pts[i].y, .0);
    }
}


void Gui::RawPoints2vtkPoints(const Point &pts, vtkSmartPointer<vtkPoints> vtkPts)
{
        vtkPts->InsertNextPoint(pts.x, pts.y, .0);
}

void Gui::Particles2vtkPoints(const ParticleVector &pts, vtkSmartPointer<vtkPoints> vtkPts)
{
    for(int i=0; i<pts.size(); i++) {
        vtkPts->InsertNextPoint(pts[i].pos.x, pts[i].pos.y, .0);
    }
}


void Gui::startInteractor()
{
    renderWindowInteractor->Start();
}

void Gui::setScanPoints(RawPoints pts)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    RawPoints2vtkPoints(pts,points);

    scanMeasurement->SetPoints(points);
    scanMeasurement->Modified();
    renderWindow->Render();
}

void Gui::setParticlePoints(ParticleVector pts)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    Particles2vtkPoints(pts,points);

    particles->SetPoints(points);
    particles->Modified();
    renderWindow->Render();
}

void Gui::clearScanPoints()
{
    scanMeasurement->SetPoints(vtkSmartPointer<vtkPoints>::New());
    scanMeasurement->Modified();
    renderWindow->Render();
}

void Gui::clearMapPoints()
{
    map->SetPoints(vtkSmartPointer<vtkPoints>::New());
//    position->SetPoints(vtkSmartPointer<vtkPoints>::New());
//    icpPosition->SetPoints(vtkSmartPointer<vtkPoints>::New());
    map->Modified();
//    position->Modified();
//    icpPosition->Modified();
    renderWindow->Render();
}


void Gui::clearPositionPoints()
{
    position->SetPoints(vtkSmartPointer<vtkPoints>::New());
    position->Modified();
    renderWindow->Render();
}

void Gui::setPosition(Point pos)
{
    vtkSmartPointer<vtkPoints> p = position->GetPoints();
    RawPoints2vtkPoints(pos,p);

    position->SetPoints(p);
    position->Modified();

    renderWindow->Render();
}


void Gui::setPointsToMap(RawPoints pts, Point pos)
{
    vtkSmartPointer<vtkPoints> p = position->GetPoints();
    RawPoints2vtkPoints(pos,p);

    position->SetPoints(p);
    position->Modified();


    vtkSmartPointer<vtkPoints> points = map->GetPoints();
    RawPoints2vtkPoints(pts,points);

    map->SetPoints(points);
    map->Modified();
    renderWindow->Render();
}

void Gui::screenshot(std::string filename) {
  renderer->ResetCamera();
  renderWindow->Render();

  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetMagnification(5); //set the resolution of the output image (3 times the current resolution of vtk render window)
  windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
  windowToImageFilter->ReadFrontBufferOff();
  windowToImageFilter->Update();

  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();

}

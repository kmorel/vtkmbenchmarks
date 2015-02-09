
#include <iostream>
#include <fstream>

static
void saveAsPly(std::vector< vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > >& vertices, std::string path)
{
  typedef vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > HandleType;
  typedef HandleType::PortalConstControl PortalConst;
  typedef HandleType::ValueType TriangleType;

  //1. count the number of total triangles
  vtkm::Id numVerts = 0;

  for(std::size_t i=0; i < vertices.size(); ++i)
      {
      numVerts += vertices[i].GetNumberOfValues();
      }

  //2. write out the file header
  std::fstream file(path.c_str(), std::fstream::out | std::fstream::trunc);
  file << "ply" << std::endl;
  file << "format ascii 1.0" << std::endl;
  file << "element vertex " << numVerts << std::endl;
  file << "property float32 x" << std::endl;
  file << "property float32 y" << std::endl;
  file << "property float32 z" << std::endl;
  file << "element face " << (numVerts / 3) << std::endl;
  file << "property list uint8 int32 vertex_indices" << std::endl;
  file << "end_header" << std::endl;

  //output the coordinates
  for(vtkm::Id i=0; i < vertices.size(); ++i)
    {
    vtkm::Id numTris = vertices[i].GetNumberOfValues();
    if( numTris > 0)
      {
      PortalConst portal = vertices[i].GetPortalConstControl();
      for(vtkm::Id cellId=0; cellId < numTris; ++cellId)
        {
        TriangleType tri = portal.Get(cellId);
        file << tri[0] << " " << tri[1] << " " << tri[2] << std::endl;
        }
      }
    }

  //output the connectivity
  vtkm::Id offset = 0;
  for(vtkm::Id i=0; i < vertices.size(); ++i)
    {
    vtkm::Id numVerts = vertices[i].GetNumberOfValues() * 3;
    if( numVerts > 0)
      {
      for(vtkm::Id vertIndex=offset; vertIndex < (offset+numVerts); vertIndex+=3)
        {
        file << "3 " << vertIndex << " " << vertIndex+1 << " " << vertIndex+2 << std::endl;
        }
      offset += numVerts;
      }
    }
  file.close();
}

static
void saveAsPly(vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> >& vertices, std::string path)
{
  typedef vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > HandleType;
  std::vector< HandleType > vec; vec.push_back(vertices);
  saveAsPly(vec, path);
}